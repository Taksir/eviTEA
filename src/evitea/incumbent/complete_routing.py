"""Optional complete routing with a partner-read gate (off by default; ``--complete-routing``).

The frozen first stage routes reads through a coverage-optimised subset of
exact TE anchor words.  This layer adds, for routing only, every other catalogue
anchor query that is not low-complexity, comes from exactly one family and has
an msBWT width in ``[1, COMPLETE_ROUTING_CAP]``.  It never changes a pair the
frozen panel already routes: only *extra* pairs are considered, i.e. pairs not
routed by the frozen panel whose extra-word hits fall on exactly one mate.

An extra pair is kept only if it is not an ordinary reference pair: the
TE-side mate (``-k 20``) has no alignment on the host-side mate's chromosome
within ``Q99`` bp of the host-side mate's primary alignment.  This is the
discordant-pair entry rule that TEMP2 and TEforest use.  Kept extra pairs join
the panel with read count 1, so the unchanged downstream code localises them
like any other one-mate pair.  Later stages keep the frozen TE library.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path

import numpy as np

from ..constants import Q99
from .anchor_census import AnchorCatalogue, AnchorQuery
from .external import run_bowtie2, write_fasta
from .pair_aware_panel import PairAwarePanel

__all__ = [
    "COMPLETE_ROUTING_CAP",
    "TE_MATE_OPTIONS",
    "HOST_MATE_OPTIONS",
    "extra_queries",
    "collect_extra_pairs",
    "concordant",
    "sam_records",
    "gate_extra_pairs",
    "augment_panel",
]

COMPLETE_ROUTING_CAP = 1024
HOST_MATE_OPTIONS = ("-k", "2", "--very-sensitive", "--no-hd")
TE_MATE_OPTIONS = ("-k", "20", "--very-sensitive", "--no-hd")


def _member(sorted_ref: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Boolean membership of ``values`` in the sorted unique array ``sorted_ref`` (O(len(values) log n))."""

    if len(sorted_ref) == 0 or len(values) == 0:
        return np.zeros(len(values), dtype=bool)
    idx = np.searchsorted(sorted_ref, values)
    idx[idx == len(sorted_ref)] = 0
    return sorted_ref[idx] == values


def extra_queries(catalogue: AnchorCatalogue, selected: Iterable[str], cap: int = COMPLETE_ROUTING_CAP) -> list[AnchorQuery]:
    """Catalogue queries outside the frozen panel that the complete routing adds (sorted by sequence)."""

    taken = set(selected)
    out = [
        q for q in catalogue.queries
        if q.query_seq not in taken and not q.low_complexity and q.source_family_count == 1 and 1 <= q.width <= cap
    ]
    out.sort(key=lambda q: q.query_seq)
    seen: set[str] = set()
    unique = []
    for q in out:
        if q.query_seq not in seen:
            seen.add(q.query_seq)
            unique.append(q)
    return unique


def collect_extra_pairs(
    queries: Sequence[AnchorQuery], frozen_pairs: Iterable[int], index, progress=None
) -> tuple[dict[int, int], dict[str, np.ndarray], dict[int, np.ndarray]]:
    """Extra one-mate pairs.

    Returns ``(te_file, pair_sets, te_dollars)``: ``te_file[pair] = file id of
    the only mate carrying extra words``; ``pair_sets[seq]`` = extra one-mate
    pairs carrying ``seq``; ``te_dollars[pair]`` = the dollar IDs (reads) of
    that mate that carry extra words (one read, kept as an array).
    """

    # one boolean mask over all pair IDs (a per-word np.isin against the frozen set costs ~1 s each)
    frozen_mask = np.zeros(int(index.pair_count), dtype=bool)
    frozen_list = np.asarray(sorted(set(int(x) for x in frozen_pairs)), dtype=np.int64)
    frozen_mask[frozen_list] = True
    per_seq: list[tuple[str, np.ndarray]] = []
    all_d, all_f, all_p = [], [], []
    for n_done, q in enumerate(queries, 1):
        if progress is not None and (n_done % 10000 == 0 or n_done == len(queries)):
            progress(f"complete routing words {n_done}/{len(queries)}")
        dollars = np.unique(np.asarray(index.dollar_ids(q.query_seq), dtype=np.uint64))
        files, pairs = index.file_and_pair(dollars)
        files = np.asarray(files, dtype=np.uint8)
        pairs = np.asarray(pairs, dtype=np.uint64)
        keep = ~frozen_mask[pairs.astype(np.int64)]
        per_seq.append((q.query_seq, np.unique(pairs[keep])))
        all_d.append(dollars[keep]); all_f.append(files[keep]); all_p.append(pairs[keep])
    if not all_d or sum(len(x) for x in all_d) == 0:
        return {}, {s: np.zeros(0, dtype=np.uint64) for s, _ in per_seq}, {}
    d = np.concatenate(all_d); f = np.concatenate(all_f).astype(np.uint64); p = np.concatenate(all_p)
    keys = np.unique(p * np.uint64(2) + f)  # distinct (pair, file) hits
    kp, kf = keys // np.uint64(2), keys % np.uint64(2)
    up, nfiles = np.unique(kp, return_counts=True)
    one = up[nfiles == 1]
    one_file = kf[_member(one, kp)]
    te_file = {int(a): int(b) for a, b in zip(one.tolist(), one_file.tolist())}
    # per-word membership by binary search (np.isin per word against ~370 k pairs cost ~0.15 s each)
    pair_sets = {s: ps[_member(one, ps)] for s, ps in per_seq}
    sel = _member(one, p)
    order = np.lexsort((d[sel], p[sel]))
    dp, dd = p[sel][order], d[sel][order]
    te_dollars: dict[int, np.ndarray] = {}
    if len(dp):
        cuts = np.flatnonzero(np.diff(dp)) + 1
        for chunk_p, chunk_d in zip(np.split(dp, cuts), np.split(dd, cuts)):
            te_dollars[int(chunk_p[0])] = np.unique(chunk_d)
    return te_file, pair_sets, te_dollars


def concordant(host_rows: Sequence[tuple[int, str, int]], te_rows: Sequence[tuple[int, str, int]], q99: int = Q99) -> bool:
    """True if the pair looks like an ordinary reference pair (then the gate drops it).

    ``host_rows``/``te_rows`` are ``(flag, chrom, pos)`` SAM rows of one read.
    The host-side mate's primary alignment is its first row that is mapped,
    not secondary (256) and not supplementary (2048).  An unmapped host mate
    never makes a pair concordant.
    """

    primary = next(((c, p) for f, c, p in host_rows if not (f & 4 or f & 256 or f & 2048)), None)
    if primary is None:
        return False
    chrom, pos = primary
    return any((not f & 4) and c == chrom and abs(p - pos) <= q99 for f, c, p in te_rows)


def sam_records(sam: Path) -> Iterator[tuple[str, list[tuple[int, str, int]]]]:
    """Stream a headerless SAM; yield ``(qname, rows)``; rows of one read are contiguous in Bowtie2 output."""

    current: str | None = None
    rows: list[tuple[int, str, int]] = []
    with Path(sam).open("r", encoding="ascii", errors="strict") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.split("\t", 5)
            if fields[0] != current:
                if current is not None:
                    yield current, rows
                current, rows = fields[0], []
            rows.append((int(fields[1]), fields[2], int(fields[3])))
    if current is not None:
        yield current, rows


def gate_extra_pairs(
    te_file: Mapping[int, int],
    index,
    *,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    threads: int,
    q99: int = Q99,
    progress=None,
) -> set[int]:
    """Kept extra pairs (not concordant reference pairs)."""

    pairs = sorted(te_file)
    if not pairs:
        return set()
    work_directory.mkdir(parents=True, exist_ok=True)
    te_d = index.mate_dollars({p: te_file[p] for p in pairs})
    host_d = index.mate_dollars({p: 1 - te_file[p] for p in pairs})
    sams = {}
    for name, dollars, options in (("gate_host_mates", host_d, HOST_MATE_OPTIONS), ("gate_te_mates", te_d, TE_MATE_OPTIONS)):
        fasta = work_directory / f"{name}.fa"
        sam = work_directory / f"{name}.sam"
        if progress is not None:
            progress(f"complete routing gate: mapping {name} ({len(pairs)} reads)")
        write_fasta(fasta, ((f"p{p}", index.recover(dollars[p])) for p in pairs))
        run_bowtie2(bowtie2, ["-x", str(host_index_prefix), "-f", "-U", str(fasta), "-S", str(sam), *options, "-p", str(max(1, threads))])
        sams[name] = sam
    host = {q: rows for q, rows in sam_records(sams["gate_host_mates"])}
    kept = set(pairs)
    for q, rows in sam_records(sams["gate_te_mates"]):
        if concordant(host.get(q, []), rows, q99):
            kept.discard(int(q[1:]))
    if progress is not None:
        progress(f"complete routing gate: kept {len(kept)} of {len(pairs)} extra pairs")
    return kept


def augment_panel(
    panel: PairAwarePanel,
    kept: Iterable[int],
    queries: Sequence[AnchorQuery],
    pair_sets: Mapping[str, np.ndarray],
    te_dollars: Mapping[int, np.ndarray],
) -> PairAwarePanel:
    """The frozen panel plus kept extra pairs (count 1); frozen pairs are untouched."""

    keep = np.asarray(sorted(set(int(x) for x in kept)), dtype=np.uint64)
    if len(keep) == 0:
        return panel
    keep_set = set(keep.tolist())
    selected = set(panel.selected_sequences)
    widths = dict(panel.sequence_widths)
    families = {s: set(v) for s, v in panel.sequence_families.items()}
    sets = dict(panel.pair_sets)
    for q in queries:
        ps = pair_sets.get(q.query_seq)
        if ps is None:
            continue
        ps = ps[_member(keep, ps)]
        if len(ps) == 0:
            continue
        selected.add(q.query_seq)
        widths[q.query_seq] = q.width
        families[q.query_seq] = set(q.source_families)
        sets[q.query_seq] = ps
    extra_reads = np.concatenate([te_dollars[p] for p in sorted(keep_set)]).astype(np.uint64)
    read_ids = np.union1d(panel.read_ids.astype(np.uint64), extra_reads)
    order = np.argsort(np.concatenate([panel.pair_ids.astype(np.uint64), keep]), kind="stable")
    pair_ids = np.concatenate([panel.pair_ids.astype(np.uint64), keep])[order]
    counts = np.concatenate([panel.pair_read_counts.astype(np.uint8), np.ones(len(keep), dtype=np.uint8)])[order]
    return PairAwarePanel(
        selected_sequences=selected,
        sequence_widths=widths,
        sequence_families=families,
        groups=list(panel.groups),
        read_ids=read_ids,
        pair_ids=pair_ids,
        pair_read_counts=counts,
        pair_sets=sets,
    )
