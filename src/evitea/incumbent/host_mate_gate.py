"""Host-mate localization gate and resident pair families (stages G1F/G1G).

Pairs with exactly one TE-positive mate send the *other* mate to Bowtie2
(``-k 2 --very-sensitive``, query ``p{pair}|h{mate}``).  A mate is
UNIQUE_HIGH when exactly one alignment is reported with MAPQ >= 20; its host
position is the primary alignment midpoint ``POS + (ref_span - 1) // 2``.
UNIQUE_HIGH midpoints are single-linkage clustered per chromosome (gap 500 bp)
and clusters with >= 2 pairs are kept.  The working candidate pairs are those
cluster pairs plus every sparse-panel pair hit through a terminal anchor.
Cluster pairs inherit the families of every selected panel string they carry.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..errors import BackendCompatibilityError
from .external import run_bowtie2, write_fasta
from .pair_aware_panel import PairAwarePanel
from .sparse_panel import CHANNEL_TERMINAL, PanelPair

__all__ = [
    "HOST_CLUSTER_GAP",
    "HOST_CLUSTER_SUPPORT",
    "UNIQUE_MAPQ",
    "HostMapping",
    "one_mate_host_requests",
    "run_host_mapping",
    "parse_primary_host_sam",
    "supported_host_cluster_pairs",
    "working_candidate_pairs",
    "resident_pair_families",
]

HOST_CLUSTER_GAP = 500
HOST_CLUSTER_SUPPORT = 2
UNIQUE_MAPQ = 20
PRIMARY_HOST_OPTIONS = ("-k", "2", "--very-sensitive", "--no-unal", "--no-hd", "--no-sq")


@dataclass(frozen=True, slots=True)
class HostMapping:
    pair_id: int
    status: str
    alignment_count: int
    chrom: str
    pos: int | None
    midpoint: int | None
    mapq: int | None
    strand: str
    cigar: str


def one_mate_host_requests(panel: PairAwarePanel, index) -> dict[int, int]:
    """``pair_id -> host file id`` for routed pairs with one TE-positive mate."""

    files, pairs = index.file_and_pair(panel.read_ids)
    one = set(panel.pair_ids[panel.pair_read_counts == 1].tolist())
    requests: dict[int, int] = {}
    for file_id, pair_id in zip(files.tolist(), pairs.tolist()):
        if pair_id in one:
            requests[pair_id] = 1 - file_id
    if len(requests) != len(one):
        raise BackendCompatibilityError("failed to recover the TE-positive mate of one-mate pairs")
    return requests


def run_host_mapping(
    records: Sequence[tuple[str, str]],
    *,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    name: str,
    options: Sequence[str],
    threads: int,
) -> Path:
    work_directory.mkdir(parents=True, exist_ok=True)
    fasta = work_directory / f"{name}.fa"
    sam = work_directory / f"{name}.sam"
    write_fasta(fasta, records)
    run_bowtie2(
        bowtie2,
        ["-x", str(host_index_prefix), "-f", "-U", str(fasta), "-S", str(sam), *options, "-p", str(max(1, threads))],
    )
    return sam


def _ref_span(cigar: str) -> int:
    if cigar == "*" or not cigar:
        return 0
    total = 0
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        if not number:
            raise BackendCompatibilityError(f"malformed CIGAR {cigar}")
        if char in "MDN=X":
            total += int(number)
        number = ""
    if number:
        raise BackendCompatibilityError(f"malformed CIGAR {cigar}")
    return total


def parse_primary_host_sam(sam: Path, expected: Iterable[int]) -> dict[int, HostMapping]:
    expected_sorted = sorted(int(x) for x in expected)
    wanted = set(expected_sorted)
    counts: dict[int, int] = defaultdict(int)
    primary: dict[int, tuple] = {}
    with Path(sam).open("r", encoding="ascii", errors="strict") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise BackendCompatibilityError("malformed SAM line")
            qname = fields[0]
            if not qname.startswith("p"):
                raise BackendCompatibilityError(f"unexpected host query {qname}")
            pair_id = int(qname.split("|", 1)[0][1:])
            if pair_id not in wanted:
                raise BackendCompatibilityError(f"unexpected pair in SAM: {pair_id}")
            flag = int(fields[1])
            if flag & 4:
                continue
            counts[pair_id] += 1
            if not flag & 256 and not flag & 2048:
                pos = int(fields[3])
                span = _ref_span(fields[5])
                primary[pair_id] = (fields[2], pos, pos + max(0, span - 1) // 2, int(fields[4]), "-" if flag & 16 else "+", fields[5])
    out: dict[int, HostMapping] = {}
    for pair_id in expected_sorted:
        n = counts.get(pair_id, 0)
        record = primary.get(pair_id)
        if n == 0 or record is None:
            out[pair_id] = HostMapping(pair_id, "UNMAPPED", 0, "", None, None, None, "", "")
            continue
        chrom, pos, midpoint, mapq, strand, cigar = record
        if n == 1 and mapq >= UNIQUE_MAPQ:
            status = "UNIQUE_HIGH"
        elif n >= 2:
            status = "MULTI"
        else:
            status = "MAPPED_LOWCONF"
        out[pair_id] = HostMapping(pair_id, status, n, chrom, pos, midpoint, mapq, strand, cigar)
    return out


def supported_host_cluster_pairs(
    mapping: Mapping[int, HostMapping], gap: int = HOST_CLUSTER_GAP, support: int = HOST_CLUSTER_SUPPORT
) -> set[int]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in mapping.values():
        if row.status == "UNIQUE_HIGH":
            by_chrom[row.chrom].append((int(row.midpoint), row.pair_id))
    kept: set[int] = set()
    for chrom in sorted(by_chrom):
        values = sorted(by_chrom[chrom])
        current: list[tuple[int, int]] = []
        last = None
        for value in values + [None]:
            if value is not None and (last is None or value[0] - last <= gap):
                current.append(value)
                last = value[0]
                continue
            pairs = {pid for _, pid in current}
            if len(pairs) >= support:
                kept.update(pairs)
            if value is not None:
                current = [value]
                last = value[0]
    return kept


def working_candidate_pairs(host_pairs: set[int], panel_pairs: Sequence[PanelPair]) -> np.ndarray:
    terminal = {pair.pair_id for pair in panel_pairs if pair.channel_mask & CHANNEL_TERMINAL}
    return np.asarray(sorted(host_pairs | terminal), dtype=np.uint64)


def resident_pair_families(host_pairs: set[int], panel: PairAwarePanel) -> dict[int, str]:
    target = np.zeros(max(host_pairs, default=-1) + 1, dtype=bool)
    if host_pairs:
        target[np.fromiter(host_pairs, dtype=np.int64)] = True
    families: dict[int, set[str]] = defaultdict(set)
    for seq in sorted(panel.selected_sequences):
        names = panel.sequence_families.get(seq, set())
        if not names:
            continue
        pairs = panel.pair_sets[seq].astype(np.int64)
        pairs = pairs[pairs < len(target)]
        for pair_id in pairs[target[pairs]].tolist():
            families[pair_id].update(names)
    return {pid: ",".join(sorted(names)) for pid, names in sorted(families.items())}
