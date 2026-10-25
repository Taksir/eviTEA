"""Detector driver: one contig at a time, then one grouping over all contigs.

Port of the research driver (``p4_genome.py piece`` / ``merge``, built
from ``p4_run.py``) without its evaluation-only targeted mode. The
per-contig scan and the grouping call the vendored ``_p4_lib`` exactly as the
research driver does, so ``places_tsv_bytes`` reproduces ``P4_PLACES.tsv`` and
``segments_fa_bytes`` reproduces ``P4_SEGMENTS.fa`` (tests/unit/test_library_free.py).
"""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence

from . import _p4_lib as P

__all__ = [
    "COLS",
    "group_and_call",
    "parse_places",
    "parse_segments",
    "places_tsv_bytes",
    "reference_multi",
    "scan_contig",
    "segments_fa_bytes",
]

COLS = ["chrom", "first", "last", "n_low", "D", "stage", "why", "bL", "bR", "delta", "tsd", "two_sided", "geno", "fracL",
        "fracR", "rep", "lstop", "rstop", "rstat", "seg_keep", "seg_reject", "dup", "cluster", "track", "fam_pair",
        "fam_vtsd", "target"]


def reference_multi(genome: Mapping[str, str]) -> P.RefMulti:
    """Multi-copy canonical 21- and 31-letter words of the whole reference (``p4_run.py refmulti``)."""

    return P.RefMulti(P.multi_codes(genome, P.REP_K), P.multi_codes(genome, P.K))


def scan_contig(counter, sequence: str, chrom: str, rm: P.RefMulti, firewall) -> list[dict]:
    """S1-S4 for every candidate of one contig, in candidate order (``p4_genome.py piece`` without targets).

    ``counter`` is a ``_p4_lib.Counter`` over the sample msBWT built with ``cache=False``. As in the research driver, its
    cache is switched on after the S1 scan and cleared every 2,000,000 entries.
    """

    ref = sequence.upper()
    pos, c = P.scan_counts(counter, ref)
    depth, nfrac = P.local_depth(c)
    cands = P.candidates(pos, c, depth, nfrac)
    counter.use = True
    places = []
    for cand in cands:
        places.append(P.process(counter, ref, cand, rm, firewall, chrom, fork=P.FORK))
        if len(getattr(counter, "cache", ())) > 2_000_000:
            counter.clear()
    return places


def group_and_call(places: list[dict]) -> list[dict]:
    """S5-S7 over every S4_pass place (``p4_run.group_and_call``); adds cluster, track, fam_pair, fam_vtsd."""

    P.dedupe(places)
    gp = [p for p in places if p["stage"] == "S4_pass" and not p.get("dup")]
    for p in gp:
        keep = p["seg_keep"].split(",")
        p["segs"] = [s for s, side in ((p.get("lseg", ""), "L"), (p.get("rseg", ""), "R")) if side in keep]
    cl, _ = P.group(gp)
    fam = P.families(gp, cl)
    for p, c in zip(gp, cl):
        p["cluster"] = c
        p["track"] = P.classify_place(p, fam[c])
        p["fam_pair"] = fam[c]["pair"]
        p["fam_vtsd"] = fam[c]["vtsd"]
    return gp


def places_tsv_bytes(places: Sequence[dict]) -> bytes:
    """``P4_PLACES.tsv``: every place except S2_through, in order (``p4_run.write_places``)."""

    out = io.StringIO(newline="")
    out.write("\t".join(COLS + ["lseg_len", "rseg_len"]) + "\n")
    for p in places:
        if p["stage"] == "S2_through":
            continue
        out.write("\t".join("" if p.get(c) is None else str(p.get(c)) for c in COLS)
                  + f"\t{len(p.get('lseg', ''))}\t{len(p.get('rseg', ''))}\n")
    return out.getvalue().encode("utf-8")


def segments_fa_bytes(grouped: Sequence[dict]) -> bytes:
    """``P4_SEGMENTS.fa``: the inserted-DNA segments of the grouped places (``p4_run.write_segments``)."""

    out = io.StringIO(newline="")
    for p in grouped:
        for side in ("lseg", "rseg"):
            s = p.get(side, "")
            if s:
                out.write(f">{p['chrom']}:{p.get('bL', p['first'])}:{side}:cl{p.get('cluster', '')}\n{s}\n")
    return out.getvalue().encode("utf-8")


def parse_places(data: bytes) -> list[dict[str, str]]:
    """Return the rows of a places table as text fields (as the research scorers read them)."""

    lines = data.decode("utf-8-sig").splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:] if line]


def parse_segments(data: bytes) -> dict[tuple[str, str, str], dict[str, str]]:
    """(chrom, bL, cluster) -> {'lseg' | 'rseg': sequence} from a segments FASTA."""

    out: dict[tuple[str, str, str], dict[str, str]] = {}
    name = None
    for line in data.decode("utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
        elif name is not None:
            chrom, bl, side, cl = name.split(":")
            entry = out.setdefault((chrom, bl, cl[2:]), {})
            if side in entry:
                raise ValueError(f"duplicate segment {name}")
            entry[side] = line
    return out
