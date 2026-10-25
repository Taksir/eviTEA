"""Bowtie2 invocations, query naming, and SAM interpretation.

Bowtie2 seeds its per-read pseudo-random choices from the read name, so the
query-name convention changes which placements a multi-mapping read reports
under ``-k``.  The convention is therefore an explicit, recorded profile:

* ``pid``            - ``P{pair_id}``: the default convention, used for
  every sample.  Every other Bowtie2 replay uses it too (full-memory finite
  state, internal 21-mer, paired LOCAL reference mapping).
* ``aka-historical`` - ``p{pair_id}|h{host_mate}`` for single-mate host
  mapping.  **Reproduction only**: it reproduces the host read names of the
  AKA017 sample, so reproducing that sample's layers byte for byte requires
  it.  It is never a default.

The paired LOCAL reference mapping names both mates ``P{pair_id}`` in both
profiles.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path

from ..constants import Q99
from ..errors import InputContractError
from ..incumbent.external import run_bowtie2
from .codebook import CodeUnit, revcomp

__all__ = [
    "HOST_READ_NAME_PROFILES",
    "PRODUCTION_HOST_READ_NAMES",
    "REPRODUCTION_ONLY_HOST_READ_NAMES",
    "HOST_OPTIONS",
    "REFERENCE_OPTIONS",
    "reference_options",
    "host_query_name",
    "write_host_fasta",
    "write_pair_fastas",
    "map_host_mates",
    "map_reference_pairs",
    "parse_host_sam",
    "parse_reference_sam",
    "cigar_ops",
    "ref_span",
    "terminal_softclips",
    "family_occurrences",
]

PRODUCTION_HOST_READ_NAMES = "pid"
REPRODUCTION_ONLY_HOST_READ_NAMES = ("aka-historical",)
HOST_READ_NAME_PROFILES = (PRODUCTION_HOST_READ_NAMES, *REPRODUCTION_ONLY_HOST_READ_NAMES)
HOST_OPTIONS = ("-k", "2", "--very-sensitive", "--no-unal", "--no-hd", "--no-sq")
REFERENCE_OPTIONS = ("--very-sensitive-local", "--no-mixed", "--no-discordant", "-I", "0", "-X", str(Q99),
                     "-k", "20", "--no-unal", "--reorder", "--no-hd", "--no-sq")


def reference_options(fragment_q99: int = Q99) -> tuple[str, ...]:
    """``REFERENCE_OPTIONS`` with the sample's fragment-size 99th percentile as the maximum fragment length."""

    options = list(REFERENCE_OPTIONS)
    options[options.index("-X") + 1] = str(fragment_q99)
    return tuple(options)
_NAME = re.compile(r"[Pp](\d+)(?:\|h[12])?")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def host_query_name(pair_id: int, host_mate: int, profile: str) -> str:
    if profile == PRODUCTION_HOST_READ_NAMES:
        return f"P{pair_id}"
    if profile == "aka-historical":
        return f"p{pair_id}|h{host_mate}"
    raise InputContractError(f"unknown host read-name profile {profile!r}")


def _pair_id(qname: str) -> int:
    match = _NAME.fullmatch(qname)
    if not match:
        raise InputContractError(f"unexpected query name {qname!r}")
    return int(match.group(1))


def write_host_fasta(path: Path, host_mates: Mapping[int, int], reads: Mapping[int, tuple[str, str]], profile: str) -> None:
    with Path(path).open("w", encoding="ascii", newline="\n") as handle:
        for pid in sorted(host_mates):
            mate = host_mates[pid]
            handle.write(f">{host_query_name(pid, mate, profile)}\n{reads[pid][mate - 1]}\n")


def write_pair_fastas(path1: Path, path2: Path, pair_ids, reads: Mapping[int, tuple[str, str]]) -> None:
    with Path(path1).open("w", encoding="ascii", newline="\n") as a, Path(path2).open("w", encoding="ascii", newline="\n") as b:
        for pid in sorted(pair_ids):
            r1, r2 = reads[pid]
            a.write(f">P{pid}\n{r1}\n")
            b.write(f">P{pid}\n{r2}\n")


def map_host_mates(bowtie2: Path, index_prefix: Path, fasta: Path, sam: Path, threads: int) -> None:
    run_bowtie2(bowtie2, ["-x", str(index_prefix), "-f", "-U", str(fasta), *HOST_OPTIONS,
                          "-S", str(sam), "-p", str(max(1, threads))])


def map_reference_pairs(bowtie2: Path, index_prefix: Path, fasta1: Path, fasta2: Path, sam: Path, threads: int,
                        fragment_q99: int = Q99) -> None:
    run_bowtie2(bowtie2, ["-x", str(index_prefix), "-f", "-1", str(fasta1), "-2", str(fasta2),
                          *reference_options(fragment_q99),
                          "-S", str(sam), "-p", str(max(1, threads))])


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    if cigar == "*":
        return []
    ops = [(int(n), op) for n, op in _CIGAR.findall(cigar)]
    if "".join(f"{n}{op}" for n, op in ops) != cigar:
        raise InputContractError(f"malformed CIGAR {cigar!r}")
    return ops


def ref_span(cigar: str) -> int:
    return sum(n for n, op in cigar_ops(cigar) if op in "MDN=X")


def _aligned_query_mask(cigar: str, qlen: int) -> list[bool]:
    mask = [False] * qlen
    q = 0
    for n, op in cigar_ops(cigar):
        if op in "M=X":
            for i in range(q, min(q + n, qlen)):
                mask[i] = True
            q += n
        elif op in "IS":
            q += n
    if q != qlen:
        raise InputContractError(f"CIGAR/query length mismatch: {cigar} consumes {q}, read {qlen}")
    return mask


def terminal_softclips(cigar: str, qlen: int) -> tuple[tuple[int, int] | None, tuple[int, int] | None]:
    """Leading / trailing query soft clips; ``S`` must be the terminal operation."""

    ops = cigar_ops(cigar)
    if sum(n for n, op in ops if op in "MIS=X") != qlen:
        raise InputContractError(f"CIGAR/query length mismatch: {cigar} vs {qlen}")
    left = (0, ops[0][0]) if ops and ops[0][1] == "S" and ops[0][0] > 0 else None
    right = (qlen - ops[-1][0], qlen) if ops and ops[-1][1] == "S" and ops[-1][0] > 0 else None
    return left, right


def family_occurrences(seq: str, units: list[CodeUnit]) -> list[tuple[int, int, int]]:
    seq = seq.upper()
    occ: set[tuple[int, int, int]] = set()
    for unit in units:
        for query in unit.sequences:
            start = 0
            while True:
                p = seq.find(query, start)
                if p < 0:
                    break
                occ.add((p, p + len(query), unit.unit_id))
                start = p + 1
    return sorted(occ)


def parse_host_sam(path: Path) -> dict[int, dict]:
    """``UNIQUE_HIGH`` = exactly one reported placement with MAPQ >= 20."""

    by: dict[int, list[dict]] = defaultdict(list)
    with Path(path).open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 11:
                continue
            flag = int(x[1])
            if flag & 0x4:
                continue
            pos0 = int(x[3]) - 1
            by[_pair_id(x[0])].append({
                "chrom": x[2], "start": pos0, "end": pos0 + ref_span(x[5]),
                "strand": "-" if flag & 0x10 else "+", "mapq": int(x[4]), "cigar": x[5],
            })
    out: dict[int, dict] = {}
    for pid, alns in by.items():
        if len(alns) == 1 and alns[0]["mapq"] >= 20:
            out[pid] = {**alns[0], "host_class": "UNIQUE_HIGH", "reported_alignments": 1}
        elif alns:
            best = max(alns, key=lambda a: a["mapq"])
            out[pid] = {**best, "host_class": "MULTI_OR_LOW", "reported_alignments": len(alns)}
    return out


def parse_reference_sam(
    path: Path,
    eligible: Mapping[int, dict],
    reads: Mapping[int, tuple[str, str]],
    by_family: Mapping[str, list[CodeUnit]],
) -> tuple[dict[int, dict], dict[int, list[dict]]]:
    """Evidence-complete reference table over LOCAL proper placement signatures.

    A pair is reference-explained if any proper placement signature (both mates
    present) has a TE-mate record whose every candidate-family safe-probe
    occurrence lies in ``M``/``=``/``X`` query bases.
    """

    grouped: dict[int, list[list[str]]] = defaultdict(list)
    with Path(path).open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) < 11:
                continue
            pid = _pair_id(x[0])
            if pid in eligible:
                grouped[pid].append(x)
    ledger: dict[int, dict] = {}
    te_records: dict[int, list[dict]] = defaultdict(list)
    for pid, meta in eligible.items():
        sigs: dict[tuple, dict[int, list[dict]]] = defaultdict(lambda: {1: [], 2: []})
        family = meta["family"]
        for x in grouped.get(pid, []):
            flag = int(x[1])
            rname = x[2]
            if rname == "*" or not (flag & 0x2) or (flag & 0x4) or (flag & 0x8) or (flag & 0x800):
                continue
            mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
            if mate == 0:
                continue
            pos, pnext, tlen = int(x[3]), int(x[7]), abs(int(x[8]))
            sig = (rname, min(pos, pnext), max(pos, pnext), tlen)
            seq = x[9]
            if seq == "*":
                seq = reads[pid][mate - 1]
                if flag & 0x10:
                    seq = revcomp(seq)
            seq = seq.upper()
            rec = {"flag": flag, "chrom": rname, "start": pos - 1, "end": pos - 1 + ref_span(x[5]),
                   "mapq": int(x[4]), "cigar": x[5], "seq": seq, "strand": "-" if flag & 0x10 else "+"}
            if mate == meta["te_mate"]:
                occ = family_occurrences(seq, by_family[family])
                mask = _aligned_query_mask(x[5], len(seq))
                bad = sum(1 for a, b, _ in occ if a < 0 or b > len(mask) or not all(mask[a:b]))
                rec.update({"probe_occurrences": len(occ), "unexplained_occurrences": bad,
                            "evidence_complete": bool(occ) and bad == 0})
            sigs[sig][mate].append(rec)
        proper = clean = max_occ = max_unexplained = 0
        records: list[dict] = []
        for d in sigs.values():
            if not d[1] or not d[2]:
                continue
            proper += 1
            trs = d[meta["te_mate"]]
            records.extend(trs)
            max_occ = max(max_occ, max((int(r.get("probe_occurrences", 0)) for r in trs), default=0))
            max_unexplained = max(max_unexplained, max((int(r.get("unexplained_occurrences", 0)) for r in trs), default=0))
            clean += int(any(r.get("evidence_complete") for r in trs))
        te_records[pid] = records
        state = ("EVIDENCE_COMPLETE_CLEAN_REFERENCE" if clean else
                 "LOCAL_PROPER_TE_UNEXPLAINED" if proper else "NO_LOCAL_PROPER")
        ledger[pid] = {**meta, "reference_state": state, "proper_signature_count": proper,
                       "evidence_complete_signature_count": clean, "evidence_complete_placement": bool(clean),
                       "candidate_family_occurrences_max": max_occ,
                       "max_unexplained_probe_occurrences": max_unexplained}
    return ledger, te_records
