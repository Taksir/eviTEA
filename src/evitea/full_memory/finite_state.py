"""Finite host / reference alternatives and the per-pair accounting.

Every safe-hit physical pair gets one accounting row with an explicit
experiment state.  For the one-mate family-resolved (``ELIGIBLE``) pairs:

* host: all distinct placements reported under Bowtie2 ``-k 2`` (deduplicated
  by (chrom, start, end, strand), ranked primary-first then MAPQ); two
  reported alignments hit the reporting cap and are *censored*, never read as
  "no third placement";
* reference: one row per proper LOCAL placement signature under ``-k 20``;
  a signature is evidence-complete when some TE-mate record has every
  candidate-family safe-probe occurrence in ``M``/``=``/``X`` bases.  Any
  complete signature -> ``CLEAN_REFERENCE``; else a hit cap (20 signatures or
  a 40-record group) -> ``CENSORED``; else ``REFERENCE_UNEXPLAINED``.

Finally ``ELIGIBLE`` pairs with no host placement become ``UNMAPPED``, and a
single low-MAPQ placement ``HOST_UNRESOLVED``.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

from ..constants import N_FOLDS
from ..errors import InputContractError
from ..anchor_pe.codebook import SafeCodebook

__all__ = [
    "HOST_K", "REF_K", "ACCOUNTING_FIELDS", "HOST_FIELDS", "REFERENCE_FIELDS",
    "family_probes", "finite_host_alternatives", "finite_reference_alternatives", "pair_accounting",
]

HOST_K = 2
REF_K = 20
ACCOUNTING_FIELDS = (
    "pair_id", "fold", "family", "families", "te_mate", "host_mate", "signature_state", "experiment_state",
    "exclusion_reason", "reference_state", "proper_signature_count", "complete_signature_count",
    "reference_cap_censored", "host_alt_count", "host_reported_alignment_count", "host_cap_censored", "host_class",
)
HOST_FIELDS = ("pair_id", "placement_rank", "cap_censored", "chrom", "start", "end", "midpoint", "strand", "mapq", "cigar", "secondary")
REFERENCE_FIELDS = (
    "pair_id", "family", "placement_rank", "reference_cap_censored", "signature_id", "chrom", "outer_start1",
    "outer_end1", "tlen", "evidence_complete", "te_start", "te_end", "te_strand", "te_mapq", "te_cigar",
    "host_start", "host_end", "host_strand", "host_mapq", "probe_occurrences", "unexplained_probe_occurrences",
    "explained_query_intervals", "unexplained_query_intervals",
)
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_QNAME = re.compile(r"(?:^P|^p)(\d+)")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _pid(qname: str) -> int:
    m = _QNAME.search(qname)
    if not m:
        raise InputContractError(f"unrecognized pair query name {qname!r}")
    return int(m.group(1))


def _ref_span(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR.findall(str(cigar or "")) if op in "MDN=X")


def _mask(cigar: str, qlen: int) -> list[bool]:
    mask, q = [False] * qlen, 0
    tokens = _CIGAR.findall(str(cigar or ""))
    if not tokens:
        return mask
    for n, op in tokens:
        n = int(n)
        if op in "M=X":
            for i in range(q, min(qlen, q + n)):
                mask[i] = True
            q += n
        elif op in "IS":
            q += n
    if q != qlen:
        raise InputContractError(f"CIGAR/query length mismatch {cigar} -> {q} != {qlen}")
    return mask


def family_probes(codebook: SafeCodebook) -> dict[str, tuple[str, ...]]:
    by: dict[str, set[str]] = defaultdict(set)
    for unit in codebook.units:
        if not unit.family:
            continue
        for seq in unit.sequences:
            if seq and not set(seq) - set("ACGT"):
                by[unit.family].add(seq)
    return {k: tuple(sorted(v, key=lambda x: (len(x), x))) for k, v in by.items()}


def _occurrences(seq: str, probes: tuple[str, ...]) -> list[tuple[int, int]]:
    out, seq = [], seq.upper()
    for probe in probes:
        i = seq.find(probe)
        while i >= 0:
            out.append((i, i + len(probe)))
            i = seq.find(probe, i + 1)
    return out


def finite_host_alternatives(sam: Path, eligible: Iterable[int]) -> tuple[list[dict], dict[int, dict]]:
    wanted = set(eligible)
    by: dict[int, list[dict]] = defaultdict(list)
    with Path(sam).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            x = line.rstrip("\r\n").split("\t")
            if len(x) < 11:
                continue
            flag = int(x[1])
            if flag & 0x4 or flag & 0x800:
                continue
            pid = _pid(x[0])
            if pid not in wanted:
                continue
            pos, span = int(x[3]), _ref_span(x[5])
            if pos <= 0 or span <= 0 or x[2] == "*":
                continue
            by[pid].append({"pair_id": pid, "chrom": x[2], "start": pos - 1, "end": pos - 1 + span,
                            "strand": "-" if flag & 0x10 else "+", "mapq": int(x[4]), "cigar": x[5],
                            "secondary": int(bool(flag & 0x100))})
    rows: list[dict] = []
    summary: dict[int, dict] = {}
    for pid in sorted(wanted):
        reported = by.get(pid, [])
        seen: set[tuple] = set()
        unique = []
        for r in sorted(reported, key=lambda z: (z["secondary"], -z["mapq"], z["chrom"], z["start"], z["end"], z["strand"])):
            key = (r["chrom"], r["start"], r["end"], r["strand"])
            if key not in seen:
                seen.add(key)
                unique.append(r)
        censored = int(len(reported) >= HOST_K or len(unique) >= HOST_K)
        for rank, r in enumerate(unique, 1):
            rows.append({**r, "placement_rank": rank, "cap_censored": censored, "midpoint": (r["start"] + r["end"]) // 2})
        summary[pid] = {
            "host_alt_count": len(unique), "host_reported_alignment_count": len(reported), "host_cap_censored": censored,
            "host_class": "UNMAPPED" if not unique else
            "UNIQUE_HIGH" if len(unique) == 1 and unique[0]["mapq"] >= 20 else "MULTI_OR_LOW",
        }
    return rows, summary


def _signature_rows(records: list[list[str]], meta: Mapping, probes: tuple[str, ...],
                    reads: Mapping[int, tuple[str, str]] | None) -> list[dict]:
    sigs: dict[tuple, dict[int, list[dict]]] = defaultdict(lambda: {1: [], 2: []})
    pid, te_mate = meta["pair_id"], meta["te_mate"]
    for x in records:
        if len(x) < 11:
            continue
        flag, rname = int(x[1]), x[2]
        if rname == "*" or not flag & 0x2 or flag & 0x4 or flag & 0x8 or flag & 0x800:
            continue
        mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
        if mate == 0:
            continue
        pos, pnext, tlen = int(x[3]), int(x[7]), abs(int(x[8] or 0))
        sig = (rname, min(pos, pnext), max(pos, pnext), tlen)
        seq = x[9]
        if seq == "*" and reads is not None:
            seq = reads[pid][mate - 1]
            if flag & 0x10:
                seq = seq.translate(_COMP)[::-1]
        rec = {"chrom": rname, "start": pos - 1, "end": pos - 1 + _ref_span(x[5]), "strand": "-" if flag & 0x10 else "+",
               "mapq": int(x[4]), "cigar": x[5]}
        if mate == te_mate:
            occ, mask = _occurrences(seq, probes), _mask(x[5], len(seq))
            unexplained = [(a, b) for a, b in occ if not all(mask[a:b])]
            rec.update({
                "probe_occurrences": len(occ), "unexplained_occurrences": len(unexplained),
                "complete": bool(occ) and not unexplained,
                "explained_intervals": ";".join(f"{a}-{b}" for a, b in occ if all(mask[a:b])),
                "unexplained_intervals": ";".join(f"{a}-{b}" for a, b in unexplained),
            })
        sigs[sig][mate].append(rec)
    out = []
    for sig, d in sigs.items():
        if not d[1] or not d[2]:
            continue
        te = d[te_mate]
        tr = max(te, key=lambda r: (int(r.get("complete", False)), r["mapq"]))
        hr = max(d[3 - te_mate], key=lambda r: r["mapq"])
        out.append({
            "pair_id": pid, "family": meta["family"], "signature_id": f"{sig[0]}:{sig[1]}:{sig[2]}:{sig[3]}",
            "chrom": sig[0], "outer_start1": sig[1], "outer_end1": sig[2], "tlen": sig[3],
            "evidence_complete": int(any(r.get("complete") for r in te)), "te_start": tr["start"], "te_end": tr["end"],
            "te_strand": tr["strand"], "te_mapq": tr["mapq"], "te_cigar": tr["cigar"], "host_start": hr["start"],
            "host_end": hr["end"], "host_strand": hr["strand"], "host_mapq": hr["mapq"],
            "probe_occurrences": tr.get("probe_occurrences", 0),
            "unexplained_probe_occurrences": tr.get("unexplained_occurrences", 0),
            "explained_query_intervals": tr.get("explained_intervals", ""),
            "unexplained_query_intervals": tr.get("unexplained_intervals", ""),
        })
    return out


def finite_reference_alternatives(
    sam: Path, meta: Mapping[int, Mapping], probes: Mapping[str, tuple[str, ...]],
    reads: Mapping[int, tuple[str, str]] | None,
    *,
    pair_probes: Mapping[int, tuple[str, ...]] | None = None,
) -> tuple[list[dict], dict[int, dict]]:
    """Stream a ``--reorder`` SAM grouped by consecutive query name.

    Evidence completeness is judged on the pair's family probes, or on
    ``pair_probes[pair_id]`` when given (INTERNAL21: the pair's own observed
    probes, both strands).
    """

    summary = {pid: {"proper_signature_count": 0, "complete_signature_count": 0,
                     "reference_state": "REFERENCE_UNEXPLAINED", "reference_cap_censored": 0} for pid in meta}
    rows: list[dict] = []

    def flush(qname, group):
        if qname is None:
            return
        pid = _pid(qname)
        if pid not in meta:
            return
        own = pair_probes[pid] if pair_probes is not None else probes.get(meta[pid]["family"], ())
        xs = _signature_rows(group, meta[pid], own, reads)
        xs.sort(key=lambda r: (r["chrom"], r["outer_start1"], r["outer_end1"], r["tlen"]))
        complete = sum(int(r["evidence_complete"]) for r in xs)
        censored = int(len(xs) >= REF_K or len(group) >= 2 * REF_K)
        state = "CLEAN_REFERENCE" if complete else ("CENSORED" if censored else "REFERENCE_UNEXPLAINED")
        summary[pid] = {"proper_signature_count": len(xs), "complete_signature_count": complete,
                        "reference_state": state, "reference_cap_censored": censored}
        for rank, r in enumerate(xs, 1):
            rows.append({**r, "placement_rank": rank, "reference_cap_censored": censored})

    current, buffer = None, []
    with Path(sam).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            x = line.rstrip("\r\n").split("\t")
            if current is None:
                current = x[0]
            if x[0] != current:
                flush(current, buffer)
                current, buffer = x[0], []
            buffer.append(x)
    flush(current, buffer)
    return rows, summary


def pair_accounting(signature_rows: Iterable[Mapping], host_summary: Mapping[int, dict],
                    reference_summary: Mapping[int, dict]) -> list[dict]:
    """One row per safe-hit pair (``anchor_pe.codebook.classify_pairs`` rows)."""

    out = []
    for row in sorted(signature_rows, key=lambda r: int(r["pair_id"])):
        pid = int(row["pair_id"])
        families = [f for f in str(row.get("families", "")).split(";") if f]
        family = str(row.get("family", "") or "")
        if not family and len(families) == 1:
            family = families[0]
        te_mate, host_mate = int(row.get("te_mate", 0) or 0), int(row.get("host_mate", 0) or 0)
        state = str(row["eligibility_state"]).upper()
        if len(families) > 1 or "FAMILY_CONFLICT" in state:
            exp, reason = "TE_SOURCE_UNRESOLVED", "FAMILY_CONFLICT"
        elif "BOTH_MATES" in state:
            exp, reason = "TE_SOURCE_UNRESOLVED", "BOTH_MATES_TE_POSITIVE"
        elif family and te_mate in (1, 2) and host_mate in (1, 2):
            exp, reason = "ELIGIBLE", ""
        else:
            exp, reason = "TE_SOURCE_UNRESOLVED", "FROZEN_ELIGIBILITY_ROLE_UNRESOLVED"
        info = {"pair_id": pid, "fold": pid % N_FOLDS, "family": family, "families": ";".join(families),
                "te_mate": te_mate, "host_mate": host_mate, "signature_state": state, "experiment_state": exp,
                "exclusion_reason": reason}
        if exp == "ELIGIBLE":
            info.update(reference_summary[pid])
            info.update(host_summary[pid])
            if info["host_alt_count"] == 0:
                info["experiment_state"], info["exclusion_reason"] = "UNMAPPED", "HOST_UNMAPPED"
            elif info["host_class"] == "MULTI_OR_LOW" and info["host_alt_count"] < 2:
                info["experiment_state"], info["exclusion_reason"] = "HOST_UNRESOLVED", "LOW_CONFIDENCE_SINGLE_HOST_PLACEMENT"
        else:
            info.update({"reference_state": "NOT_APPLICABLE", "proper_signature_count": 0,
                         "complete_signature_count": 0, "reference_cap_censored": 0})
        out.append(info)
    return out
