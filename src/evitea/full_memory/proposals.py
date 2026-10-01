"""Non-certifying terminal proposal memory (frozen safe-pair candidate memory).

Population: uncensored ``ELIGIBLE`` + ``REFERENCE_UNEXPLAINED`` +
``UNIQUE_HIGH`` pairs with a family, placed by their rank-1 finite host
alternative.  Each contributes the q99 interval ``[end, end + 805]`` (``+``) or
``[max(0, start - 805), start]`` (``-``) — not clipped at contig ends.  On the
primary chromosomes only, a sweep keeps every elementary segment covered by
>= 2 pairs over >= 2 folds; adjacent/overlapping segments merge into regions.
A rotation's regions come from its three discovery folds and are *confirmed*
when each held-out fold has >= 1 overlapping interval.  Proposals never
certify an event.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from ..constants import N_FOLDS, Q99
from ..errors import InputContractError
from .common import PRIMARY_CHROMS, ROTATIONS, as_bool, as_int, normfam, rotation_label

__all__ = ["ProposalRecord", "load_proposal_records", "discover_regions", "full_data_proposals", "rotation_proposals"]


@dataclass(frozen=True, slots=True)
class ProposalRecord:
    pair_id: int
    fold: int
    family: str
    chrom: str
    candidate_lo: int
    candidate_hi: int


def load_proposal_records(accounting: Iterable[Mapping[str, str]],
                          host_alternatives: Iterable[Mapping[str, str]]) -> list[ProposalRecord]:
    selected: dict[int, tuple[int, str]] = {}
    for row in accounting:
        pid, fold = as_int(row.get("pair_id")), as_int(row.get("fold"))
        if pid < 0 or fold != pid % N_FOLDS:
            raise InputContractError(f"pair accounting fold contract violated for {pid}")
        if (row.get("experiment_state") == "ELIGIBLE" and row.get("reference_state") == "REFERENCE_UNEXPLAINED"
                and row.get("host_class") == "UNIQUE_HIGH" and not as_bool(row.get("reference_cap_censored"))
                and not as_bool(row.get("host_cap_censored"))):
            family = normfam(row.get("family", ""))
            if not family:
                continue
            if pid in selected:
                raise InputContractError(f"duplicate pair accounting row {pid}")
            selected[pid] = (fold, family)
    rank1: dict[int, tuple[str, int, int, str]] = {}
    for row in host_alternatives:
        pid = as_int(row.get("pair_id"))
        if pid not in selected or as_int(row.get("placement_rank")) != 1:
            continue
        if as_bool(row.get("cap_censored")):
            raise InputContractError(f"UNIQUE_HIGH pair {pid} has censored host alternatives")
        if pid in rank1:
            raise InputContractError(f"more than one rank-1 host alternative for {pid}")
        strand, start, end = str(row.get("strand", "")), as_int(row.get("start")), as_int(row.get("end"))
        if strand not in {"+", "-"} or start < 0 or end <= start:
            raise InputContractError(f"invalid host alternative for {pid}")
        rank1[pid] = (str(row.get("chrom", "")), start, end, strand)
    missing = set(selected) - set(rank1)
    if missing:
        raise InputContractError(f"{len(missing)} UNIQUE_HIGH pairs lack a rank-1 host alternative")
    records = []
    for pid in sorted(selected):
        fold, family = selected[pid]
        chrom, start, end, strand = rank1[pid]
        lo, hi = (end, end + Q99) if strand == "+" else (max(0, start - Q99), start)
        records.append(ProposalRecord(pid, fold, family, chrom, lo, hi))
    return records


def discover_regions(records: Iterable[ProposalRecord], discovery_folds: set[int]) -> list[dict]:
    by: dict[tuple[str, str], list[ProposalRecord]] = defaultdict(list)
    for r in records:
        if r.chrom in PRIMARY_CHROMS and r.fold in discovery_folds:
            by[(r.family, r.chrom)].append(r)
    out: list[dict] = []
    for (family, chrom), values in sorted(by.items()):
        events: dict[int, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: {"start": [], "end": []})
        for r in values:
            events[r.candidate_lo]["start"].append((r.pair_id, r.fold))
            events[r.candidate_hi + 1]["end"].append((r.pair_id, r.fold))
        xs = sorted(events)
        active: dict[int, int] = {}
        segments = []
        for i, x in enumerate(xs):
            for pid, _ in events[x]["end"]:
                active.pop(pid, None)
            for pid, fold in events[x]["start"]:
                active[pid] = fold
            if i + 1 >= len(xs):
                break
            nx = xs[i + 1]
            if nx <= x:
                continue
            folds = set(active.values())
            if len(active) >= 2 and len(folds) >= 2:
                segments.append((x, nx - 1, set(active), folds))
        current = None
        for lo, hi, pids, folds in segments:
            if current is None or lo > current["hi"] + 1:
                if current is not None:
                    out.append(current)
                current = {"family": family, "chrom": chrom, "lo": lo, "hi": hi, "pair_ids": set(pids), "folds": set(folds)}
            else:
                current["hi"] = max(current["hi"], hi)
                current["pair_ids"] |= pids
                current["folds"] |= folds
        if current is not None:
            out.append(current)
    out.sort(key=lambda c: (c["chrom"], c["lo"], c["hi"], c["family"]))
    return out


def full_data_proposals(sample: str, records: list[ProposalRecord]) -> list[dict]:
    rows = []
    for j, c in enumerate(discover_regions(records, set(range(N_FOLDS))), 1):
        rows.append({
            "sample": sample, "candidate_id": f"{sample}_FULL_{j:06d}", "family": c["family"], "chrom": c["chrom"],
            "candidate_lo": c["lo"], "candidate_hi": c["hi"], "candidate_width": c["hi"] - c["lo"] + 1,
            "support_pair_count": len(c["pair_ids"]), "support_fold_count": len(c["folds"]),
            "support_pair_ids": ";".join(map(str, sorted(c["pair_ids"]))),
            "support_folds": ";".join(map(str, sorted(c["folds"]))),
        })
    return rows


def rotation_proposals(sample: str, records: list[ProposalRecord]) -> list[dict]:
    rows = []
    for held in ROTATIONS:
        discovery = set(range(N_FOLDS)) - set(held)
        index: dict[tuple[str, str, int], list[tuple[int, int, int]]] = defaultdict(list)
        for r in records:
            if r.chrom in PRIMARY_CHROMS and r.fold in held:
                index[(r.family, r.chrom, r.fold)].append((r.candidate_lo, r.candidate_hi, r.pair_id))
        for key in index:
            index[key].sort()
        label = rotation_label(held)
        for j, c in enumerate(discover_regions(records, discovery), 1):
            support = {}
            for fold in held:
                hits = []
                for a, b, pid in index.get((c["family"], c["chrom"], fold), []):
                    if a > c["hi"]:
                        break
                    if b >= c["lo"]:
                        hits.append(pid)
                support[fold] = hits
            confirmed = all(support[f] for f in held)
            # The frozen confirmation stops at the first empty held-out fold.
            b_ids = support[held[1]] if support[held[0]] else []
            rows.append({
                "sample": sample, "rotation": label, "candidate_id": f"{sample}_{label}_{j:06d}",
                "family": c["family"], "chrom": c["chrom"], "candidate_lo": c["lo"], "candidate_hi": c["hi"],
                "candidate_width": c["hi"] - c["lo"] + 1, "discovery_pair_count": len(c["pair_ids"]),
                "discovery_fold_count": len(c["folds"]), "discovery_pair_ids": ";".join(map(str, sorted(c["pair_ids"]))),
                "discovery_folds": ";".join(map(str, sorted(c["folds"]))), "heldout_fold_a": held[0], "heldout_fold_b": held[1],
                "heldout_a_pair_ids": ";".join(map(str, support[held[0]])), "heldout_b_pair_ids": ";".join(map(str, b_ids)),
                "heldout_confirmed": int(confirmed),
            })
    return rows
