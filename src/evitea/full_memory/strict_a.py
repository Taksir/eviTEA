"""Global exact junction evidence ("A") and the full-data event table.

Witness: an uncensored ``ELIGIBLE`` + ``REFERENCE_UNEXPLAINED`` pair with a
finite, uncensored, *not* evidence-complete proper reference alternative whose
TE-mate CIGAR has a terminal soft clip overlapping an unexplained candidate-
family query interval (half-open overlap); the clip edge is the exact
interbase breakpoint (a leading clip gives RIGHT at ``te_start``; a trailing
clip gives LEFT at ``te_end``).  Exact anchors need >= 2 pairs over >= 2 folds; per
``(family, chrom, breakpoint)`` one side is kept (most pairs, most folds,
LEFT first).  A full-data event is an anchor contained in >= 1 full-data
proposal; ``proposal_only_A`` marks exact keys absent from the paired-end
anchor layer's witness set.  The event table's paired-end anchor support state ("B") excludes
the anchor's own A pairs.
"""

from __future__ import annotations

import bisect
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping

from ..constants import MIN_FOLDS, MIN_PAIRS, N_FOLDS, Q99
from ..errors import InputContractError
from .common import as_bool, as_int, normfam

__all__ = [
    "BIndex", "allowed_pairs", "derive_witnesses", "exact_anchors", "a_groups",
    "build_b_index", "score_b", "full_data_events", "historical_keys", "EVENT_LEDGER_FIELDS",
]

EVENT_LEDGER_FIELDS = (
    "sample", "family", "chrom", "breakpoint0", "anchor_host_side", "A_pairs", "A_folds", "A_pair_ids", "B_state",
    "left_B_pairs", "left_B_folds", "left_B_pair_ids", "right_B_pairs", "right_B_folds", "right_B_pair_ids",
    "candidate_ids", "proposal_only_A", "validated_anchorpe",
)
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _cigar(cigar: str) -> list[tuple[int, str]]:
    out, pos = [], 0
    for m in _CIGAR.finditer(str(cigar or "")):
        if m.start() != pos:
            return []
        out.append((int(m.group(1)), m.group(2)))
        pos = m.end()
    return out if pos == len(str(cigar or "")) and out else []


def _intervals(text: str) -> list[tuple[int, int]]:
    out = []
    for token in re.split(r"[;,]", str(text or "")):
        m = re.fullmatch(r"(\d+)-(\d+)", token.strip())
        if m and int(m.group(2)) > int(m.group(1)):
            out.append((int(m.group(1)), int(m.group(2))))
    return out


def allowed_pairs(accounting: Iterable[Mapping[str, str]]) -> dict[int, str]:
    """Map each uncensored ELIGIBLE + REFERENCE_UNEXPLAINED pair to its family."""

    allowed: dict[int, str] = {}
    for row in accounting:
        pid = as_int(row.get("pair_id"))
        if pid < 0 or as_int(row.get("fold")) != pid % N_FOLDS:
            raise InputContractError(f"pair accounting row for pair {pid} has an invalid pair ID or fold")
        if (row.get("experiment_state") != "ELIGIBLE" or row.get("reference_state") != "REFERENCE_UNEXPLAINED"
                or as_bool(row.get("reference_cap_censored"))):
            continue
        family = normfam(row.get("family", ""))
        if family:
            allowed[pid] = family
    return allowed


def derive_witnesses(allowed: Mapping[int, str], reference_alternatives: Iterable[Mapping[str, str]]) -> list[dict]:
    rows: list[dict] = []
    seen: set[tuple] = set()
    for r in reference_alternatives:
        pid = as_int(r.get("pair_id"))
        family = allowed.get(pid)
        if family is None or as_bool(r.get("reference_cap_censored")) or as_bool(r.get("evidence_complete")):
            continue
        if normfam(r.get("family", "")) != family:
            raise InputContractError(f"family mismatch for pair {pid}")
        ops, unexplained = _cigar(r.get("te_cigar", "")), _intervals(r.get("unexplained_query_intervals", ""))
        if not ops or not unexplained:
            continue
        qlen = sum(n for n, op in ops if op in "MIS=X")
        cases = []
        if ops[0][1] == "S" and ops[0][0] > 0:
            cases.append((0, ops[0][0], "RIGHT", as_int(r.get("te_start"))))
        if ops[-1][1] == "S" and ops[-1][0] > 0:
            cases.append((qlen - ops[-1][0], qlen, "LEFT", as_int(r.get("te_end"))))
        chrom = str(r.get("chrom", ""))
        for a, b, side, bp in cases:
            if bp < 0 or not any(max(a, c) < min(b, d) for c, d in unexplained):
                continue
            key = (pid, family, chrom, bp, side)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"pair_id": pid, "fold": pid % N_FOLDS, "family": family, "chrom": chrom,
                         "breakpoint0": bp, "anchor_host_side": side})
    return rows


def a_groups(witnesses: Iterable[Mapping]) -> dict[tuple[str, str, int, str], set[int]]:
    groups: dict[tuple[str, str, int, str], set[int]] = defaultdict(set)
    for w in witnesses:
        groups[(normfam(w["family"]), str(w["chrom"]), int(w["breakpoint0"]), str(w["anchor_host_side"]))].add(int(w["pair_id"]))
    return groups


def exact_anchors(groups: Mapping[tuple[str, str, int, str], set[int]], folds: set[int]) -> list[dict]:
    by_coord: dict[tuple[str, str, int], list[dict]] = defaultdict(list)
    for (family, chrom, bp, side), pids in groups.items():
        kept = sorted(p for p in pids if p % N_FOLDS in folds)
        fold_ids = sorted({p % N_FOLDS for p in kept})
        if len(kept) >= MIN_PAIRS and len(fold_ids) >= MIN_FOLDS:
            by_coord[(family, chrom, bp)].append({"family": family, "chrom": chrom, "breakpoint0": bp,
                                                   "anchor_host_side": side, "A_pairs": len(kept),
                                                   "A_folds": len(fold_ids), "A_pair_ids": kept})
    out = []
    for key in sorted(by_coord):
        rows = sorted(by_coord[key], key=lambda r: (-r["A_pairs"], -r["A_folds"], 0 if r["anchor_host_side"] == "LEFT" else 1))
        out.append(rows[0])
    return out


class BIndex:
    """Directional q99 intervals keyed by (family, chrom, side)."""

    def __init__(self, q99: int = Q99) -> None:
        self.q99 = q99
        self._raw: dict[tuple[str, str, str], list[tuple[int, int, int, int]]] = defaultdict(list)
        self._idx: dict[tuple[str, str, str], tuple[list[int], list[tuple[int, int, int, int]]]] = {}

    def add(self, family: str, chrom: str, side: str, lo: int, hi: int, pid: int, fold: int) -> None:
        self._raw[(family, chrom, side)].append((int(lo), int(hi), int(pid), int(fold)))

    def freeze(self) -> "BIndex":
        for key, xs in self._raw.items():
            xs.sort()
            self._idx[key] = ([x[0] for x in xs], xs)
        self._raw.clear()
        return self

    def support(self, family: str, chrom: str, side: str, bp: int, folds: set[int] | None = None,
                excluded: set[int] = frozenset()) -> list[tuple[int, int, int, int]]:
        item = self._idx.get((family, chrom, side))
        if item is None:
            return []
        los, xs = item
        a, z = bisect.bisect_left(los, bp - self.q99), bisect.bisect_right(los, bp)
        return [x for x in xs[a:z] if (folds is None or x[3] in folds) and x[2] not in excluded and x[0] <= bp <= x[1]]


def build_b_index(b_universe: Iterable[Mapping], q99: int = Q99) -> tuple[BIndex, set[int]]:
    """Build the index from paired-end anchor support universe rows (host strand geometry); q99 is the fragment-size 99th percentile."""

    index, pids = BIndex(q99), set()
    for r in b_universe:
        pid, fold = int(r["pair_id"]), int(r["fold"])
        if fold != pid % N_FOLDS:
            raise InputContractError("B fold drift")
        strand, start, end = str(r["strand"]), int(r["start"]), int(r["end"])
        if strand == "+":
            side, lo, hi = "LEFT", end, end + q99
        elif strand == "-":
            side, lo, hi = "RIGHT", max(0, start - q99), start
        else:
            continue
        index.add(normfam(r["family"]), str(r["chrom"]), side, lo, hi, pid, fold)
        pids.add(pid)
    return index.freeze(), pids


def score_b(anchor: Mapping, index: BIndex, folds: set[int] | None, *, exclude_a: bool) -> dict:
    excluded = set(anchor["A_pair_ids"]) if exclude_a else set()
    sides = {}
    for side in ("LEFT", "RIGHT"):
        xs = index.support(anchor["family"], anchor["chrom"], side, int(anchor["breakpoint0"]), folds, excluded)
        sides[side] = (sorted({x[2] for x in xs}), sorted({x[3] for x in xs}))
    left = len(sides["LEFT"][0]) >= MIN_PAIRS and len(sides["LEFT"][1]) >= MIN_FOLDS
    right = len(sides["RIGHT"][0]) >= MIN_PAIRS and len(sides["RIGHT"][1]) >= MIN_FOLDS
    state = "BILATERAL" if left and right else "LEFT_ONLY" if left else "RIGHT_ONLY" if right else "NEITHER"
    return {"B_state": state, "left_B_pairs": len(sides["LEFT"][0]), "left_B_folds": len(sides["LEFT"][1]),
            "left_B_pair_ids": sides["LEFT"][0], "right_B_pairs": len(sides["RIGHT"][0]),
            "right_B_folds": len(sides["RIGHT"][1]), "right_B_pair_ids": sides["RIGHT"][0]}


def historical_keys(witnesses: Iterable[Mapping]) -> set[tuple[str, str, int, str]]:
    """Exact keys of the paired-end anchor layer's (corrected-cluster) exact-junction witnesses."""

    return {(normfam(w["family"]), str(w["chrom"]), int(w["breakpoint0"]), str(w["anchor_host_side"]))
            for w in witnesses if int(w["pair_id"]) >= 0 and int(w["breakpoint0"]) >= 0
            and w["anchor_host_side"] in {"LEFT", "RIGHT"}}


def full_data_events(sample: str, full_proposals: list[Mapping], groups, index: BIndex,
                     old_keys: set[tuple[str, str, int, str]]) -> list[dict]:
    by_key: dict[tuple, dict] = {}
    for anchor in exact_anchors(groups, set(range(N_FOLDS))):
        hits = [p["candidate_id"] for p in full_proposals
                if normfam(p.get("family", "")) == anchor["family"] and str(p.get("chrom", "")) == anchor["chrom"]
                and as_int(p.get("candidate_lo")) <= anchor["breakpoint0"] <= as_int(p.get("candidate_hi"))]
        if not hits:
            continue
        key = (anchor["family"], anchor["chrom"], anchor["breakpoint0"], anchor["anchor_host_side"])
        b = score_b(anchor, index, set(range(N_FOLDS)), exclude_a=True)
        by_key[key] = {"sample": sample, **anchor, **b, "candidate_ids": hits,
                       "proposal_only_A": int(key not in old_keys), "validated_anchorpe": int(b["B_state"] != "NEITHER")}
    return list(by_key.values())
