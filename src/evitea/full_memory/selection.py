"""A-only proposal selection and full-data Channel-B support.

For every full-data proposal, among the ledger's *proposal-only* A events that
it contains, one anchor is chosen by A evidence only (most pairs, most folds,
lower breakpoint, LEFT first) — B is not consulted.  The chosen anchor is then
scored against the full-data Channel-B universe without excluding its A pairs
(the frozen hostile-audit scoring used by the production freeze).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping

from ..shadow_layers import FullMemoryEvent
from .common import as_int, normfam
from .strict_a import BIndex, score_b

__all__ = ["select_anchors", "full_memory_events"]


def _ids(value) -> list[int]:
    if isinstance(value, (list, tuple, set)):
        return sorted(int(x) for x in value)
    return [int(float(x)) for x in str(value or "").split(";") if x.strip()]


def select_anchors(event_ledger: Iterable[Mapping], sample: str) -> list[dict]:
    by: dict[str, list[Mapping]] = defaultdict(list)
    for row in event_ledger:
        if row["sample"] != sample or as_int(row.get("proposal_only_A")) != 1:
            continue
        candidates = row["candidate_ids"]
        for cid in (candidates if isinstance(candidates, (list, tuple)) else str(candidates).split(";")):
            if cid:
                by[cid].append(row)
    out = []
    for cid, rows in sorted(by.items()):
        best = min(rows, key=lambda r: (-as_int(r["A_pairs"]), -as_int(r["A_folds"]), as_int(r["breakpoint0"]),
                                        0 if r["anchor_host_side"] == "LEFT" else 1))
        out.append({"candidate_id": cid, "family": normfam(best["family"]), "chrom": best["chrom"],
                    "breakpoint0": as_int(best["breakpoint0"]), "anchor_host_side": best["anchor_host_side"],
                    "A_pairs": as_int(best["A_pairs"]), "A_folds": as_int(best["A_folds"]),
                    "A_pair_ids": _ids(best["A_pair_ids"])})
    return out


def full_memory_events(anchors: Iterable[Mapping], index: BIndex) -> tuple[FullMemoryEvent, ...]:
    events = []
    for a in anchors:
        b = score_b(a, index, None, exclude_a=False)
        events.append(FullMemoryEvent(
            candidate_id=a["candidate_id"], family=a["family"], chrom=a["chrom"], breakpoint0=a["breakpoint0"],
            anchor_host_side=a["anchor_host_side"], a_pair_ids=tuple(a["A_pair_ids"]),
            left_b_pair_ids=tuple(b["left_B_pair_ids"]), right_b_pair_ids=tuple(b["right_B_pair_ids"]),
            b_state=b["B_state"],
        ))
    return tuple(events)
