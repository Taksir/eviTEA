"""Weak-site pass of the library layer (opt-in, ``--weak-site-pass``; never run by default).

The default library layer calls a locus only when a stage finds two fragments on each side
(``mate_rescue.MIN_SIDE_FRAGMENTS``), an exact split junction of score 4 or more, or the repeat-aware
2-LEFT + 2-RIGHT rule.  At low depth a true insertion often has one fragment on a side.  This pass takes the evidence the
stages already saved and calls a site that no call covers when the reads of the two sides, pooled over the stages, are at
least one on each side, or at least two split reads sit on one side.

Evidence is grouped as the desk replay of 3 Nov 2026 grouped it (``e1/replay_floors.py``, function ``sites``):

* split events (``SplitJunctionResult.all_events``): one-sided left or right events give their fragment ids to that side of
  the site; a two-sided event gives none but widens the site;
* read pairs (``MateRescueResult.ledger``): the left and right pair ids of each locus;
* junction-group reads (``JunctionGroupsResult.ledger``): the strict left and right ids of each locus;
* items whose starts lie within ``SITE_GAP`` (100) bases of the running end of the site join one site.

Per side the distinct fragment ids of the split, pair and strict items are pooled.  A site with a call of the layer within
``NEAR_BP`` (500) bases becomes no call.  The family of a new call is the family with the most supporting ids at the site
(ties: the family name that sorts first); its coordinates are the span of that family's items; ``left_support`` and
``right_support`` are the pooled counts per side and ``fragment_pair_ids`` the pooled ids, so the support of a call is the
count the rule used.  The event id is ``WEAK_<n>``.  Nothing here reads a reference annotation: the new rows go through the
same own-family copy filter and breakpoint refinement as every other library call.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from .events import IncumbentCall, IncumbentCallClass
from .junction_groups import JunctionGroupsResult
from .mate_rescue import MateRescueResult
from .split_junction import SplitJunctionResult

__all__ = [
    "EVIDENCE_TYPE", "EVENT_PREFIX", "SITE_GAP", "NEAR_BP", "POOLED_PER_SIDE", "SPLIT_ONE_SIDE",
    "WeakSiteResult", "evidence_sites", "weak_site_calls",
]

SITE_GAP = 100          # evidence items closer than this join one site
NEAR_BP = 500           # a site with a call of the layer this close becomes no call
POOLED_PER_SIDE = 1     # pooled reads needed on each side
SPLIT_ONE_SIDE = 2      # or split reads needed on one side
EVIDENCE_TYPE = "weak_site_pooled_floor"
EVENT_PREFIX = "WEAK_"

_COUNTED = ("split", "pair", "strict")


@dataclass(frozen=True, slots=True)
class WeakSiteResult:
    calls: tuple[IncumbentCall, ...]
    sites: int              # evidence sites
    uncovered: int          # sites with no call of the layer within NEAR_BP


def _items(split: SplitJunctionResult, mates: MateRescueResult, groups: JunctionGroupsResult):
    """(chrom, start, end, kind, side, family, ids) of every evidence item."""

    items: list[tuple[str, int, int, str, str, str, frozenset[int]]] = []
    for event in split.all_events:
        ids = frozenset(event.fragment_pair_ids)
        if event.evidence_type == "one_sided_left":
            items.append((event.chrom, event.start, event.end, "split", "L", event.family, ids))
        elif event.evidence_type == "one_sided_right":
            items.append((event.chrom, event.start, event.end, "split", "R", event.family, ids))
        else:
            items.append((event.chrom, event.start, event.end, "split2", "B", event.family, ids))
    for locus in mates.ledger:
        start = int(min(locus.candidate_start, locus.candidate_end))
        end = int(max(locus.candidate_start, locus.candidate_end)) + 1
        if locus.left_pair_ids:
            items.append((locus.chrom, start, end, "pair", "L", locus.family, frozenset(locus.left_pair_ids)))
        if locus.right_pair_ids:
            items.append((locus.chrom, start, end, "pair", "R", locus.family, frozenset(locus.right_pair_ids)))
    for locus in groups.ledger:
        start = int(min(locus.candidate_start, locus.candidate_end))
        end = int(max(locus.candidate_start, locus.candidate_end)) + 1
        if locus.strict_left_ids:
            items.append((locus.chrom, start, end, "strict", "L", locus.family, frozenset(locus.strict_left_ids)))
        if locus.strict_right_ids:
            items.append((locus.chrom, start, end, "strict", "R", locus.family, frozenset(locus.strict_right_ids)))
    return items


def evidence_sites(split: SplitJunctionResult, mates: MateRescueResult, groups: JunctionGroupsResult) -> list[dict]:
    """Evidence sites in (chrom, start) order; each holds ``chrom``, ``start``, ``end``, ``items`` and, per side, the pooled
    (``left``/``right``) and split-only (``split_left``/``split_right``) distinct fragment ids."""

    items = sorted(_items(split, mates, groups), key=lambda item: (item[0], item[1]))
    sites: list[dict] = []
    current: dict | None = None
    for item in items:
        if current is not None and item[0] == current["chrom"] and item[1] <= current["end"] + SITE_GAP:
            current["end"] = max(current["end"], item[2])
            current["items"].append(item)
        else:
            current = {"chrom": item[0], "start": item[1], "end": item[2], "items": [item]}
            sites.append(current)
    for site in sites:
        left: set[int] = set()
        right: set[int] = set()
        split_left: set[int] = set()
        split_right: set[int] = set()
        for _chrom, _start, _end, kind, side, _family, ids in site["items"]:
            if kind not in _COUNTED:
                continue
            (left if side == "L" else right).update(ids)
            if kind == "split":
                (split_left if side == "L" else split_right).update(ids)
        site.update(left=left, right=right, split_left=split_left, split_right=split_right)
    return sites


def passes(site: dict) -> bool:
    """The rule: pooled reads on each side, or split reads on one side."""

    return (min(len(site["left"]), len(site["right"])) >= POOLED_PER_SIDE
            or max(len(site["split_left"]), len(site["split_right"])) >= SPLIT_ONE_SIDE)


def _covered(site: dict, by_chrom: dict[str, list[tuple[int, int]]]) -> bool:
    return any(start <= site["end"] + NEAR_BP and end >= site["start"] - NEAR_BP
               for start, end in by_chrom.get(site["chrom"], ()))


def weak_site_calls(
    sample_id: str,
    split: SplitJunctionResult,
    mates: MateRescueResult,
    groups: JunctionGroupsResult,
    existing_calls: Iterable[IncumbentCall],
) -> WeakSiteResult:
    """Calls for the evidence sites that pass the rule and have no call of the layer within ``NEAR_BP``."""

    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for call in existing_calls:
        by_chrom[call.chrom].append((call.start, call.end))
    sites = evidence_sites(split, mates, groups)
    calls: list[IncumbentCall] = []
    uncovered = 0
    for site in sites:
        if _covered(site, by_chrom):
            continue
        uncovered += 1
        if not passes(site):
            continue
        votes: dict[str, set[int]] = defaultdict(set)
        spans: dict[str, list[int]] = {}
        for _chrom, start, end, kind, _side, family, ids in site["items"]:
            if kind in _COUNTED:
                votes[family].update(ids)
            span = spans.setdefault(family, [start, end])
            span[0] = min(span[0], start)
            span[1] = max(span[1], end)
        family = min(votes, key=lambda name: (-len(votes[name]), name))
        start, end = spans[family]
        union = tuple(sorted(site["left"] | site["right"]))
        calls.append(IncumbentCall(
            sample_id, f"{EVENT_PREFIX}{len(calls) + 1:06d}", IncumbentCallClass.WEAK_SITE, site["chrom"],
            start, max(end, start + 1), family, None, None, None, len(site["left"]), len(site["right"]),
            EVIDENCE_TYPE, float(len(union)), union,
        ))
    return WeakSiteResult(tuple(calls), len(sites), uncovered)
