"""Reselect the coverage panel by newly routed pairs.

The B = 1024 coverage panel fixes, per group, a coverage and
terminal target.  Starting from the sparse-panel and Hamming-1 21-mer strings,
groups (hardest first) re-add its exact actions.
Strings already paid for elsewhere are free.  Terminal sides are restored
first.  Then coverage is raised to 99.5% of the panel's coverage by choosing
the action with the most new windows per *newly routed physical pair*.  Routed
reads and pairs are the union over selected strings.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..errors import BackendCompatibilityError
from .anchor_census import AnchorCatalogue
from .coverage_panel import (
    READ_WINDOW,
    RESCUE_PANEL_CAP,
    Group,
    RescuePanel,
    ShortLocus,
    baseline_anchor_loci,
    coverage_interval,
    merge_intervals,
    terminal_flags,
)
from .sparse_panel import PanelQuery

__all__ = ["FROZEN_RESCUE_CAP", "PairAwareGroup", "PairAwarePanel", "select_pair_aware_panel"]

FROZEN_RESCUE_CAP = RESCUE_PANEL_CAP
TARGET_RETAIN_FRACTION = 0.995


@dataclass(frozen=True, slots=True)
class PairAwareGroup:
    te_id: str
    family: str
    orientation: str
    window_count: int
    selected_coverage: float


@dataclass(slots=True)
class PairAwarePanel:
    selected_sequences: set[str]
    sequence_widths: dict[str, int]
    sequence_families: dict[str, set[str]]
    groups: list[PairAwareGroup]
    read_ids: np.ndarray
    pair_ids: np.ndarray
    pair_read_counts: np.ndarray
    pair_sets: dict[str, np.ndarray]


def _mean(mask: np.ndarray) -> float:
    return float(mask.mean()) if len(mask) else 1.0


def _apply(mask: np.ndarray, intervals) -> None:
    for lo, hi in intervals:
        mask[lo : hi + 1] = True


def _gain(mask: np.ndarray, intervals) -> int:
    return sum(int(np.count_nonzero(~mask[lo : hi + 1])) for lo, hi in intervals)


def select_pair_aware_panel(
    catalogue: AnchorCatalogue,
    panel: Sequence[PanelQuery],
    short: Sequence[ShortLocus],
    rescue: RescuePanel,
    index,
    read_window: int = READ_WINDOW,
) -> PairAwarePanel:
    if rescue.cap != FROZEN_RESCUE_CAP:
        raise BackendCompatibilityError(
            f"pair-aware reselection accepts only the B={FROZEN_RESCUE_CAP} coverage panel; got B={rescue.cap}"
        )
    records = catalogue.records
    family_by_te = catalogue.family_by_te()
    lengths: dict[Group, int] = {}
    families: dict[Group, str] = {}
    for te_id, sequence in sorted(records.items()):
        for orientation in ("F", "RC"):
            lengths[(te_id, orientation)] = len(sequence)
            families[(te_id, orientation)] = family_by_te.get(te_id, te_id)

    baseline = baseline_anchor_loci(catalogue, panel)
    masks: dict[Group, np.ndarray] = {}
    terms: dict[Group, tuple[bool, bool]] = {}
    for key, length in lengths.items():
        mask = np.zeros(max(0, length - read_window + 1), dtype=bool)
        loci = [(s, k) for s, k, _, _ in baseline.get(key, ())]
        for s, k in loci:
            iv = coverage_interval(length, s, k, read_window)
            if iv is not None:
                mask[iv[0] : iv[1] + 1] = True
        masks[key] = mask
        terms[key] = terminal_flags(length, loci)

    wanted: dict[tuple[str, str, str], dict] = {}
    h1_fixed: dict[Group, list[dict]] = defaultdict(list)
    for action in rescue.selected_actions:
        key = (action.te_id, action.orientation)
        if action.action_type == "SHORT_EXACT":
            seq = action.action_id.split(":", 3)[3]
            wanted[(action.te_id, action.orientation, seq)] = {"seq": seq, "family": action.family, "loci": []}
        elif action.action_type == "H1_COMPLETE":
            start = int(action.action_id.split(":")[3])
            length = lengths[key]
            t5, t3 = terminal_flags(length, [(start, action.anchor_k)])
            h1_fixed[key].append({
                "intervals": merge_intervals([coverage_interval(length, start, action.anchor_k, read_window)]),
                "terminal5": t5,
                "terminal3": t3,
            })
    for locus in short:
        entry = wanted.get((locus.te_id, locus.orientation, locus.query_seq))
        if entry is not None:
            entry["loci"].append((locus.start, locus.k))
    exact_actions: dict[Group, list[dict]] = defaultdict(list)
    for (te_id, orientation, _), entry in wanted.items():
        if not entry["loci"]:
            raise BackendCompatibilityError("could not reconstruct a selected exact action")
        length = lengths[(te_id, orientation)]
        entry["intervals"] = merge_intervals(coverage_interval(length, s, k, read_window) for s, k in entry["loci"])
        entry["terminal5"], entry["terminal3"] = terminal_flags(length, entry["loci"])
        exact_actions[(te_id, orientation)].append(entry)

    widths = dict(rescue.selected_sequences)
    sources = rescue.selected_sources
    read_sets: dict[str, np.ndarray] = {}
    pair_sets: dict[str, np.ndarray] = {}
    for seq in sorted(widths):
        dollars = index.dollar_ids(seq)
        if len(dollars) != widths[seq]:
            raise BackendCompatibilityError(f"interval width changed for {seq}")
        unique = np.unique(dollars)
        read_sets[seq] = unique
        pair_sets[seq] = np.unique(index.file_and_pair(unique)[1])

    selected = {s for s, src in sources.items() if "G1B_BASELINE" in src or "H1_COMPLETE" in src}
    routed_pair = np.zeros(index.pair_count, dtype=bool)
    routed_read = np.zeros(index.record_count, dtype=bool)
    for seq in sorted(selected):
        routed_read[read_sets[seq]] = True
        routed_pair[pair_sets[seq]] = True

    targets = rescue.group_targets
    order = sorted(targets, key=lambda k: (targets[k].selected_coverage, -targets[k].window_count, k))
    selected_family_rows: list[tuple[str, str]] = []
    groups: list[PairAwareGroup] = []
    for key in order:
        mask = masks[key].copy()
        cur5, cur3 = terms[key]
        for fixed in h1_fixed.get(key, ()):
            _apply(mask, fixed["intervals"])
            cur5 = cur5 or fixed["terminal5"]
            cur3 = cur3 or fixed["terminal3"]
        remaining = []
        for action in exact_actions.get(key, ()):
            if action["seq"] in selected:
                _apply(mask, action["intervals"])
                cur5 = cur5 or action["terminal5"]
                cur3 = cur3 or action["terminal3"]
            else:
                remaining.append(action)
        target = targets[key]
        target_cov = max(_mean(mask), target.selected_coverage * TARGET_RETAIN_FRACTION)

        def cost_of(action) -> int:
            return int(np.count_nonzero(~routed_pair[pair_sets[action["seq"]]]))

        def take(action) -> None:
            nonlocal cur5, cur3
            selected.add(action["seq"])
            routed_read[read_sets[action["seq"]]] = True
            routed_pair[pair_sets[action["seq"]]] = True
            _apply(mask, action["intervals"])
            cur5 = cur5 or action["terminal5"]
            cur3 = cur3 or action["terminal3"]
            selected_family_rows.append((action["seq"], action["family"]))
            remaining.remove(action)

        need5 = target.selected_terminal5 and not cur5
        need3 = target.selected_terminal3 and not cur3
        while need5 or need3:
            best = best_key = None
            for action in remaining:
                covered = int(need5 and action["terminal5"]) + int(need3 and action["terminal3"])
                if not covered:
                    continue
                cost = cost_of(action)
                rank = (cost / covered, cost, -covered, action["seq"])
                if best_key is None or rank < best_key:
                    best_key, best = rank, action
            if best is None:
                raise BackendCompatibilityError(f"cannot preserve terminal target for {key}")
            take(best)
            need5 = target.selected_terminal5 and not cur5
            need3 = target.selected_terminal3 and not cur3
        while _mean(mask) + 1e-15 < target_cov:
            best = best_key = None
            for action in remaining:
                gain = _gain(mask, action["intervals"])
                if gain <= 0:
                    continue
                cost = cost_of(action)
                score = math.inf if cost == 0 else gain / cost
                rank = (-score, -gain, cost, action["seq"])
                if best_key is None or rank < best_key:
                    best_key, best = rank, action
            if best is None:
                raise BackendCompatibilityError(f"coverage target unreachable for {key}")
            take(best)
        groups.append(PairAwareGroup(key[0], families[key], key[1], len(mask), _mean(mask)))

    read_ids = np.flatnonzero(routed_read).astype(np.uint64)
    pair_ids = np.flatnonzero(routed_pair).astype(np.uint64)
    read_pairs = index.file_and_pair(read_ids)[1]
    check_pairs, counts = np.unique(read_pairs, return_counts=True)
    if not np.array_equal(check_pairs, pair_ids) or np.any((counts < 1) | (counts > 2)):
        raise BackendCompatibilityError("pair-aware routed reads and pairs are inconsistent")

    sequence_families: dict[str, set[str]] = defaultdict(set)
    for entry in panel:
        sequence_families[entry.query_seq].update(entry.source_families)
    for seq, family in selected_family_rows:
        if seq and family:
            sequence_families[seq].add(family)
    return PairAwarePanel(
        selected_sequences=selected,
        sequence_widths={s: widths[s] for s in selected},
        sequence_families={s: set(sequence_families.get(s, set())) for s in selected},
        groups=groups,
        read_ids=read_ids,
        pair_ids=pair_ids,
        pair_read_counts=counts.astype(np.uint8),
        pair_sets={s: pair_sets[s] for s in selected},
    )
