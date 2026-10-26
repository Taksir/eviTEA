"""Coverage panel: exact short anchors, Hamming-1 loci, optimized panel.

A TE record orientation is a *group*. A 125 bp TE-only read window is covered
when it wholly contains an admitted exact anchor. A terminal side is covered
when a 35 bp terminal segment contains one.

A targeted exact census covers 19/17/15-mers (terminal starts at stride 1,
interior at stride 10). Where the sparse panel plus exact short
anchors with width <= 1024 leave a terminal side or read windows
uncovered, it adds 21-mer Hamming-1 (H1) loci.

Actions are one per distinct exact string per group, and one per H1 locus
whose source and one-substitution variants all occur with width <= B. For B in
256/512/1024, actions are selected greedily per group (hardest groups first): terminal sides first, then lazy weighted maximum
coverage up to min(0.99, 0.99 * attainable). Among the panels under the
5,000,000-row materialization cap, it chooses the one with the best
(p10, mean, terminal fraction, -rows).
"""

from __future__ import annotations

import heapq
import math
import statistics
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from ..errors import BackendCompatibilityError
from .anchor_census import AnchorCatalogue
from .sequence import revcomp
from .sparse_panel import PanelQuery

__all__ = [
    "READ_WINDOW",
    "TERMINAL_SEGMENT",
    "ShortLocus",
    "H1Locus",
    "GroupTarget",
    "PanelAction",
    "RescuePanel",
    "baseline_anchor_loci",
    "short_exact_census",
    "hard_gap_h1_loci",
    "optimize_rescue_panel",
    "RESCUE_PANEL_CAP",
    "coverage_interval",
    "merge_intervals",
    "terminal_flags",
]

# Default read length (bp) of the window geometry; a run passes the sample's own read length.
READ_WINDOW = 125
TERMINAL_SEGMENT = 35
SHORT_K_VALUES = (19, 17, 15)
SHORT_INTERNAL_STRIDE = 10
SHORT_TERMINAL_SPAN = 120
H1_K = 21
H1_INTERNAL_LOCUS_STRIDE = 50
H1_GAP_CAP = 1024
H1_MAX_DISTINCT_VARIANTS = 1_500_000
PANEL_CAPS = (256, 512, 1024)
# The coverage panel that the pair-aware reselection is built for (the widest cap). Only this cap is built: the narrower caps
# can only tie with it (when every TE word is rare) or lose, and the reselection refuses the 256-row panel.
RESCUE_PANEL_CAP = PANEL_CAPS[-1]
TARGET_ABSOLUTE_COVERAGE = 0.99
RETAIN_ATTAINABLE_FRACTION = 0.99
MATERIALIZE_ROW_HARD_CAP = 5_000_000
_DNA = frozenset("ACGT")

Group = tuple[str, str]


@dataclass(frozen=True, slots=True)
class ShortLocus:
    te_id: str
    family: str
    orientation: str
    start: int
    k: int
    query_seq: str
    width: int


@dataclass(frozen=True, slots=True)
class H1Locus:
    te_id: str
    family: str
    orientation: str
    start: int
    source_21mer: str


@dataclass(frozen=True, slots=True)
class PanelAction:
    action_id: str
    action_type: str
    te_id: str
    family: str
    orientation: str
    intervals: tuple[tuple[int, int], ...]
    terminal5: bool
    terminal3: bool
    members: tuple[tuple[str, int], ...]
    anchor_k: int


@dataclass(frozen=True, slots=True)
class GroupTarget:
    selected_coverage: float
    selected_terminal5: bool
    selected_terminal3: bool
    window_count: int


@dataclass(slots=True)
class RescuePanel:
    cap: int
    group_targets: dict[Group, GroupTarget]
    selected_actions: list[PanelAction]
    selected_sequences: dict[str, int]
    selected_sources: dict[str, set[str]] = field(default_factory=dict)


# -- geometry ------------------------------------------------------------
def coverage_interval(length: int, start: int, k: int, read_window: int = READ_WINDOW) -> tuple[int, int] | None:
    """Return read windows (start offsets of a read of ``read_window`` bp) that contain the k-mer at ``start``."""

    n = max(0, length - read_window + 1)
    if n <= 0:
        return None
    lo = max(0, start + k - read_window)
    hi = min(start, n - 1)
    return (lo, hi) if lo <= hi else None


def merge_intervals(intervals) -> tuple[tuple[int, int], ...]:
    values = sorted(x for x in intervals if x is not None)
    if not values:
        return ()
    out = [list(values[0])]
    for lo, hi in values[1:]:
        if lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return tuple((a, b) for a, b in out)


def terminal_flags(length: int, loci) -> tuple[bool, bool]:
    t5 = t3 = False
    for start, k in loci:
        if start >= 0 and start + k <= min(TERMINAL_SEGMENT, length):
            t5 = True
        if start >= max(0, length - TERMINAL_SEGMENT) and start + k <= length:
            t3 = True
    return t5, t3


def _apply(mask: np.ndarray, intervals) -> None:
    for lo, hi in intervals:
        mask[lo : hi + 1] = True


def _gain(prefix_uncovered: np.ndarray, intervals) -> int:
    return sum(int(prefix_uncovered[hi + 1] - prefix_uncovered[lo]) for lo, hi in intervals)


def _mean(mask: np.ndarray) -> float:
    return float(mask.mean()) if len(mask) else 1.0


def _window_mask(length: int, loci, read_window: int = READ_WINDOW) -> np.ndarray:
    mask = np.zeros(max(0, length - read_window + 1), dtype=bool)
    for start, k in loci:
        iv = coverage_interval(length, start, k, read_window)
        if iv is not None:
            mask[iv[0] : iv[1] + 1] = True
    return mask


# -- exact short anchors and Hamming-1 loci ------------------------------
def baseline_anchor_loci(catalogue: AnchorCatalogue, panel: Sequence[PanelQuery]):
    """Return every queryable anchor locus whose query is in the sparse panel, by group."""

    widths = {entry.query_id: entry.width for entry in panel}
    loci: dict[Group, list[tuple[int, int, str, int]]] = defaultdict(list)
    for source in catalogue.sources:
        if source.query_id in widths and source.queryable:
            loci[(source.te_id, source.orientation)].append(
                (source.start, source.k, source.query_seq, widths[source.query_id])
            )
    return loci


def _short_positions(length: int, k: int) -> list[int]:
    if length < k:
        return []
    last = length - k
    starts = set(range(0, min(SHORT_TERMINAL_SPAN - k, last) + 1))
    starts.update(range(max(0, length - SHORT_TERMINAL_SPAN), last + 1))
    starts.update(range(0, last + 1, SHORT_INTERNAL_STRIDE))
    return sorted(starts)


def short_exact_census(records: Mapping[str, str], family_by_te: Mapping[str, str], index) -> list[ShortLocus]:
    raw: list[tuple[str, str, str, int, int, str]] = []
    for te_id, sequence in sorted(records.items()):
        family = family_by_te.get(te_id, te_id)
        for k in SHORT_K_VALUES:
            for start in _short_positions(len(sequence), k):
                source = sequence[start : start + k]
                if len(source) != k or any(base not in _DNA for base in source):
                    continue
                raw.append((te_id, family, "F", start, k, source))
                raw.append((te_id, family, "RC", start, k, revcomp(source)))
    widths = {query: index.width(query) for query in sorted({row[5] for row in raw})}
    return [ShortLocus(te, fam, ori, start, k, query, widths[query]) for te, fam, ori, start, k, query in raw]


def _terminal_hit(length: int, loci, side: str) -> bool:
    if length <= 0:
        return False
    lo, hi = (0, min(TERMINAL_SEGMENT, length)) if side == "5p" else (max(0, length - TERMINAL_SEGMENT), length)
    return any(start >= lo and start + k <= hi for start, k in loci)


def _uncovered_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    runs = []
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            i += 1
            continue
        j = i + 1
        while j < n and not mask[j]:
            j += 1
        runs.append((i, j))
        i = j
    return runs


def hard_gap_h1_loci(
    records: Mapping[str, str],
    family_by_te: Mapping[str, str],
    baseline: Mapping[Group, list[tuple[int, int, str, int]]],
    short: Sequence[ShortLocus],
    read_window: int = READ_WINDOW,
) -> list[H1Locus]:
    admitted: dict[Group, list[tuple[int, int]]] = defaultdict(list)
    for key, loci in baseline.items():
        admitted[key].extend((start, k) for start, k, _, _ in loci)
    for locus in short:
        if 0 < locus.width <= H1_GAP_CAP:
            admitted[(locus.te_id, locus.orientation)].append((locus.start, locus.k))
    chosen: set[tuple[str, str, str, int]] = set()
    for te_id, sequence in sorted(records.items()):
        family = family_by_te.get(te_id, te_id)
        length = len(sequence)
        if length < H1_K:
            continue
        for orientation in ("F", "RC"):
            current = admitted.get((te_id, orientation), [])
            if not _terminal_hit(length, current, "5p"):
                for start in range(0, min(TERMINAL_SEGMENT - H1_K, length - H1_K) + 1):
                    chosen.add((te_id, family, orientation, start))
            if not _terminal_hit(length, current, "3p"):
                for start in range(max(0, length - TERMINAL_SEGMENT), length - H1_K + 1):
                    chosen.add((te_id, family, orientation, start))
            mask = _window_mask(length, current, read_window)
            for run_start, run_end in _uncovered_runs(mask):
                for window_start in range(run_start, run_end, H1_INTERNAL_LOCUS_STRIDE):
                    start = max(0, min(length - H1_K, window_start + (read_window - H1_K) // 2))
                    chosen.add((te_id, family, orientation, start))
    out: list[H1Locus] = []
    variants: set[str] = set()
    for te_id, family, orientation, start in sorted(chosen):
        source = records[te_id][start : start + H1_K]
        query = source if orientation == "F" else revcomp(source)
        if len(source) == H1_K and not any(base not in _DNA for base in source):
            variants.update(_h1_family(query)[1:])
        out.append(H1Locus(te_id, family, orientation, start, query if len(source) == H1_K and set(source) <= _DNA else source))
    if len(variants) > H1_MAX_DISTINCT_VARIANTS:
        raise BackendCompatibilityError("H1 variant count exceeds the 1,500,000-query cap")
    return out


def _h1_family(source: str) -> list[str]:
    if len(source) != H1_K or any(base not in _DNA for base in source):
        return []
    out = [source]
    for i, old in enumerate(source):
        for base in "ACGT":
            if base != old:
                out.append(source[:i] + base + source[i + 1 :])
    return out


# -- action selection ----------------------------------------------------
def _short_actions(
    short: Sequence[ShortLocus], lengths, family_by_group, read_window: int = READ_WINDOW
) -> dict[Group, list[PanelAction]]:
    grouped: dict[Group, dict[str, dict]] = defaultdict(dict)
    for locus in short:
        if locus.width <= 0:
            continue
        key = (locus.te_id, locus.orientation)
        entry = grouped[key].setdefault(locus.query_seq, {"width": locus.width, "k": locus.k, "loci": []})
        if entry["width"] != locus.width:
            raise BackendCompatibilityError(f"inconsistent interval width for {locus.query_seq}")
        entry["loci"].append((locus.start, locus.k))
    actions: dict[Group, list[PanelAction]] = defaultdict(list)
    for key, by_seq in grouped.items():
        length = lengths[key]
        family = family_by_group.get(key, key[0])
        for seq, entry in by_seq.items():
            loci = entry["loci"]
            t5, t3 = terminal_flags(length, loci)
            actions[key].append(
                PanelAction(
                    f"EXACT:{key[0]}:{key[1]}:{seq}", "SHORT_EXACT", key[0], family, key[1],
                    merge_intervals(coverage_interval(length, s, k, read_window) for s, k in loci),
                    t5, t3, ((seq, int(entry["width"])),), int(entry["k"]),
                )
            )
    return actions


def _h1_actions(
    h1: Sequence[H1Locus], lengths, family_by_group, index, read_window: int = READ_WINDOW
) -> dict[int, dict[Group, list[PanelAction]]]:
    members_by_locus = [(locus, _h1_family(locus.source_21mer)) for locus in h1]
    widths = {seq: index.width(seq) for seq in sorted({m for _, ms in members_by_locus for m in ms})}
    by_cap: dict[int, dict[Group, list[PanelAction]]] = {cap: defaultdict(list) for cap in PANEL_CAPS}
    for locus, members in members_by_locus:
        key = (locus.te_id, locus.orientation)
        length = lengths[key]
        family = family_by_group.get(key, locus.family)
        nonempty = tuple((seq, widths[seq]) for seq in members if widths[seq] > 0)
        widest = max((w for _, w in nonempty), default=0)
        for cap in PANEL_CAPS:
            if nonempty and widest <= cap:
                t5, t3 = terminal_flags(length, [(locus.start, H1_K)])
                by_cap[cap][key].append(
                    PanelAction(
                        f"H1:{locus.te_id}:{locus.orientation}:{locus.start}", "H1_COMPLETE",
                        locus.te_id, family, locus.orientation,
                        merge_intervals([coverage_interval(length, locus.start, H1_K, read_window)]),
                        t5, t3, nonempty, H1_K,
                    )
                )
    return by_cap


def _cost(action: PanelAction, selected: set[str]) -> int:
    return sum(w for seq, w in action.members if seq not in selected)


def _optimize_group(
    key, length, baseline, actions, selected_sequences, selected_sources, selected_actions, read_window=READ_WINDOW
):
    base_loci = baseline.get(key, ())
    base_mask = _window_mask(length, [(s, k) for s, k, _, _ in base_loci], read_window)
    base_t5, base_t3 = terminal_flags(length, [(s, k) for s, k, _, _ in base_loci])
    full_mask = base_mask.copy()
    full_t5, full_t3 = base_t5, base_t3
    for action in actions:
        _apply(full_mask, action.intervals)
        full_t5 = full_t5 or action.terminal5
        full_t3 = full_t3 or action.terminal3
    full_cov = _mean(full_mask)
    base_cov = _mean(base_mask)
    target = max(base_cov, min(TARGET_ABSOLUTE_COVERAGE, full_cov * RETAIN_ATTAINABLE_FRACTION))

    def add(action: PanelAction) -> None:
        _apply(mask, action.intervals)
        for seq, width in action.members:
            old = selected_sequences.get(seq)
            if old is not None and old != width:
                raise BackendCompatibilityError(f"width mismatch for selected sequence {seq}")
            selected_sequences[seq] = width
            selected_sources.setdefault(seq, set()).add(action.action_type)
        selected_actions.append(action)

    mask = base_mask.copy()
    cur_t5, cur_t3 = base_t5, base_t3
    remaining = list(actions)
    while (full_t5 and not cur_t5) or (full_t3 and not cur_t3):
        need5 = full_t5 and not cur_t5
        need3 = full_t3 and not cur_t3
        chosen = None
        chosen_key = None
        already = set(selected_sequences)
        for action in remaining:
            covered = int(need5 and action.terminal5) + int(need3 and action.terminal3)
            if covered == 0:
                continue
            cost = _cost(action, already)
            rank = (cost / covered if covered else math.inf, cost, -covered, action.action_type, action.action_id)
            if chosen_key is None or rank < chosen_key:
                chosen_key, chosen = rank, action
        if chosen is None:
            break
        add(chosen)
        cur_t5 = cur_t5 or chosen.terminal5
        cur_t3 = cur_t3 or chosen.terminal3
        remaining.remove(chosen)

    already = set(selected_sequences)
    static_cost = {action.action_id: _cost(action, already) for action in remaining}

    def prefix() -> np.ndarray:
        values = np.empty(len(mask) + 1, dtype=np.int64)
        values[0] = 0
        if len(mask):
            np.cumsum(~mask, dtype=np.int64, out=values[1:])
        return values

    heap: list[tuple[float, str, int]] = []
    by_id: dict[str, PanelAction] = {}
    uncovered = prefix()
    for action in remaining:
        gain = _gain(uncovered, action.intervals)
        if gain <= 0:
            continue
        cost = static_cost[action.action_id]
        score = float("inf") if cost == 0 else gain / cost
        heapq.heappush(heap, (-score, action.action_id, gain))
        by_id[action.action_id] = action
    while _mean(mask) + 1e-15 < target and heap:
        _, action_id, _ = heapq.heappop(heap)
        action = by_id[action_id]
        gain = _gain(prefix(), action.intervals)
        if gain <= 0:
            continue
        cost = static_cost[action_id]
        score = float("inf") if cost == 0 else gain / cost
        next_upper = -heap[0][0] if heap else -1.0
        if score + 1e-18 >= next_upper:
            add(action)
        else:
            heapq.heappush(heap, (-score, action_id, gain))
    return GroupTarget(_mean(mask), cur_t5, cur_t3, len(mask)), full_t5, full_t3


def _optimize_panel(cap, lengths, baseline, baseline_widths, short_actions, h1_actions, read_window=READ_WINDOW):
    selected_sequences = dict(baseline_widths)
    selected_sources = {seq: {"G1B_BASELINE"} for seq in baseline_widths}
    selected_actions: list[PanelAction] = []
    group_actions: dict[Group, list[PanelAction]] = {}
    ordering = []
    for key, length in lengths.items():
        actions = [a for a in short_actions.get(key, ()) if a.members[0][1] <= cap]
        actions.extend(h1_actions.get(key, ()))
        group_actions[key] = actions
        full = _window_mask(length, [(s, k) for s, k, _, _ in baseline.get(key, ())], read_window)
        for action in actions:
            _apply(full, action.intervals)
        ordering.append((_mean(full), -len(full), key))
    targets: dict[Group, GroupTarget] = {}
    terminal_selected = 0
    for _, _, key in sorted(ordering):
        target, _, _ = _optimize_group(
            key, lengths[key], baseline, group_actions[key], selected_sequences, selected_sources, selected_actions,
            read_window,
        )
        targets[key] = target
        terminal_selected += int(target.selected_terminal5) + int(target.selected_terminal3)
    coverages = [t.selected_coverage for t in targets.values()]
    summary = (
        float(np.quantile(np.asarray(coverages), 0.10)),
        statistics.mean(coverages),
        terminal_selected / (2 * len(targets)),
        -sum(selected_sequences.values()),
    )
    return RescuePanel(cap, targets, selected_actions, selected_sequences, selected_sources), summary


def optimize_rescue_panel(
    catalogue: AnchorCatalogue,
    panel: Sequence[PanelQuery],
    short: Sequence[ShortLocus],
    h1: Sequence[H1Locus],
    index,
    read_window: int = READ_WINDOW,
) -> RescuePanel:
    records = catalogue.records
    family_by_te = catalogue.family_by_te()
    baseline = baseline_anchor_loci(catalogue, panel)
    baseline_widths = {entry.query_seq: entry.width for entry in panel}
    lengths: dict[Group, int] = {}
    family_by_group: dict[Group, str] = {}
    for te_id, sequence in sorted(records.items()):
        for orientation in ("F", "RC"):
            lengths[(te_id, orientation)] = len(sequence)
            family_by_group[(te_id, orientation)] = family_by_te.get(te_id, te_id)
    short_actions = _short_actions(short, lengths, family_by_group, read_window)
    h1_actions = _h1_actions(h1, lengths, family_by_group, index, read_window)
    rescue, summary = _optimize_panel(
        RESCUE_PANEL_CAP, lengths, baseline, baseline_widths, short_actions, h1_actions[RESCUE_PANEL_CAP], read_window
    )
    rows = -summary[3]
    if rows > MATERIALIZE_ROW_HARD_CAP:
        raise BackendCompatibilityError(
            f"the {RESCUE_PANEL_CAP}-row coverage panel needs {rows:,} msBWT interval rows, more than the "
            f"{MATERIALIZE_ROW_HARD_CAP:,}-row materialization limit; this TE library has too many frequent words for the "
            "library layer"
        )
    return rescue
