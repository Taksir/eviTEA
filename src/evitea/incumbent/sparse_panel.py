"""Select the sparse exact-anchor panel and route reads and pairs through it.

For every anchor locus (record, start, region, orientation) the narrowest
nonempty exact query of width <= 64 is chosen (ties: not low-complexity,
longer k, query id).  Families observed in the sample but without such a query
receive up to four fallback queries of width 65..1024.  Every msBWT row of the
panel is enumerated.  A read inherits the channel mask (terminal, internal or
fallback bits) and the families of all panel queries it carries.  A pair
inherits the union over its mates.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..errors import BackendCompatibilityError
from .anchor_census import AnchorCatalogue, AnchorQuery

__all__ = [
    "CHANNEL_TERMINAL",
    "CHANNEL_INTERNAL",
    "CHANNEL_FALLBACK",
    "PanelQuery",
    "PanelPair",
    "select_sparse_panel",
    "route_sparse_panel",
]

CORE_WIDTH_CAP = 64
FALLBACK_WIDTH_CAP = 1024
FALLBACK_PER_FAMILY = 4
PANEL_ROW_HARD_CAP = 2_000_000
CHANNEL_TERMINAL = 1
CHANNEL_INTERNAL = 2
CHANNEL_FALLBACK = 4


@dataclass(slots=True)
class PanelQuery:
    query_id: str
    query_seq: str
    k: int
    width: int
    low_complexity: bool
    channel_mask: int
    source_families: set[str]


@dataclass(frozen=True, slots=True)
class PanelPair:
    pair_id: int
    mate1_hit: bool
    mate2_hit: bool
    channel_mask: int
    families: tuple[str, ...]


def _region_mask(region: str) -> int:
    if region in ("5p_terminal", "3p_terminal"):
        return CHANNEL_TERMINAL
    if region == "internal":
        return CHANNEL_INTERNAL
    return 0


def select_sparse_panel(catalogue: AnchorCatalogue) -> list[PanelQuery]:
    queries = catalogue.by_id()
    all_families: set[str] = set()
    observed: set[str] = set()
    for query in catalogue.queries:
        all_families.update(query.source_families)
        if query.width > 0:
            observed.update(query.source_families)

    best: dict[tuple[str, int, str, str], tuple[tuple, object, AnchorQuery]] = {}
    for source in catalogue.sources:
        if not source.queryable:
            continue
        query = queries.get(source.query_id)
        if query is None or query.width <= 0 or query.width > CORE_WIDTH_CAP:
            continue
        key = (source.te_id, source.start, source.region, source.orientation)
        score = (query.width, int(query.low_complexity), -query.k, query.query_id)
        old = best.get(key)
        if old is None or score < old[0]:
            best[key] = (score, source, query)

    selected: dict[str, PanelQuery] = {}

    def item(query: AnchorQuery) -> PanelQuery:
        return selected.setdefault(
            query.query_id,
            PanelQuery(query.query_id, query.query_seq, query.k, query.width, query.low_complexity, 0, set()),
        )

    for _, source, query in best.values():
        entry = item(query)
        entry.source_families.add(source.family)
        entry.channel_mask |= _region_mask(source.region)

    core_families = set()
    for entry in selected.values():
        core_families.update(entry.source_families)
    for family in sorted(observed - core_families):
        candidates = sorted(
            (
                (int(q.low_complexity), q.width, -q.k, q.query_id), q)
            for q in catalogue.queries
            if family in q.source_families and CORE_WIDTH_CAP < q.width <= FALLBACK_WIDTH_CAP
        )
        for _, query in candidates[:FALLBACK_PER_FAMILY]:
            entry = item(query)
            entry.channel_mask |= CHANNEL_FALLBACK
            entry.source_families.add(family)
            if query.has_terminal:
                entry.channel_mask |= CHANNEL_TERMINAL
            if query.has_internal:
                entry.channel_mask |= CHANNEL_INTERNAL

    panel = [selected[qid] for qid in sorted(selected)]
    if sum(entry.width for entry in panel) > PANEL_ROW_HARD_CAP:
        raise BackendCompatibilityError("sparse panel exceeds the 2,000,000-row enumeration cap")
    return panel


def route_sparse_panel(panel: Sequence[PanelQuery], index) -> list[PanelPair]:
    read_mask: dict[int, int] = {}
    read_families: dict[int, set[str]] = {}
    for entry in panel:
        for dollar in index.dollar_ids(entry.query_seq).tolist():
            read_mask[dollar] = read_mask.get(dollar, 0) | entry.channel_mask
            read_families.setdefault(dollar, set()).update(entry.source_families)
    dollars = np.asarray(sorted(read_mask), dtype=np.uint64)
    files, pairs = index.file_and_pair(dollars)
    pair_mask: dict[int, int] = {}
    pair_mates: dict[int, int] = {}
    pair_families: dict[int, set[str]] = {}
    for dollar, file_id, pair_id in zip(dollars.tolist(), files.tolist(), pairs.tolist()):
        pair_mask[pair_id] = pair_mask.get(pair_id, 0) | read_mask[dollar]
        pair_mates[pair_id] = pair_mates.get(pair_id, 0) | (1 << file_id)
        pair_families.setdefault(pair_id, set()).update(read_families[dollar])
    return [
        PanelPair(
            pair_id,
            bool(pair_mates[pair_id] & 1),
            bool(pair_mates[pair_id] & 2),
            pair_mask[pair_id],
            tuple(sorted(pair_families[pair_id])),
        )
        for pair_id in sorted(pair_mask)
    ]
