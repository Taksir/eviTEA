"""Exact TE anchor catalogue and its sample interval census (lineage stage G1A).

Every TE consensus record contributes k-mers (k = 21/25/31/35) from its two
120 bp terminal windows (stride 5, last start always kept), its interior
(stride 50 starting at 120, last interior start always kept) or, for records
too short to have an interior, every remaining start.  Each source k-mer is
queried in both orientations.  Distinct ``(k, query)`` strings are numbered
``Q0000001...`` in sorted order and counted exactly in the sample msBWT.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from .sequence import revcomp

__all__ = [
    "ANCHOR_K_VALUES",
    "AnchorQuery",
    "AnchorSource",
    "AnchorCatalogue",
    "parse_te_consensus",
    "load_te_family_map",
    "build_anchor_catalogue",
    "census_anchor_intervals",
]

ANCHOR_K_VALUES = (21, 25, 31, 35)
TERMINAL_WINDOW = 120
TERMINAL_STRIDE = 5
INTERNAL_STRIDE = 50
_DNA = frozenset("ACGT")


@dataclass(frozen=True, slots=True)
class AnchorSource:
    """One (record, k, start, orientation) exact anchor; ``query_id`` empty if unqueryable."""

    query_id: str
    te_id: str
    family: str
    k: int
    start: int
    region: str
    orientation: str
    query_seq: str
    queryable: bool


@dataclass(slots=True)
class AnchorQuery:
    query_id: str
    k: int
    query_seq: str
    source_te_ids: tuple[str, ...]
    source_families: tuple[str, ...]
    source_regions: tuple[str, ...]
    has_terminal: bool
    has_internal: bool
    low_complexity: bool
    interval_lo: int = -1
    interval_hi: int = -1

    @property
    def width(self) -> int:
        return self.interval_hi - self.interval_lo

    @property
    def source_family_count(self) -> int:
        return len(self.source_families)


@dataclass(frozen=True, slots=True)
class AnchorCatalogue:
    records: Mapping[str, str]
    queries: tuple[AnchorQuery, ...]
    sources: tuple[AnchorSource, ...]

    def by_id(self) -> dict[str, AnchorQuery]:
        return {query.query_id: query for query in self.queries}

    def family_by_te(self) -> dict[str, str]:
        return {source.te_id: source.family for source in self.sources}


def parse_te_consensus(path: str | Path) -> dict[str, str]:
    """TE consensus FASTA keyed by the first token of each header (before ``|``)."""

    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with Path(path).open("rt", encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper()
                name = line[1:].split()[0].split("|")[0]
                if not name:
                    raise ValueError("empty FASTA identifier")
                if name in records:
                    raise ValueError(f"duplicate FASTA identifier: {name}")
                parts = []
            else:
                if name is None:
                    raise ValueError("sequence before FASTA header")
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    if not records:
        raise ValueError("empty TE consensus FASTA")
    return records


def load_te_family_map(path: str | Path) -> dict[str, str]:
    """Raw ``FBte_ID -> FlyBase_Family`` labels (no alias normalization)."""

    import csv

    mapping: dict[str, str] = {}
    with Path(path).open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "FBte_ID" not in reader.fieldnames or "FlyBase_Family" not in reader.fieldnames:
            raise ValueError("family map must contain FBte_ID and FlyBase_Family")
        for row in reader:
            te_id = (row.get("FBte_ID") or "").strip()
            family = (row.get("FlyBase_Family") or "").strip()
            if te_id and family:
                mapping[te_id] = family
    return mapping


def _terminal_positions(length: int, k: int, side: str) -> list[int]:
    if length < k:
        return []
    if side == "5p":
        last = min(length, TERMINAL_WINDOW) - k
        if last < 0:
            return []
        positions = set(range(0, last + 1, TERMINAL_STRIDE))
        positions.add(0)
    else:
        window_start = max(0, length - TERMINAL_WINDOW)
        last = length - k
        if last < window_start:
            return []
        positions = set(range(window_start, last + 1, TERMINAL_STRIDE))
        positions.add(length - k)
    return sorted(positions)


def _source_positions(sequence: str, k: int) -> list[tuple[int, str]]:
    length = len(sequence)
    if length < k:
        return []
    five = _terminal_positions(length, k, "5p")
    three = _terminal_positions(length, k, "3p")
    five_set = set(five)
    terminal = five_set | set(three)
    out = [(p, "5p_terminal") for p in five] + [(p, "3p_terminal") for p in three if p not in five_set]
    lo, hi = TERMINAL_WINDOW, length - TERMINAL_WINDOW - k
    if lo <= hi:
        interior = set(range(lo, hi + 1, INTERNAL_STRIDE))
        interior.add(hi)
        out.extend((p, "internal") for p in sorted(interior) if p not in terminal)
    else:
        out.extend((p, "short_internal_fallback") for p in range(0, length - k + 1) if p not in terminal)
    seen: set[int] = set()
    unique: list[tuple[int, str]] = []
    for position, region in out:
        if position not in seen:
            seen.add(position)
            unique.append((position, region))
    return unique


def _low_complexity(sequence: str) -> bool:
    counts = Counter(sequence)
    n = len(sequence)
    entropy = -sum((v / n) * math.log2(v / n) for v in counts.values()) if n else 0.0
    best = run = 1 if sequence else 0
    for a, b in zip(sequence, sequence[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return entropy < 1.30 or best * 2 >= n or len(counts) <= 2


def build_anchor_catalogue(records: Mapping[str, str], family_map: Mapping[str, str]) -> AnchorCatalogue:
    raw_sources: list[tuple[str, str, int, int, str, str, str, bool]] = []
    for te_id in sorted(records):
        sequence = records[te_id]
        family = family_map.get(te_id, te_id)
        for k in ANCHOR_K_VALUES:
            for start, region in _source_positions(sequence, k):
                source = sequence[start : start + k]
                queryable = set(source) <= _DNA
                raw_sources.append((te_id, family, k, start, region, "F", source if queryable else source, queryable))
                raw_sources.append((te_id, family, k, start, region, "RC", revcomp(source) if queryable else "", queryable))
    provenance: dict[tuple[int, str], list[int]] = defaultdict(list)
    for index, (_, _, k, _, _, _, query, queryable) in enumerate(raw_sources):
        if queryable:
            provenance[(k, query)].append(index)
    queries: list[AnchorQuery] = []
    query_ids: dict[tuple[int, str], str] = {}
    for number, key in enumerate(sorted(provenance), 1):
        k, query = key
        rows = [raw_sources[i] for i in provenance[key]]
        query_id = f"Q{number:07d}"
        query_ids[key] = query_id
        queries.append(
            AnchorQuery(
                query_id=query_id,
                k=k,
                query_seq=query,
                source_te_ids=tuple(sorted({row[0] for row in rows})),
                source_families=tuple(sorted({row[1] for row in rows})),
                source_regions=tuple(sorted({row[4] for row in rows})),
                has_terminal=any("terminal" in row[4] for row in rows),
                has_internal=any("internal" in row[4] for row in rows),
                low_complexity=_low_complexity(query),
            )
        )
    sources = tuple(
        AnchorSource(
            query_ids.get((k, query), "") if queryable else "",
            te_id, family, k, start, region, orientation, query, queryable,
        )
        for te_id, family, k, start, region, orientation, query, queryable in raw_sources
    )
    return AnchorCatalogue(dict(records), tuple(queries), sources)


def census_anchor_intervals(catalogue: AnchorCatalogue, index) -> None:
    """Attach exact sample msBWT intervals to every catalogue query (in place)."""

    for query in sorted(catalogue.queries, key=lambda q: (q.k, q.query_seq)):
        query.interval_lo, query.interval_hi = index.interval(query.query_seq)
