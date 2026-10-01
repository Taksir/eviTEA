"""Mapper-observation classification for the frozen candidate contract."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from ..contracts import HostMappingState, ReadMate
from ..errors import InputContractError
from .protocols import RawAlignment


@dataclass(frozen=True, slots=True)
class HostMappingSummary:
    pair_id: int
    mate: ReadMate
    state: HostMappingState
    reported_placements: int
    alignments: tuple[RawAlignment, ...]


def classify_host_mapping(
    alignments: Iterable[RawAlignment],
    *,
    min_mapq: int = 20,
) -> HostMappingSummary:
    """Classify one read's complete ``-k 2`` SAM observation set."""

    if isinstance(min_mapq, bool) or not isinstance(min_mapq, int) or min_mapq < 0:
        raise InputContractError("min_mapq must be a non-negative integer")
    rows = tuple(alignments)
    if not rows:
        raise InputContractError("host mapping classification requires SAM observations")
    pair_mates = {(row.pair_id, row.mate) for row in rows}
    if len(pair_mates) != 1:
        raise InputContractError("host mapping rows mix physical read identities")
    pair_id, mate = next(iter(pair_mates))
    placements = tuple(
        row
        for row in rows
        if row.is_mapped and not row.is_supplementary
    )
    if not placements:
        state = HostMappingState.UNMAPPED
    elif len(placements) == 1 and not placements[0].is_secondary and placements[0].mapq >= min_mapq:
        state = HostMappingState.UNIQUE_HIGH
    else:
        state = HostMappingState.NONUNIQUE_OR_LOW
    return HostMappingSummary(
        pair_id=pair_id,
        mate=mate,
        state=state,
        reported_placements=len(placements),
        alignments=rows,
    )
