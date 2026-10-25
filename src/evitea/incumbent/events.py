"""Typed library layer calls and their combined-call-table rows."""

from __future__ import annotations

import csv
import io
from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum

from ..errors import ArchitectureInvariantError

LEGACY_EVENT_FIELDS = (
    "event_id", "chr", "start", "end", "family", "families",
    "left_coord", "right_coord", "delta", "left_support",
    "right_support", "total_support", "evidence_type", "score",
    "fragment_ids",
)

__all__ = ["IncumbentCallClass", "IncumbentCall", "LEGACY_EVENT_FIELDS", "legacy_rows", "legacy_tsv_bytes"]


class IncumbentCallClass(StrEnum):
    SPLIT_JUNCTION = "SPLIT_JUNCTION"
    MATE_RESCUE = "MATE_RESCUE"
    REPEAT_AWARE_RESCUE = "REPEAT_AWARE_RESCUE"


@dataclass(frozen=True, slots=True)
class IncumbentCall:
    """One call of the library layer.

    ``start``/``end`` are the interbase envelope. Split-junction calls carry
    exact side coordinates (``-1`` = side absent); range classes carry none.
    ``family`` is the emitted label (raw, not alias-normalized).
    """

    sample_id: str
    legacy_event_id: str
    call_class: IncumbentCallClass
    chrom: str
    start: int
    end: int
    family: str
    left_coord: int | None
    right_coord: int | None
    delta: int | None
    left_support: int
    right_support: int
    evidence_type: str
    score: float
    fragment_pair_ids: tuple[int, ...]

    def __post_init__(self) -> None:
        if tuple(sorted(set(self.fragment_pair_ids))) != self.fragment_pair_ids:
            raise ArchitectureInvariantError("fragment pair IDs must be sorted and unique")
        if self.end < self.start:
            raise ArchitectureInvariantError("envelope end precedes start")
        exact = self.call_class is IncumbentCallClass.SPLIT_JUNCTION
        if exact != (self.left_coord is not None and self.right_coord is not None and self.delta is not None):
            raise ArchitectureInvariantError("only split-junction calls carry exact side coordinates")

    @property
    def total_support(self) -> int:
        return len(self.fragment_pair_ids)

    def legacy_row(self) -> dict[str, str]:
        def text(value: object) -> str:
            return "" if value is None else str(value)

        return {
            "event_id": self.legacy_event_id,
            "chr": self.chrom,
            "start": str(self.start),
            "end": str(self.end),
            "family": self.family,
            "families": repr([self.family]) if self.call_class is IncumbentCallClass.SPLIT_JUNCTION else str([self.family]),
            "left_coord": text(self.left_coord),
            "right_coord": text(self.right_coord),
            "delta": text(self.delta),
            "left_support": str(self.left_support),
            "right_support": str(self.right_support),
            "total_support": str(self.total_support),
            "evidence_type": self.evidence_type,
            "score": str(float(self.score)),
            "fragment_ids": ",".join(str(value) for value in self.fragment_pair_ids),
        }


def legacy_rows(calls: Iterable[IncumbentCall]) -> list[dict[str, str]]:
    return [call.legacy_row() for call in calls]


def legacy_tsv_bytes(calls: Iterable[IncumbentCall]) -> bytes:
    """Return the ``*_COMBINED_CALLS.tsv`` content as bytes (LF line endings)."""

    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=list(LEGACY_EVENT_FIELDS), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in legacy_rows(calls):
        writer.writerow(row)
    return handle.getvalue().encode("utf-8")
