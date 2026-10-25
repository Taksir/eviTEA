"""Backend observation models and structural protocols.

Backends expose raw observations only. Interpretation stays in the ``evidence``
and ``adjudication`` packages.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

from ..contracts import ReadMate
from ..errors import InputContractError


@dataclass(frozen=True, slots=True, order=True)
class ReadIdentity:
    """Immutable EviTEA identity assigned before backend indexing."""

    backend_read_id: int
    pair_id: int
    mate: ReadMate

    def __post_init__(self) -> None:
        for value, field in (
            (self.backend_read_id, "backend_read_id"),
            (self.pair_id, "pair_id"),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise InputContractError(f"{field} must be a non-negative integer")
        if not isinstance(self.mate, ReadMate):
            raise InputContractError("mate must be ReadMate.READ1 or READ2")


@dataclass(frozen=True, slots=True)
class RawReadHit:
    """One exact-query occurrence, including multiplicity and provenance."""

    sequence: str
    backend_row: int
    backend_read_id: int
    pair_id: int
    mate: ReadMate
    occurrence_ordinal0: int


@dataclass(frozen=True, slots=True)
class MapperRead:
    """One query sequence provided to a host mapper."""

    query_name: str
    sequence: str
    pair_id: int
    mate: ReadMate
    quality: str | None = None


@dataclass(frozen=True, slots=True)
class RawAlignment:
    """A single parsed SAM observation with zero-based host position."""

    query_name: str
    pair_id: int
    mate: ReadMate
    flag: int
    reference_name: str | None
    position0: int | None
    mapq: int
    cigar: str | None
    sequence: str
    tags: tuple[str, ...] = ()
    mate_reference_name: str | None = None
    mate_position0: int | None = None
    template_length: int = 0

    @property
    def is_mapped(self) -> bool:
        return not bool(self.flag & 0x4)

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & 0x10)

    @property
    def is_secondary(self) -> bool:
        return bool(self.flag & 0x100)

    @property
    def is_supplementary(self) -> bool:
        return bool(self.flag & 0x800)


@runtime_checkable
class ExactQueryBackend(Protocol):
    """Required exact substring-query surface."""

    def build(self, reads: Iterable[MapperRead], output_dir: str | Path) -> None: ...

    def count_exact(self, sequence: str) -> int: ...

    def locate_reads(self, sequence: str) -> Iterable[RawReadHit]: ...


@runtime_checkable
class HostMapperBackend(Protocol):
    """Required host-mapper observation surface."""

    def map_reads(
        self,
        reads: Sequence[MapperRead],
        prepared_reference: object,
        *,
        max_reported: int = 2,
    ) -> Iterable[RawAlignment]: ...
