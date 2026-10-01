"""Exact projection of safe probes into immutable physical reads."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from types import MappingProxyType

from ..backends.protocols import ExactQueryBackend, RawReadHit
from ..contracts import BaseSpan0, ReadMate, TEEventOrientation, TEProbe, TEProbeHit
from ..errors import BackendCompatibilityError, InputContractError
from .codebook import reverse_complement


@dataclass(frozen=True, slots=True)
class SampleReadIndex:
    sample_id: str
    sequences: Mapping[tuple[int, ReadMate], str]

    def __post_init__(self) -> None:
        if not self.sample_id:
            raise InputContractError("sample_id must be non-empty text")
        normalized: dict[tuple[int, ReadMate], str] = {}
        for key, sequence in self.sequences.items():
            if (
                not isinstance(key, tuple)
                or len(key) != 2
                or isinstance(key[0], bool)
                or not isinstance(key[0], int)
                or key[0] < 0
                or not isinstance(key[1], ReadMate)
            ):
                raise InputContractError("sample read key must be (pair_id, ReadMate)")
            value = sequence.upper()
            if not value:
                raise InputContractError("sample read sequence must be non-empty")
            normalized[key] = value
        object.__setattr__(self, "sequences", MappingProxyType(normalized))


def _occurrences(sequence: str, query: str) -> tuple[BaseSpan0, ...]:
    spans: list[BaseSpan0] = []
    start = sequence.find(query)
    while start >= 0:
        spans.append(BaseSpan0(start, start + len(query)))
        start = sequence.find(query, start + 1)
    return tuple(spans)


def _project_orientation(
    sample_index: SampleReadIndex,
    probe: TEProbe,
    query: str,
    orientation: TEEventOrientation,
    backend: ExactQueryBackend,
) -> Iterable[TEProbeHit]:
    hits = tuple(backend.locate_reads(query))
    if backend.count_exact(query) != len(hits):
        raise BackendCompatibilityError("exact-query count/locate multiplicity mismatch")
    grouped: dict[tuple[int, ReadMate], list[RawReadHit]] = defaultdict(list)
    for hit in hits:
        if hit.sequence.upper() != query:
            raise BackendCompatibilityError("exact-query backend returned another query")
        grouped[(hit.pair_id, hit.mate)].append(hit)
    for identity in sorted(grouped, key=lambda item: (item[0], int(item[1]))):
        try:
            read_sequence = sample_index.sequences[identity]
        except KeyError as exc:
            raise BackendCompatibilityError(
                f"exact-query hit has no immutable input read: {identity}"
            ) from exc
        spans = _occurrences(read_sequence, query)
        if len(spans) != len(grouped[identity]):
            raise BackendCompatibilityError(
                "exact-query occurrence multiplicity disagrees with input read"
            )
        for span in spans:
            yield TEProbeHit(
                sample_id=sample_index.sample_id,
                pair_id=identity[0],
                mate_no=identity[1],
                probe_id=probe.probe_id,
                query_span0=span,
                hit_orientation=orientation,
            )


def project_probes(
    sample_index: SampleReadIndex,
    probes: Iterable[TEProbe],
    backend: ExactQueryBackend,
) -> Iterable[TEProbeHit]:
    """Project both orientations and preserve occurrence, not support, rows."""

    for probe in sorted(probes, key=lambda item: item.probe_id):
        forward = probe.sequence.upper()
        yield from _project_orientation(
            sample_index, probe, forward, TEEventOrientation.FORWARD, backend
        )
        reverse = reverse_complement(forward)
        if reverse != forward:
            yield from _project_orientation(
                sample_index, probe, reverse, TEEventOrientation.REVERSE, backend
            )
