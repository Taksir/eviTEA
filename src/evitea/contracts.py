"""Immutable scientific value objects and deterministic serialization."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum, IntEnum, StrEnum
from pathlib import Path
from types import MappingProxyType
from typing import Any

from .constants import (
    CERTIFYING_PAIRED_END_STATES,
    MIN_FOLDS,
    MIN_PAIRS,
    N_FOLDS,
)
from .errors import ArchitectureInvariantError, InputContractError


def _nonnegative_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ArchitectureInvariantError(f"{name} must be a non-negative integer")


def _positive_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ArchitectureInvariantError(f"{name} must be a positive integer")


def _nonempty_text(value: str, name: str) -> None:
    if not isinstance(value, str) or not value:
        raise ArchitectureInvariantError(f"{name} must be non-empty text")


def _canonical_int_tuple(values: Sequence[int], name: str) -> tuple[int, ...]:
    result: set[int] = set()
    for value in values:
        _nonnegative_int(value, name)
        result.add(value)
    return tuple(sorted(result))


def _canonical_text_tuple(values: Sequence[str], name: str) -> tuple[str, ...]:
    result: set[str] = set()
    for value in values:
        _nonempty_text(value, name)
        result.add(value)
    return tuple(sorted(result))


@dataclass(frozen=True, slots=True, order=True)
class BaseSpan0:
    """Zero-based half-open base or query span ``[start0, end0)``."""

    start0: int
    end0: int

    def __post_init__(self) -> None:
        _nonnegative_int(self.start0, "start0")
        _nonnegative_int(self.end0, "end0")
        if self.end0 < self.start0:
            raise ArchitectureInvariantError("BaseSpan0 end0 must be >= start0")

    @property
    def length(self) -> int:
        return self.end0 - self.start0

    def overlap_length(self, other: BaseSpan0) -> int:
        return max(0, min(self.end0, other.end0) - max(self.start0, other.start0))


@dataclass(frozen=True, slots=True, order=True)
class InterbasePoint0:
    """Zero-based interbase position."""

    position0: int

    def __post_init__(self) -> None:
        _nonnegative_int(self.position0, "position0")

    def require_within(self, contig_length: int) -> InterbasePoint0:
        _nonnegative_int(contig_length, "contig_length")
        if self.position0 > contig_length:
            raise ArchitectureInvariantError("interbase point lies beyond contig length")
        return self


@dataclass(frozen=True, slots=True, order=True)
class ClosedBreakpointInterval0:
    """Closed inclusive interval over interbase positions ``[lo0, hi0]``."""

    lo0: int
    hi0: int

    def __post_init__(self) -> None:
        _nonnegative_int(self.lo0, "lo0")
        _nonnegative_int(self.hi0, "hi0")
        if self.hi0 < self.lo0:
            raise ArchitectureInvariantError(
                "ClosedBreakpointInterval0 hi0 must be >= lo0"
            )

    def contains(self, point: int | InterbasePoint0) -> bool:
        value = point.position0 if isinstance(point, InterbasePoint0) else point
        _nonnegative_int(value, "point")
        return self.lo0 <= value <= self.hi0

    def distance_to(
        self, other: ClosedBreakpointInterval0 | InterbasePoint0
    ) -> int:
        other_lo0 = other.position0 if isinstance(other, InterbasePoint0) else other.lo0
        other_hi0 = other.position0 if isinstance(other, InterbasePoint0) else other.hi0
        if self.hi0 < other_lo0:
            return other_lo0 - self.hi0
        if other_hi0 < self.lo0:
            return self.lo0 - other_hi0
        return 0


class PairIdentityMode(StrEnum):
    EVITEA_ORDINAL0 = "EVITEA_ORDINAL0"
    LEGACY_FROZEN = "LEGACY_FROZEN"


class ReadMate(IntEnum):
    READ1 = 1
    READ2 = 2


class FamilyResolutionState(StrEnum):
    NO_SAFE_EVIDENCE = "NO_SAFE_EVIDENCE"
    RESOLVED_SINGLETON = "RESOLVED_SINGLETON"
    FAMILY_CONFLICT = "FAMILY_CONFLICT"


class TEEligibilityState(StrEnum):
    ELIGIBLE = "ELIGIBLE"
    NO_SAFE_EVIDENCE = "NO_SAFE_EVIDENCE"
    UNSUPPORTED_SAFE_TE_RECORD = "UNSUPPORTED_SAFE_TE_RECORD"
    INSUFFICIENT_INDEPENDENT_BLOCKS = "INSUFFICIENT_INDEPENDENT_BLOCKS"
    BOTH_MATES_TE_POSITIVE = "BOTH_MATES_TE_POSITIVE"
    FAMILY_CONFLICT = "FAMILY_CONFLICT"


class HostMappingState(StrEnum):
    UNASSESSED = "UNASSESSED"
    UNMAPPED = "UNMAPPED"
    NONUNIQUE_OR_LOW = "NONUNIQUE_OR_LOW"
    UNIQUE_HIGH = "UNIQUE_HIGH"


class ReferenceExplanationCode(StrEnum):
    UNASSESSED = "UNASSESSED"
    EVIDENCE_COMPLETE = "EVIDENCE_COMPLETE"
    REFERENCE_UNEXPLAINED = "REFERENCE_UNEXPLAINED"
    NO_LOCAL_PROPER = "NO_LOCAL_PROPER"
    INDETERMINATE_CENSORED = "INDETERMINATE_CENSORED"


class CensorState(StrEnum):
    COMPLETE = "COMPLETE"
    HOST_RESOURCE_CENSORED = "HOST_RESOURCE_CENSORED"
    REFERENCE_RESOURCE_CENSORED = "REFERENCE_RESOURCE_CENSORED"
    HOST_AND_REFERENCE_RESOURCE_CENSORED = "HOST_AND_REFERENCE_RESOURCE_CENSORED"
    SOURCE_INCOMPLETE = "SOURCE_INCOMPLETE"
    MODEL_INCOMPLETE = "MODEL_INCOMPLETE"


class BreakpointSide(StrEnum):
    LEFT = "LEFT"
    RIGHT = "RIGHT"


class PairedEndSupportState(StrEnum):
    BILATERAL = "BILATERAL"
    LEFT_ONLY = "LEFT_ONLY"
    RIGHT_ONLY = "RIGHT_ONLY"
    NEITHER = "NEITHER"


class CoordinateKind(StrEnum):
    INTERBASE_POINT = "INTERBASE_POINT"
    INTERBASE_CLOSED_RANGE = "INTERBASE_CLOSED_RANGE"


class TEEventOrientation(StrEnum):
    FORWARD = "FORWARD"
    REVERSE = "REVERSE"
    UNKNOWN = "UNKNOWN"


class ProbeClass(StrEnum):
    TERMINAL_21MER = "TERMINAL_21MER"
    INTERNAL_21MER = "INTERNAL_21MER"


class SafeProbeRecordStatus(StrEnum):
    MULTI_INDEPENDENT_BLOCK = "MULTI_INDEPENDENT_BLOCK"
    SINGLE_INDEPENDENT_BLOCK = "SINGLE_INDEPENDENT_BLOCK"
    UNSUPPORTED_TERMINAL_SAFE = "UNSUPPORTED_TERMINAL_SAFE"


@dataclass(frozen=True, slots=True)
class ReadRecord:
    header: str
    pair_name: str
    sequence: str
    quality: str
    mate: ReadMate

    def __post_init__(self) -> None:
        _nonempty_text(self.header, "header")
        _nonempty_text(self.pair_name, "pair_name")
        _nonempty_text(self.sequence, "sequence")
        if len(self.sequence) != len(self.quality):
            raise InputContractError("FASTQ sequence and quality lengths differ")


@dataclass(frozen=True, slots=True)
class RawReadPair:
    pair_name: str
    read1: ReadRecord
    read2: ReadRecord

    def __post_init__(self) -> None:
        _nonempty_text(self.pair_name, "pair_name")
        if self.read1.mate is not ReadMate.READ1 or self.read2.mate is not ReadMate.READ2:
            raise ArchitectureInvariantError("RawReadPair mate roles are invalid")
        if self.read1.pair_name != self.pair_name or self.read2.pair_name != self.pair_name:
            raise ArchitectureInvariantError("RawReadPair names are inconsistent")


@dataclass(frozen=True, slots=True)
class PhysicalPair:
    sample_id: str
    raw_pair_ordinal0: int
    pair_id: int
    fold_id: int
    pair_name: str
    read1: ReadRecord
    read2: ReadRecord
    identity_mode: PairIdentityMode = PairIdentityMode.EVITEA_ORDINAL0
    legacy_pair_id: int | None = None

    def __post_init__(self) -> None:
        _nonempty_text(self.sample_id, "sample_id")
        _nonnegative_int(self.raw_pair_ordinal0, "raw_pair_ordinal0")
        _nonnegative_int(self.pair_id, "pair_id")
        _nonnegative_int(self.fold_id, "fold_id")
        if self.fold_id >= N_FOLDS:
            raise ArchitectureInvariantError("fold_id is outside the frozen fold range")
        if self.identity_mode is PairIdentityMode.EVITEA_ORDINAL0:
            if self.pair_id != self.raw_pair_ordinal0:
                raise ArchitectureInvariantError("ordinal-zero pair identity was renumbered")
            if self.legacy_pair_id is not None:
                raise ArchitectureInvariantError("new-data identity cannot carry legacy_pair_id")
        if self.fold_id != self.pair_id % N_FOLDS:
            raise ArchitectureInvariantError("fold_id must equal pair_id % 5")
        if (
            self.read1.mate is not ReadMate.READ1
            or self.read2.mate is not ReadMate.READ2
            or self.read1.pair_name != self.pair_name
            or self.read2.pair_name != self.pair_name
        ):
            raise ArchitectureInvariantError("PhysicalPair read identities are inconsistent")


@dataclass(frozen=True, slots=True)
class FastaRecord:
    record_id: str
    description: str
    sequence: str
    rank: int

    def __post_init__(self) -> None:
        _nonempty_text(self.record_id, "record_id")
        _nonempty_text(self.sequence, "sequence")
        _nonnegative_int(self.rank, "rank")


@dataclass(frozen=True, slots=True)
class ReferenceInput:
    path: Path
    records: tuple[FastaRecord, ...]
    raw_sha256: str
    semantic_sha256: str
    contig_lengths: Mapping[str, int]

    def __post_init__(self) -> None:
        if not self.records:
            raise InputContractError("reference FASTA is empty")
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(
            self,
            "contig_lengths",
            MappingProxyType(dict(self.contig_lengths)),
        )


@dataclass(frozen=True, slots=True)
class TELibraryInput:
    path: Path
    records: tuple[FastaRecord, ...]
    raw_sha256: str
    semantic_sha256: str

    def __post_init__(self) -> None:
        if not self.records:
            raise InputContractError("TE FASTA is empty")
        object.__setattr__(self, "path", Path(self.path))


@dataclass(frozen=True, slots=True)
class RepeatRecord:
    chrom: str
    span0: BaseSpan0
    name: str
    score: str
    strand: str
    source_row0: int


@dataclass(frozen=True, slots=True)
class RepeatContext:
    path: Path
    records: tuple[RepeatRecord, ...]
    raw_sha256: str
    reference_semantic_sha256: str


@dataclass(frozen=True, slots=True)
class TEProbe:
    probe_id: str
    probe_class: ProbeClass
    sequence: str
    family: str
    te_record_id: str
    source_span0: BaseSpan0
    probe_source_strand: str
    independent_block_id: str
    record_safe_status: SafeProbeRecordStatus


@dataclass(frozen=True, slots=True)
class TEProbeHit:
    sample_id: str
    pair_id: int
    mate_no: ReadMate
    probe_id: str
    query_span0: BaseSpan0
    hit_orientation: TEEventOrientation


@dataclass(frozen=True, slots=True)
class TEFamilyEvidence:
    sample_id: str
    pair_id: int
    family_state: FamilyResolutionState
    eligibility_state: TEEligibilityState
    resolved_family: str | None = None
    te_positive_mate: ReadMate | None = None
    families_observed: tuple[str, ...] = ()
    probe_ids: tuple[str, ...] = ()
    independent_block_ids: tuple[str, ...] = ()
    contributing_te_record_ids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class HostPlacement:
    sample_id: str
    pair_id: int
    mate_no: ReadMate
    placement_rank: int
    chrom: str
    span0: BaseSpan0
    host_strand: str
    mapq: int
    cigar: str
    mapping_state: HostMappingState


@dataclass(frozen=True, slots=True)
class ReferenceAlternative:
    sample_id: str
    pair_id: int
    alternative_id: str
    placement_rank: int
    family: str
    chrom: str
    te_alignment_span0: BaseSpan0
    te_cigar: str
    is_local_proper: bool
    evidence_complete: bool
    unexplained_query_spans0: tuple[BaseSpan0, ...] = ()
    censor_state: CensorState = CensorState.COMPLETE
    signature_id: str | None = None
    mate_no: ReadMate | None = None


@dataclass(frozen=True, slots=True)
class ReferenceExplanationState:
    sample_id: str
    pair_id: int
    code: ReferenceExplanationCode
    censor_state: CensorState = CensorState.COMPLETE
    finite_alternative_ids: tuple[str, ...] = ()
    evidence_complete_alternative_ids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class ProposalInterval:
    sample_id: str
    proposal_id: str
    family: str
    chrom: str
    interval0: ClosedBreakpointInterval0
    support_pair_ids: tuple[int, ...]
    support_fold_ids: tuple[int, ...]
    rotation_id: str | None = None
    discovery_fold_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        _nonempty_text(self.sample_id, "proposal sample_id")
        _nonempty_text(self.proposal_id, "proposal_id")
        _nonempty_text(self.family, "proposal family")
        _nonempty_text(self.chrom, "proposal chrom")
        object.__setattr__(
            self, "support_pair_ids", _canonical_int_tuple(self.support_pair_ids, "pair_id")
        )
        object.__setattr__(
            self, "support_fold_ids", _canonical_int_tuple(self.support_fold_ids, "fold_id")
        )
        object.__setattr__(
            self,
            "discovery_fold_ids",
            _canonical_int_tuple(self.discovery_fold_ids, "discovery fold_id"),
        )
        if len(self.support_pair_ids) < MIN_PAIRS or len(self.support_fold_ids) < MIN_FOLDS:
            raise ArchitectureInvariantError(
                "proposal lacks frozen pair/fold support"
            )
        if any(fold_id >= N_FOLDS for fold_id in self.support_fold_ids):
            raise ArchitectureInvariantError("proposal fold_id lies outside frozen folds")
        if any(fold_id >= N_FOLDS for fold_id in self.discovery_fold_ids):
            raise ArchitectureInvariantError("discovery fold_id lies outside frozen folds")
        expected_folds = tuple(
            sorted({pair_id % N_FOLDS for pair_id in self.support_pair_ids})
        )
        if self.support_fold_ids != expected_folds:
            raise ArchitectureInvariantError(
                "proposal folds disagree with physical-pair identities"
            )
        if self.discovery_fold_ids and not set(self.support_fold_ids).issubset(
            self.discovery_fold_ids
        ):
            raise ArchitectureInvariantError(
                "proposal support includes a non-discovery fold"
            )


@dataclass(frozen=True, slots=True)
class ExactJunctionWitness:
    sample_id: str
    pair_id: int
    fold_id: int
    family: str
    chrom: str
    breakpoint0: int
    host_anchor_side: BreakpointSide
    te_clip_side: BreakpointSide
    te_cigar: str
    reference_alternative_id: str
    unexplained_query_spans0: tuple[BaseSpan0, ...]
    evidence_origin: str
    signature_id: str | None = None

    def __post_init__(self) -> None:
        _nonnegative_int(self.pair_id, "junction pair_id")
        _nonnegative_int(self.fold_id, "junction fold_id")
        _nonnegative_int(self.breakpoint0, "junction breakpoint0")
        if self.fold_id != self.pair_id % N_FOLDS:
            raise ArchitectureInvariantError(
                "junction fold_id must equal pair_id modulo the frozen fold count"
            )
        if not self.unexplained_query_spans0:
            raise ArchitectureInvariantError(
                "junction witness requires an unexplained query span"
            )


@dataclass(frozen=True, slots=True)
class ExactJunctionEvent:
    sample_id: str
    family: str
    chrom: str
    breakpoint0: int
    host_anchor_side: BreakpointSide
    support_pair_ids: tuple[int, ...]
    support_fold_ids: tuple[int, ...]
    evidence_origins: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _nonnegative_int(self.breakpoint0, "breakpoint0")
        pair_ids = _canonical_int_tuple(self.support_pair_ids, "junction pair_id")
        fold_ids = _canonical_int_tuple(self.support_fold_ids, "junction fold_id")
        if len(pair_ids) < MIN_PAIRS or len(fold_ids) < MIN_FOLDS:
            raise ArchitectureInvariantError(
                "exact junction event lacks frozen pair/fold support"
            )
        if fold_ids != tuple(sorted({pair_id % N_FOLDS for pair_id in pair_ids})):
            raise ArchitectureInvariantError(
                "exact junction folds disagree with physical-pair identities"
            )
        object.__setattr__(self, "support_pair_ids", pair_ids)
        object.__setattr__(self, "support_fold_ids", fold_ids)
        object.__setattr__(
            self,
            "evidence_origins",
            _canonical_text_tuple(self.evidence_origins, "junction evidence origin"),
        )


@dataclass(frozen=True, slots=True)
class PairedEndAnchorInterval:
    sample_id: str
    pair_id: int
    fold_id: int
    family: str
    chrom: str
    anchor_side: BreakpointSide
    interval0: ClosedBreakpointInterval0
    host_placement_rank: int
    host_strand: str
    reference_code: ReferenceExplanationCode

    def __post_init__(self) -> None:
        _nonnegative_int(self.pair_id, "paired-end pair_id")
        _nonnegative_int(self.fold_id, "paired-end fold_id")
        if self.fold_id != self.pair_id % N_FOLDS:
            raise ArchitectureInvariantError(
                "paired-end fold_id must equal pair_id modulo the frozen fold count"
            )
        if self.host_strand not in {"+", "-"}:
            raise ArchitectureInvariantError("paired-end host strand must be + or -")


@dataclass(frozen=True, slots=True)
class PairedEndAnchorSupport:
    sample_id: str
    family: str
    chrom: str
    breakpoint0: int
    excluded_exact_pair_ids: tuple[int, ...]
    left_pair_ids: tuple[int, ...]
    left_fold_ids: tuple[int, ...]
    right_pair_ids: tuple[int, ...]
    right_fold_ids: tuple[int, ...]
    support_state: PairedEndSupportState

    def __post_init__(self) -> None:
        _nonnegative_int(self.breakpoint0, "paired-end breakpoint0")
        for field_name in (
            "excluded_exact_pair_ids",
            "left_pair_ids",
            "left_fold_ids",
            "right_pair_ids",
            "right_fold_ids",
        ):
            object.__setattr__(
                self,
                field_name,
                _canonical_int_tuple(getattr(self, field_name), field_name),
            )
        if set(self.excluded_exact_pair_ids) & (
            set(self.left_pair_ids) | set(self.right_pair_ids)
        ):
            raise ArchitectureInvariantError(
                "exact-junction pairs cannot contribute paired-end support"
            )
        if self.left_fold_ids != tuple(
            sorted({pair_id % N_FOLDS for pair_id in self.left_pair_ids})
        ) or self.right_fold_ids != tuple(
            sorted({pair_id % N_FOLDS for pair_id in self.right_pair_ids})
        ):
            raise ArchitectureInvariantError(
                "paired-end folds disagree with physical-pair identities"
            )
        left_passes = (
            len(self.left_pair_ids) >= MIN_PAIRS
            and len(self.left_fold_ids) >= MIN_FOLDS
        )
        right_passes = (
            len(self.right_pair_ids) >= MIN_PAIRS
            and len(self.right_fold_ids) >= MIN_FOLDS
        )
        expected_state = (
            PairedEndSupportState.BILATERAL
            if left_passes and right_passes
            else PairedEndSupportState.LEFT_ONLY
            if left_passes
            else PairedEndSupportState.RIGHT_ONLY
            if right_passes
            else PairedEndSupportState.NEITHER
        )
        if self.support_state is not expected_state:
            raise ArchitectureInvariantError(
                "paired-end support state disagrees with frozen pair/fold thresholds"
            )


@dataclass(frozen=True, slots=True)
class CertifiedEvent:
    sample_id: str
    family: str
    chrom: str
    breakpoint0: int
    host_anchor_side: BreakpointSide
    exact_junction: ExactJunctionEvent
    paired_end_support: PairedEndAnchorSupport
    proposal_ids: tuple[str, ...]
    certification_rotation_ids: tuple[str, ...] = ()
    call_origin: str = "EVIDENCE_EXTENSION"

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "proposal_ids", _canonical_text_tuple(self.proposal_ids, "proposal_id")
        )
        object.__setattr__(
            self,
            "certification_rotation_ids",
            _canonical_text_tuple(
                self.certification_rotation_ids, "certification rotation_id"
            ),
        )
        if not self.proposal_ids:
            raise ArchitectureInvariantError(
                "certified event requires proposal provenance"
            )
        if (
            self.paired_end_support.support_state.value
            not in CERTIFYING_PAIRED_END_STATES
        ):
            raise ArchitectureInvariantError(
                "certified event requires positive paired-end anchor support"
            )
        if (
            self.exact_junction.sample_id != self.sample_id
            or self.exact_junction.family != self.family
            or self.exact_junction.chrom != self.chrom
            or self.exact_junction.breakpoint0 != self.breakpoint0
            or self.exact_junction.host_anchor_side is not self.host_anchor_side
            or self.paired_end_support.sample_id != self.sample_id
            or self.paired_end_support.family != self.family
            or self.paired_end_support.chrom != self.chrom
            or self.paired_end_support.breakpoint0 != self.breakpoint0
        ):
            raise ArchitectureInvariantError(
                "support evidence mutated the certified exact key"
            )
        if set(self.exact_junction.support_pair_ids) & (
            set(self.paired_end_support.left_pair_ids)
            | set(self.paired_end_support.right_pair_ids)
        ):
            raise ArchitectureInvariantError(
                "certified exact and paired-end evidence share a physical pair"
            )
        if self.paired_end_support.excluded_exact_pair_ids != (
            self.exact_junction.support_pair_ids
        ):
            raise ArchitectureInvariantError(
                "paired-end exclusion set must equal exact-junction support pairs"
            )


@dataclass(frozen=True, slots=True)
class RepeatAwareEvent:
    sample_id: str
    legacy_event_id: str
    family: str
    chrom: str
    interval0: ClosedBreakpointInterval0
    left_coord0: int | None
    right_coord0: int | None
    left_support_pairs: int
    right_support_pairs: int
    total_support_pairs: int
    fragment_pair_ids: tuple[int, ...]
    evidence_type: str
    call_origin: str = "REPEAT_AWARE_INCUMBENT"

    def __post_init__(self) -> None:
        for value, field in (
            (self.sample_id, "sample_id"),
            (self.legacy_event_id, "legacy_event_id"),
            (self.family, "family"),
            (self.chrom, "chrom"),
            (self.evidence_type, "evidence_type"),
            (self.call_origin, "call_origin"),
        ):
            _nonempty_text(value, field)
        if not isinstance(self.interval0, ClosedBreakpointInterval0):
            raise ArchitectureInvariantError(
                "repeat-aware interval must be a closed interbase interval"
            )
        for value, field in (
            (self.left_coord0, "left_coord0"),
            (self.right_coord0, "right_coord0"),
        ):
            if value is not None:
                _nonnegative_int(value, field)
        for value, field in (
            (self.left_support_pairs, "left_support_pairs"),
            (self.right_support_pairs, "right_support_pairs"),
            (self.total_support_pairs, "total_support_pairs"),
        ):
            _nonnegative_int(value, field)
        pair_ids = _canonical_int_tuple(self.fragment_pair_ids, "fragment pair_id")
        if pair_ids != self.fragment_pair_ids:
            raise ArchitectureInvariantError(
                "repeat-aware fragment pair IDs must be unique and sorted"
            )
        if self.total_support_pairs != len(pair_ids):
            raise ArchitectureInvariantError(
                "repeat-aware total support must count distinct fragment pair IDs"
            )
        if (
            self.left_support_pairs < MIN_PAIRS
            or self.right_support_pairs < MIN_PAIRS
            or self.left_support_pairs > self.total_support_pairs
            or self.right_support_pairs > self.total_support_pairs
            or self.total_support_pairs
            > self.left_support_pairs + self.right_support_pairs
        ):
            raise ArchitectureInvariantError(
                "repeat-aware event violates bilateral physical-pair support"
            )


@dataclass(frozen=True, slots=True)
class UnresolvedCompetition:
    sample_id: str
    unresolved_id: str
    stage: str
    reason: str
    chrom: str | None = None
    interval0: ClosedBreakpointInterval0 | None = None
    families: tuple[str, ...] = ()
    pair_ids: tuple[int, ...] = ()
    censor_state: CensorState = CensorState.COMPLETE
    participant_ids: tuple[str, ...] = ()
    detail: Mapping[str, Any] = MappingProxyType({})

    def __post_init__(self) -> None:
        for value, name in (
            (self.sample_id, "unresolved sample_id"),
            (self.unresolved_id, "unresolved_id"),
            (self.stage, "unresolved stage"),
            (self.reason, "unresolved reason"),
        ):
            _nonempty_text(value, name)
        object.__setattr__(
            self, "families", _canonical_text_tuple(self.families, "unresolved family")
        )
        object.__setattr__(
            self, "pair_ids", _canonical_int_tuple(self.pair_ids, "unresolved pair_id")
        )
        object.__setattr__(
            self,
            "participant_ids",
            _canonical_text_tuple(
                self.participant_ids, "unresolved participant_id"
            ),
        )
        object.__setattr__(self, "detail", MappingProxyType(dict(self.detail)))


@dataclass(frozen=True, slots=True)
class FinalCall:
    sample_id: str
    call_id: str
    chrom: str
    family: str
    coordinate_kind: CoordinateKind
    interval0: ClosedBreakpointInterval0
    exact_breakpoint0: int | None
    host_anchor_side: BreakpointSide | None
    te_event_orientation: TEEventOrientation
    call_origin: str
    paired_end_support_state: PairedEndSupportState | None = None
    exact_support_pair_ids: tuple[int, ...] = ()
    left_anchor_pair_ids: tuple[int, ...] = ()
    right_anchor_pair_ids: tuple[int, ...] = ()
    all_support_pair_ids: tuple[int, ...] = ()
    legacy_event_id: str | None = None

    def __post_init__(self) -> None:
        for value, name in (
            (self.sample_id, "final call sample_id"),
            (self.call_id, "final call_id"),
            (self.chrom, "final call chrom"),
            (self.family, "final call family"),
            (self.call_origin, "final call origin"),
        ):
            _nonempty_text(value, name)
        for field_name in (
            "exact_support_pair_ids",
            "left_anchor_pair_ids",
            "right_anchor_pair_ids",
            "all_support_pair_ids",
        ):
            object.__setattr__(
                self,
                field_name,
                _canonical_int_tuple(getattr(self, field_name), field_name),
            )
        if self.coordinate_kind is CoordinateKind.INTERBASE_POINT:
            if (
                self.exact_breakpoint0 is None
                or self.interval0.lo0 != self.interval0.hi0
                or self.interval0.lo0 != self.exact_breakpoint0
            ):
                raise ArchitectureInvariantError("exact call coordinates are inconsistent")
        elif self.exact_breakpoint0 is not None:
            raise ArchitectureInvariantError("range call cannot invent an exact breakpoint")


@dataclass(frozen=True, slots=True)
class ProvenanceRecord:
    schema_version: str
    architecture_id: str
    values: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", MappingProxyType(dict(self.values)))


def stable_to_primitive(value: Any) -> Any:
    """Convert supported typed values into deterministic JSON primitives."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Enum):
        return stable_to_primitive(value.value)
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: stable_to_primitive(getattr(value, field.name))
            for field in fields(value)
        }
    if isinstance(value, Mapping):
        result: dict[str, Any] = {}
        for key in sorted(value):
            if not isinstance(key, str):
                raise ArchitectureInvariantError("stable mappings require string keys")
            result[key] = stable_to_primitive(value[key])
        return result
    if isinstance(value, (set, frozenset)):
        converted = [stable_to_primitive(item) for item in value]
        return sorted(
            converted,
            key=lambda item: json.dumps(item, sort_keys=True, separators=(",", ":")),
        )
    if isinstance(value, (tuple, list)):
        return [stable_to_primitive(item) for item in value]
    raise ArchitectureInvariantError(
        f"unsupported stable-serialization type: {type(value).__name__}"
    )


def stable_json_dumps(value: Any) -> str:
    """Serialize a supported value as compact, key-sorted UTF-8 JSON text."""

    return json.dumps(
        stable_to_primitive(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
