"""Host-mate mapping classification from raw mapper alignments."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from ..assets.reference import PreparedReference
from ..backends.host_mapper import classify_host_mapping as _classify_summary
from ..backends.protocols import HostMapperBackend, MapperRead, RawAlignment
from ..contracts import (
    BaseSpan0,
    CensorState,
    HostMappingState,
    HostPlacement,
    PhysicalPair,
    ReadMate,
    TEEligibilityState,
    TEFamilyEvidence,
)
from ..errors import InputContractError
from .reference import cigar_reference_length


@dataclass(frozen=True, slots=True)
class HostMappingResult:
    sample_id: str
    pair_id: int
    mate_no: ReadMate
    state: HostMappingState
    placements: tuple[HostPlacement, ...]
    censor_state: CensorState = CensorState.COMPLETE


def classify_host_mapping(
    raw_alignments: Iterable[RawAlignment],
    *,
    min_mapq: int = 20,
) -> HostMappingState:
    return _classify_summary(raw_alignments, min_mapq=min_mapq).state


def adjudicate_host_alignments(
    sample_id: str,
    raw_alignments: Iterable[RawAlignment],
    *,
    min_mapq: int = 20,
    censor_state: CensorState = CensorState.COMPLETE,
) -> HostMappingResult:
    if not sample_id:
        raise InputContractError("sample_id must be non-empty text")
    if not isinstance(censor_state, CensorState):
        raise InputContractError("censor_state must be a CensorState")
    rows = tuple(raw_alignments)
    summary = _classify_summary(rows, min_mapq=min_mapq)
    placements: list[HostPlacement] = []
    placement_rank = 1
    for row in rows:
        if not row.is_mapped or row.is_supplementary:
            continue
        if (
            row.reference_name is None
            or row.position0 is None
            or row.cigar is None
        ):
            raise InputContractError(
                "mapped host alignment requires reference, position, and CIGAR"
            )
        reference_length = cigar_reference_length(row.cigar)
        if reference_length == 0:
            raise InputContractError("mapped host CIGAR consumes no reference bases")
        span0 = BaseSpan0(
            row.position0,
            row.position0 + reference_length,
        )
        placements.append(
            HostPlacement(
                sample_id=sample_id,
                pair_id=row.pair_id,
                mate_no=row.mate,
                placement_rank=placement_rank,
                chrom=row.reference_name,
                span0=span0,
                host_strand="-" if row.is_reverse else "+",
                mapq=row.mapq,
                cigar=row.cigar,
                mapping_state=summary.state,
            )
        )
        placement_rank += 1
    return HostMappingResult(
        sample_id=sample_id,
        pair_id=summary.pair_id,
        mate_no=summary.mate,
        state=summary.state,
        placements=tuple(placements),
        censor_state=censor_state,
    )


def map_host_mate(
    pair: PhysicalPair,
    family_evidence: TEFamilyEvidence,
    mapper: HostMapperBackend,
    prepared_reference: PreparedReference,
    *,
    max_reported: int = 2,
    min_mapq: int = 20,
) -> HostMappingResult:
    """Map only the mate opposite the eligible TE-positive mate."""

    if (
        family_evidence.sample_id != pair.sample_id
        or family_evidence.pair_id != pair.pair_id
        or family_evidence.eligibility_state is not TEEligibilityState.ELIGIBLE
        or family_evidence.te_positive_mate is None
    ):
        raise InputContractError("host mapping requires matching eligible pair evidence")
    host_mate = (
        ReadMate.READ2
        if family_evidence.te_positive_mate is ReadMate.READ1
        else ReadMate.READ1
    )
    read = pair.read2 if host_mate is ReadMate.READ2 else pair.read1
    mapper_read = MapperRead(
        query_name=f"evitea:{pair.pair_id}:{int(host_mate)}",
        sequence=read.sequence,
        quality=read.quality,
        pair_id=pair.pair_id,
        mate=host_mate,
    )
    rows = tuple(
        mapper.map_reads(
            [mapper_read], prepared_reference, max_reported=max_reported
        )
    )
    return adjudicate_host_alignments(
        pair.sample_id, rows, min_mapq=min_mapq
    )
