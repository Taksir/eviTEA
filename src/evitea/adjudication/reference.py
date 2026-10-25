"""Evidence-complete finite reference alternatives."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from ..backends.protocols import RawAlignment
from ..contracts import (
    BaseSpan0,
    CensorState,
    ReadMate,
    ReferenceAlternative,
    ReferenceExplanationCode,
    ReferenceExplanationState,
    TEFamilyEvidence,
)
from ..errors import InputContractError

_CIGAR_TOKEN = re.compile(r"([1-9][0-9]*)([MIDNSHP=X])")
_QUERY_CONSUMING = frozenset("MIS=X")
_REFERENCE_CONSUMING = frozenset("MDN=X")
_EXPLAINING = frozenset("M=X")


def _cigar_tokens(cigar: str) -> tuple[tuple[int, str], ...]:
    if not isinstance(cigar, str) or not cigar or cigar == "*":
        raise InputContractError("CIGAR must be non-empty mapped SAM text")
    tokens: list[tuple[int, str]] = []
    cursor = 0
    for match in _CIGAR_TOKEN.finditer(cigar):
        if match.start() != cursor:
            raise InputContractError(f"invalid CIGAR: {cigar}")
        tokens.append((int(match.group(1)), match.group(2)))
        cursor = match.end()
    if cursor != len(cigar) or not tokens:
        raise InputContractError(f"invalid CIGAR: {cigar}")
    return tuple(tokens)


def aligned_query_spans(cigar: str) -> tuple[BaseSpan0, ...]:
    """Return query bases explained by CIGAR ``M``, ``=``, or ``X``."""

    query0 = 0
    spans: list[BaseSpan0] = []
    for length, operation in _cigar_tokens(cigar):
        if operation in _EXPLAINING:
            span = BaseSpan0(query0, query0 + length)
            if spans and spans[-1].end0 == span.start0:
                spans[-1] = BaseSpan0(spans[-1].start0, span.end0)
            else:
                spans.append(span)
        if operation in _QUERY_CONSUMING:
            query0 += length
    return tuple(spans)


def cigar_reference_length(cigar: str) -> int:
    return sum(
        length for length, operation in _cigar_tokens(cigar)
        if operation in _REFERENCE_CONSUMING
    )


def _subtract_span(
    required: BaseSpan0,
    explained: Iterable[BaseSpan0],
) -> tuple[BaseSpan0, ...]:
    remaining = [required]
    for coverage in sorted(explained):
        next_remaining: list[BaseSpan0] = []
        for span in remaining:
            if coverage.end0 <= span.start0 or coverage.start0 >= span.end0:
                next_remaining.append(span)
                continue
            if coverage.start0 > span.start0:
                next_remaining.append(BaseSpan0(span.start0, coverage.start0))
            if coverage.end0 < span.end0:
                next_remaining.append(BaseSpan0(coverage.end0, span.end0))
        remaining = next_remaining
    return tuple(remaining)


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentObservation:
    alternative_id: str
    placement_rank: int
    chrom: str
    te_alignment_span0: BaseSpan0
    cigar: str
    is_local_proper: bool
    censor_state: CensorState = CensorState.COMPLETE
    signature_id: str | None = None

    def __post_init__(self) -> None:
        if not self.alternative_id or not self.chrom:
            raise InputContractError("reference alternative identity must be non-empty")
        if (
            isinstance(self.placement_rank, bool)
            or not isinstance(self.placement_rank, int)
            or self.placement_rank < 0
        ):
            raise InputContractError("placement_rank must be non-negative")
        if not isinstance(self.is_local_proper, bool):
            raise InputContractError("is_local_proper must be boolean")
        if not isinstance(self.censor_state, CensorState):
            raise InputContractError("censor_state must be a CensorState")
        reference_length = cigar_reference_length(self.cigar)
        if self.te_alignment_span0.length != reference_length:
            raise InputContractError(
                "TE alignment span length disagrees with reference-consuming CIGAR"
            )


@dataclass(frozen=True, slots=True)
class ReferenceObservationResult:
    """Finite reference alternatives and their reduction to one explanation state.

    A censored search gives an indeterminate state.
    """

    observations: tuple[ReferenceAlignmentObservation, ...]
    alternatives: tuple[ReferenceAlternative, ...]
    state: ReferenceExplanationState


def _proper_signature(
    row: RawAlignment,
) -> tuple[str, int, int, int] | None:
    """Return the LOCAL proper-pair signature for one SAM row."""

    if (
        not row.is_mapped
        or not (row.flag & 0x2)
        or bool(row.flag & 0x8)
        or row.is_supplementary
    ):
        return None
    if (
        row.reference_name is None
        or row.position0 is None
        or row.mate_position0 is None
    ):
        raise InputContractError(
            "proper reference alignment lacks pair geometry"
        )
    return (
        row.reference_name,
        min(row.position0, row.mate_position0),
        max(row.position0, row.mate_position0),
        abs(row.template_length),
    )


def reference_observations_from_alignments(
    pair_evidence: TEFamilyEvidence,
    raw_alignments: Iterable[RawAlignment],
    *,
    censor_state: CensorState = CensorState.COMPLETE,
) -> tuple[ReferenceAlignmentObservation, ...]:
    """Produce observations from every complete LOCAL proper signature.

    Secondary rows are retained.  A signature is finite only when both
    physical mates are represented.
    """

    if pair_evidence.te_positive_mate is None:
        raise InputContractError(
            "building reference observations requires a TE-positive mate"
        )
    if not isinstance(censor_state, CensorState):
        raise InputContractError("censor_state must be a CensorState")
    grouped: dict[
        tuple[str, int, int, int],
        dict[ReadMate, list[tuple[int, RawAlignment]]],
    ] = defaultdict(lambda: {ReadMate.READ1: [], ReadMate.READ2: []})
    for ordinal0, row in enumerate(raw_alignments):
        if row.pair_id != pair_evidence.pair_id:
            raise InputContractError("reference alignments mix physical pairs")
        signature = _proper_signature(row)
        if signature is None:
            continue
        grouped[signature][row.mate].append((ordinal0, row))

    observations: list[ReferenceAlignmentObservation] = []
    placement_rank = 1
    for signature_ordinal, signature in enumerate(sorted(grouped), start=1):
        by_mate = grouped[signature]
        if not by_mate[ReadMate.READ1] or not by_mate[ReadMate.READ2]:
            continue
        chrom, left0, right0, template_length = signature
        signature_id = (
            f"{chrom}:{left0 + 1}:{right0 + 1}:{template_length}"
        )
        te_rows = by_mate[pair_evidence.te_positive_mate]
        for row_ordinal, (_, row) in enumerate(te_rows, start=1):
            if row.cigar is None or row.position0 is None:
                raise InputContractError(
                    "mapped TE-mate alignment lacks position or CIGAR"
                )
            reference_length = cigar_reference_length(row.cigar)
            if reference_length == 0:
                raise InputContractError(
                    "mapped TE-mate CIGAR consumes no reference bases"
                )
            observations.append(
                ReferenceAlignmentObservation(
                    alternative_id=(
                        f"pair-{pair_evidence.pair_id}:signature-"
                        f"{signature_ordinal}:te-row-{row_ordinal}"
                    ),
                    placement_rank=placement_rank,
                    chrom=chrom,
                    te_alignment_span0=BaseSpan0(
                        row.position0, row.position0 + reference_length
                    ),
                    cigar=row.cigar,
                    is_local_proper=True,
                    censor_state=censor_state,
                    signature_id=signature_id,
                )
            )
            placement_rank += 1
    return tuple(observations)


def produce_reference_observation_state(
    pair_evidence: TEFamilyEvidence,
    relevant_query_spans0: Iterable[BaseSpan0],
    raw_alignments: Iterable[RawAlignment],
    *,
    search_censor_state: CensorState = CensorState.COMPLETE,
) -> ReferenceObservationResult:
    """Reduce alignments to observations, alternatives and one reference explanation state."""

    observations = reference_observations_from_alignments(
        pair_evidence,
        raw_alignments,
        censor_state=search_censor_state,
    )
    alternatives = build_reference_alternatives(
        pair_evidence,
        relevant_query_spans0,
        observations,
    )
    state = evaluate_reference_explanation(
        pair_evidence,
        alternatives,
        search_censor_state=search_censor_state,
    )
    return ReferenceObservationResult(observations, alternatives, state)


def build_reference_alternatives(
    pair_evidence: TEFamilyEvidence,
    relevant_query_spans0: Iterable[BaseSpan0],
    observations: Iterable[ReferenceAlignmentObservation],
) -> tuple[ReferenceAlternative, ...]:
    """Evaluate every provided finite placement; never select the first or nearest."""

    if pair_evidence.resolved_family is None:
        raise InputContractError("reference competition requires resolved family evidence")
    if pair_evidence.te_positive_mate is None:
        raise InputContractError(
            "reference competition requires an identified TE-positive mate"
        )
    required = tuple(sorted(set(relevant_query_spans0)))
    if not required:
        raise InputContractError("reference competition requires safe-probe occurrences")
    alternatives: list[ReferenceAlternative] = []
    seen_ids: set[str] = set()
    for observation in observations:
        if observation.alternative_id in seen_ids:
            raise InputContractError("duplicate reference alternative_id")
        seen_ids.add(observation.alternative_id)
        explained = aligned_query_spans(observation.cigar)
        unexplained = tuple(
            fragment
            for span in required
            for fragment in _subtract_span(span, explained)
        )
        alternatives.append(
            ReferenceAlternative(
                sample_id=pair_evidence.sample_id,
                pair_id=pair_evidence.pair_id,
                alternative_id=observation.alternative_id,
                placement_rank=observation.placement_rank,
                family=pair_evidence.resolved_family,
                chrom=observation.chrom,
                te_alignment_span0=observation.te_alignment_span0,
                te_cigar=observation.cigar,
                is_local_proper=observation.is_local_proper,
                evidence_complete=observation.is_local_proper and not unexplained,
                unexplained_query_spans0=unexplained,
                censor_state=observation.censor_state,
                signature_id=observation.signature_id,
                mate_no=pair_evidence.te_positive_mate,
            )
        )
    return tuple(sorted(alternatives, key=lambda item: (item.placement_rank, item.alternative_id)))


def combine_censor_states(states: Iterable[CensorState]) -> CensorState:
    values = {state for state in states if state is not CensorState.COMPLETE}
    if not values:
        return CensorState.COMPLETE
    if CensorState.HOST_AND_REFERENCE_RESOURCE_CENSORED in values or {
        CensorState.HOST_RESOURCE_CENSORED,
        CensorState.REFERENCE_RESOURCE_CENSORED,
    }.issubset(values):
        return CensorState.HOST_AND_REFERENCE_RESOURCE_CENSORED
    if CensorState.SOURCE_INCOMPLETE in values:
        return CensorState.SOURCE_INCOMPLETE
    if CensorState.MODEL_INCOMPLETE in values:
        return CensorState.MODEL_INCOMPLETE
    return sorted(values, key=lambda state: state.value)[0]


def evaluate_reference_explanation(
    pair_evidence: TEFamilyEvidence,
    alternatives: Iterable[ReferenceAlternative],
    *,
    search_censor_state: CensorState = CensorState.COMPLETE,
) -> ReferenceExplanationState:
    """Reduce all finite alternatives while preserving ambiguity and censor."""

    rows = tuple(
        sorted(alternatives, key=lambda item: (item.placement_rank, item.alternative_id))
    )
    if len({row.alternative_id for row in rows}) != len(rows):
        raise InputContractError("duplicate reference alternative_id")
    family = pair_evidence.resolved_family
    for alternative in rows:
        if (
            alternative.sample_id != pair_evidence.sample_id
            or alternative.pair_id != pair_evidence.pair_id
            or family is None
            or alternative.family != family
        ):
            raise InputContractError("reference alternatives mix pair/family identity")
    censor = combine_censor_states(
        [search_censor_state, *(row.censor_state for row in rows)]
    )
    finite_ids = tuple(row.alternative_id for row in rows)
    complete_ids = tuple(
        row.alternative_id for row in rows if row.evidence_complete
    )
    if censor is not CensorState.COMPLETE:
        code = ReferenceExplanationCode.INDETERMINATE_CENSORED
    elif complete_ids:
        code = ReferenceExplanationCode.EVIDENCE_COMPLETE
    elif not any(row.is_local_proper for row in rows):
        code = ReferenceExplanationCode.NO_LOCAL_PROPER
    else:
        code = ReferenceExplanationCode.REFERENCE_UNEXPLAINED
    return ReferenceExplanationState(
        sample_id=pair_evidence.sample_id,
        pair_id=pair_evidence.pair_id,
        code=code,
        censor_state=censor,
        finite_alternative_ids=finite_ids,
        evidence_complete_alternative_ids=complete_ids,
    )
