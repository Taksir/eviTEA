"""Host and evidence-complete reference adjudication."""

from .host import HostMappingResult, adjudicate_host_alignments
from .reference import (
    ReferenceAlignmentObservation,
    ReferenceObservationResult,
    build_reference_alternatives,
    evaluate_reference_explanation,
    produce_reference_observation_state,
    reference_observations_from_alignments,
)

__all__ = [
    "HostMappingResult",
    "ReferenceAlignmentObservation",
    "ReferenceObservationResult",
    "adjudicate_host_alignments",
    "build_reference_alternatives",
    "evaluate_reference_explanation",
    "produce_reference_observation_state",
    "reference_observations_from_alignments",
]
