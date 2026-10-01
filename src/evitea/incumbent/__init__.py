"""Frozen incumbent (repeat-aware R9E) caller branch.

``run_incumbent_caller`` derives the incumbent callset from the prepared
sample msBWT, the TE consensus library and family map, the host reference and
its Bowtie2 index.  The typed repeat-aware core (``run_repeat_aware_incumbent``)
remains available for observation-level use.
"""

from .acquisition import CandidateAcquisition, acquire_candidate_loci
from .candidates import CandidateLocus, CandidatePair
from .engine import IncumbentInputs, IncumbentResult, run_incumbent_caller
from .events import IncumbentCall, IncumbentCallClass, legacy_tsv_bytes
from .msbwt_queries import MsbwtSample
from .repeat_aware import (
    RepeatAwareGroupDecision,
    RepeatAwareGroupObservation,
    RepeatAwareLocusDecision,
    RepeatAwareLocusObservation,
    RepeatAwareObservations,
    RepeatAwarePreparedAssets,
    RepeatAwareResult,
    RepeatAwareSampleContext,
    RepeatAwareVariant,
    alignment_variant_ops,
    equivalent_indel_positions,
    indel_haplotype_admissibility,
    run_repeat_aware_incumbent,
)

__all__ = [
    "CandidateAcquisition",
    "CandidateLocus",
    "CandidatePair",
    "IncumbentCall",
    "IncumbentCallClass",
    "IncumbentInputs",
    "IncumbentResult",
    "MsbwtSample",
    "RepeatAwareGroupDecision",
    "RepeatAwareGroupObservation",
    "RepeatAwareLocusDecision",
    "RepeatAwareLocusObservation",
    "RepeatAwareObservations",
    "RepeatAwarePreparedAssets",
    "RepeatAwareResult",
    "RepeatAwareSampleContext",
    "RepeatAwareVariant",
    "acquire_candidate_loci",
    "alignment_variant_ops",
    "equivalent_indel_positions",
    "indel_haplotype_admissibility",
    "legacy_tsv_bytes",
    "run_incumbent_caller",
    "run_repeat_aware_incumbent",
]
