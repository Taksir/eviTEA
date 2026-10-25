"""Extension engine: library-layer callset bytes to the four extension layers.

Runs the evidence producers in dependency order on one prepared sample and
composes the four extension layers (``shadow_layers``)::

    anchor_pe (candidates)  -> bilateral, asymmetric
    full_memory (events)    -> full-memory
    internal21 (loci)       -> internal 21-mer
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .anchor_pe.engine import AnchorPEInputs, AnchorPEResult, run_anchor_pe
from .full_memory.engine import FullMemoryResult, run_full_memory
from .full_memory.strict_a import build_b_index
from .internal21.acquisition import Internal21Acquisition, run_internal21_acquisition
from .internal21.observability import IncrementalLocus, novel_unconfounded_loci, observe_increment
from .shadow_layers import (
    ShadowCallset,
    ShadowLabels,
    ShadowLayerResult,
    append_asymmetric_layer,
    append_bilateral_layer,
    append_full_memory_layer,
    append_internal21_layer,
)

__all__ = ["ExtensionRun", "run_extension"]


@dataclass(slots=True)
class ExtensionRun:
    anchor_pe: AnchorPEResult
    full_memory: FullMemoryResult
    internal21: Internal21Acquisition
    internal21_events: list[dict]
    internal21_loci: list[IncrementalLocus]
    layers: tuple[ShadowLayerResult, ...]


def run_extension(
    inputs: AnchorPEInputs,
    internal_codebook: Path,
    sample,
    incumbent: bytes,
    labels: ShadowLabels,
    work_directory: Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> ExtensionRun:
    work = Path(work_directory)
    anchor_pe = run_anchor_pe(inputs, sample, work, progress=progress)
    full_memory = run_full_memory(inputs, anchor_pe, work / "full_memory", progress=progress)
    internal = run_internal21_acquisition(inputs, internal_codebook, sample, work / "internal21", progress=progress)
    index, _ = build_b_index(anchor_pe.b_universe, inputs.fragment_q99)
    events, loci = observe_increment(inputs.sample_id, full_memory.rotation_proposals, full_memory.witnesses,
                                     internal.witnesses, index, full_memory.event_ledger)
    base = ShadowCallset.from_bytes(incumbent)
    bilateral = append_bilateral_layer(base, anchor_pe.candidates, labels)
    asymmetric = append_asymmetric_layer(bilateral.callset, anchor_pe.candidates, labels)
    full = append_full_memory_layer(asymmetric.callset, full_memory.events, labels, sample_id=inputs.sample_id)
    internal21 = append_internal21_layer(full.callset, novel_unconfounded_loci(loci), labels)
    return ExtensionRun(anchor_pe, full_memory, internal, events, loci, (bilateral, asymmetric, full, internal21))
