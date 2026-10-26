"""Full-memory layer engine: finite state, proposals, exact junction evidence, events."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from ..anchor_pe.alignments import (
    PRODUCTION_HOST_READ_NAMES,
    map_host_mates,
    map_reference_pairs,
    write_host_fasta,
    write_pair_fastas,
)
from ..anchor_pe.codebook import load_safe_codebook
from ..anchor_pe.engine import AnchorPEInputs, AnchorPEResult, _fingerprint
from ..incumbent.engine import StageStore
from ..shadow_layers import FullMemoryEvent
from .finite_state import family_probes, finite_host_alternatives, finite_reference_alternatives, pair_accounting
from .proposals import full_data_proposals, load_proposal_records, rotation_proposals
from .selection import full_memory_events, select_anchors
from .strict_a import a_groups, allowed_pairs, build_b_index, derive_witnesses, full_data_events, historical_keys

__all__ = ["FullMemoryResult", "run_full_memory"]

# The finite-state replay names every host query P{pid}, whatever the
# paired-end anchor host read-name profile.
# The AKA017 profile therefore maps host mates twice.
FINITE_STATE_HOST_PROFILE = PRODUCTION_HOST_READ_NAMES


@dataclass(slots=True)
class FullMemoryResult:
    events: tuple[FullMemoryEvent, ...]
    accounting: list[dict]
    host_alternatives: list[dict]
    full_proposals: list[dict]
    rotation_proposals: list[dict]
    event_ledger: list[dict]
    witnesses: list[dict] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)


def run_full_memory(
    inputs: AnchorPEInputs,
    anchor_pe: AnchorPEResult,
    work_directory: Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> FullMemoryResult:
    def say(message: str) -> None:
        if progress is not None:
            progress(message)

    work = Path(work_directory)
    work.mkdir(parents=True, exist_ok=True)
    store = StageStore(work / "checkpoints", {**_fingerprint(inputs), "schema": "evitea.full_memory.checkpoint.v1"})
    codebook = load_safe_codebook(inputs.safe_codebook)
    eligible, reads = anchor_pe.eligible, anchor_pe.reads

    host_sam = anchor_pe.host_sams.get(FINITE_STATE_HOST_PROFILE)
    if host_sam is None or not Path(host_sam).is_file():
        host_sam = work / f"host_mates_{FINITE_STATE_HOST_PROFILE}.sam"
        if not host_sam.is_file():
            fasta = work / f"host_mates_{FINITE_STATE_HOST_PROFILE}.fa"
            write_host_fasta(fasta, {p: m["host_mate"] for p, m in eligible.items()}, reads, FINITE_STATE_HOST_PROFILE)
            map_host_mates(inputs.bowtie2, inputs.host_index_prefix, fasta, host_sam, inputs.threads)
    reference_sam = anchor_pe.reference_sam
    if reference_sam is None or not Path(reference_sam).is_file():
        reference_sam = work / "local_reference.sam"
        if not reference_sam.is_file():
            r1, r2 = work / "pairs_R1.fa", work / "pairs_R2.fa"
            write_pair_fastas(r1, r2, eligible, reads)
            map_reference_pairs(inputs.bowtie2, inputs.host_index_prefix, r1, r2, reference_sam, inputs.threads,
                                inputs.fragment_q99)

    def finite_state():
        host_rows, host_summary = finite_host_alternatives(Path(host_sam), eligible)
        meta = {p: {"pair_id": p, "family": m["family"], "te_mate": m["te_mate"], "host_mate": m["host_mate"]}
                for p, m in eligible.items()}
        ref_rows, ref_summary = finite_reference_alternatives(Path(reference_sam), meta, family_probes(codebook), reads)
        accounting = pair_accounting(anchor_pe.signature_rows, host_summary, ref_summary)
        # Only uncensored, not evidence-complete alternatives can witness exact junction evidence.
        ref_rows = [r for r in ref_rows if not r["evidence_complete"] and not r["reference_cap_censored"]]
        return accounting, host_rows, ref_rows

    accounting, host_rows, ref_rows = store.get("finite_state", finite_state)
    say(f"finite state: {len(accounting):,} accounted pairs, {len(host_rows):,} host alternatives")
    records = load_proposal_records(accounting, host_rows, inputs.fragment_q99)
    primary = frozenset(inputs.primary_contigs)
    full = full_data_proposals(inputs.sample_id, records, primary)
    rotations = rotation_proposals(inputs.sample_id, records, primary)
    witnesses = derive_witnesses(allowed_pairs(accounting), ref_rows)
    index, _ = build_b_index(anchor_pe.b_universe, inputs.fragment_q99)
    ledger = full_data_events(inputs.sample_id, full, a_groups(witnesses), index, historical_keys(anchor_pe.witnesses))
    events = full_memory_events(select_anchors(ledger, inputs.sample_id), index)
    counts = {
        "proposal_records": len(records), "full_proposals": len(full), "rotation_proposals": len(rotations),
        "strict_a_witnesses": len(witnesses), "full_data_events": len(ledger),
        "proposal_only_events": sum(e["proposal_only_A"] for e in ledger), "selected_anchors": len(events),
        "certified_events": sum(e.b_state != "NEITHER" for e in events),
    }
    say(f"full-memory counts: {counts}")
    return FullMemoryResult(events, accounting, host_rows, full, rotations, ledger, witnesses, counts)
