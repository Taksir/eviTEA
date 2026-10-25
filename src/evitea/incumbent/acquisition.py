"""Candidate-locus acquisition for the incumbent branch (stages G1A..G1H-R).

This composes the anchor census, sparse panel, coverage rescue, pair-aware
reselection, host-mate gate, resident pair families and recall repair into
the candidate physical pairs and loci consumed by the split-junction layer.
Every quantity comes from the prepared sample index, the TE consensus
library, its family map and the host Bowtie2 index; no historical artifact is
read.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..constants import Q99
from .anchor_census import build_anchor_catalogue, census_anchor_intervals
from .candidates import (
    CandidateLocus,
    CandidatePair,
    HostMateRecord,
    PairPanelEvidence,
    cluster_candidate_loci,
    reconstruct_candidate_pairs,
)
from .coverage_panel import READ_WINDOW, hard_gap_h1_loci, optimize_rescue_panel, short_exact_census
from .host_mate_gate import (
    PRIMARY_HOST_OPTIONS,
    HostMapping,
    one_mate_host_requests,
    parse_primary_host_sam,
    resident_pair_families,
    run_host_mapping,
    supported_host_cluster_pairs,
    working_candidate_pairs,
)
from .pair_aware_panel import select_pair_aware_panel
from .recall_repair import (
    RESCUE_HOST_OPTIONS,
    rescue_query_name,
    enumerate_rescue_pairs,
    parse_rescue_host_sam,
    MAX_RESCUE_INTERVAL_MASS,
    plan_rescue_anchors,
    rescue_host_requests,
    supported_rescue_pairs,
)
from .sparse_panel import route_sparse_panel, select_sparse_panel

__all__ = ["CandidateAcquisition", "acquire_candidate_loci"]


@dataclass(frozen=True, slots=True)
class CandidateAcquisition:
    loci: tuple[CandidateLocus, ...]
    candidates: tuple[CandidatePair, ...]
    counts: Mapping[str, int] = field(default_factory=dict)


def _host_record(row: HostMapping) -> HostMateRecord:
    return HostMateRecord(row.pair_id, row.status, row.chrom, "" if row.midpoint is None else str(row.midpoint), row.strand)


def acquire_candidate_loci(
    *,
    te_records: Mapping[str, str],
    te_family_map: Mapping[str, str],
    index,
    bowtie2: str | Path,
    host_index_prefix: str | Path,
    work_directory: Path,
    threads: int = 1,
    progress: Callable[[str], None] | None = None,
    rescue_read_names: str = "pair",
    complete_routing: bool = False,
    max_rescue_interval_mass: int = MAX_RESCUE_INTERVAL_MASS,
    read_length: int = READ_WINDOW,
    fragment_q99: int = Q99,
) -> CandidateAcquisition:
    def say(message: str) -> None:
        if progress is not None:
            progress(message)

    counts: dict[str, int] = {}
    index.verify_paired_layout()
    say("anchor census")
    catalogue = build_anchor_catalogue(te_records, te_family_map)
    census_anchor_intervals(catalogue, index)
    counts["anchor_queries"] = len(catalogue.queries)

    say("sparse panel routing")
    panel = select_sparse_panel(catalogue)
    panel_pairs = route_sparse_panel(panel, index)
    counts["sparse_panel_queries"] = len(panel)
    counts["sparse_panel_pairs"] = len(panel_pairs)

    say("coverage rescue census")
    family_by_te = catalogue.family_by_te()
    short = short_exact_census(catalogue.records, family_by_te, index)
    from .coverage_panel import baseline_anchor_loci

    h1 = hard_gap_h1_loci(catalogue.records, family_by_te, baseline_anchor_loci(catalogue, panel), short, read_length)
    say("rescue panel optimization")
    rescue = optimize_rescue_panel(catalogue, panel, short, h1, index, read_length)
    counts["rescue_panel_sequences"] = len(rescue.selected_sequences)

    say("pair-aware reselection")
    pair_panel = select_pair_aware_panel(catalogue, panel, short, rescue, index, read_length)
    counts["pair_aware_pairs"] = int(len(pair_panel.pair_ids))
    if complete_routing:
        from .complete_routing import augment_panel, collect_extra_pairs, extra_queries, gate_extra_pairs

        say("complete routing: extra words and pairs")
        extra = extra_queries(catalogue, pair_panel.selected_sequences)
        te_file, extra_sets, te_dollars = collect_extra_pairs(extra, pair_panel.pair_ids.tolist(), index, progress=say)
        say(f"complete routing: {len(extra)} words, {len(te_file)} extra one-mate pairs; partner-read gate")
        kept = gate_extra_pairs(
            te_file, index, bowtie2=bowtie2, host_index_prefix=host_index_prefix,
            work_directory=work_directory / "complete_routing", threads=threads, q99=fragment_q99, progress=say,
        )
        pair_panel = augment_panel(pair_panel, kept, extra, extra_sets, te_dollars)
        counts["complete_routing_words"] = len(extra)
        counts["complete_routing_extra_pairs"] = len(te_file)
        counts["complete_routing_kept_pairs"] = len(kept)

    say("host-mate localization")
    requests = one_mate_host_requests(pair_panel, index)
    dollars = index.mate_dollars(requests)
    records = [(f"p{pid}|h{requests[pid] + 1}", index.recover(dollars[pid])) for pid in sorted(requests)]
    sam = run_host_mapping(
        records, bowtie2=bowtie2, host_index_prefix=host_index_prefix, work_directory=work_directory,
        name="host_mates", options=PRIMARY_HOST_OPTIONS, threads=threads,
    )
    primary = parse_primary_host_sam(sam, requests)
    host_pairs = supported_host_cluster_pairs(primary)
    working = working_candidate_pairs(host_pairs, panel_pairs)
    families = resident_pair_families(host_pairs, pair_panel)
    counts["one_mate_pairs"] = len(requests)
    counts["working_candidate_pairs"] = int(len(working))

    say("candidate recall repair")
    anchors = plan_rescue_anchors(
        catalogue, pair_panel.groups, max_interval_mass=max_rescue_interval_mass, read_window=read_length
    )
    evidence = enumerate_rescue_pairs(anchors, pair_panel.pair_ids, index)
    rescue_requests = rescue_host_requests(evidence)
    rescue_dollars = index.mate_dollars(rescue_requests)
    rescue_records = [
        (rescue_query_name(pid, rescue_read_names), index.recover(rescue_dollars[pid]))
        for pid in sorted(rescue_requests)
    ]
    rescue_sam = run_host_mapping(
        rescue_records, bowtie2=bowtie2, host_index_prefix=host_index_prefix, work_directory=work_directory,
        name="rescue_host_mates", options=RESCUE_HOST_OPTIONS, threads=threads,
    )
    rescue_mapping = parse_rescue_host_sam(rescue_sam, rescue_requests)
    repaired = np.union1d(working, np.asarray(sorted(supported_rescue_pairs(rescue_mapping)), dtype=np.uint64))
    counts["rescue_anchor_count"] = len(anchors)
    counts["repaired_candidate_pairs"] = int(len(repaired))

    candidates, meta = reconstruct_candidate_pairs(
        final_pair_ids=repaired.tolist(),
        working_pair_ids=working.tolist(),
        primary_host={pid: _host_record(row) for pid, row in primary.items()},
        rescue_host={pid: _host_record(row) for pid, row in rescue_mapping.items()},
        resident_families=families,
        rescue_evidence={
            pid: PairPanelEvidence(pid, ",".join(row.families), row.mate1_hit, row.mate2_hit)
            for pid, row in evidence.items()
        },
        panel_evidence={
            row.pair_id: PairPanelEvidence(row.pair_id, ",".join(row.families), row.mate1_hit, row.mate2_hit)
            for row in panel_pairs
        },
    )
    loci = cluster_candidate_loci(candidates)
    counts.update({key: int(value) for key, value in meta.items()})
    counts["candidate_loci"] = len(loci)
    return CandidateAcquisition(loci, candidates, counts)
