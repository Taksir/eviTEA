"""Library layer caller, run end to end over prepared inputs.

Steps:

1. candidate acquisition (anchor census to the wide-anchor step) giving
   loci;
2. split-junction reference competition, giving calls with
   score >= 4;
3. orientation-corrected two-sided mate-anchored search, giving
   calls and the locus record;
4. one-sided targets, exact host-seed routing and reference-discriminative
   strict pair confirmation;
5. insertion-specific junction discovery, exact-group location, group
   competition and locus linking;
6. full TE-segment family confirmation of linked groups;
7. repeat-aware divergence step (retains a group whose edits recur in other
   groups; leave-one-group-out) and 2-LEFT + 2-RIGHT reaggregation.

The call set lists split-junction calls, then mate-anchored calls, then
repeat-aware calls.
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from ..constants import Q99
from ..errors import ResumeStateError
from .acquisition import CandidateAcquisition, acquire_candidate_loci
from .anchor_census import load_te_family_map, parse_te_consensus
from .events import IncumbentCall, legacy_tsv_bytes
from .candidates import CANDIDATE_MIN_PAIRS
from .external import run_bowtie2_build
from .host_mate_gate import HOST_CLUSTER_SUPPORT
from .junction_groups import JunctionGroupsResult, run_junction_groups
from .mate_rescue import (
    FRAGMENT_Q95,
    MateRescueResult,
    RoutingResult,
    StrictPairResult,
    certify_strict_pairs,
    finalize_mate_rescue,
    guard_candidate_loci,
    one_sided_targets,
    route_one_sided_targets,
)
from .parallel import map_locus_batches
from .recall_repair import MAX_RESCUE_INTERVAL_MASS, READ_WINDOW
from .repeat_rescue import crossvalidate, reaggregate, rescue_calls
from .sequence import fasta_records_full_header, load_reference_contigs
from .split_junction import SplitJunctionResult, finalize_split_junctions
from .te_certification import align_segments, certify_groups, full_te_segment, used_groups
from .te_library import (
    LocalTEIndex,
    load_family_kmer_panels,
    load_heuristic_family_map,
    load_record_family_map,
    load_te_records,
    normalize_family,
)

__all__ = ["IncumbentInputs", "IncumbentResult", "run_incumbent_caller"]

FAMILY_KMER = 21
FAMILY_KMER_MIN_ENTROPY = 1.10


@dataclass(frozen=True, slots=True)
class IncumbentInputs:
    sample_id: str
    te_fasta: Path
    te_family_map: Path
    reference_fasta: Path
    host_index_prefix: Path
    bowtie2: Path
    bowtie2_build: Path | None = None
    te_index_prefix: Path | None = None
    threads: int = 1
    workers: int = 1
    # Query names of the wide-anchor host mapping (see
    # ``recall_repair.RESCUE_READ_NAME_PROFILES``); the default is ``pair``.
    rescue_read_names: str = "pair"
    # Optional complete routing with a partner-read filter (off = default panel only).
    complete_routing: bool = False
    # Work limit of the wide-anchor step (default ``recall_repair.MAX_RESCUE_INTERVAL_MASS``).
    max_rescue_interval_mass: int = 1_500_000
    # Read length of the sample (bp): the window of the coverage panels and of the wide-anchor step.
    read_length: int = 125
    # 95th and 99th percentiles of the sample's fragment size (bp); see ``sample_geometry``.
    fragment_q95: float = FRAGMENT_Q95
    fragment_q99: int = Q99
    # Smallest numbers of the candidate step (``host_mate_gate.HOST_CLUSTER_SUPPORT``, ``candidates.CANDIDATE_MIN_PAIRS``):
    # host-mapped mates in a cluster, and pairs in a candidate place. The defaults are the research rule.
    host_cluster_support: int = HOST_CLUSTER_SUPPORT
    candidate_min_pairs: int = CANDIDATE_MIN_PAIRS
    # Side of a host mate in the mate-anchored search (``mate_rescue.MATE_SIDE_RULES``): ``midpoint`` (the research rule) or
    # ``strand`` (split point from the strands of the place's host anchors).
    mate_side_rule: str = "midpoint"
    # Side test of the missing-side searches (strict pairs, junction groups; ``mate_rescue.JUNCTION_SIDE_RULES``): ``midpoint`` or
    # ``split`` (the per-locus point of the mate-anchored search's host anchors).
    junction_side_rule: str = "midpoint"

    def __post_init__(self) -> None:
        from .recall_repair import RESCUE_READ_NAME_PROFILES
        from .mate_rescue import JUNCTION_SIDE_RULES, MATE_SIDE_RULES

        if self.mate_side_rule not in MATE_SIDE_RULES:
            raise ValueError(f"mate_side_rule must be one of {MATE_SIDE_RULES}")
        if self.junction_side_rule not in JUNCTION_SIDE_RULES:
            raise ValueError(f"junction_side_rule must be one of {JUNCTION_SIDE_RULES}")

        if self.rescue_read_names not in RESCUE_READ_NAME_PROFILES:
            raise ValueError(f"rescue_read_names must be one of {RESCUE_READ_NAME_PROFILES}")
        value = self.max_rescue_interval_mass
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError("max_rescue_interval_mass must be a positive integer")
        for name in ("host_cluster_support", "candidate_min_pairs"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")


@dataclass(slots=True)
class IncumbentResult:
    calls: tuple[IncumbentCall, ...]
    acquisition: CandidateAcquisition
    split_junction: SplitJunctionResult
    mate_rescue: MateRescueResult
    routing: RoutingResult
    strict_pairs: StrictPairResult
    junction_groups: JunctionGroupsResult
    rescue_calls: tuple[IncumbentCall, ...]
    counts: dict[str, int] = field(default_factory=dict)

    def legacy_tsv(self) -> bytes:
        return legacy_tsv_bytes(self.calls)


def _fingerprint(inputs: IncumbentInputs) -> dict[str, object]:
    """Return the identity of a checkpoint's inputs."""

    def stat(path: Path | None) -> object:
        if path is None:
            return None
        target = Path(path)
        return [str(target), target.stat().st_size, int(target.stat().st_mtime)] if target.is_file() else str(target)

    identity = {
        "schema": "evitea.incumbent.checkpoint.v1",
        "sample_id": inputs.sample_id,
        "te_fasta": stat(inputs.te_fasta),
        "te_family_map": stat(inputs.te_family_map),
        "reference_fasta": stat(inputs.reference_fasta),
        "host_index_prefix": str(inputs.host_index_prefix),
    }
    # Recorded only when not the default, so default runs keep the same
    # fingerprint and stay resumable.
    if inputs.rescue_read_names != "pair":
        identity["rescue_read_names"] = inputs.rescue_read_names
    if inputs.complete_routing:
        identity["complete_routing"] = True
    if inputs.max_rescue_interval_mass != MAX_RESCUE_INTERVAL_MASS:
        identity["max_rescue_interval_mass"] = inputs.max_rescue_interval_mass
    if inputs.read_length != READ_WINDOW:
        identity["read_length"] = inputs.read_length
    if inputs.fragment_q95 != FRAGMENT_Q95:
        identity["fragment_q95"] = inputs.fragment_q95
    if inputs.fragment_q99 != Q99:
        identity["fragment_q99"] = inputs.fragment_q99
    if inputs.host_cluster_support != HOST_CLUSTER_SUPPORT:
        identity["host_cluster_support"] = inputs.host_cluster_support
    if inputs.candidate_min_pairs != CANDIDATE_MIN_PAIRS:
        identity["candidate_min_pairs"] = inputs.candidate_min_pairs
    return identity


class StageStore:
    """Pickle checkpoints of completed steps, marker written last.

    A checkpoint is reused only when its input fingerprint equals the
    current one. Files live in the run's private work directory.
    """

    def __init__(self, root: Path, fingerprint: dict[str, object], later: dict[str, object] | None = None,
                 later_stages: Sequence[str] = ()):
        self.root = Path(root)
        self.fingerprint = fingerprint
        # extra identity of the stages that read an option which the earlier stages do not (empty by default: every stage then has
        # the one fingerprint); the earlier stages' checkpoints stay reusable when only that option changes
        self.later = dict(later or {})
        self.later_stages = frozenset(later_stages)
        self.root.mkdir(parents=True, exist_ok=True)

    def fingerprint_of(self, stage: str) -> dict[str, object]:
        return {**self.fingerprint, **self.later} if self.later and stage in self.later_stages else self.fingerprint

    def get(self, stage: str, compute):
        data = self.root / f"{stage}.pkl"
        marker = self.root / f"{stage}.complete.json"
        fingerprint = self.fingerprint_of(stage)
        if data.is_file() and marker.is_file():
            recorded = json.loads(marker.read_text(encoding="utf-8"))
            if recorded.get("fingerprint") != fingerprint:
                raise ResumeStateError(f"checkpoint {stage} was computed from different inputs")
            if hashlib.sha256(data.read_bytes()).hexdigest() != recorded.get("sha256"):
                raise ResumeStateError(f"checkpoint {stage} is corrupt")
            with data.open("rb") as handle:
                return pickle.load(handle)
        value = compute()
        payload = pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)
        temporary = data.with_suffix(".tmp")
        temporary.write_bytes(payload)
        os.replace(temporary, data)
        marker.write_text(
            json.dumps({"fingerprint": fingerprint, "sha256": hashlib.sha256(payload).hexdigest()}, indent=2),
            encoding="utf-8",
        )
        return value


def _te_index_prefix(inputs: IncumbentInputs, work: Path) -> Path:
    if inputs.te_index_prefix is not None:
        return inputs.te_index_prefix
    if inputs.bowtie2_build is None:
        raise ValueError("a TE-library Bowtie2 index prefix or a bowtie2-build executable is required")
    prefix = work / "te_library_index" / "te_library"
    if not Path(f"{prefix}.1.bt2").is_file():
        run_bowtie2_build(inputs.bowtie2_build, inputs.te_fasta, prefix)
    return prefix


SIDE_TEST_COUNTS = (
    # (name in the manifest counts, count, where): the pairs each side test rejects
    ("side_test_rejects_main_loop", "mate_rescue", "PAIR_GEOMETRY_INCONSISTENT"),
    ("side_test_rejects_local_geometry", "strict_pairs", "PAIR_GEOMETRY_INCONSISTENT_LOCAL"),
    ("side_test_rejects_local_not_missing_side", "strict_pairs", "LOCAL_HOST_NOT_MISSING_SIDE"),
    ("side_test_rejects_global_geometry", "strict_pairs", "REJECT_GLOBAL_HOST_GEOMETRY_INCONSISTENT"),
    ("side_test_rejects_global_not_missing_side", "strict_pairs", "REJECT_GLOBAL_HOST_NOT_MISSING_SIDE"),
    ("side_test_rejects_junction_unresolved", "junction_groups", "REJECT_JUNCTION_SIDE_UNRESOLVED"),
    ("side_test_rejects_junction_not_missing_side", "junction_groups", "REJECT_JUNCTION_NOT_MISSING_SIDE"),
)


def side_test_counts(mates, strict, groups) -> dict[str, int]:
    """Pairs rejected by the side tests of the mate-anchored search (main loop), the two missing-side tests of the strict-pair
    step and the side test of the junction groups; the same under either rule, so two runs can be compared."""

    sources = {"mate_rescue": mates.counts, "strict_pairs": strict.counts, "junction_groups": groups.counts}
    return {name: int(sources[where].get(key, 0)) for name, where, key in SIDE_TEST_COUNTS}


def _mate_side_counts(mates, strict, groups) -> dict[str, int]:
    out = {"mate_side_" + key.lower(): int(value) for key, value in sorted(mates.counts.items()) if key.startswith("STRAND_RULE_")}
    out.update(side_test_counts(mates, strict, groups))
    return out


def run_incumbent_caller(
    inputs: IncumbentInputs,
    index,
    work_directory: Path,
    *,
    progress: Callable[[str], None] | None = None,
) -> IncumbentResult:
    def say(message: str) -> None:
        if progress is not None:
            progress(message)

    work = Path(work_directory)
    work.mkdir(parents=True, exist_ok=True)
    # the mate-side rule is read from the mate-anchored search on; acquisition and split-junction do not depend on it
    store = StageStore(
        work / "checkpoints", _fingerprint(inputs),
        {**({"mate_side_rule": inputs.mate_side_rule} if inputs.mate_side_rule != "midpoint" else {}),
         **({"junction_side_rule": inputs.junction_side_rule} if inputs.junction_side_rule != "midpoint" else {})} or None,
        ("mate_rescue", "host_seed_routing", "strict_pairs", "junction_groups"),
    )
    te_records = parse_te_consensus(inputs.te_fasta)
    acquisition = store.get("candidate_acquisition", lambda: acquire_candidate_loci(
        te_records=te_records,
        te_family_map=load_te_family_map(inputs.te_family_map),
        index=index,
        bowtie2=inputs.bowtie2,
        host_index_prefix=inputs.host_index_prefix,
        work_directory=work / "acquisition",
        threads=inputs.threads,
        progress=progress,
        rescue_read_names=inputs.rescue_read_names,
        complete_routing=inputs.complete_routing,
        max_rescue_interval_mass=inputs.max_rescue_interval_mass,
        read_length=inputs.read_length,
        fragment_q99=inputs.fragment_q99,
        host_cluster_support=inputs.host_cluster_support,
        candidate_min_pairs=inputs.candidate_min_pairs,
    ))
    loci = acquisition.loci
    families = {locus.family for locus in loci}
    family_kmers, _ = load_family_kmer_panels(
        inputs.te_fasta, inputs.te_family_map, families, FAMILY_KMER, FAMILY_KMER_MIN_ENTROPY
    )
    te_index = LocalTEIndex(load_te_records(inputs.te_fasta, load_record_family_map(inputs.te_family_map)))
    reference = load_reference_contigs(inputs.reference_fasta, {locus.chrom for locus in loci})

    def batch_progress(label: str):
        return None if progress is None else (lambda done, total: progress(f"{label} batch {done}/{total}"))

    reads: dict[int, tuple[str, str]] | None = None

    def locus_reads() -> dict[int, tuple[str, str]]:
        nonlocal reads
        if reads is None:
            reads = index.pair_reads(pid for locus in loci for pid in locus.pair_ids)
        return reads

    state = lambda: {  # noqa: E731
        "reads": locus_reads(), "reference": reference, "family_kmers": family_kmers, "te_index": te_index,
    }

    def mate_state() -> dict:
        # the option is passed only when it is not the default, so a default run calls the search as before
        return {**state(), **({"side_rule": inputs.mate_side_rule} if inputs.mate_side_rule != "midpoint" else {}),
                **({"record_split": True} if inputs.junction_side_rule != "midpoint" else {})}

    def run_split() -> SplitJunctionResult:
        say("split-junction reference competition")
        rows, counts = map_locus_batches(
            "split_junction", loci, state(), workers=inputs.workers, progress=batch_progress("split-junction"),
        )
        return finalize_split_junctions(inputs.sample_id, rows, counts)

    split = store.get("split_junction", run_split)

    def run_mates() -> MateRescueResult:
        say("two-sided mate-anchored search")
        eligible, guarded = guard_candidate_loci(loci, split.calls)
        rows, counts = map_locus_batches(
            "mate_rescue", eligible, mate_state(), workers=inputs.workers, progress=batch_progress("mate-rescue"),
        )
        return finalize_mate_rescue(inputs.sample_id, eligible, guarded, rows, counts)

    mates = store.get("mate_rescue", run_mates)
    reads = None

    say("one-sided host-seed routing")
    targets = one_sided_targets(mates.ledger)
    reserved = {pid for call in (*split.calls, *mates.calls) for pid in call.fragment_pair_ids}
    routing = store.get("host_seed_routing", lambda: route_one_sided_targets(
        targets=targets, reserved_pair_ids=reserved, reference=reference, locator=index,
    ))
    routed_reads = index.pair_reads(routing.pair_to_loci)
    # the points of --junction-side-rule split (nothing at the default, so the later steps are called as before)
    side_args: dict = {}
    if inputs.junction_side_rule != "midpoint":
        side_args = {"side_points": {key: row.junction_point for key, row in getattr(mates, "splits", {}).items()}}

    def run_strict() -> StrictPairResult:
        say("reference-discriminative strict pairs")
        return certify_strict_pairs(
            routing=routing, reads=routed_reads, reference=reference, family_kmers=family_kmers,
            te_index=te_index, bowtie2=inputs.bowtie2, host_index_prefix=inputs.host_index_prefix,
            work_directory=work / "strict_pairs", fragment_q95=inputs.fragment_q95, **side_args,
        )

    strict = store.get("strict_pairs", run_strict)

    def run_groups() -> JunctionGroupsResult:
        say("insertion-specific junction groups")
        return run_junction_groups(
            targets=targets, routing=routing, strict_winners=strict.winners, reserved_call_pair_ids=reserved,
            reads=routed_reads, reference=reference, family_kmers=family_kmers, locator=index,
            bowtie2=inputs.bowtie2, host_index_prefix=inputs.host_index_prefix,
            work_directory=work / "junction_groups", pair_reader=index.pair_reads, **side_args,
        )

    groups = store.get("junction_groups", run_groups)

    say("full TE-segment family confirmation")
    used = used_groups(groups.ledger, groups.winners, groups.accepted_proposals)
    source_reads = {pid: routed_reads[pid] for pid in {g.source_pair_id for g in used} if pid in routed_reads}
    missing = {g.source_pair_id for g in used} - set(source_reads)
    if missing:
        source_reads.update(index.pair_reads(missing))
    segments = {
        g.group_hash: full_te_segment(source_reads[g.source_pair_id][g.source_mate - 1], g.read_view,
                                      g.junction_semantic, g.split)
        for g in used
    }
    alignments = align_segments(
        segments, bowtie2=inputs.bowtie2, te_index_prefix=_te_index_prefix(inputs, work),
        work_directory=work / "te_segments",
    )
    record_families = {
        key: normalize_family(value) for key, value in load_heuristic_family_map(inputs.te_family_map).items()
    }
    certifications = certify_groups(used, segments, alignments, record_families)

    say("repeat-aware step")
    te_sequences = {header.split()[0]: sequence for header, sequence in fasta_records_full_header(inputs.te_fasta)}
    decisions = crossvalidate(certifications, te_sequences, record_families)
    group_pass = {
        row.group.group_hash: row.passes or decisions[row.group.group_hash].crossvalidated_pass
        for row in certifications
    }
    rescued = rescue_calls(inputs.sample_id, groups.ledger, reaggregate(groups.ledger, group_pass))
    calls = (*split.calls, *mates.calls, *rescued)
    counts = dict(acquisition.counts)
    counts.update({
        "split_junction_calls": len(split.calls),
        "mate_rescue_calls": len(mates.calls),
        "one_sided_targets": len(targets),
        "strict_pair_winners": len(strict.winners),
        "linked_junction_groups": len(used),
        "repeat_aware_rescue_calls": len(rescued),
        "incumbent_calls": len(calls),
    })
    if inputs.mate_side_rule != "midpoint" or inputs.junction_side_rule != "midpoint":    # only for an opt-in rule: default counts keep their bytes
        counts.update(_mate_side_counts(mates, strict, groups))
    return IncumbentResult(calls, acquisition, split, mates, routing, strict, groups, rescued, counts)
