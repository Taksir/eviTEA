"""Fail-closed orchestration, marker-last checkpoints, and finalization."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import uuid
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .assets.reference import PreparedReference, verify_reference_asset
from .assets.te import PreparedTEAssets, verify_te_assets
from .backends.exact_query import (
    PreparedSample,
    load_prepared_sample,
    verify_prepared_sample,
)
from .backends.msbwt import build_prepared_msbwt_sample
from .config import EviTEAConfig
from .sample_geometry import resolve_read_length
from .constants import ARCHITECTURE_ID
from .contracts import ProvenanceRecord, stable_json_dumps
from .errors import (
    EviTEAError,
    InputContractError,
    OutputFinalizationError,
    OutputStateError,
    ResumeStateError,
)
from .callset import CallsetState, calls_bed_bytes, calls_tsv_bytes, unresolved_tsv_bytes
from .outputs import (
    write_bytes,
    write_run_log,
    write_run_manifest,
    write_run_summary,
)
from .provenance import (
    RunManifestContext,
    build_run_manifest,
    finalize_output_hashes,
    hash_file,
)

CANONICAL_OUTPUTS = (
    "calls.native.tsv",
    "calls.tsv",
    "calls.bed",
    "unresolved.tsv",
    "run_manifest.json",
    "run_summary.json",
    "run.log",
)
BREAKPOINTS_OUTPUT = "breakpoints.tsv"  # written only when refinement.breakpoints is on
ANNOTATIONS_OUTPUT = "call_annotations.tsv"  # written only when refinement.call_annotations is on
CHECKPOINT_SCHEMA = "evitea.checkpoint.v1"
RAW_SAMPLE_DIRECTORY = "prepared-sample"
RAW_SAMPLE_STAGING_PREFIX = ".sample-build-"

RawSampleBuilder = Callable[..., PreparedSample]


@dataclass(frozen=True, slots=True)
class PreparedAssets:
    reference: PreparedReference
    te: PreparedTEAssets
    sample: PreparedSample | None = None


@dataclass(frozen=True, slots=True)
class CheckpointIdentity:
    relevant_config: Mapping[str, Any]
    upstream_hashes: Mapping[str, str]
    backend_versions: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class RunResult:
    output_directory: Path
    state: CallsetState
    manifest: ProvenanceRecord
    resumed_stages: tuple[str, ...] = ()


def _canonical_bytes(value: Any) -> bytes:
    return (stable_json_dumps(value) + "\n").encode("utf-8")


def _atomic_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as handle:
            temporary_name = handle.name
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except OSError:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
        raise


def prepare_output_directory(path: str | Path, *, resume: bool) -> Path:
    """Apply absent/empty/incomplete/completed/arbitrary directory policy."""

    output = Path(path).expanduser().resolve()
    if not output.exists():
        output.mkdir(parents=True)
        return output
    if not output.is_dir():
        raise OutputStateError(f"output path is not a directory: {output}")
    entries = tuple(output.iterdir())
    if not entries:
        return output
    names = {entry.name for entry in entries}
    if "run_manifest.json" in names or set(CANONICAL_OUTPUTS).issubset(names):
        raise OutputStateError("completed canonical output cannot be overwritten")
    if names == {".evitea-work"} or names <= {".evitea-work", "run.log"}:
        if not resume:
            raise ResumeStateError(
                "incomplete EviTEA work state requires explicit resume"
            )
        return output
    raise OutputStateError("refusing arbitrary nonempty output directory")


class CheckpointManager:
    """Transactional JSON checkpoints with completion marker written last."""

    def __init__(self, work_directory: str | Path, *, resume: bool) -> None:
        self.root = Path(work_directory).expanduser().resolve()
        self.resume = resume
        self.committed_stages: list[str] = []
        self.resumed_stages: list[str] = []
        self.root.mkdir(parents=True, exist_ok=True)

    def _stage_directory(self, stage: str) -> Path:
        if not stage or any(character in stage for character in "/\\"):
            raise InputContractError("checkpoint stage must be a simple non-empty name")
        return self.root / stage

    def commit(
        self,
        stage: str,
        payload: Mapping[str, Any],
        identity: CheckpointIdentity,
    ) -> Path:
        stage_directory = self._stage_directory(stage)
        if stage_directory.exists() and any(stage_directory.iterdir()):
            raise ResumeStateError(f"checkpoint stage already contains state: {stage}")
        stage_directory.mkdir(parents=True, exist_ok=True)
        data = _canonical_bytes(payload)
        data_hash = hashlib.sha256(data).hexdigest()
        manifest = {
            "schema_version": CHECKPOINT_SCHEMA,
            "architecture_id": ARCHITECTURE_ID,
            "stage": stage,
            "relevant_config": identity.relevant_config,
            "upstream_hashes": identity.upstream_hashes,
            "backend_versions": identity.backend_versions,
            "counts": {"payload_keys": len(payload)},
            "outputs": {"stage_data.json": data_hash},
            "completion_status": "COMPLETE",
        }
        _atomic_bytes(stage_directory / "stage_data.json", data)
        _atomic_bytes(
            stage_directory / "stage_manifest.json", _canonical_bytes(manifest)
        )
        _atomic_bytes(stage_directory / "complete", b"complete\n")
        self.committed_stages.append(stage)
        return stage_directory

    def load(
        self,
        stage: str,
        identity: CheckpointIdentity,
    ) -> dict[str, Any] | None:
        stage_directory = self._stage_directory(stage)
        if not stage_directory.exists():
            return None
        if not self.resume:
            raise ResumeStateError(
                f"checkpoint exists but resume was not requested: {stage}"
            )
        marker = stage_directory / "complete"
        manifest_path = stage_directory / "stage_manifest.json"
        data_path = stage_directory / "stage_data.json"
        if not marker.is_file() or not manifest_path.is_file() or not data_path.is_file():
            raise ResumeStateError(f"checkpoint is partial or uncommitted: {stage}")
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            payload = json.loads(data_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ResumeStateError(f"checkpoint is unreadable: {stage}") from exc
        expected_identity = {
            "schema_version": CHECKPOINT_SCHEMA,
            "architecture_id": ARCHITECTURE_ID,
            "stage": stage,
            "relevant_config": identity.relevant_config,
            "upstream_hashes": identity.upstream_hashes,
            "backend_versions": identity.backend_versions,
        }
        for key, expected in expected_identity.items():
            if manifest.get(key) != expected:
                raise ResumeStateError(
                    f"checkpoint identity mismatch for {stage}: {key}"
                )
        actual_hash = hashlib.sha256(data_path.read_bytes()).hexdigest()
        if manifest.get("outputs", {}).get("stage_data.json") != actual_hash:
            raise ResumeStateError(f"checkpoint output hash mismatch: {stage}")
        if manifest.get("completion_status") != "COMPLETE":
            raise ResumeStateError(f"checkpoint is not complete: {stage}")
        if not isinstance(payload, dict):
            raise ResumeStateError(f"checkpoint payload must be an object: {stage}")
        self.resumed_stages.append(stage)
        return payload


def prepare_assets(config: EviTEAConfig) -> PreparedAssets:
    """Load and cross-check prepared assets without performing scientific work."""

    if config.assets.reference_assets is None or config.assets.te_assets is None:
        raise InputContractError("reference_assets and te_assets are required")
    from .assets.reference import load_prepared_reference
    from .assets.te import load_prepared_te_assets

    reference = load_prepared_reference(config.assets.reference_assets)
    te = load_prepared_te_assets(config.assets.te_assets)
    verify_reference_asset(reference)
    verify_te_assets(te, reference)
    sample = None
    if config.inputs.sample_assets is not None:
        from .backends.exact_query import load_prepared_sample

        sample = load_prepared_sample(config.inputs.sample_assets)
        verify_prepared_sample(sample)
    return PreparedAssets(reference, te, sample)


def derive_sample_id(read1: str | Path, sample_id: str | None = None) -> str:
    """Use the same explicit-or-R1-basename rule for every sample preparation path."""

    if sample_id:
        return sample_id
    name = Path(read1).name.removesuffix(".gz")
    return Path(name).stem or "sample"


def _verify_raw_sample_identity(
    sample: PreparedSample,
    *,
    sample_id: str,
    read1_sha256: str,
    read2_sha256: str,
) -> None:
    if sample.sample_id != sample_id:
        raise ResumeStateError("prepared sample sample_id does not match current input")
    if sample.read1_raw_sha256 != read1_sha256:
        raise ResumeStateError("prepared sample read1 hash does not match current input")
    if sample.read2_raw_sha256 != read2_sha256:
        raise ResumeStateError("prepared sample read2 hash does not match current input")


def _prepare_raw_sample(
    config: EviTEAConfig,
    work_directory: Path,
    *,
    builder: RawSampleBuilder | None,
) -> PreparedSample:
    read1 = config.inputs.read1
    read2 = config.inputs.read2
    if read1 is None or read2 is None:
        raise InputContractError("raw sample preparation requires paired read paths")

    sample_id = derive_sample_id(read1, config.inputs.sample_id)
    read1_sha256 = hash_file(read1)
    read2_sha256 = hash_file(read2)
    sample_directory = work_directory / RAW_SAMPLE_DIRECTORY
    partial_builds = tuple(work_directory.glob(f"{RAW_SAMPLE_STAGING_PREFIX}*"))

    if config.runtime.resume:
        if partial_builds:
            raise ResumeStateError("prepared sample build is partial or uncommitted")
        if not sample_directory.exists():
            raise ResumeStateError("resume requires a completed prepared sample")
        try:
            sample = load_prepared_sample(sample_directory)
            verify_prepared_sample(sample)
        except (EviTEAError, OSError) as exc:
            raise ResumeStateError(
                "prepared sample is partial, unreadable, or incompatible"
            ) from exc
        _verify_raw_sample_identity(
            sample,
            sample_id=sample_id,
            read1_sha256=read1_sha256,
            read2_sha256=read2_sha256,
        )
        return sample

    if sample_directory.exists() or partial_builds:
        raise ResumeStateError(
            "prepared sample state already exists; explicit resume is required"
        )

    staging_directory = work_directory / (
        f"{RAW_SAMPLE_STAGING_PREFIX}{uuid.uuid4().hex[:12]}"
    )
    staging_directory.mkdir()
    sample_builder = builder or build_prepared_msbwt_sample
    try:
        sample = sample_builder(
            read1,
            read2,
            staging_directory,
            sample_id=sample_id,
            threads=config.runtime.resolved_threads,
        )
        if not isinstance(sample, PreparedSample):
            raise InputContractError("sample builder must return PreparedSample")
        if sample.root.resolve() != staging_directory.resolve():
            raise InputContractError("sample builder returned an unexpected root")
        verify_prepared_sample(sample)
        _verify_raw_sample_identity(
            sample,
            sample_id=sample_id,
            read1_sha256=read1_sha256,
            read2_sha256=read2_sha256,
        )
        if hash_file(read1) != read1_sha256 or hash_file(read2) != read2_sha256:
            raise InputContractError("raw FASTQ inputs changed during sample build")
        os.replace(staging_directory, sample_directory)
        completed_sample = load_prepared_sample(sample_directory)
        verify_prepared_sample(completed_sample)
        _verify_raw_sample_identity(
            completed_sample,
            sample_id=sample_id,
            read1_sha256=read1_sha256,
            read2_sha256=read2_sha256,
        )
        return completed_sample
    except Exception:
        if staging_directory.exists():
            shutil.rmtree(staging_directory)
        raise


def run_sample(
    config: EviTEAConfig,
    prepared_assets: PreparedAssets,
    call_engine: Callable[[EviTEAConfig, PreparedAssets, CheckpointManager], CallsetState],
    *,
    manifest_context: RunManifestContext,
    sample_builder: RawSampleBuilder | None = None,
    breakpoint_refiner: Callable[[CallsetState], Sequence[Any]] | None = None,
) -> RunResult:
    """Orchestrate a supplied scientific engine and finalize canonical outputs.

    ``breakpoint_refiner`` is required exactly when ``refinement.breakpoints``
    is on.  It never changes ``calls.native.tsv``; refined ranges go to the
    ``calls.tsv`` / ``calls.bed`` views and to ``breakpoints.tsv``.
    """

    if config.output.directory is None:
        raise InputContractError("output directory is required")
    if config.refinement.breakpoints != (breakpoint_refiner is not None):
        raise InputContractError(
            "refinement.breakpoints and the breakpoint refiner must be given together"
        )
    output = prepare_output_directory(
        config.output.directory, resume=config.runtime.resume
    )
    work = config.runtime.work_dir or (output / ".evitea-work")
    checkpoints = CheckpointManager(work, resume=config.runtime.resume)
    engine_assets = prepared_assets
    finalized_manifest_context = manifest_context
    if config.inputs.sample_assets is None and config.inputs.read1 is not None:
        sample = _prepare_raw_sample(
            config,
            checkpoints.root,
            builder=sample_builder,
        )
        engine_assets = replace(prepared_assets, sample=sample)
        prepared_identities = dict(manifest_context.prepared_identities)
        prepared_identities["sample"] = sample.identity_sha256
        finalized_manifest_context = replace(
            manifest_context,
            prepared_identities=prepared_identities,
        )
    started_at = manifest_context.started_at or datetime.now(UTC).isoformat()
    state = call_engine(config, engine_assets, checkpoints)
    if not isinstance(state, CallsetState):
        raise InputContractError("call engine must return CallsetState")
    breakpoint_rows = None
    breakpoint_pins = None
    annotations = None
    refined = None
    if breakpoint_refiner is not None:
        from .breakpoints.refine import RefinementResult, refined_coordinates

        result = RefinementResult.of(breakpoint_refiner(state))
        if (result.pins is not None) != config.refinement.pin_one_sided:
            raise InputContractError("refinement.pin_one_sided and the refiner's pins must be given together")
        if (result.annotations is not None) != config.refinement.call_annotations:
            raise InputContractError("refinement.call_annotations and the refiner's annotations must be given together")
        breakpoint_rows, breakpoint_pins, annotations = result.rows, result.pins, result.annotations
        refined = refined_coordinates(breakpoint_rows)
        if annotations is not None:
            from .shadow_layers import ShadowCallset

            native_ids = [row["event_id"] for row in ShadowCallset.from_bytes(state.native).rows]
            if [row.event_id for row in annotations] != native_ids:
                raise InputContractError("call annotations must cover every native call in native order")
    ended_at = manifest_context.ended_at or datetime.now(UTC).isoformat()

    stage_history = (
        *manifest_context.stage_history,
        *(
            {"stage": stage, "status": "RESUMED"}
            for stage in checkpoints.resumed_stages
        ),
        *(
            {"stage": stage, "status": "COMMITTED"}
            for stage in checkpoints.committed_stages
        ),
    )
    finalized_context = replace(
        finalized_manifest_context,
        run_parameters={**finalized_manifest_context.run_parameters, **state.run_parameters},
        stage_history=stage_history,
        started_at=started_at,
        ended_at=ended_at,
    )

    try:
        write_bytes(output / "calls.native.tsv", state.native)
        write_bytes(output / "calls.tsv", calls_tsv_bytes(state, refined))
        write_bytes(
            output / "calls.bed",
            calls_bed_bytes(state, dict(prepared_assets.reference.contigs), refined),
        )
        write_bytes(output / "unresolved.tsv", unresolved_tsv_bytes(state))
        summary = state.summary()
        hashed = tuple(name for name in CANONICAL_OUTPUTS if name != "run_manifest.json")
        log_lines = (
            "run_completed=true",
            f"profile={state.profile}",
            f"final_call_count={summary['final_call_count']}",
            f"native_callset_sha256={summary['native_callset_sha256']}",
        )
        if breakpoint_rows is not None:
            from .breakpoints.refine import breakpoints_tsv_bytes, refinement_summary

            write_bytes(output / BREAKPOINTS_OUTPUT, breakpoints_tsv_bytes(breakpoint_rows, pins=breakpoint_pins))
            summary["breakpoint_refinement"] = refinement_summary(breakpoint_rows)
            hashed = (*hashed, BREAKPOINTS_OUTPUT)
            log_lines = (*log_lines, f"breakpoint_refined_calls={len(refined)}")
            if breakpoint_pins is not None:
                summary["breakpoint_refinement"]["pinned_calls"] = len(breakpoint_pins)
                log_lines = (*log_lines, f"breakpoint_pinned_calls={len(breakpoint_pins)}")
        if annotations is not None:
            from .breakpoints.annotate import call_annotations_tsv_bytes

            write_bytes(output / ANNOTATIONS_OUTPUT, call_annotations_tsv_bytes(annotations))
            hashed = (*hashed, ANNOTATIONS_OUTPUT)
            log_lines = (*log_lines, f"annotated_calls={len(annotations)}")
        write_run_summary(output / "run_summary.json", summary)
        write_run_log(output / "run.log", log_lines)
        manifest = build_run_manifest(finalized_context)
        manifest = finalize_output_hashes(
            manifest,
            output,
            hashed,
            output_counts={
                "calls": summary["final_call_count"],
                "unresolved": summary["unresolved_count"],
            },
        )
        write_run_manifest(output / "run_manifest.json", manifest)
    except (InputContractError, OSError) as exc:
        raise OutputFinalizationError(
            "failed to finalize canonical outputs and provenance"
        ) from exc
    return RunResult(output, state, manifest, tuple(checkpoints.resumed_stages))


# ---------------------------------------------------------------------------
# Default scientific engine assembly (incumbent + frozen v3 extension)
# ---------------------------------------------------------------------------

ENGINE_STAGE = "engine"


def engine_inputs(config: EviTEAConfig, assets: PreparedAssets) -> "EngineInputs":
    """Resolve prepared assets and runtime tools into the engine's inputs."""

    from .engine import EngineInputs

    if config.inputs.msbwt_index is not None:
        msbwt_index = config.inputs.msbwt_index
        sample_id = config.inputs.sample_id
    elif assets.sample is not None:
        profile = assets.sample.backend_profile
        msbwt_index = assets.sample.root / str(profile.get("index_directory", "index"))
        sample_id = config.inputs.sample_id or assets.sample.sample_id
        if sample_id != assets.sample.sample_id:
            raise InputContractError("sample_id differs from the prepared sample's id")
    else:
        raise InputContractError("the engine needs a prepared sample or an msBWT index")
    if not (Path(msbwt_index) / "about.npy").is_file():
        raise InputContractError(f"msBWT index has no about.npy: {msbwt_index}")
    bowtie2 = config.runtime.bowtie2 or (
        Path(found) if (found := shutil.which("bowtie2")) else None
    )
    if bowtie2 is None:
        from .errors import OperationalResourceError

        raise OperationalResourceError("bowtie2 executable is unavailable")
    reference, te = assets.reference, assets.te
    read_length = resolve_read_length(
        config.runtime.read_length,
        assets.sample.backend_profile if assets.sample is not None and config.inputs.msbwt_index is None else None,
        msbwt_index,
    )
    return EngineInputs(
        sample_id=str(sample_id),
        msbwt_index=Path(msbwt_index),
        te_fasta=te.source_path("te_library"),
        te_family_map=te.source_path("family_map"),
        reference_fasta=reference.reference_fasta_path,
        repeatmasker_bed=reference.repeat_annotation_path,
        host_index_prefix=reference.bowtie2_index_prefix,
        terminal_codebook=te.language_path("terminal_final"),
        internal_codebook=te.language_path("internal_only"),
        bowtie2=Path(bowtie2),
        bowtie2_build=config.runtime.bowtie2_build,
        profile=config.inputs.profile,
        threads=config.runtime.resolved_threads,
        workers=config.runtime.resolved_workers,
        edlib_python=config.runtime.edlib_python,
        complete_routing=config.refinement.complete_routing,
        max_rescue_interval_mass=config.runtime.max_rescue_interval_mass,
        read_length=read_length.length,
        read_length_source=read_length.source,
        fragment_q95=config.runtime.fragment_q95,
        fragment_q99=config.runtime.fragment_q99,
    )


def default_call_engine(
    config: EviTEAConfig,
    assets: PreparedAssets,
    checkpoints: CheckpointManager,
    *,
    progress: Callable[[str], None] | None = None,
) -> CallsetState:
    """Run the clean incumbent and the v3 extension; return the native callset.

    Stage checkpoints live under ``<work>/engine``; an existing engine work
    directory is reused only with explicit resume.
    """

    from .engine import run_engine

    inputs = engine_inputs(config, assets)
    engine_work = checkpoints.root / ENGINE_STAGE / "work"
    if engine_work.exists() and not checkpoints.resume:
        raise ResumeStateError("engine work state already exists; explicit resume is required")
    identity = CheckpointIdentity(
        relevant_config=inputs.identity(),
        upstream_hashes={
            "reference": assets.reference.identity_sha256,
            "te": assets.te.identity_sha256,
            **({"sample": assets.sample.identity_sha256} if assets.sample is not None else {}),
        },
        backend_versions={},
    )
    result = run_engine(inputs, engine_work, progress=progress)
    state = CallsetState(
        sample_id=inputs.sample_id,
        native=result.final.data,
        unresolved=result.unresolved,
        layer_rows=result.layer_rows,
        layer_sha256=result.layer_sha256,
        counts=result.counts,
        profile=result.profile.name,
        run_parameters={**inputs.run_parameters(), **(result.fragment.record() if result.fragment else {})},
    )
    record = {
        "native_callset_sha256": state.native_sha256,
        "layer_rows": dict(state.layer_rows),
        "layer_sha256": dict(state.layer_sha256),
    }
    previous = checkpoints.load("engine-summary", identity) if checkpoints.resume else None
    if previous is None:
        checkpoints.commit("engine-summary", record, identity)
    elif previous != record:
        raise ResumeStateError("resumed engine reproduced a different callset")
    return state
