from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

import evitea.cli as cli
from evitea.backends.exact_query import (
    make_read_identity_map,
    prepare_sample_manifest,
    write_read_identity_map,
)
from evitea.backends.protocols import ReadIdentity
from evitea.cli import _validate_run
from evitea.config import EviTEAConfig, InputConfig, OutputConfig, RuntimeConfig
from evitea.constants import ARCHITECTURE_ID
from evitea.contracts import ReadMate
from evitea.errors import InputContractError, ResumeStateError
from evitea.callset import CallsetState
from evitea.pipeline import (
    CheckpointIdentity,
    CheckpointManager,
    PreparedAssets,
    prepare_output_directory,
    run_sample,
)
from evitea.provenance import RunManifestContext
from evitea.provenance import build_run_manifest, finalize_output_hashes


def _identity(*, config_value: str = "same") -> CheckpointIdentity:
    return CheckpointIdentity(
        relevant_config={"mode": config_value},
        upstream_hashes={"input": "a" * 64},
        backend_versions={"exact_query": "fixture"},
    )


NATIVE = (
    "event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\t"
    "left_support\tright_support\ttotal_support\tevidence_type\tscore\tfragment_ids\n"
    "V2FCRC1_EVENT_000001\t2L\t10\t10\troo\t['roo']\t10\t-1\t0\t2\t0\t2\t"
    "one_sided_left\t2.0\t6,1\n"
).encode()


def _state() -> CallsetState:
    from evitea.shadow_layers import UnresolvedShadowEvent

    unresolved = UnresolvedShadowEvent(
        "INTERNAL21", "roo:2L:40-45", "UNRESOLVED_NEARBY_DIFFERENT_FAMILY_CALL", "2L", "roo", 40, 45, ("X",)
    )
    return CallsetState("sample", NATIVE, (unresolved,), {"incumbent": 1})


def _context(config: EviTEAConfig) -> RunManifestContext:
    return RunManifestContext(
        normalized_config=config.normalized(),
        command_line=("evitea", "call"),
        source_identity={"architecture": ARCHITECTURE_ID},
        input_assets={},
        prepared_identities={"reference": "r", "te": "t"},
        backend_versions={"exact_query": "fixture"},
        resolved_threads=config.runtime.resolved_threads,
    )


def _assets() -> PreparedAssets:
    return PreparedAssets(
        reference=SimpleNamespace(contigs=(("2L", 100),)),
        te=SimpleNamespace(),
    )


def test_checkpoint_marker_last_identity_and_hash_enforcement(tmp_path: Path) -> None:
    writer = CheckpointManager(tmp_path / "work", resume=False)
    stage = writer.commit("projection", {"pairs": [1, 2]}, _identity())
    assert (stage / "complete").read_bytes() == b"complete\n"
    reader = CheckpointManager(tmp_path / "work", resume=True)
    assert reader.load("projection", _identity()) == {"pairs": [1, 2]}
    with pytest.raises(ResumeStateError, match="identity mismatch"):
        reader.load("projection", _identity(config_value="changed"))
    (stage / "stage_data.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(ResumeStateError, match="hash mismatch"):
        reader.load("projection", _identity())


def test_partial_checkpoint_is_never_reusable(tmp_path: Path) -> None:
    stage = tmp_path / "work" / "projection"
    stage.mkdir(parents=True)
    (stage / "stage_data.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(ResumeStateError, match="partial"):
        CheckpointManager(tmp_path / "work", resume=True).load(
            "projection", _identity()
        )


def test_output_directory_policy_fails_closed(tmp_path: Path) -> None:
    fresh = prepare_output_directory(tmp_path / "fresh", resume=False)
    assert fresh.is_dir()
    arbitrary = tmp_path / "arbitrary"
    arbitrary.mkdir()
    (arbitrary / "notes.txt").write_text("mine", encoding="utf-8")
    with pytest.raises(InputContractError, match="nonempty"):
        prepare_output_directory(arbitrary, resume=False)
    incomplete = tmp_path / "incomplete"
    (incomplete / ".evitea-work").mkdir(parents=True)
    with pytest.raises(ResumeStateError, match="explicit resume"):
        prepare_output_directory(incomplete, resume=False)
    assert prepare_output_directory(incomplete, resume=True) == incomplete.resolve()
    completed = tmp_path / "completed"
    completed.mkdir()
    (completed / "run_manifest.json").write_text("{}\n", encoding="utf-8")
    with pytest.raises(InputContractError, match="cannot be overwritten"):
        prepare_output_directory(completed, resume=True)


def test_manifest_hashes_outputs_but_never_itself(tmp_path: Path) -> None:
    config = EviTEAConfig(output=OutputConfig(tmp_path))
    manifest = build_run_manifest(_context(config))
    (tmp_path / "calls.tsv").write_bytes(b"stable\n")
    finalized = finalize_output_hashes(manifest, tmp_path, ("calls.tsv",))
    assert tuple(finalized.values["output_hashes"]) == ("calls.tsv",)
    with pytest.raises(InputContractError, match="recursively"):
        finalize_output_hashes(manifest, tmp_path, ("run_manifest.json",))


def test_resume_and_thread_count_preserve_scientific_output_bytes(
    tmp_path: Path,
) -> None:
    work = tmp_path / "shared-work"
    first_config = EviTEAConfig(
        runtime=RuntimeConfig(threads=1, work_dir=work),
        output=OutputConfig(tmp_path / "first"),
    )

    def first_engine(_config, _assets, checkpoints):
        checkpoints.commit("integration", {"complete": True}, _identity())
        return _state()

    first_result = run_sample(
        first_config,
        _assets(),
        first_engine,
        manifest_context=_context(first_config),
    )

    resumed_config = EviTEAConfig(
        runtime=RuntimeConfig(threads=4, work_dir=work, resume=True),
        output=OutputConfig(tmp_path / "resumed"),
    )

    def resumed_engine(_config, _assets, checkpoints):
        assert checkpoints.load("integration", _identity()) == {"complete": True}
        return _state()

    resumed_result = run_sample(
        resumed_config,
        _assets(),
        resumed_engine,
        manifest_context=_context(resumed_config),
    )
    assert first_result.manifest.values["started_at"]
    assert first_result.manifest.values["ended_at"]
    assert first_result.manifest.values["stage_history"] == [
        {"stage": "integration", "status": "COMMITTED"}
    ]
    assert resumed_result.resumed_stages == ("integration",)
    assert resumed_result.manifest.values["stage_history"] == [
        {"stage": "integration", "status": "RESUMED"}
    ]
    _validate_run(tmp_path / "first")
    _validate_run(tmp_path / "resumed")
    for filename in (
        "calls.native.tsv",
        "calls.tsv",
        "calls.bed",
        "unresolved.tsv",
        "run_summary.json",
        "run.log",
    ):
        assert (tmp_path / "first" / filename).read_bytes() == (
            tmp_path / "resumed" / filename
        ).read_bytes()


def _raw_sample_builder(calls: list[tuple[Path, str, int]]):
    def build(read1, read2, output_dir, *, sample_id, threads):
        output = Path(output_dir)
        calls.append((output, sample_id, threads))
        (output / "identity").mkdir()
        (output / "index").mkdir()
        write_read_identity_map(
            output / "identity" / "reads.json",
            make_read_identity_map(
                (
                    ReadIdentity(0, 0, ReadMate.READ1),
                    ReadIdentity(1, 0, ReadMate.READ2),
                )
            ),
        )
        (output / "index" / "sample.msbwt").write_bytes(b"test index")
        return prepare_sample_manifest(
            output,
            sample_id=sample_id,
            pair_count=1,
            read1_path=read1,
            read2_path=read2,
            identity_map_relative_path="identity/reads.json",
            backend_profile={"builder": "unit-test"},
            index_relative_paths=("index/sample.msbwt",),
        )

    return build


def _raw_config(tmp_path: Path, *, resume: bool = False, sample_id=None):
    tmp_path = tmp_path.resolve()
    reads = tmp_path / "reads"
    reads.mkdir(exist_ok=True)
    read1 = reads / "specimen_R1.fastq.gz"
    read2 = reads / "specimen_R2.fastq.gz"
    if not read1.exists():
        read1.write_bytes(b"R1 fixture")
        read2.write_bytes(b"R2 fixture")
    return EviTEAConfig(
        inputs=InputConfig(sample_id=sample_id, read1=read1, read2=read2),
        runtime=RuntimeConfig(
            threads=3,
            work_dir=tmp_path / "shared-work",
            resume=resume,
        ),
        output=OutputConfig(tmp_path / ("resumed-output" if resume else "output")),
    )


def test_cli_raw_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _raw_config(tmp_path)
    call_assets = PreparedAssets(
        reference=SimpleNamespace(contigs=(("2L", 100),), identity_sha256="r"),
        te=SimpleNamespace(identity_sha256="t"),
    )
    builder_calls: list[tuple[Path, str, int]] = []
    engine_samples = []
    monkeypatch.setattr(cli, "_call_config", lambda _args: config)
    monkeypatch.setattr(cli, "prepare_assets", lambda _config: call_assets)
    monkeypatch.setattr(cli, "_source_identity", lambda: {"source": "fixture"})
    monkeypatch.setattr(cli, "_distribution_version", lambda *_names: "fixture")
    monkeypatch.setattr(cli, "_executable_version", lambda _executable: "fixture")

    def engine(_config, assets, _checkpoints):
        engine_samples.append(assets.sample)
        assert assets.sample is not None
        assert assets.sample.root == config.runtime.work_dir / "prepared-sample"
        return _state()

    result = cli._call(
        SimpleNamespace(command_line=("evitea", "call")),
        engine,
        _raw_sample_builder(builder_calls),
    )

    assert len(builder_calls) == 1
    staging_path, sample_id, threads = builder_calls[0]
    assert staging_path.parent == config.runtime.work_dir
    assert staging_path.name.startswith(".sample-build-")
    assert sample_id == "specimen_R1"
    assert threads == 3
    assert len(engine_samples) == 1
    assert engine_samples[0] is not None
    assert engine_samples[0].root == config.runtime.work_dir / "prepared-sample"
    assert result.manifest.values["prepared_identities"]["sample"] == (
        engine_samples[0].identity_sha256
    )
    assert not staging_path.exists()
    assert (config.runtime.work_dir / "prepared-sample" / "sample_manifest.json").is_file()


@pytest.mark.parametrize("mismatch", ("read1", "read2", "sample_id"))
def test_resume_raw_mismatch(
    tmp_path: Path, mismatch: str
) -> None:
    first_config = _raw_config(tmp_path)
    builder_calls: list[tuple[Path, str, int]] = []
    run_sample(
        first_config,
        _assets(),
        lambda _config, _assets, _checkpoints: _state(),
        manifest_context=_context(first_config),
        sample_builder=_raw_sample_builder(builder_calls),
    )
    resumed_config = _raw_config(tmp_path, resume=True)
    if mismatch == "read1":
        resumed_config.inputs.read1.write_bytes(b"changed R1")
    elif mismatch == "read2":
        resumed_config.inputs.read2.write_bytes(b"changed R2")
    else:
        resumed_config = EviTEAConfig(
            inputs=InputConfig(
                sample_id="different-sample",
                read1=resumed_config.inputs.read1,
                read2=resumed_config.inputs.read2,
            ),
            runtime=resumed_config.runtime,
            output=resumed_config.output,
        )

    with pytest.raises(ResumeStateError):
        run_sample(
            resumed_config,
            _assets(),
            lambda _config, _assets, _checkpoints: _state(),
            manifest_context=_context(resumed_config),
            sample_builder=lambda *_args, **_kwargs: pytest.fail(
                "resume must load the completed prepared sample"
            ),
        )
    assert len(builder_calls) == 1


def test_resume_raw_partial(tmp_path: Path) -> None:
    config = _raw_config(tmp_path, resume=True)
    partial = config.runtime.work_dir / ".sample-build-interrupted"
    partial.mkdir(parents=True)
    with pytest.raises(ResumeStateError, match="partial or uncommitted"):
        run_sample(
            config,
            _assets(),
            lambda _config, _assets, _checkpoints: pytest.fail(
                "engine must not receive a partial sample"
            ),
            manifest_context=_context(config),
            sample_builder=lambda *_args, **_kwargs: pytest.fail(
                "resume must not rebuild a partial sample"
            ),
        )


def test_raw_build_failure(tmp_path: Path) -> None:
    config = _raw_config(tmp_path)

    def fail_after_partial_write(_read1, _read2, output_dir, **_kwargs):
        Path(output_dir, "partial.txt").write_text("incomplete", encoding="utf-8")
        raise RuntimeError("builder interrupted")

    with pytest.raises(RuntimeError, match="builder interrupted"):
        run_sample(
            config,
            _assets(),
            lambda _config, _assets, _checkpoints: pytest.fail(
                "engine must not receive an incomplete sample"
            ),
            manifest_context=_context(config),
            sample_builder=fail_after_partial_write,
        )
    assert not (config.runtime.work_dir / "prepared-sample").exists()
    assert not tuple(config.runtime.work_dir.glob(".sample-build-*"))

