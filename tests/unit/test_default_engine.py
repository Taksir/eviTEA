"""Default `evitea call` assembly: profiles, wiring, languages, outputs."""

from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import pytest

import evitea.engine as engine_module
import evitea.pipeline as pipeline_module
from evitea.anchor_pe.alignments import PRODUCTION_HOST_READ_NAMES
from evitea.callset import CallsetState
from evitea.config import EviTEAConfig, InputConfig, OutputConfig, RuntimeConfig
from evitea.engine import EngineInputs, resolve_profile, run_engine
from evitea.sample_geometry import FragmentSizes
from evitea.errors import AssetCompatibilityError, InputContractError, ResumeStateError
from evitea.pipeline import CheckpointManager, PreparedAssets, default_call_engine, engine_inputs
from evitea.shadow_layers import ShadowCallset, ShadowLabels

HEADER = (
    "event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\t"
    "left_support\tright_support\ttotal_support\tevidence_type\tscore\tfragment_ids"
)
INCUMBENT = (HEADER + "\nV2FCRC1_EVENT_000001\tc1\t10\t10\troo\t['roo']\t10\t-1\t0\t2\t0\t2\t"
             "one_sided_left\t2.0\t1,2\n").encode()
TERMINAL = (
    "family\tte_id\tprobe_rank_within_te_record\tterminal_side\tk\tcanonical_strand_probe\t"
    "query_sequences_both_strands\tsource_start\tsource_end\tpositionally_interpretable\n"
    "roo\tFBte1\t1\t3p\t21\tACGTACGTTGCAACGTTGCAA\tACGTACGTTGCAACGTTGCAA;TTGCAACGTTGCAACGTACGT\t1\t22\t1\n"
)
INTERNAL = (
    "canonical_strand_probe\tte_id\tfamily\tacquisition_safe\tinternal_only_query_unit\n"
    "ACGTACGTTGCAACGTTGCTT\tFBte1\troo\t1\t1\n"
)


def test_production_profile_uses_pid_names_and_sample_labels():
    profile = resolve_profile(None, "my-sample.2")
    assert profile.name == "production"
    assert profile.host_read_names == PRODUCTION_HOST_READ_NAMES == "pid"
    assert profile.labels == ShadowLabels.for_sample("my-sample.2") == ShadowLabels("MYSAMPLE2_", "MYSAMPLE2", "MYSAMPLE2", "MYSAMPLE2")
    assert not profile.reproduction_only


def test_reproduction_profiles_are_bound_to_their_samples():
    aka = resolve_profile("AKA017", "AKA017")
    assert aka.host_read_names == "aka-historical"
    assert aka.labels == ShadowLabels("", "AKA", "AKA", "AKA")
    iso = resolve_profile("ISO1", "ISO1")
    assert iso.host_read_names == "pid"
    assert iso.labels == ShadowLabels("ISO1_", "ISO1", "ISO", "ISO")
    with pytest.raises(InputContractError, match="reproduces only"):
        resolve_profile("AKA017", "NEW1")
    with pytest.raises(InputContractError, match="unknown run profile"):
        resolve_profile("aka-historical", "AKA017")


def _inputs(tmp_path: Path, profile: str = "production", sample_id: str = "S1") -> EngineInputs:
    return EngineInputs(
        sample_id=sample_id, msbwt_index=tmp_path / "msbwt", te_fasta=tmp_path / "te.fa",
        te_family_map=tmp_path / "map.tsv", reference_fasta=tmp_path / "ref.fa",
        repeatmasker_bed=tmp_path / "rm.bed", host_index_prefix=tmp_path / "idx",
        terminal_codebook=tmp_path / "terminal.tsv", internal_codebook=tmp_path / "internal.tsv",
        bowtie2=tmp_path / "bowtie2", bowtie2_build=tmp_path / "bowtie2-build", profile=profile,
    )


def _fake_extension(incumbent: bytes, labels) -> SimpleNamespace:
    base = ShadowCallset.from_bytes(incumbent)
    layer = SimpleNamespace(layer="INTERNAL21", callset=base, unresolved=())
    return SimpleNamespace(
        layers=(layer,),
        anchor_pe=SimpleNamespace(counts={"a": 1}),
        full_memory=SimpleNamespace(counts={"f": 2}),
        internal21=SimpleNamespace(counts={"i": 3}),
        labels=labels,
    )


@pytest.mark.parametrize(("profile", "sample", "names"), (("production", "S1", "pid"), ("AKA017", "AKA017", "aka-historical")))
def test_run_engine_wires_incumbent_bytes_codebooks_and_profile(tmp_path, monkeypatch, profile, sample, names):
    seen = {}

    def incumbent(inputs, sample_index, work, progress=None):
        seen["incumbent"] = (inputs, work)
        return SimpleNamespace(calls=(1,), counts={"incumbent_calls": 1}, legacy_tsv=lambda: INCUMBENT)

    def extension(anchor_inputs, internal_codebook, sample_index, incumbent_bytes, labels, work, progress=None):
        seen["extension"] = (anchor_inputs, internal_codebook, incumbent_bytes, labels, work)
        return _fake_extension(incumbent_bytes, labels)

    monkeypatch.setattr(engine_module, "MsbwtSample", lambda path: ("msbwt", path))
    monkeypatch.setattr(engine_module, "run_incumbent_caller", incumbent)
    monkeypatch.setattr(engine_module, "run_extension", extension)
    monkeypatch.setattr(
        engine_module, "estimate_fragment_sizes", lambda *a, **k: FragmentSizes(300.0, 420, 5000, "estimate")
    )
    inputs = _inputs(tmp_path, profile, sample)
    result = run_engine(inputs, tmp_path / "work")
    anchor_inputs, internal_codebook, incumbent_bytes, labels, work = seen["extension"]
    # the sample's fragment-size estimate reaches both layers
    assert (seen["incumbent"][0].fragment_q95, seen["incumbent"][0].fragment_q99) == (300.0, 420)
    assert anchor_inputs.fragment_q99 == 420
    assert result.fragment.record() == {
        "fragment_q95": 300.0, "fragment_q99": 420, "fragment_pairs": 5000, "fragment_source": "estimate",
    }
    assert anchor_inputs.host_read_names == names
    assert anchor_inputs.safe_codebook == inputs.terminal_codebook
    assert internal_codebook == inputs.internal_codebook
    assert incumbent_bytes == INCUMBENT
    assert labels == resolve_profile(profile, sample).labels
    assert work == tmp_path / "work" / "extension"
    assert seen["incumbent"][1] == tmp_path / "work" / "incumbent"
    assert result.final.data == INCUMBENT
    assert result.layer_rows == {"incumbent": 1, "internal21": 1}
    digest = hashlib.sha256(INCUMBENT).hexdigest()
    assert result.layer_sha256 == {"incumbent": digest, "internal21": digest}


def _assets(tmp_path: Path) -> PreparedAssets:
    (tmp_path / "msbwt").mkdir(exist_ok=True)
    (tmp_path / "msbwt" / "about.npy").write_bytes(b"x")
    reference = SimpleNamespace(
        identity_sha256="r", contigs=(("c1", 100),), reference_fasta_path=tmp_path / "ref.fa",
        repeat_annotation_path=tmp_path / "rm.bed", bowtie2_index_prefix=tmp_path / "idx",
    )
    te = SimpleNamespace(
        identity_sha256="t",
        source_path=lambda prefix: tmp_path / f"{prefix}.src",
        language_path=lambda key: tmp_path / f"{key}.tsv",
    )
    return PreparedAssets(reference, te)


def _config(tmp_path: Path, *, resume=False, profile="production") -> EviTEAConfig:
    return EviTEAConfig(
        inputs=InputConfig(sample_id="S1", msbwt_index=tmp_path / "msbwt", profile=profile),
        runtime=RuntimeConfig(threads=2, workers=1, bowtie2=tmp_path / "bowtie2", resume=resume),
        output=OutputConfig(tmp_path / "out"),
    )


def test_engine_inputs_resolve_prepared_assets(tmp_path):
    inputs = engine_inputs(_config(tmp_path), _assets(tmp_path))
    assert inputs.msbwt_index == tmp_path / "msbwt"
    assert inputs.terminal_codebook == tmp_path / "terminal_final.tsv"
    assert inputs.internal_codebook == tmp_path / "internal_only.tsv"
    assert inputs.te_fasta == tmp_path / "te_library.src"
    assert inputs.te_family_map == tmp_path / "family_map.src"
    assert inputs.repeatmasker_bed == tmp_path / "rm.bed"
    assert (inputs.threads, inputs.workers) == (2, 1)
    with pytest.raises(InputContractError, match="reproduces only"):
        engine_inputs(_config(tmp_path, profile="ISO1"), _assets(tmp_path))


def test_msbwt_index_mode_requires_sample_id(tmp_path):
    with pytest.raises(InputContractError, match="sample_id"):
        InputConfig(msbwt_index=tmp_path)
    with pytest.raises(InputContractError, match="not several"):
        InputConfig(sample_id="S", msbwt_index=tmp_path, sample_assets=tmp_path)


def test_default_call_engine_records_and_rechecks_the_native_callset(tmp_path, monkeypatch):
    runs = []

    def fake_run(inputs, work, progress=None):
        runs.append(work)
        Path(work).mkdir(parents=True, exist_ok=True)
        return SimpleNamespace(
            final=ShadowCallset.from_bytes(INCUMBENT), unresolved=(), layer_rows={"incumbent": 1},
            layer_sha256={"incumbent": hashlib.sha256(INCUMBENT).hexdigest()},
            counts={"incumbent": {"x": 1}}, profile=SimpleNamespace(name="production"),
            fragment=FragmentSizes(494.0, 805, 0, "option"),
        )

    monkeypatch.setattr(engine_module, "run_engine", fake_run)
    assets = _assets(tmp_path)
    state = default_call_engine(_config(tmp_path), assets, CheckpointManager(tmp_path / "w", resume=False))
    assert isinstance(state, CallsetState)
    assert state.native_sha256 == hashlib.sha256(INCUMBENT).hexdigest()
    assert runs == [tmp_path / "w" / "engine" / "work"]
    assert state.run_parameters["fragment_q99"] == 805
    assert state.run_parameters["fragment_source"] == "option"
    with pytest.raises(ResumeStateError, match="explicit resume"):
        default_call_engine(_config(tmp_path), assets, CheckpointManager(tmp_path / "w", resume=False))
    resumed = default_call_engine(_config(tmp_path, resume=True), assets, CheckpointManager(tmp_path / "w", resume=True))
    assert resumed.native == state.native


def test_prepared_te_languages_are_bound_and_tamper_evident(tmp_path):
    from evitea.assets.te import load_prepared_te_assets, prepare_te_assets, verify_te_assets
    from evitea.inputs import load_family_map, validate_te_fasta
    from test_phase3_asset_manifests import FIXTURES, _prepare_reference

    reference = _prepare_reference(tmp_path)
    (tmp_path / "terminal.tsv").write_text(TERMINAL, encoding="utf-8")
    (tmp_path / "internal.tsv").write_text(INTERNAL, encoding="utf-8")
    te = prepare_te_assets(
        validate_te_fasta(FIXTURES / "te.fa"), load_family_map(FIXTURES / "te_families.tsv"),
        FIXTURES / "te_families.tsv", reference, tmp_path / "t",
        terminal_codebook_path=tmp_path / "terminal.tsv", internal_codebook_path=tmp_path / "internal.tsv",
    )
    assert set(te.languages) == {"terminal_final", "internal_only"}
    assert te.languages["terminal_final"]["source"] == "supplied"
    assert te.language_path("internal_only").read_text(encoding="utf-8") == INTERNAL
    verify_te_assets(load_prepared_te_assets(tmp_path / "t"), reference)
    te.language_path("terminal_final").write_text(TERMINAL.replace("roo\tFBte1", "gypsy\tFBte1"), encoding="utf-8")
    with pytest.raises(AssetCompatibilityError):
        verify_te_assets(load_prepared_te_assets(tmp_path / "t"), reference)
    with pytest.raises(InputContractError, match="together"):
        prepare_te_assets(
            validate_te_fasta(FIXTURES / "te.fa"), load_family_map(FIXTURES / "te_families.tsv"),
            FIXTURES / "te_families.tsv", reference, tmp_path / "t2",
            terminal_codebook_path=tmp_path / "terminal.tsv",
        )


def test_missing_languages_fail_closed_with_guidance(tmp_path):
    from test_phase3_asset_manifests import _prepare_te_assets

    te, _ = _prepare_te_assets(tmp_path)
    with pytest.raises(AssetCompatibilityError, match="--terminal-codebook"):
        te.language_path("terminal_final")


def test_flybase_family_table_is_accepted(tmp_path):
    from evitea.inputs import is_flybase_family_table, load_family_map

    table = tmp_path / "fbte.tsv"
    table.write_bytes(b"FBte_ID\tFlyBase_Family\tStructural_Class\r\nFBte1\trooA\tLTR\r\nFBte2\tjockey\tLINE\r\n")
    assert is_flybase_family_table(table)
    assert dict(load_family_map(table).items()) == {"FBte1": "roo", "FBte2": "jockey"}


def test_recall_repair_read_names_are_a_profile():
    from evitea.engine import REPRODUCTION_PROFILES
    from evitea.incumbent.engine import IncumbentInputs, _fingerprint
    from evitea.incumbent.recall_repair import rescue_query_name

    assert rescue_query_name(42) == "pair_42"
    assert rescue_query_name(42, "numeric") == "42"
    assert resolve_profile(None, "S").rescue_read_names == "pair"
    assert REPRODUCTION_PROFILES["AKA017"].rescue_read_names == "pair"
    assert REPRODUCTION_PROFILES["ISO1"].rescue_read_names == "numeric"
    base = dict(sample_id="S", te_fasta=Path("t"), te_family_map=Path("m"), reference_fasta=Path("r"),
                host_index_prefix=Path("i"), bowtie2=Path("b"))
    assert "rescue_read_names" not in _fingerprint(IncumbentInputs(**base))  # old work dirs resume
    assert _fingerprint(IncumbentInputs(**base, rescue_read_names="numeric"))["rescue_read_names"] == "numeric"
    with pytest.raises(ValueError):
        IncumbentInputs(**base, rescue_read_names="P")
