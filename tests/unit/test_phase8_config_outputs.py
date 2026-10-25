from __future__ import annotations

import json
from pathlib import Path

import pytest

from evitea.cli import (
    _call_config,
    _validate_run,
    build_parser,
    exit_code_for_error,
    main,
)
from evitea.config import load_config
from evitea.constants import ARCHITECTURE_ID
from evitea.errors import (
    ArchitectureInvariantError,
    AssetCompatibilityError,
    BackendCompatibilityError,
    InputContractError,
    OperationalResourceError,
    OutputFinalizationError,
    OutputStateError,
    PairSynchronizationError,
    ResumeStateError,
)
from evitea.callset import CallsetState, calls_bed_bytes, calls_tsv_bytes, unresolved_tsv_bytes
from evitea.outputs import write_run_summary
from evitea.shadow_layers import UnresolvedShadowEvent
from evitea.provenance import RUN_MANIFEST_SCHEMA


HEADER = (
    "event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\t"
    "left_support\tright_support\ttotal_support\tevidence_type\tscore\tfragment_ids"
)


def _row(event_id, start, end, family="roo", evidence="one_sided_left", ids="3,1,3"):
    return (
        f"{event_id}\t2L\t{start}\t{end}\t{family}\t['{family}']\t{start}\t-1\t0\t2\t0\t2\t"
        f"{evidence}\t2.0\t{ids}"
    )


def _state(*rows: str, unresolved=()) -> CallsetState:
    return CallsetState("sample", ("\n".join((HEADER, *rows)) + "\n").encode(), tuple(unresolved))


def test_toml_paths_and_explicit_cli_precedence(tmp_path: Path) -> None:
    config_path = tmp_path / "sample.toml"
    config_path.write_text(
        """
[inputs]
sample_id = "configured"
read1 = "reads/R1.fastq.gz"
read2 = "reads/R2.fastq.gz"
reference_assets = "reference"
te_assets = "te"

[runtime]
threads = 3
resume = true

[output]
directory = "results"
log_level = "DEBUG"
""",
        encoding="utf-8",
    )
    loaded = load_config(config_path)
    assert loaded.inputs.read1 == (tmp_path / "reads/R1.fastq.gz").resolve()
    assert loaded.assets.reference_assets == (tmp_path / "reference").resolve()

    args = build_parser().parse_args(
        ["call", "--config", str(config_path), "--threads", "7"]
    )
    merged = _call_config(args)
    assert merged.runtime.threads == 7
    assert merged.runtime.resume is True
    assert merged.output.log_level == "DEBUG"
    assert merged.output.directory == (tmp_path / "results").resolve()

    prepared = build_parser().parse_args(
        [
            "call",
            "--config",
            str(config_path),
            "--sample-assets",
            "prepared-sample",
        ]
    )
    prepared_config = _call_config(prepared)
    assert prepared_config.inputs.read1 is None
    assert prepared_config.inputs.read2 is None
    assert prepared_config.inputs.sample_assets == Path("prepared-sample").resolve()


def test_config_rejects_science_keys_and_cli_has_no_science_switch(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "forbidden.toml"
    config_path.write_text("[runtime]\nq99 = 12\n", encoding="utf-8")
    with pytest.raises(InputContractError, match="scientific"):
        load_config(config_path)
    with pytest.raises(SystemExit):
        build_parser().parse_args(["call", "--q99", "12"])


def test_config_rejects_non_boolean_runtime_flags(tmp_path: Path) -> None:
    config_path = tmp_path / "wrong-type.toml"
    config_path.write_text('[runtime]\nresume = "false"\n', encoding="utf-8")
    with pytest.raises(InputContractError, match="boolean"):
        load_config(config_path)


def test_outputs_are_order_independent_utf8_lf_and_canonical(tmp_path: Path) -> None:
    exact = _row("z", 100, 100, ids="5,3,1,3")
    ranged = _row("range", 20, 25, family="gypsy", evidence="two_sided_mate_cluster_rescue", ids="")
    unresolved = UnresolvedShadowEvent(
        "INTERNAL21", "k", "UNRESOLVED_NEW_CROSS_FAMILY_COMPETITION", "2L", "roo", 7, 9, ("b", "a")
    )
    left = _state(exact, ranged, unresolved=(unresolved,))
    right = _state(ranged, exact, unresolved=(unresolved,))
    assert left.native != right.native  # native order is the producers' fixed order
    assert calls_tsv_bytes(left) == calls_tsv_bytes(right)
    assert unresolved_tsv_bytes(left) == unresolved_tsv_bytes(right)
    text = calls_tsv_bytes(left)
    assert b"\r" not in text
    assert b"[1,3,5]" in text
    assert b"INTERBASE_POINT" in text and b"RANGE" in text
    assert b"INCUMBENT_MATE_RESCUE" in text
    assert b'["a","b"]' in unresolved_tsv_bytes(left)


def test_unknown_native_evidence_type_fails_closed() -> None:
    state = _state(_row("x", 5, 5, evidence="made_up"))
    with pytest.raises(InputContractError, match="unknown evidence_type"):
        calls_tsv_bytes(state)


def test_bed_proxy_boundary_rules(tmp_path: Path) -> None:
    state = _state(
        _row("range", 20, 25, evidence="two_sided_mate_cluster_rescue"),
        _row("middle", 10, 10),
        _row("end", 100, 100),
    )
    rows = {
        row.split("\t")[3]: row.split("\t")
        for row in calls_bed_bytes(state, {"2L": 100}).decode().splitlines()[1:]
    }
    assert rows["range"][1:3] == ["20", "25"]
    assert rows["middle"][1:3] == ["10", "11"]
    assert rows["end"][1:3] == ["99", "100"]
    with pytest.raises(InputContractError, match="outside"):
        calls_bed_bytes(_state(_row("far", 101, 101)), {"2L": 100})


def test_summary_truth_firewall_and_version_json(tmp_path: Path, capsys) -> None:
    with pytest.raises(InputContractError, match="benchmark"):
        write_run_summary(tmp_path / "summary.json", {"truth_score": 1})
    assert main(["version", "--json"]) == 0
    document = json.loads(capsys.readouterr().out)
    assert document["architecture_id"]
    assert set(document) >= {
        "bowtie2_version",
        "exact_query_backend",
        "native_backend",
        "numpy_version",
        "python_version",
    }


def test_validate_run_requires_every_canonical_output_hash(tmp_path: Path) -> None:
    (tmp_path / "run_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": RUN_MANIFEST_SCHEMA,
                "architecture_id": ARCHITECTURE_ID,
                "values": {"output_hashes": {}},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(AssetCompatibilityError, match="every canonical"):
        _validate_run(tmp_path)


@pytest.mark.parametrize(
    ("error", "code"),
    (
        (InputContractError("usage"), 2),
        (PairSynchronizationError("pairing"), 3),
        (AssetCompatibilityError("asset"), 4),
        (BackendCompatibilityError("backend"), 5),
        (ResumeStateError("resume"), 6),
        (OutputStateError("output state"), 7),
        (OSError("filesystem"), 7),
        (OperationalResourceError("resource"), 8),
        (ArchitectureInvariantError("invariant"), 9),
        (OutputFinalizationError("output"), 10),
    ),
)
def test_stable_exit_codes(error: BaseException, code: int) -> None:
    assert exit_code_for_error(error) == code



def test_check_install_lists_versions_and_missing_programs(monkeypatch, capsys) -> None:
    import json as _json
    import shutil as _shutil

    from evitea import cli
    from evitea.cli import main

    assert main(["check-install", "--json"]) in (0, 8)
    rows = {row["name"]: row for row in _json.loads(capsys.readouterr().out)}
    assert {"python", "pymsbwt", "numpy", "bowtie2", "bowtie2-build", "edlib", "mappy", "parasail", "makeblastdb",
            "blastn"} == set(rows)
    assert rows["numpy"]["version"] and rows["bowtie2"]["required"] and not rows["blastn"]["required"]
    monkeypatch.setattr(cli.shutil, "which", lambda name: None)
    assert _shutil.which is cli.shutil.which
    assert main(["check-install"]) == 8
    text = capsys.readouterr().out
    assert "bowtie2        MISSING (required)" in text and "blastn         missing" in text
