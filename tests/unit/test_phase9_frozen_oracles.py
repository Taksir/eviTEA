"""Tests for the validation-only verifier of fixed reference outputs (row counts and SHA-256)."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path


_VERIFIER_PATH = (
    Path(__file__).parents[2]
    / "validation"
    / "frozen_regression"
    / "verify_oracles.py"
)
_VERIFIER_SPEC = importlib.util.spec_from_file_location(
    "phase9_frozen_oracle_verifier", _VERIFIER_PATH
)
assert _VERIFIER_SPEC is not None and _VERIFIER_SPEC.loader is not None
verifier = importlib.util.module_from_spec(_VERIFIER_SPEC)
_VERIFIER_SPEC.loader.exec_module(verifier)


def _check(relative_file: str, payload: bytes, *, sha256: str | None = None):
    return {
        "layer": "test-layer",
        "dataset": "AKA017",
        "route_root": "route_b_v2",
        "file": relative_file,
        "data_rows": 2,
        "sha256": sha256 or hashlib.sha256(payload).hexdigest(),
    }


def _write(root: Path, relative_file: str, payload: bytes) -> None:
    target = root / relative_file
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)


def test_verifier_passes_for_matching_file(tmp_path):
    root_v2 = tmp_path / "route-v2"
    root_v3 = tmp_path / "route-v3"
    payload = b"header\nrow-1\nrow-2\n"
    _write(root_v2, "frozen/calls.tsv", payload)

    report = verifier.verify_oracles(
        root_v2, root_v3, checks=[_check("frozen/calls.tsv", payload)]
    )

    assert report["status"] == "PASS"
    assert report["scope"] == "FROZEN_ORACLE_INTEGRITY_ONLY"
    assert report["method_equivalence"] == "NOT_EVALUATED"
    assert len(report["checks"]) == 1
    assert report["checks"][0]["actual_data_rows"] == 2
    assert report["checks"][0]["actual_sha256"] == hashlib.sha256(payload).hexdigest()


def test_verifier_stops_at_earliest_hash_mismatch(tmp_path):
    root_v2 = tmp_path / "route-v2"
    payload = b"header\nrow-1\nrow-2\n"
    _write(root_v2, "first.tsv", payload)
    checks = [
        _check("first.tsv", payload, sha256="0" * 64),
        _check("later-missing.tsv", payload),
    ]

    report = verifier.verify_oracles(root_v2, tmp_path / "route-v3", checks=checks)

    assert report["status"] == "FAIL"
    assert len(report["checks"]) == 1
    assert report["checks"][0]["file"] == "first.tsv"
    assert report["checks"][0]["reason"] == "sha256_mismatch"


def test_verifier_reports_missing_file(tmp_path):
    report = verifier.verify_oracles(
        tmp_path / "route-v2",
        tmp_path / "route-v3",
        checks=[_check("missing.tsv", b"header\nrow-1\nrow-2\n")],
    )

    assert report["status"] == "FAIL"
    assert report["checks"][0]["reason"] == "missing_file"
    assert report["checks"][0]["actual_data_rows"] is None


def test_cli_output_is_deterministic_and_saved(tmp_path, monkeypatch, capsys):
    root_v2 = tmp_path / "route-v2"
    root_v3 = tmp_path / "route-v3"
    payload = b"header\nrow-1\nrow-2\n"
    _write(root_v2, "frozen/calls.tsv", payload)
    checks = [_check("frozen/calls.tsv", payload)]
    monkeypatch.setattr(verifier, "_load_oracles", lambda: checks)
    output_path = tmp_path / "report.json"
    arguments = [
        "--route-b-v2",
        str(root_v2),
        "--route-b-v3",
        str(root_v3),
        "--output",
        str(output_path),
    ]

    assert verifier.main(arguments) == 0
    first_stdout = capsys.readouterr().out
    first_file = output_path.read_text(encoding="utf-8")
    assert verifier.main(arguments) == 0
    second_stdout = capsys.readouterr().out

    assert first_stdout == second_stdout == first_file
    assert json.loads(first_file)["status"] == "PASS"
