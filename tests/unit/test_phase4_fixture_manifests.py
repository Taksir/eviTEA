"""Fixture-level provenance and cardinality checks for the Phase 4 batch."""

import hashlib
import json
from pathlib import Path

FIXTURES = Path(__file__).parents[1] / "fixtures" / "phase4"


def _load(fixture_id):
    root = FIXTURES / fixture_id
    metadata = json.loads((root / "fixture.json").read_text(encoding="utf-8"))
    input_doc = json.loads((root / metadata["input"]).read_text(encoding="utf-8"))
    expected_doc = json.loads((root / metadata["expected"]).read_text(encoding="utf-8"))
    return root, metadata, input_doc, expected_doc


def test_fixture_manifests_bind_synthetic_input_and_expected_bytes():
    for fixture_id in ("F-TE-001", "F-HOST-001", "F-REF-001"):
        root, metadata, _input_doc, _expected_doc = _load(fixture_id)
        assert metadata["fixture_id"] == fixture_id
        assert metadata["synthetic"] is True
        referenced_files = {metadata["input"], metadata["expected"]}
        if "projection" in metadata:
            referenced_files.add(metadata["projection"])
        assert set(metadata["sha256"]) == referenced_files
        for relative_path, expected_hash in metadata["sha256"].items():
            actual_hash = hashlib.sha256((root / relative_path).read_bytes()).hexdigest()
            assert actual_hash == expected_hash


def test_fixture_row_counts_and_pair_id_sets_are_explicit():
    _root, te_meta, te_input, te_expected = _load("F-TE-001")
    assert te_meta["pair_rows"] == len(te_input["pair_ids"]) == 7
    assert te_meta["expected_pair_ids"] == te_input["pair_ids"]
    assert sorted(map(int, te_expected)) == te_input["pair_ids"]

    _root, host_meta, host_input, host_expected = _load("F-HOST-001")
    host_pair_ids = sorted(
        rows[0]["pair_id"] for rows in host_input["alignments"].values()
    )
    assert host_meta["case_rows"] == len(host_input["alignments"]) == 4
    assert host_meta["expected_pair_ids"] == host_pair_ids
    assert len(host_expected) == host_meta["case_rows"]

    _root, ref_meta, ref_input, ref_expected = _load("F-REF-001")
    alternative_rows = sum(
        len(ref_input[key])
        for key in ("pair0_alternatives", "pair1_alternatives", "pair2_alternatives")
    )
    assert ref_meta["alternative_rows"] == alternative_rows == 5
    assert ref_meta["expected_pair_ids"] == [0, 1, 2]
    assert set(ref_expected) == {f"pair{pair_id}" for pair_id in ref_meta["expected_pair_ids"]}

