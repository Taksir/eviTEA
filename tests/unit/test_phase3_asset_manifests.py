"""Synthetic Phase 3 checks for prepared reference and TE asset manifests."""

import json
import shutil
from pathlib import Path

import pytest

from evitea.assets.context import canonical_json_bytes
from evitea.assets.reference import (
    load_prepared_reference,
    prepare_reference,
    verify_reference_asset,
)
from evitea.assets.te import (
    load_prepared_te_assets,
    prepare_te_assets,
    verify_te_assets,
)
from evitea.errors import AssetCompatibilityError
from evitea.inputs import (
    load_family_map,
    load_repeat_annotation,
    validate_reference_fasta,
    validate_te_fasta,
)

FIXTURES = Path(__file__).parents[1] / "fixtures" / "phase3"
CANDIDATE_UNAPPROVED_PROFILE = json.loads(
    (FIXTURES / "bowtie2_candidate_profile.json").read_text(encoding="utf-8")
)


def _prepare_reference(tmp_path: Path, *, include_repeat: bool = True):
    root = tmp_path / "r"
    root.mkdir()
    index_paths = ["index/reference.1.bt2", "index/reference.2.bt2"]
    (root / "index").mkdir()
    for relative, fixture in zip(
        index_paths,
        ("reference_index.1.bt2", "reference_index.2.bt2"),
        strict=True,
    ):
        shutil.copyfile(FIXTURES / fixture, root / relative)

    reference_input = validate_reference_fasta(FIXTURES / "reference.fa")
    repeat_context = (
        load_repeat_annotation(FIXTURES / "repeats.gff3", reference_input)
        if include_repeat
        else None
    )
    prepared = prepare_reference(
        reference_input,
        root,
        mapper_profile=CANDIDATE_UNAPPROVED_PROFILE,
        index_relative_paths=index_paths,
        repeat_context=repeat_context,
    )
    return prepared


def _prepare_te_assets(tmp_path: Path):
    reference = _prepare_reference(tmp_path)
    root = tmp_path / "t"
    root.mkdir()
    (root / "index").mkdir()
    shutil.copyfile(FIXTURES / "te_assets.fa", root / "index" / "te_assets.fa")
    shutil.copyfile(FIXTURES / "te_assets.tsv", root / "index" / "te_assets.tsv")
    te_library = validate_te_fasta(FIXTURES / "te.fa")
    family_map_path = FIXTURES / "te_families.tsv"
    prepared = prepare_te_assets(
        te_library,
        load_family_map(family_map_path),
        family_map_path,
        reference,
        root,
        asset_relative_paths=["index/te_assets.fa", "index/te_assets.tsv"],
    )
    return prepared, reference


def test_canon():
    left = {"z": [3, {"b": 2, "a": 1}], "a": True}
    right = {"a": True, "z": [3, {"a": 1, "b": 2}]}

    encoded = canonical_json_bytes(left)

    assert encoded == canonical_json_bytes(right)
    assert encoded == b'{"a":true,"z":[3,{"a":1,"b":2}]}\n'


def test_ref_manifest(tmp_path):
    prepared = _prepare_reference(tmp_path)

    assert prepared.mapper_profile["version"] == "2.5.5"
    assert prepared.mapper_profile["approved"] is False
    assert prepared.mapper_profile["options"]["k"] == 2
    assert prepared.contigs == (("chrToy", 16), ("chrToy2", 16))
    assert prepared.repeat_raw_sha256
    assert set(prepared.file_hashes) == {
        "index/reference.1.bt2",
        "index/reference.2.bt2",
        "sources/reference.fa",
        "sources/repeat_annotation.gff3",
    }
    verify_reference_asset(prepared)

    reloaded = load_prepared_reference(prepared.root)
    assert reloaded.identity_sha256 == prepared.identity_sha256
    manifest_path = prepared.root / "reference_manifest.json"
    assert manifest_path.read_bytes() == canonical_json_bytes(json.loads(manifest_path.read_text()))


@pytest.mark.parametrize(
    "relative_path",
    [
        pytest.param("sources/reference.fa", id="src"),
        pytest.param("sources/repeat_annotation.gff3", id="rep"),
        pytest.param("index/reference.1.bt2", id="idx1"),
        pytest.param("index/reference.2.bt2", id="idx2"),
    ],
)
def test_ref_change(tmp_path, relative_path):
    prepared = _prepare_reference(tmp_path)
    bound_file = prepared.root / relative_path
    bound_file.write_bytes(bound_file.read_bytes() + b"\nchanged\n")

    with pytest.raises(AssetCompatibilityError):
        verify_reference_asset(prepared)


def test_te_manifest(tmp_path):
    prepared, reference = _prepare_te_assets(tmp_path)

    assert prepared.reference_identity_sha256 == reference.identity_sha256
    assert dict(prepared.family_assignments) == {"TEalpha": "FamilyA", "TEbeta": "FamilyB"}
    assert set(prepared.file_hashes) == {
        "index/te_assets.fa",
        "index/te_assets.tsv",
        "sources/family_map.tsv",
        "sources/te_library.fa",
    }
    verify_te_assets(prepared, reference)

    reloaded = load_prepared_te_assets(prepared.root)
    assert reloaded.identity_sha256 == prepared.identity_sha256
    document = reloaded.to_document()
    assert "safe_language" not in document
    assert "phase4" not in " ".join(document.keys()).casefold()


@pytest.mark.parametrize(
    "relative_path",
    [
        pytest.param("sources/te_library.fa", id="src"),
        pytest.param("sources/family_map.tsv", id="map"),
        pytest.param("index/te_assets.fa", id="fa"),
        pytest.param("index/te_assets.tsv", id="tsv"),
    ],
)
def test_te_change(tmp_path, relative_path):
    prepared, reference = _prepare_te_assets(tmp_path)
    bound_file = prepared.root / relative_path
    bound_file.write_bytes(bound_file.read_bytes() + b"\nchanged\n")

    with pytest.raises(AssetCompatibilityError):
        verify_te_assets(prepared, reference)

