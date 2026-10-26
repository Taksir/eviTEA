"""Truth-blind checks for candidate codebook preparation."""

import csv
import hashlib
import io
import json
import random
import zipfile
from pathlib import Path

import pytest

from evitea.contracts import (
    FastaRecord,
    ReferenceInput,
    RepeatContext,
    RepeatRecord,
    TELibraryInput,
    BaseSpan0,
)
from evitea.errors import AssetCompatibilityError
from evitea.evidence import frozen_candidates
from evitea.evidence.frozen_candidates import build_frozen_candidate_codebooks


def _mer(seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(21))


def _tsv(rows, fields):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _resource_inputs(tmp_path, monkeypatch, *, overlap_candidates=False):
    alpha, internal, cross_family = _mer(10), _mer(20), _mer(30)
    catalogue_fields = [
        "query_id", "k", "query_seq", "canonical_strand_key",
        "has_terminal_provenance", "has_internal_provenance", "low_complexity",
    ]
    provenance_fields = [
        "query_id", "te_id", "family", "k", "start", "region",
        "source_orientation", "source_seq", "query_seq",
    ]
    catalogue = [
        {
            "query_id": "Q_TERM", "k": "21", "query_seq": alpha,
            "canonical_strand_key": min(alpha, alpha.translate(str.maketrans("ACGT", "TGCA"))[::-1]),
            "has_terminal_provenance": "1", "has_internal_provenance": "1", "low_complexity": "0",
        },
        {
            "query_id": "Q_INTERNAL", "k": "21", "query_seq": internal,
            "canonical_strand_key": min(internal, internal.translate(str.maketrans("ACGT", "TGCA"))[::-1]),
            "has_terminal_provenance": "0", "has_internal_provenance": "1", "low_complexity": "0",
        },
        {
            "query_id": "Q_CROSS_FAMILY", "k": "21", "query_seq": cross_family,
            "canonical_strand_key": min(cross_family, cross_family.translate(str.maketrans("ACGT", "TGCA"))[::-1]),
            "has_terminal_provenance": "0", "has_internal_provenance": "1", "low_complexity": "0",
        },
    ]
    overlap_suffix = "GTCAC" if overlap_candidates else ""
    te_alpha = alpha + overlap_suffix + "NNN" + internal + "NNN" + cross_family
    internal_start = len(alpha) + len(overlap_suffix) + 3
    cross_family_start = internal_start + 21 + 3
    if overlap_candidates:
        overlapping = te_alpha[5:26]
        catalogue.append(
            {
                "query_id": "Q_OVERLAP", "k": "21", "query_seq": overlapping,
                "canonical_strand_key": min(overlapping, overlapping.translate(str.maketrans("ACGT", "TGCA"))[::-1]),
                "has_terminal_provenance": "1", "has_internal_provenance": "0", "low_complexity": "0",
            }
        )
    provenance = [
        {
            "query_id": "Q_TERM", "te_id": "TEalpha", "family": "FamilyA", "k": "21",
            "start": "0", "region": "5p_terminal", "source_orientation": "F",
            "source_seq": alpha, "query_seq": alpha,
        },
        {
            "query_id": "Q_TERM", "te_id": "TEgamma", "family": "FamilyA", "k": "21",
            "start": "5", "region": "internal", "source_orientation": "F",
            "source_seq": alpha, "query_seq": alpha,
        },
        {
            "query_id": "Q_INTERNAL", "te_id": "TEalpha", "family": "FamilyA", "k": "21",
            "start": str(internal_start), "region": "internal", "source_orientation": "F",
            "source_seq": internal, "query_seq": internal,
        },
        {
            "query_id": "Q_CROSS_FAMILY", "te_id": "TEalpha", "family": "FamilyA", "k": "21",
            "start": str(cross_family_start), "region": "internal", "source_orientation": "F",
            "source_seq": cross_family, "query_seq": cross_family,
        },
    ]
    if overlap_candidates:
        overlapping = te_alpha[5:26]
        provenance.append(
            {
                "query_id": "Q_OVERLAP", "te_id": "TEalpha", "family": "FamilyA", "k": "21",
                "start": "5", "region": "5p_terminal", "source_orientation": "F",
                "source_seq": overlapping, "query_seq": overlapping,
            }
        )
    catalogue_bytes = _tsv(catalogue, catalogue_fields)
    provenance_bytes = _tsv(provenance, provenance_fields)
    archive_path = tmp_path / "G1A_FROZEN_RESULTS_FOR_CODEBOOK.zip"
    members = {
        "G1A_ANCHOR_CATALOGUE.tsv": catalogue_bytes,
        "G1A_ANCHOR_PROVENANCE.tsv": provenance_bytes,
    }
    manifest = {
        "data_root": "frozen-test-data",
        "files": [
            {
                "bytes": len(content),
                "file": name,
                "sha256": hashlib.sha256(content).hexdigest(),
            }
            for name, content in sorted(members.items())
        ],
        "generated_at": "frozen-test",
        "project_root": "frozen-test-project",
        "python": "test",
    }
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("G1A_RUN_MANIFEST.json", json.dumps(manifest, sort_keys=True))
        for name, content in members.items():
            archive.writestr(name, content)
    crosswalk_path = tmp_path / "fbte_to_flybase_family_FROZEN.tsv"
    crosswalk_path.write_text(
        "FBte_ID\tFlyBase_Family\nFBte0000001\tFamilyA\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(frozen_candidates, "G1A_ARCHIVE_SHA256", frozen_candidates.sha256_file(archive_path))
    monkeypatch.setattr(
        frozen_candidates,
        "FAMILY_CROSSWALK_SHA256",
        frozen_candidates.sha256_file(crosswalk_path),
    )

    te_records = (
        FastaRecord("TEalpha", "", te_alpha, 0),
        FastaRecord("TEbeta", "", "NNNNN" + cross_family, 1),
        FastaRecord("TEgamma", "", "NNNNN" + alpha, 2),
    )
    te_library = TELibraryInput(archive_path, te_records, "te-raw", "te-semantic")
    ref_record = FastaRecord("chrToy", "", alpha, 0)
    reference = ReferenceInput(
        tmp_path / "reference.fa",
        (ref_record,),
        "ref-raw",
        "ref-semantic",
        {"chrToy": len(alpha)},
    )
    repeat_context = RepeatContext(
        tmp_path / "repeats.bed",
        (RepeatRecord("chrToy", BaseSpan0(0, 21), "FBte0000001", "0", "+", 0),),
        "repeat-raw",
        "ref-semantic",
    )
    return archive_path, crosswalk_path, te_library, reference, repeat_context, alpha, internal


def _build(tmp_path, monkeypatch, *, overlap_candidates=False):
    archive, crosswalk, te_library, reference, repeats, alpha, internal = _resource_inputs(
        tmp_path, monkeypatch, overlap_candidates=overlap_candidates
    )
    built = build_frozen_candidate_codebooks(
        archive,
        crosswalk,
        te_library,
        {"TEalpha": "FamilyA", "TEbeta": "FamilyB", "TEgamma": "FamilyA"},
        reference,
        repeats,
    )
    terminal = list(csv.DictReader(io.StringIO(built.terminal_tsv.decode()), delimiter="\t"))
    internal_rows = list(csv.DictReader(io.StringIO(built.internal_tsv.decode()), delimiter="\t"))
    internal_only_rows = list(
        csv.DictReader(io.StringIO(built.internal_only_tsv.decode()), delimiter="\t")
    )
    return built, terminal, internal_rows, internal_only_rows, alpha, internal


def test_safe_languages_keep_candidate_scope_and_repeatmasker_context(tmp_path, monkeypatch):
    built, terminal, internal_rows, internal_only_rows, alpha, internal = _build(tmp_path, monkeypatch)

    assert len(terminal) == 1
    assert terminal[0]["sequence"] == min(alpha, alpha.translate(str.maketrans("ACGT", "TGCA"))[::-1])
    assert terminal[0]["reference_occurrence_count"] == "1"
    assert terminal[0]["reference_class"] == "TARGET_CONTEXT_ANNOTATED_TE_ONLY"
    assert terminal[0]["positionally_interpretable"] == "1"
    assert len(internal_rows) == 2
    assert {row["sequence"] for row in internal_rows} == {
        min(alpha, alpha.translate(str.maketrans("ACGT", "TGCA"))[::-1]),
        min(internal, internal.translate(str.maketrans("ACGT", "TGCA"))[::-1]),
    }
    assert {row["sequence"] for row in internal_only_rows} == {
        min(internal, internal.translate(str.maketrans("ACGT", "TGCA"))[::-1])
    }
    assert built.candidate_counts["terminal_candidate_sequences"] == 1
    assert built.candidate_counts["internal_candidate_sequences"] == 3
    assert all("Q_CROSS_FAMILY" not in row["probe_id"] for row in internal_rows)


def test_terminal_selection_is_deterministic_and_source_nonoverlapping(tmp_path, monkeypatch):
    _, terminal, _, _, _, _ = _build(tmp_path, monkeypatch, overlap_candidates=True)

    assert len(terminal) == 1
    assert terminal[0]["source_start"] == "5"
    assert terminal[0]["probe_rank_within_te_record"] == "1"


def test_frozen_archive_hash_drift_fails_closed(tmp_path, monkeypatch):
    archive, crosswalk, te_library, reference, repeats, *_ = _resource_inputs(tmp_path, monkeypatch)
    archive.write_bytes(archive.read_bytes() + b"drift")

    with pytest.raises(AssetCompatibilityError, match="archive SHA-256 does"):
        build_frozen_candidate_codebooks(
            archive,
            crosswalk,
            te_library,
            {"TEalpha": "FamilyA", "TEbeta": "FamilyB", "TEgamma": "FamilyA"},
            reference,
            repeats,
        )


def test_missing_frozen_catalogue_member_fails_closed(tmp_path, monkeypatch):
    archive, crosswalk, te_library, reference, repeats, *_ = _resource_inputs(tmp_path, monkeypatch)
    with zipfile.ZipFile(archive) as source:
        members = {
            name: source.read(name)
            for name in source.namelist()
            if name != "G1A_ANCHOR_CATALOGUE.tsv"
        }
    rewritten = tmp_path / "missing-member.zip"
    manifest = json.loads(members["G1A_RUN_MANIFEST.json"])
    manifest["files"] = [entry for entry in manifest["files"] if entry["file"] != "G1A_ANCHOR_CATALOGUE.tsv"]
    members["G1A_RUN_MANIFEST.json"] = json.dumps(manifest, sort_keys=True).encode()
    with zipfile.ZipFile(rewritten, "w", compression=zipfile.ZIP_DEFLATED) as target:
        for name, content in members.items():
            target.writestr(name, content)
    monkeypatch.setattr(frozen_candidates, "G1A_ARCHIVE_SHA256", frozen_candidates.sha256_file(rewritten))

    with pytest.raises(AssetCompatibilityError, match="lacks catalogue/provenance"):
        build_frozen_candidate_codebooks(
            rewritten,
            crosswalk,
            te_library,
            {"TEalpha": "FamilyA", "TEbeta": "FamilyB", "TEgamma": "FamilyA"},
            reference,
            repeats,
        )


def test_frozen_crosswalk_hash_drift_fails_closed(tmp_path, monkeypatch):
    archive, crosswalk, te_library, reference, repeats, *_ = _resource_inputs(tmp_path, monkeypatch)
    crosswalk.write_text(crosswalk.read_text(encoding="utf-8") + "FBte0000002\tFamilyB\n", encoding="utf-8")

    with pytest.raises(AssetCompatibilityError, match="crosswalk SHA-256 does"):
        build_frozen_candidate_codebooks(
            archive,
            crosswalk,
            te_library,
            {"TEalpha": "FamilyA", "TEbeta": "FamilyB", "TEgamma": "FamilyA"},
            reference,
            repeats,
        )
