"""Tiny truth-blind parity checks for Phase 3 backend plumbing."""

from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import pytest

from evitea.backends.bowtie2 import (
    CANDIDATE_STATUS,
    build_bowtie2_command,
    candidate_bowtie2_2_5_5,
    iter_sam,
    map_paired_reference_reads,
)
from evitea.backends.exact_query import (
    load_read_identity_map,
    make_read_identity_map,
    prepare_sample_manifest,
    verify_prepared_sample,
    write_read_identity_map,
)
from evitea.backends.host_mapper import classify_host_mapping
from evitea.backends.msbwt import MsbwtExactQueryBackend
from evitea.backends.protocols import ReadIdentity
from evitea.contracts import HostMappingState, ReadMate
from evitea.assets.reference import prepare_reference
from evitea.evidence.projection import SampleReadIndex
from evitea.errors import AssetCompatibilityError, BackendCompatibilityError
from evitea.inputs import validate_reference_fasta

FIXTURES = Path(__file__).parents[1] / "fixtures" / "phase3"


class _TinyBwt:
    def __init__(self) -> None:
        self._dollar_ids = (10, 10, 21, 11)

    def findIndicesOfStr(self, sequence: str) -> tuple[int, int]:
        return (0, 4) if sequence == "ACGT" else (0, 0)

    def getSequenceDollarID(self, row: int) -> int:
        return self._dollar_ids[row]

    def recoverString(self, dollar_id: int) -> bytes:
        return {
            10: b"$ACGTACGT",
            11: b"$TTTTACGT",
            20: b"$CCCCCCCC",
            21: b"$GGGGACGT",
        }[dollar_id]


def test_msbwt_tiny_parity_preserves_count_identity_and_multiplicity():
    identity_map = load_read_identity_map(FIXTURES / "identity_map.json")
    backend = MsbwtExactQueryBackend(
        FIXTURES,
        identity_map,
        native_bwt=_TinyBwt(),
    )

    hits = tuple(backend.locate_reads("acgt"))

    assert backend.count_exact("ACGT") == 4
    assert [hit.backend_read_id for hit in hits] == [10, 10, 21, 11]
    assert [(hit.pair_id, hit.mate) for hit in hits] == [
        (0, ReadMate.READ1),
        (0, ReadMate.READ1),
        (1, ReadMate.READ2),
        (0, ReadMate.READ2),
    ]
    assert [hit.occurrence_ordinal0 for hit in hits] == [0, 1, 2, 3]


def test_msbwt_missing_explicit_identity_fails_closed():
    class _UnknownBwt(_TinyBwt):
        def getSequenceDollarID(self, row: int) -> int:
            return 999

    backend = MsbwtExactQueryBackend(
        FIXTURES,
        load_read_identity_map(FIXTURES / "identity_map.json"),
        native_bwt=_UnknownBwt(),
    )

    with pytest.raises(AssetCompatibilityError):
        tuple(backend.locate_reads("ACGT"))


def test_msbwt_recovers_complete_sample_index_from_bound_dollar_ids():
    backend = MsbwtExactQueryBackend(
        FIXTURES,
        load_read_identity_map(FIXTURES / "identity_map.json"),
        native_bwt=_TinyBwt(),
    )

    sample_index = backend.recover_sample_read_index("toy")

    assert sample_index.sample_id == "toy"
    assert sample_index.sequences[(0, ReadMate.READ1)] == "ACGTACGT"
    assert sample_index.sequences[(1, ReadMate.READ2)] == "GGGGACGT"


def test_msbwt_read_recovery_fails_closed_on_missing_sentinel():
    class _BadRecovery(_TinyBwt):
        def recoverString(self, dollar_id: int) -> bytes:
            return b"ACGT"

    backend = MsbwtExactQueryBackend(
        FIXTURES,
        load_read_identity_map(FIXTURES / "identity_map.json"),
        native_bwt=_BadRecovery(),
    )

    with pytest.raises(BackendCompatibilityError):
        backend.recover_read(10)


def _sam_identities() -> dict[str, ReadIdentity]:
    return {
        "toy-pair-0/1": ReadIdentity(10, 0, ReadMate.READ1),
        "toy-pair-0/2": ReadIdentity(11, 0, ReadMate.READ2),
        "toy-pair-1/1": ReadIdentity(20, 1, ReadMate.READ1),
        "toy-pair-1/2": ReadIdentity(21, 1, ReadMate.READ2),
        "toy-unmapped/1": ReadIdentity(30, 2, ReadMate.READ1),
        "toy-unmapped/2": ReadIdentity(31, 2, ReadMate.READ2),
    }


def test_bowtie2_candidate_fixture_reports_frozen_semantic_states():
    profile = candidate_bowtie2_2_5_5()
    assert profile.status == CANDIDATE_STATUS
    assert profile.max_reported == 2 and profile.min_mapq == 20

    with (FIXTURES / "bowtie2_candidate.sam").open(encoding="utf-8") as handle:
        rows = tuple(iter_sam(handle, _sam_identities()))
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.pair_id, row.mate)].append(row)
    states = {
        key: classify_host_mapping(value, min_mapq=profile.min_mapq).state
        for key, value in grouped.items()
    }

    assert states == {
        (0, ReadMate.READ1): HostMappingState.NONUNIQUE_OR_LOW,
        (0, ReadMate.READ2): HostMappingState.UNIQUE_HIGH,
        (1, ReadMate.READ1): HostMappingState.NONUNIQUE_OR_LOW,
        (1, ReadMate.READ2): HostMappingState.UNIQUE_HIGH,
        (2, ReadMate.READ1): HostMappingState.UNMAPPED,
        (2, ReadMate.READ2): HostMappingState.UNMAPPED,
    }
    assert rows[0].position0 == 0
    assert rows[1].is_secondary


def test_bowtie2_command_is_argument_vector_with_frozen_k2():
    command = build_bowtie2_command(
        "bowtie2",
        "path with spaces/reference",
        candidate_bowtie2_2_5_5(options=("--very-sensitive",)),
    )
    assert command[0] == "bowtie2"
    assert command[command.index("-x") + 1] == "path with spaces/reference"
    assert command[command.index("-k") + 1] == "2"
    assert "--very-sensitive" in command

    with pytest.raises(BackendCompatibilityError):
        build_bowtie2_command(
            "bowtie2",
            "reference",
            candidate_bowtie2_2_5_5(),
            max_reported=3,
        )


def test_paired_reference_runner_binds_shared_qname_to_both_mates(
    tmp_path: Path, monkeypatch
):
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(">chrToy\n" + "A" * 100 + "\n", encoding="ascii")
    prepared_root = tmp_path / "prepared-reference"
    (prepared_root / "bowtie2").mkdir(parents=True)
    (prepared_root / "bowtie2" / "reference.1.bt2").write_bytes(b"index")
    prepared = prepare_reference(
        validate_reference_fasta(reference_path),
        prepared_root,
        mapper_profile={"index_prefix": "bowtie2/reference"},
        index_relative_paths=("bowtie2/reference.1.bt2",),
    )

    def fake_run(command, **_kwargs):
        sam_path = Path(command[command.index("-S") + 1])
        sam_path.write_text(
            "P0\t99\tchrToy\t1\t42\t8M\t=\t11\t18\tACGTACGT\t*\n"
            "P0\t147\tchrToy\t11\t42\t8M\t=\t1\t-18\tTGCATGCA\t*\n",
            encoding="ascii",
        )
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr("evitea.backends.bowtie2.subprocess.run", fake_run)
    batch = map_paired_reference_reads(
        "bowtie2",
        SampleReadIndex(
            "toy",
            {
                (0, ReadMate.READ1): "ACGTACGT",
                (0, ReadMate.READ2): "TGCATGCA",
            },
        ),
        (0,),
        prepared,
        tmp_path / "work",
        threads=2,
    )

    assert [(row.pair_id, row.mate) for row in batch.alignments] == [
        (0, ReadMate.READ1),
        (0, ReadMate.READ2),
    ]
    assert batch.alignments[0].mate_position0 == 10
    assert batch.capped_pair_ids == ()


def test_prepared_sample_binds_sources_identity_map_and_index(tmp_path: Path):
    root = tmp_path / "sample"
    (root / "identity").mkdir(parents=True)
    (root / "index").mkdir()
    shutil.copyfile(FIXTURES / "identity_map.json", root / "identity" / "reads.json")
    shutil.copyfile(FIXTURES / "reference_index.1.bt2", root / "index" / "sample.msbwt")

    prepared = prepare_sample_manifest(
        root,
        sample_id="toy",
        pair_count=2,
        read1_path=FIXTURES / "sample_R1.fastq.txt",
        read2_path=FIXTURES / "sample_R2.fastq.txt",
        identity_map_relative_path="identity/reads.json",
        backend_profile={"backend": "pymsbwt", "version": "0.5.0"},
        index_relative_paths=["index/sample.msbwt"],
    )

    identity_map = verify_prepared_sample(prepared)
    assert len(identity_map.identities) == 4
    assert set(prepared.file_hashes) == {
        "identity/reads.json",
        "index/sample.msbwt",
    }
    assert len(prepared.read1_raw_sha256) == len(prepared.read2_raw_sha256) == 64

    (root / "index" / "sample.msbwt").write_bytes(b"changed")
    with pytest.raises(AssetCompatibilityError):
        verify_prepared_sample(prepared)


def test_prepared_sample_rejects_incomplete_pair_mate_identity(tmp_path: Path):
    root = tmp_path / "sample"
    (root / "identity").mkdir(parents=True)
    (root / "index").mkdir()
    identity_map = make_read_identity_map(
        [
            ReadIdentity(10, 0, ReadMate.READ1),
            ReadIdentity(11, 0, ReadMate.READ1),
        ]
    )
    write_read_identity_map(root / "identity" / "reads.json", identity_map)
    (root / "index" / "sample.msbwt").write_bytes(b"index")

    with pytest.raises(AssetCompatibilityError):
        prepare_sample_manifest(
            root,
            sample_id="toy",
            pair_count=1,
            read1_path=FIXTURES / "sample_R1.fastq.txt",
            read2_path=FIXTURES / "sample_R2.fastq.txt",
            identity_map_relative_path="identity/reads.json",
            backend_profile={"backend": "pymsbwt", "version": "0.5.0"},
            index_relative_paths=["index/sample.msbwt"],
        )
