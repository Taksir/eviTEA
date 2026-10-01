"""Phase 2 checks for raw FASTQ and user-supplied reference assets."""

import gzip
import shutil
from pathlib import Path

import pytest

from evitea.contracts import ReadMate
from evitea.errors import InputContractError, PairSynchronizationError
from evitea.inputs import (
    iter_paired_fastq,
    load_family_map,
    load_repeat_annotation,
    validate_reference_fasta,
    validate_te_fasta,
)

FIXTURES = Path(__file__).parents[1] / "fixtures"


def _write_fastq(path: Path, records):
    path.write_text(
        "".join(
            f"@{header}\n{sequence}\n+\n{quality}\n"
            for header, sequence, quality in records
        ),
        encoding="ascii",
    )
    return path


def _write_fasta(path: Path, text: str):
    path.write_text(text, encoding="ascii")
    return path


def _copy_six_pair_fixture(tmp_path):
    read1 = tmp_path / "phase2_pairs_R1.fastq"
    read2 = tmp_path / "phase2_pairs_R2.fastq"
    shutil.copyfile(FIXTURES / "phase2_pairs_R1.fastq.txt", read1)
    shutil.copyfile(FIXTURES / "phase2_pairs_R2.fastq.txt", read2)
    return read1, read2


def test_tiny_fastq_fixture_covers_header_forms_and_uppercases_sequences(tmp_path):
    read1, read2 = _copy_six_pair_fixture(tmp_path)
    pairs = list(iter_paired_fastq(read1, read2))

    assert len(pairs) == 6
    assert [pair.pair_name for pair in pairs] == [
        "pair0",
        "repeat",
        "casava",
        "pair3",
        "pair4",
        "repeat",
    ]
    assert pairs[0].read1.sequence == "ACGT"
    assert pairs[0].read1.mate is ReadMate.READ1
    assert pairs[0].read2.mate is ReadMate.READ2


def test_fastq_mismatched_names_fail_without_resynchronization(tmp_path):
    read1 = _write_fastq(tmp_path / "r1.fastq", [("alpha/1", "ACGT", "IIII")])
    read2 = _write_fastq(tmp_path / "r2.fastq", [("beta/2", "TGCA", "IIII")])

    with pytest.raises(PairSynchronizationError):
        list(iter_paired_fastq(read1, read2))


def test_fastq_wrong_mate_marker_fails_closed(tmp_path):
    read1 = _write_fastq(tmp_path / "r1.fastq", [("alpha/2", "ACGT", "IIII")])
    read2 = _write_fastq(tmp_path / "r2.fastq", [("alpha/1", "TGCA", "IIII")])

    with pytest.raises(InputContractError):
        list(iter_paired_fastq(read1, read2))


def test_fastq_casava_mate_marker_must_match_stream(tmp_path):
    read1 = _write_fastq(tmp_path / "r1.fastq", [("alpha 1:N:0:1", "ACGT", "IIII")])
    read2 = _write_fastq(tmp_path / "r2.fastq", [("alpha 1:N:0:1", "TGCA", "IIII")])

    with pytest.raises(InputContractError):
        list(iter_paired_fastq(read1, read2))


def test_fastq_unequal_sequence_and_quality_lengths_fail_closed(tmp_path):
    read1 = _write_fastq(tmp_path / "r1.fastq", [("alpha/1", "ACGT", "III")])
    read2 = _write_fastq(tmp_path / "r2.fastq", [("alpha/2", "TGCA", "IIII")])

    with pytest.raises(InputContractError):
        list(iter_paired_fastq(read1, read2))


def test_fastq_short_r2_stream_fails_closed(tmp_path):
    read1 = _write_fastq(
        tmp_path / "r1.fastq",
        [("alpha/1", "ACGT", "IIII"), ("beta/1", "TGCA", "IIII")],
    )
    read2 = _write_fastq(tmp_path / "r2.fastq", [("alpha/2", "TGCA", "IIII")])

    with pytest.raises(PairSynchronizationError):
        list(iter_paired_fastq(read1, read2))


def test_truncated_fastq_record_fails_closed(tmp_path):
    read1 = tmp_path / "r1.fastq"
    read1.write_text("@alpha/1\nACGT\n+\n", encoding="ascii")
    read2 = _write_fastq(tmp_path / "r2.fastq", [("alpha/2", "TGCA", "IIII")])

    with pytest.raises(InputContractError):
        list(iter_paired_fastq(read1, read2))


def test_fastq_plain_and_gzip_streams_are_both_accepted(tmp_path):
    plain = _write_fastq(tmp_path / "r1.fastq", [("alpha/1", "ACGT", "IIII")])
    zipped = tmp_path / "r2.fq.gz"
    with gzip.open(zipped, "wt", encoding="ascii") as stream:
        stream.write("@alpha/2\nTGCA\n+\nIIII\n")

    pair = next(iter_paired_fastq(plain, zipped))
    assert pair.pair_name == "alpha"


def test_reference_fasta_uses_first_tokens_order_and_uppercase_sequences(tmp_path):
    fasta = _write_fasta(
        tmp_path / "reference.fa",
        ">chrA primary contig\nacgt\n>chrB secondary contig\nTTaa\n",
    )

    reference = validate_reference_fasta(fasta)

    assert [(record.record_id, record.sequence, record.rank) for record in reference.records] == [
        ("chrA", "ACGT", 0),
        ("chrB", "TTAA", 1),
    ]
    assert dict(reference.contig_lengths) == {"chrA": 4, "chrB": 4}


def test_reference_fasta_rejects_duplicate_first_token_ids(tmp_path):
    fasta = _write_fasta(
        tmp_path / "duplicate.fa",
        ">chrA first\nACGT\n>chrA second\nTGCA\n",
    )

    with pytest.raises(InputContractError):
        validate_reference_fasta(fasta)


def test_te_fasta_and_family_map_preserve_only_frozen_aliases(tmp_path):
    fasta = _write_fasta(
        tmp_path / "te.fa",
        ">te_roo description\nACGT\n>te_rt1 description\ntgca\n",
    )
    te_library = validate_te_fasta(fasta)
    assert [record.record_id for record in te_library.records] == ["te_roo", "te_rt1"]
    assert [record.sequence for record in te_library.records] == ["ACGT", "TGCA"]

    family_path = tmp_path / "families.tsv"
    family_path.write_text("te_record_id\tfamily\nte_roo\trooA\nte_rt1\tRt1a\n", encoding="utf-8")
    family_map = load_family_map(family_path)
    assert dict(family_map) == {"te_roo": "roo", "te_rt1": "Rt1b"}


def test_family_map_header_may_follow_comments(tmp_path):
    family_path = tmp_path / "families.tsv"
    family_path.write_text(
        "# generated fixture\nte_record_id\tfamily\nte_roo\trooA\n",
        encoding="ascii",
    )

    assert dict(load_family_map(family_path)) == {"te_roo": "roo"}


def test_te_fasta_rejects_duplicate_ids_and_family_map_rejects_conflicts(tmp_path):
    fasta = _write_fasta(
        tmp_path / "duplicate_te.fa",
        ">te1 first\nACGT\n>te1 second\nTGCA\n",
    )
    with pytest.raises(InputContractError):
        validate_te_fasta(fasta)

    family_path = tmp_path / "conflict.tsv"
    family_path.write_text("te_record_id\tfamily\nte1\troo\nte1\tRt1b\n", encoding="utf-8")
    with pytest.raises(InputContractError):
        load_family_map(family_path)


def test_repeat_annotation_is_bound_to_reference_and_checks_bed6_bounds(tmp_path):
    reference_path = _write_fasta(tmp_path / "reference.fa", ">chr1\nACGTACGTAA\n")
    reference = validate_reference_fasta(reference_path)
    bed = tmp_path / "repeats.bed"
    bed.write_text("chr1\t2\t6\trepeatA\t900\t+\n", encoding="utf-8")

    repeats = load_repeat_annotation(bed, reference)

    assert len(repeats.records) == 1
    assert repeats.reference_semantic_sha256 == reference.semantic_sha256
    row = repeats.records[0]
    assert (row.chrom, row.span0.start0, row.span0.end0) == ("chr1", 2, 6)
    assert (row.name, row.score, row.strand, row.source_row0) == ("repeatA", "900", "+", 0)

    bed.write_text("chr1\t8\t11\trepeatB\t0\t-\n", encoding="utf-8")
    with pytest.raises(InputContractError):
        load_repeat_annotation(bed, reference)
