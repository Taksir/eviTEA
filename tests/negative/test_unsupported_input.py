"""Unsupported input stops at the start of a run with a plain message (one test per check)."""

from __future__ import annotations

from pathlib import Path

import pytest

from evitea.cli import main


def _fastq(path: Path, mate: int, reads: list[str]) -> Path:
    path.write_text("".join(f"@pair{i}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads)),
                    encoding="ascii")
    return path


def test_single_end_input_is_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    reads = _fastq(tmp_path / "reads.fastq", 1, ["ACGT" * 25] * 3)
    code = main(["prepare-sample", "--read1", str(reads), "--read2", str(reads), "--output", str(tmp_path / "out")])
    assert code == 3
    assert "same file" in capsys.readouterr().err
    other = tmp_path / "other.fastq"
    other.write_text("@different/2\n" + "ACGT" * 25 + "\n+\n" + "I" * 100 + "\n", encoding="ascii")
    assert main(["prepare-sample", "--read1", str(reads), "--read2", str(other), "--output", str(tmp_path / "o2")]) == 3
    assert "do not pair up" in capsys.readouterr().err


def test_reads_of_several_lengths_are_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    read1 = _fastq(tmp_path / "r1.fastq", 1, ["ACGT" * 25, "ACGT" * 30])
    read2 = _fastq(tmp_path / "r2.fastq", 2, ["TTGC" * 25, "TTGC" * 30])
    code = main(["prepare-sample", "--read1", str(read1), "--read2", str(read2), "--output", str(tmp_path / "out")])
    assert code == 3
    message = capsys.readouterr().err
    assert "more than one length (100, 120 bp)" in message and "one length" in message
    assert not (tmp_path / "out").exists()


def test_contig_names_that_the_reference_lacks_are_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr2L\n" + "ACGT" * 100 + "\n", encoding="ascii")
    annotation = tmp_path / "repeats.bed"
    annotation.write_text("2L\t10\t50\tFBte0000001\t0\t+\n", encoding="utf-8")
    code = main(["prepare-reference", "--reference", str(reference), "--repeat-annotation", str(annotation),
                 "--output", str(tmp_path / "out")])
    assert code == 2
    assert "the reference names it 'chr2L'" in capsys.readouterr().err


def test_te_records_missing_from_the_family_map_are_refused(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    library = tmp_path / "te.fa"
    library.write_text(">B1_Mm\n" + "ACGT" * 50 + "\n>Ty1\n" + "TTGCA" * 40 + "\n", encoding="ascii")
    flybase = tmp_path / "map.tsv"
    flybase.write_text("FBte_ID\tFlyBase_Family\nFBte0000001\troo\n", encoding="utf-8")
    code = main(["prepare-te", "--te-library", str(library), "--family-map", str(flybase), "--reference-assets",
                 str(tmp_path / "missing"), "--output", str(tmp_path / "out")])
    assert code == 2
    assert "no TE library record is named in the family map" in capsys.readouterr().err
