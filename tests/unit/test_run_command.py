from __future__ import annotations

from pathlib import Path

import pytest

from evitea.cli import build_parser
from evitea.errors import InputContractError
from evitea.run_command import run


def _args(tmp_path: Path, *extra: str):
    return build_parser().parse_args([
        "run", "--reference", "ref.fa", "--repeat-annotation", "rm.bed", "--read1", "s_R1.fastq", "--read2", "s_R2.fastq",
        "--sample-id", "s", "--output", str(tmp_path / "out"), "--threads", "2", *extra,
    ])


def test_library_run_chains_the_separate_commands(tmp_path: Path) -> None:
    seen: list[list[str]] = []
    args = _args(tmp_path, "--te-library", "te.fa", "--family-map", "map.tsv", "--library-free", "--contigs", "chr1",
                 "--read-length", "150")
    assert run(args, lambda words: seen.append(list(words)) or 0, lambda m: None) == 0
    assert [w[0] for w in seen] == ["prepare-reference", "prepare-te", "prepare-sample", "call", "library-free"]
    out = tmp_path / "out"
    call = seen[3]
    assert call[call.index("--te-assets") + 1] == str(out / "te") and call[call.index("--read-length") + 1] == "150"
    assert seen[4][seen[4].index("--contigs") + 1] == "chr1"
    # an interrupted run continues: done steps are skipped
    seen.clear()
    (out / "steps" / "call.done").unlink()
    assert run(args, lambda words: seen.append(list(words)) or 0, lambda m: None) == 0
    assert [w[0] for w in seen] == ["call"]
    assert (out / "run_steps.json").is_file()


def test_no_library_run_follows_the_slurm_chain(tmp_path: Path) -> None:
    seen: list[list[str]] = []
    args = _args(tmp_path, "--no-library")
    assert run(args, lambda words: seen.append(list(words)) or 0, lambda m: None) == 0
    assert [w[0] for w in seen] == ["prepare-reference", "prepare-sample", "build-repeat-map", "library-free",
                                    "build-family-records", "prepare-te", "call", "add-record-calls"]
    out = tmp_path / "out"
    stub = (out / "nolib_stub" / "calls.native.tsv").read_text().splitlines()
    assert stub == ["event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\tleft_support\t"
                    "right_support\ttotal_support\tevidence_type\tscore\tfragment_ids"]
    detect = seen[3]
    assert detect[detect.index("--run-dir") + 1] == str(out / "nolib_stub")
    assert "--rule-c" in detect and detect[detect.index("--repeat-map") + 1].endswith("repeat_map.bed")
    records = seen[4]
    assert "--second-segments" not in records and records[records.index("--flag-rule") + 1] == "kept"


def test_a_failed_step_stops_the_run(tmp_path: Path) -> None:
    args = _args(tmp_path, "--te-library", "te.fa", "--family-map", "map.tsv")
    assert run(args, lambda words: 4 if words[0] == "prepare-te" else 0, lambda m: None) == 4
    assert not (tmp_path / "out" / "steps" / "te.done").exists()
    with pytest.raises(InputContractError):
        run(_args(tmp_path, "--te-library", "te.fa"), lambda words: 0, lambda m: None)
