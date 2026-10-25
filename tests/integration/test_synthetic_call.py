"""End-to-end `evitea call` on the synthetic data set of tools/synthetic/make_call_data.py.

Opt-in (a few minutes per case): EVITEA_RUN_SYNTHETIC_CALL=1, with Bowtie2 on PATH and pymsbwt installed.
The non-fly case uses ordinary TE record names with a two-column family map, contigs chr1..chr5 and an ``MT``
mitochondrial contig, and every value taken from the sample and the reference (none of the development values is given).
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from evitea.cli import main

ROOT = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(
    os.environ.get("EVITEA_RUN_SYNTHETIC_CALL") != "1" or shutil.which("bowtie2") is None,
    reason="set EVITEA_RUN_SYNTHETIC_CALL=1 (needs bowtie2 on PATH and pymsbwt) for the end-to-end synthetic call",
)
CASES = {
    "fly-names": ["--naming", "flybase"],
    "plain-names": ["--naming", "plain", "--contigs", "chr1,chr2,chr3,chr4,chr5", "--mito", "MT"],
}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_synthetic_call(tmp: Path, generator_args: list[str], call_args: list[str] = ()) -> Path:
    data = tmp / "data"
    subprocess.run([sys.executable, str(ROOT / "tools" / "synthetic" / "make_call_data.py"), str(data), *generator_args],
                   check=True)
    assert main(["prepare-reference", "--reference", str(data / "reference.fa"), "--repeat-annotation",
                 str(data / "repeat_annotation.bed"), "--output", str(tmp / "ref")]) == 0
    assert main(["prepare-te", "--te-library", str(data / "te_library.fa"), "--family-map", str(data / "family_map.tsv"),
                 "--build-languages", "--reference-assets", str(tmp / "ref"), "--output", str(tmp / "te")]) == 0
    assert main(["prepare-sample", "--read1", str(data / "sample_R1.fastq"), "--read2", str(data / "sample_R2.fastq"),
                 "--output", str(tmp / "sample"), "--sample-id", "synth", "--threads", "2"]) == 0
    assert main(["call", "--sample-assets", str(tmp / "sample"), "--reference-assets", str(tmp / "ref"),
                 "--te-assets", str(tmp / "te"), "--output", str(tmp / "call"), "--threads", "2", "--workers", "2",
                 *call_args]) == 0
    return tmp


@pytest.mark.parametrize("case", sorted(CASES))
def test_synthetic_insertions_are_called_with_their_families(tmp_path: Path, case: str) -> None:
    run = run_synthetic_call(tmp_path, CASES[case])
    truth = _rows(run / "data" / "truth_insertions.tsv")
    calls = _rows(run / "call" / "calls.tsv")
    found = {
        (t["contig"], int(t["position0"]))
        for t in truth
        for c in calls
        if c["chrom"] == t["contig"] and c["family"] == t["family"]
        and int(c["start"]) - 20 <= int(t["position0"]) <= int(c["end"]) + 20
    }
    assert len(found) == len(truth), sorted({(t["contig"], int(t["position0"])) for t in truth} - found)
    assert len(calls) == len(truth)
    values = json.loads((run / "call" / "run_manifest.json").read_text(encoding="utf-8"))["values"]
    parameters = values["run_parameters"]
    assert parameters["read_length"] == 125 and parameters["read_length_source"] == "prepared-sample"
    assert parameters["fragment_source"] == "estimate" and 300 < parameters["fragment_q95"] < parameters["fragment_q99"]
    contigs = [line.split("\t")[0] for line in (run / "data" / "reference.fa").read_text().splitlines()
               if line.startswith(">")]
    assert parameters["primary_contigs"] == [c[1:] for c in contigs if c not in (">chrM", ">MT")]
    assert "pymsbwt" in values["backend_versions"]["exact_query"]
