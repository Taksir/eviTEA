from __future__ import annotations

import os
from pathlib import Path

import pytest

from evitea.backends.exact_query import load_prepared_sample, verify_prepared_sample
from evitea.backends.msbwt import MsbwtExactQueryBackend
from evitea.cli import main
from evitea.contracts import ReadMate

pytestmark = pytest.mark.skipif(
    os.environ.get("EVITEA_RUN_MSBWT_LIVE") != "1",
    reason="set EVITEA_RUN_MSBWT_LIVE=1 for compiled msBWT parity",
)


def _write_fastq(path: Path, mate: int, sequences: tuple[str, ...]) -> None:
    records = []
    for pair_id, sequence in enumerate(sequences):
        records.extend(
            (
                f"@pair{pair_id}/{mate}",
                sequence,
                "+",
                "I" * len(sequence),
            )
        )
    path.write_text("\n".join(records) + "\n", encoding="ascii", newline="\n")


def test_compiled_builder_preserves_duplicate_pair_and_mate_identity(
    tmp_path: Path,
) -> None:
    read1 = tmp_path / "R1.fastq"
    read2 = tmp_path / "R2.fastq"
    _write_fastq(read1, 1, ("ACGT", "TTTT", "ACGT"))
    _write_fastq(read2, 2, ("CCCC", "ACGT", "GGGG"))
    output = tmp_path / "prepared"
    assert (
        main(
            [
                "prepare-sample",
                "--read1",
                str(read1),
                "--read2",
                str(read2),
                "--output",
                str(output),
                "--sample-id",
                "fixture",
                "--threads",
                "1",
            ]
        )
        == 0
    )
    asset = load_prepared_sample(output)
    identity_map = verify_prepared_sample(asset)
    assert set(asset.file_hashes) == {
        "read_identity_map.json",
        "index/about.npy",
        "index/fmIndex.npy",
        "index/msbwt.npy",
        "index/totalCounts.npy",
    }
    backend = MsbwtExactQueryBackend(
        asset.root / str(asset.backend_profile["index_directory"]), identity_map
    )

    hits = tuple(backend.locate_reads("ACGT"))
    assert backend.count_exact("ACGT") == 3
    assert {(hit.pair_id, hit.mate) for hit in hits} == {
        (0, ReadMate.READ1),
        (2, ReadMate.READ1),
        (1, ReadMate.READ2),
    }
    sample_index = backend.recover_sample_read_index("fixture")
    assert sample_index.sequences[(0, ReadMate.READ1)] == "ACGT"
    assert sample_index.sequences[(1, ReadMate.READ2)] == "ACGT"
    reloaded = MsbwtExactQueryBackend(
        asset.root / str(asset.backend_profile["index_directory"]), identity_map
    )
    assert reloaded.count_exact("ACGT") == 3

