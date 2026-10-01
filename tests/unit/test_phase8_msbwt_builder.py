from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from evitea.backends.exact_query import verify_prepared_sample
from evitea.backends.msbwt import build_prepared_msbwt_sample
from evitea.contracts import ReadMate
from evitea.errors import AssetCompatibilityError, BackendCompatibilityError


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


def _fake_builder_with_about(order: tuple[tuple[int, int], ...]):
    def builder(fastqs, output_dir, threads, are_uniform, _logger):
        assert len(fastqs) == 2
        assert threads == 2
        assert are_uniform is True
        root = Path(output_dir)
        dtype = np.dtype([("file_id", "<u1"), ("read_id", "<u8")])
        about = np.empty(len(order), dtype=dtype)
        about["file_id"] = [row[0] for row in order]
        about["read_id"] = [row[1] for row in order]
        np.save(root / "about.npy", about)
        np.save(root / "msbwt.npy", np.asarray([0, 1, 2], dtype=np.uint8))

    return builder


def _fake_initializer(index_dir: Path):
    np.save(index_dir / "fmIndex.npy", np.asarray([0, 1], dtype=np.uint64))
    np.save(index_dir / "totalCounts.npy", np.asarray([1, 2], dtype=np.uint64))
    return object()


def test_uniform_builder_freezes_pair_and_mate_identity(tmp_path: Path) -> None:
    read1 = tmp_path / "R1.fastq"
    read2 = tmp_path / "R2.fastq"
    _write_fastq(read1, 1, ("AAAA", "CCCC"))
    _write_fastq(read2, 2, ("TTTT", "GGGG"))
    module = tmp_path / "builder.pyd"
    module.write_bytes(b"fixture-builder")
    output = tmp_path / "prepared"
    asset = build_prepared_msbwt_sample(
        read1,
        read2,
        output,
        sample_id="sample",
        threads=2,
        builder=_fake_builder_with_about(((0, 0), (1, 1), (0, 1), (1, 0))),
        initializer=_fake_initializer,
        backend_version="0.3.0",
        builder_module_path=module,
    )

    identity_map = verify_prepared_sample(asset)
    assert identity_map.resolve(0).pair_id == 0
    assert identity_map.resolve(0).mate is ReadMate.READ1
    assert identity_map.resolve(1).pair_id == 1
    assert identity_map.resolve(1).mate is ReadMate.READ2
    assert identity_map.resolve(2).pair_id == 1
    assert identity_map.resolve(2).mate is ReadMate.READ1
    assert identity_map.resolve(3).pair_id == 0
    assert identity_map.resolve(3).mate is ReadMate.READ2
    assert asset.backend_profile["identity_source"] == "index/about.npy"
    assert asset.backend_profile["input_file_order"] == ["read1", "read2"]
    assert not (output / ".builder-input").exists()


def test_nonuniform_reads_fail_before_builder(tmp_path: Path) -> None:
    read1 = tmp_path / "R1.fastq"
    read2 = tmp_path / "R2.fastq"
    _write_fastq(read1, 1, ("AAAA", "CCCCC"))
    _write_fastq(read2, 2, ("TTTT", "GGGG"))
    called = False

    def builder(*_args):
        nonlocal called
        called = True

    with pytest.raises(BackendCompatibilityError, match="uniform read length"):
        build_prepared_msbwt_sample(
            read1,
            read2,
            tmp_path / "prepared",
            sample_id="sample",
            threads=1,
            builder=builder,
            initializer=_fake_initializer,
        )
    assert called is False


def test_about_map_must_cover_both_mates_for_every_pair(tmp_path: Path) -> None:
    read1 = tmp_path / "R1.fastq"
    read2 = tmp_path / "R2.fastq"
    _write_fastq(read1, 1, ("AAAA", "CCCC"))
    _write_fastq(read2, 2, ("TTTT", "GGGG"))
    with pytest.raises(AssetCompatibilityError, match="every physical pair"):
        build_prepared_msbwt_sample(
            read1,
            read2,
            tmp_path / "prepared",
            sample_id="sample",
            threads=2,
            builder=_fake_builder_with_about(((0, 0), (1, 0), (0, 1))),
            initializer=_fake_initializer,
        )

