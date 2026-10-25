from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from evitea.engine import EngineInputs
from evitea.errors import InputContractError
from evitea.incumbent import coverage_panel, recall_repair
from evitea.sample_geometry import (
    DEFAULT_READ_LENGTH,
    read_length_from_index,
    resolve_read_length,
    validate_read_length,
)


def test_window_geometry_follows_the_read_length() -> None:
    # 125 bp (the default) and the same call with the default left out agree
    assert coverage_panel.coverage_interval(300, 100, 21) == coverage_panel.coverage_interval(300, 100, 21, 125)
    assert coverage_panel.coverage_interval(300, 100, 21, 125) == (0, 100)
    # 100 bp reads: a read must start at or after 21 to hold the word at 100..121, and 201 windows exist
    assert coverage_panel.coverage_interval(300, 100, 21, 100) == (21, 100)
    assert coverage_panel.coverage_interval(300, 250, 21, 100) == (171, 200)
    # 150 bp reads on a 140 bp record: no window fits
    assert coverage_panel.coverage_interval(140, 10, 21, 150) is None
    assert len(coverage_panel._window_mask(300, [(100, 21)], 150)) == 151
    assert coverage_panel._window_mask(300, [(100, 21)], 100).sum() == 80
    assert recall_repair._window_interval(100, 21, 201, 100) == (21, 100)
    assert recall_repair._window_interval(100, 21, 176) == (0, 100)


def test_h1_loci_centre_on_the_read_length() -> None:
    records = {"TE1": "ACGT" * 100}
    family = {"TE1": "fam"}
    short = []
    loci_125 = {(l.orientation, l.start) for l in coverage_panel.hard_gap_h1_loci(records, family, {}, short)}
    loci_100 = {(l.orientation, l.start) for l in coverage_panel.hard_gap_h1_loci(records, family, {}, short, 100)}
    # the first internal window starts at 0; its 21-mer is centred in the read: (125 - 21) // 2 = 52 vs (100 - 21) // 2 = 39
    assert ("F", 52) in loci_125 and ("F", 52) not in loci_100
    assert ("F", 39) in loci_100


def _counts(path: Path, values) -> Path:
    path.mkdir()
    np.save(path / "totalCounts.npy", np.asarray(values, dtype=np.uint64))
    return path


def test_read_length_from_the_index_symbol_counts(tmp_path: Path) -> None:
    # 10 reads of 100 bp: 10 end markers and 1,000 bases
    assert read_length_from_index(_counts(tmp_path / "a", [10, 246, 239, 262, 0, 253])) == 100
    with pytest.raises(InputContractError, match="one length"):
        read_length_from_index(_counts(tmp_path / "b", [10, 246, 239, 262, 0, 254]))


def test_read_length_resolution_order(tmp_path: Path) -> None:
    index = _counts(tmp_path / "idx", [4, 150, 150, 150, 0, 150])  # 4 reads of 150 bp
    assert resolve_read_length(None, None, None).length == DEFAULT_READ_LENGTH
    assert resolve_read_length(None, None, index) == resolve_read_length(None, {}, index)
    assert resolve_read_length(None, None, index).length == 150
    assert resolve_read_length(None, None, index).source == "msbwt-index"
    assert resolve_read_length(None, {"uniform_length": 151}, index).length == 151
    assert resolve_read_length(None, {"uniform_length": 151}, index).source == "prepared-sample"
    assert resolve_read_length(100, {"uniform_length": 151}, index).source == "option"
    with pytest.raises(InputContractError):
        validate_read_length(20)
    with pytest.raises(InputContractError):
        validate_read_length(True)


def test_engine_identity_records_only_a_non_default_read_length(tmp_path: Path) -> None:
    paths = {name: tmp_path / name for name in (
        "msbwt_index", "te_fasta", "te_family_map", "reference_fasta", "repeatmasker_bed",
        "host_index_prefix", "terminal_codebook", "internal_codebook", "bowtie2",
    )}
    default = EngineInputs(sample_id="s", **paths)
    assert "read_length" not in default.identity()
    longer = EngineInputs(sample_id="s", read_length=151, read_length_source="prepared-sample", **paths)
    assert longer.identity()["read_length"] == 151
    assert longer.run_parameters() == {"read_length": 151, "read_length_source": "prepared-sample"}
