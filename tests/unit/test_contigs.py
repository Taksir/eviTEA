from __future__ import annotations

from pathlib import Path

import pytest

from evitea.config import load_config
from evitea.contigs import (
    DEVELOPMENT_PRIMARY_CONTIGS,
    default_primary_contigs,
    parse_contig_list,
    resolve_primary_contigs,
)
from evitea.errors import InputContractError
from evitea.full_memory.proposals import ProposalRecord, full_data_proposals

MB = 1_000_000

DM6 = {"chr2L": 23_513_712, "chr2R": 25_286_936, "chr3L": 28_110_227, "chr3R": 32_079_331, "chr4": 1_348_131,
       "chrM": 19_524, "chrX": 23_542_271, "chrY": 3_667_352, "chrX_KI270757v1_alt": 120_000,
       "chrUn_DS483562v1": 30_000, "chr2L_random": 200_000}
MOUSE = {**{f"chr{i}": (200 - 7 * i) * MB for i in range(1, 20)}, "chrX": 171 * MB, "chrY": 91 * MB, "chrM": 16_299,
         "chr1_GL456210_random": 169_725, "chrUn_GL456239": 40_056, "JH584299.1": 953_012}
YEAST = {name: length for name, length in zip(
    ("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito"),
    (230_218, 813_184, 316_620, 1_531_933, 576_874, 270_161, 1_090_940, 562_643, 439_888, 745_751, 666_816,
     1_078_177, 924_431, 784_333, 1_091_291, 948_066, 85_779))}


def test_default_primary_contigs_of_non_fly_references() -> None:
    assert default_primary_contigs(MOUSE) == tuple(f"chr{i}" for i in range(1, 20)) + ("chrX", "chrY")
    assert default_primary_contigs(YEAST) == ("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII",
                                              "XIII", "XIV", "XV", "XVI")
    # dm6: the rule keeps chrY (11 % of chr3R) and leaves out chr4 (4 %); the paper runs pass the five arms
    assert default_primary_contigs(DM6) == ("chr2L", "chr2R", "chr3L", "chr3R", "chrX", "chrY")
    assert default_primary_contigs({"ChrC": 154_478, "Chr1": 30 * MB, "ChrM": 367_808}) == ("Chr1",)
    assert default_primary_contigs({}) == ()


def test_primary_contig_option_wins_and_is_checked() -> None:
    names, source = resolve_primary_contigs(DEVELOPMENT_PRIMARY_CONTIGS, DM6)
    assert (names, source) == (("chr2L", "chr2R", "chr3L", "chr3R", "chrX"), "option")
    assert resolve_primary_contigs(None, YEAST)[1] == "reference-rule"
    with pytest.raises(InputContractError, match="not in the reference"):
        resolve_primary_contigs(("chr2L", "2L"), DM6)


def test_contig_list_from_text_or_file(tmp_path: Path) -> None:
    assert parse_contig_list("chr2L, chr2R,chrX") == ("chr2L", "chr2R", "chrX")
    listed = tmp_path / "primary.txt"
    listed.write_text("# fly arms\nchr2L\nchr2R\n\nchrX  # last\n", encoding="utf-8")
    assert parse_contig_list(str(listed)) == ("chr2L", "chr2R", "chrX")
    with pytest.raises(InputContractError):
        parse_contig_list("chr2L,chr2L")
    config = tmp_path / "run.toml"
    config.write_text('[runtime]\nprimary_contigs = ["I", "II"]\n', encoding="utf-8")
    assert load_config(config).runtime.primary_contigs == ("I", "II")


def test_proposals_look_only_at_the_primary_contigs() -> None:
    # two pairs over two folds on chr1 make one region when chr1 is primary, none with the fly arms
    records = [ProposalRecord(pid, pid % 5, "fam", "chr1", 1_000, 1_800) for pid in (1, 2)]
    assert full_data_proposals("s", records) == []
    rows = full_data_proposals("s", records, frozenset({"chr1"}))
    assert [(r["chrom"], r["candidate_lo"], r["candidate_hi"], r["support_pair_count"]) for r in rows] == [
        ("chr1", 1_000, 1_800, 2)
    ]
