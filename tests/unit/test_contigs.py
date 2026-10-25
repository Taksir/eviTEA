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


# -- bare names and the mitochondrial contig --------------------------------------------------------------------


def test_norm_chrom_adds_chr_only_when_the_reference_has_it() -> None:
    from evitea.resident_filter import CopyIndex, call_contigs, norm_chrom

    assert norm_chrom("1", {"chr1", "chr2"}) == "chr1"
    assert norm_chrom("1", {"1", "chr1"}) == "1"
    assert norm_chrom("1", {"chr2"}) == "1"
    assert norm_chrom(" 2L ", {"chr2L"}) == "chr2L"
    assert norm_chrom("I", {"I", "II"}) == "I"
    index = CopyIndex([("1", 100, 200, "L1MdA_I", "5.0"), ("chrUn_x", 1, 5, "B1_Mm", "1.0")], {"chr1", "chrUn_x"})
    assert index.match("chr1", 150, 160, frozenset({"L1MdA_I"})) == ("chr1", 100, 200, "L1MdA_I", "5.0")
    calls = b"event_id\tchr\tstart\tend\tfamily\nE1\tchr1\t1\t2\tB1\nE2\tchr7\t1\t2\tB1\n"
    assert call_contigs(calls) == {"chr1", "chr7"}


def test_library_free_copy_index_uses_the_same_rule(tmp_path: Path) -> None:
    from evitea.library_free.resident import CopyIndex

    bed = tmp_path / "copies.bed"
    bed.write_text("1\t10\t20\tB1\n2L\t5\t9\tx\nchrX\t1\t3\ty\n", encoding="utf-8")
    index = CopyIndex.from_file(bed, contigs={"chr1", "chr2L", "chrX"})
    assert index.inside("chr1", 15) and index.inside("chr2L", 6) and index.inside("chrX", 2)
    assert not CopyIndex.from_file(bed, contigs={"1", "chr1"}).inside("chr1", 15)


def test_mitochondrial_contig_lookup() -> None:
    from evitea.contigs import find_mito_contig

    assert find_mito_contig(["1", "2", "MT"]) == "MT"
    assert find_mito_contig(["chr1", "chrM", "MT"]) == "chrM"
    assert find_mito_contig(["I", "II", "M"]) == "M"
    assert find_mito_contig(["chr1", "chrMT"]) == "chrMT"
    assert find_mito_contig(["chr1", "chr2"]) is None
    assert find_mito_contig(["chr1", "Mito"], "Mito") == "Mito"
    with pytest.raises(InputContractError, match="not a contig"):
        find_mito_contig(["chr1"], "chrM")


def test_library_free_detection_without_a_mitochondrial_contig_runs_and_says_so(tmp_path: Path) -> None:
    pytest.importorskip("mappy")
    from evitea.library_free import _p4_lib as P
    from evitea.library_free.options import Options
    from evitea.library_free.runtime import _detect

    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + "ACGTTGCAAGGCTTAC" * 200 + "\n", encoding="ascii")
    firewall = P.MappyFirewall(str(reference), None)
    assert firewall.numt("ACGTTGCAAGGCTTAC" * 5) is False
    messages: list[str] = []
    tables, info = _detect(reference_fasta=reference, genome={"chr1": "ACGT" * 800}, msbwt_index=tmp_path / "none",
                           contigs=(), work=tmp_path / "work", progress=messages.append, opts=Options())
    assert info["mito_contig"] is None
    assert any("no mitochondrial contig found" in m for m in messages)
    messages.clear()
    _, info = _detect(reference_fasta=reference, genome={"chr1": "ACGT" * 800, "MT": "ACGT" * 100},
                      msbwt_index=tmp_path / "none", contigs=(), work=tmp_path / "work2", progress=messages.append,
                      opts=Options())
    assert info["mito_contig"] == "MT"
    assert any("mitochondrial contig MT" in m for m in messages)
