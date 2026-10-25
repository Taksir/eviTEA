from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from evitea.engine import EngineInputs
from evitea.errors import InputContractError
from evitea.incumbent import coverage_panel, recall_repair
from evitea.sample_geometry import (
    DEFAULT_FRAGMENT_Q95,
    DEFAULT_FRAGMENT_Q99,
    DEFAULT_READ_LENGTH,
    FragmentSizes,
    estimate_fragment_sizes,
    fragment_sizes_from_spans,
    nearest_rank_percentile,
    properly_paired_spans,
    resolve_fragment_sizes,
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


# -- fragment size ----------------------------------------------------------------------------------------------


def _sam_pair(name: str, span: int, mapq: int = 42, flag1: int = 99, flag2: int = 147) -> list[str]:
    return [
        f"{name}\t{flag1}\tchr1\t100\t{mapq}\t50M\t=\t{100 + span - 50}\t{span}\tACGT\tIIII",
        f"{name}\t{flag2}\tchr1\t{100 + span - 50}\t{mapq}\t50M\t=\t100\t{-span}\tACGT\tIIII",
    ]


def test_properly_paired_spans_keep_only_proper_primary_confident_pairs() -> None:
    lines = ["@HD\tVN:1.0", *_sam_pair("a", 300), *_sam_pair("b", 410), *_sam_pair("low", 350, mapq=3),
             *_sam_pair("notproper", 360, flag1=97, flag2=145), *_sam_pair("secondary", 370, flag1=355, flag2=403)]
    assert sorted(properly_paired_spans(lines)) == [300, 410]


def test_fragment_percentiles_from_known_spans() -> None:
    # 1,000 pairs with spans 301..1300: nearest-rank 95th = 1250, 99th = 1290
    spans = list(range(301, 1301))
    sizes = fragment_sizes_from_spans(reversed(spans))
    assert (sizes.q95, sizes.q99, sizes.pairs, sizes.source) == (1250.0, 1290, 1000, "estimate")
    assert nearest_rank_percentile([1, 2, 3, 4], 50) == 2
    assert nearest_rank_percentile([7], 99) == 7
    # outliers far above median + 10 MAD (here spans of 50 kb) are dropped before the percentiles
    with_outliers = spans + [50_000] * 30
    assert fragment_sizes_from_spans(with_outliers) == sizes


def test_fragment_fallback_with_too_few_pairs() -> None:
    sizes = fragment_sizes_from_spans([400] * 999)
    assert (sizes.q95, sizes.q99, sizes.pairs, sizes.source) == (DEFAULT_FRAGMENT_Q95, DEFAULT_FRAGMENT_Q99, 999, "default")
    assert (DEFAULT_FRAGMENT_Q95, DEFAULT_FRAGMENT_Q99) == (494.0, 805)
    warnings: list[str] = []
    resolved = resolve_fragment_sizes(None, None, lambda: sizes, warn=warnings.append)
    assert resolved == sizes
    assert len(warnings) == 1 and "494" in warnings[0] and "805" in warnings[0]


def test_fragment_options_override_the_estimate() -> None:
    def never():
        raise AssertionError("both options given: no estimate")

    assert resolve_fragment_sizes(494, 805, never) == FragmentSizes(494.0, 805, 0, "option")
    estimate = FragmentSizes(380.0, 520, 4000, "estimate")
    assert resolve_fragment_sizes(None, 900, lambda: estimate) == FragmentSizes(380.0, 900, 4000, "estimate+option")
    assert resolve_fragment_sizes(450.5, None, lambda: estimate) == FragmentSizes(450.5, 520, 4000, "estimate+option")


class _Pairs:
    """A sample of ``pair_count`` pairs of known span: read 1 at a fragment's start, read 2 its far end."""

    def __init__(self, genome: str, spans: list[int], length: int = 50):
        self.genome, self.spans, self.length = genome, spans, length
        self.pair_count = len(spans)

    def pair_reads(self, pair_ids):
        from evitea.incumbent.sequence import revcomp

        out = {}
        for pid in pair_ids:
            start = (pid * 37) % (len(self.genome) - 2_000)
            fragment = self.genome[start:start + self.spans[pid]]
            out[pid] = (fragment[: self.length], revcomp(fragment[-self.length:]))
        return out


@pytest.mark.skipif(shutil.which("bowtie2") is None, reason="needs bowtie2 on PATH")
def test_estimate_maps_pairs_and_is_deterministic(tmp_path: Path) -> None:
    import random
    import subprocess

    rng = random.Random(3)
    genome = "".join(rng.choice("ACGT") for _ in range(60_000))
    (tmp_path / "ref.fa").write_text(f">chr1\n{genome}\n", encoding="ascii")
    subprocess.run(["bowtie2-build", "-q", str(tmp_path / "ref.fa"), str(tmp_path / "ref")], check=True)
    spans = [rng.randrange(200, 600) for _ in range(3_000)]
    sample = _Pairs(genome, spans)
    first = estimate_fragment_sizes(sample, bowtie2="bowtie2", host_index_prefix=tmp_path / "ref",
                                    work_directory=tmp_path / "a", sample_pairs=2_000)
    again = estimate_fragment_sizes(sample, bowtie2="bowtie2", host_index_prefix=tmp_path / "ref",
                                    work_directory=tmp_path / "b", sample_pairs=2_000, threads=2)
    used = sorted(spans[i * 3_000 // 2_000] for i in range(2_000))
    assert first == again
    assert first.source == "estimate" and first.pairs == 2_000
    assert (first.q95, first.q99) == (float(nearest_rank_percentile(used, 95)), nearest_rank_percentile(used, 99))
    few = estimate_fragment_sizes(_Pairs(genome, spans[:50]), bowtie2="bowtie2", host_index_prefix=tmp_path / "ref",
                                  work_directory=tmp_path / "c")
    assert few.source == "default" and (few.q95, few.q99) == (494.0, 805)


def test_reference_mappings_take_the_fragment_q99() -> None:
    from evitea.anchor_pe.alignments import REFERENCE_OPTIONS, reference_options
    from evitea.backends.bowtie2 import build_paired_reference_command

    assert reference_options() == REFERENCE_OPTIONS
    assert reference_options(1100)[reference_options(1100).index("-X") + 1] == "1100"
    default = build_paired_reference_command("bowtie2", "ref", "1.fa", "2.fa", "o.sam")
    assert default[default.index("-X") + 1] == "805"
    wider = build_paired_reference_command("bowtie2", "ref", "1.fa", "2.fa", "o.sam", max_insert=1100)
    assert wider[wider.index("-X") + 1] == "1100"


def test_q99_windows_follow_the_fragment_q99() -> None:
    from evitea.anchor_pe.anchors import score_candidates
    from evitea.full_memory.strict_a import build_b_index

    rows = [{"pair_id": p, "fold": p % 5, "family": "fam", "chrom": "chr1", "start": 1_000, "end": 1_100,
             "strand": "+", "mapq": 42} for p in (1, 2)]
    anchor = {"cluster_id": 1, "chrom": "chr1", "breakpoint0": 1_100 + 900, "family": "fam",
              "anchor_host_side": "LEFT", "channel_A_pairs": 2, "channel_A_folds": 2,
              "channel_A_pair_ids": "7;8", "channel_A_fold_ids": "2;3"}
    # the breakpoint is 900 bp past the mates: outside the default 805 window, inside a 1000 bp one
    assert score_candidates([anchor], rows)[0][0]["left_B_pairs"] == 0
    assert score_candidates([anchor], rows, 1_000)[0][0]["left_B_pairs"] == 2
    default, _ = build_b_index(rows)
    wide, _ = build_b_index(rows, 1_000)
    assert default.support("fam", "chr1", "LEFT", 2_000) == []
    assert len(wide.support("fam", "chr1", "LEFT", 2_000)) == 2
