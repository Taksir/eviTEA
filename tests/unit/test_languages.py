"""v3 language builders: scanner parity and stage mechanics on synthetic data.

Byte equality with the frozen dm6 tables is checked by
validation/language_checks/run_language_equivalence.py.
"""

from __future__ import annotations

import random

import pytest

from evitea.languages.common import (
    build_interval_index,
    encode,
    oriented_code_map,
    overlap_query,
    reference_hit_class,
    revcomp,
    scan_exact,
    tsv_bytes,
)
from evitea.languages.terminal import _choose_record_codebook, _select_best_nonoverlapping, finalize_terminal


def _naive(seq, qmap, k):
    """The research scanner (rolling 2-bit code, reset on non-ACGT)."""

    enc = {c: v for c, v in zip("ACGTacgt", (0, 1, 2, 3) * 2)}
    mask, code, valid = (1 << (2 * k)) - 1, 0, 0
    for i, c in enumerate(seq):
        v = enc.get(c, -1)
        if v < 0:
            code, valid = 0, 0
            continue
        code = ((code << 2) | v) & mask
        valid += 1
        if valid >= k and code in qmap:
            yield i - k + 1, i + 1, qmap[code]


@pytest.mark.parametrize("k", (21, 25, 31, 35))
def test_vectorized_scan_equals_research_scanner(k):
    rng = random.Random(k)
    seq = "".join(rng.choice("ACGTN" if i % 97 == 0 else "ACGT") for i in range(20_000))
    probes = {min(seq[i:i + k], revcomp(seq[i:i + k])) for i in range(0, 19_000, 377) if "N" not in seq[i:i + k]}
    probes.add("A" * k)
    qmap = oriented_code_map(probes)
    seq = seq + "A" * (k + 3) + "n" + revcomp(seq[500:500 + k])
    assert list(scan_exact(seq, qmap, k)) == list(_naive(seq, qmap, k))


def test_encode_and_collision_guard():
    assert encode("ACGT") == 0b00011011
    with pytest.raises(ValueError):
        encode("ACGN")


def test_reference_hit_classes():
    rm = build_interval_index({"c": [(0, 100, "FBte1", "FBte1", "roo"), (150, 160, "FBte2", "FBte2", "gypsy")]})
    assert reference_hit_class("roo", "FBte1", 10, 31, overlap_query(rm, "c", 10, 31)) == "TARGET_CONTEXT_ANNOTATED_TE"
    assert reference_hit_class("gypsy", "FBte9", 10, 31, overlap_query(rm, "c", 10, 31)) == "OTHER_FAMILY_ANNOTATED_TE"
    assert reference_hit_class("roo", "FBte1", 90, 111, overlap_query(rm, "c", 90, 111)) == "TE_BOUNDARY_OR_PARTIAL"
    assert reference_hit_class("roo", "FBte1", 200, 221, overlap_query(rm, "c", 200, 221)) == "NON_TE_OR_UNANNOTATED_HOST"


def test_tsv_styles_follow_the_research_writers():
    rows = [{"a": True, "b": ["x", "y"], "c": None}]
    assert tsv_bytes(rows, ("a", "b", "c")) == b"a\tb\tc\n1\tx;y\t\n"
    assert tsv_bytes(rows, ("a", "b"), style="internal") == b"a\tb\n1\tx;y\n"
    assert tsv_bytes([], style="internal") == b""


def _cand(side, start, probe, k=21, positional=1, ref=0, distance=0):
    return {"terminal_side": side, "source_start": start, "source_end": start + k, "k": k,
            "canonical_strand_probe": probe, "positionally_interpretable": positional,
            "reference_exact_occurrence_count": ref, "terminal_distance_bp": distance}


def test_record_selection_prefers_two_per_side_without_overlap():
    cands = [_cand("5p", 0, "A1"), _cand("5p", 10, "A2"), _cand("5p", 40, "A3"), _cand("3p", 900, "B1"),
             _cand("3p", 950, "B2", positional=0)]
    chosen = _choose_record_codebook(cands)
    assert [c["canonical_strand_probe"] for c in chosen] == ["B1", "B2", "A1", "A3"]
    assert [c["probe_rank_within_te_record"] for c in chosen] == [1, 2, 3, 4]


def test_gap_fill_maximizes_nonoverlapping_count_then_balance():
    cands = [_cand("5p", 0, "L1", k=35), _cand("5p", 20, "L2"), _cand("3p", 100, "R1"), _cand("3p", 110, "R2", k=25)]
    chosen = _select_best_nonoverlapping(cands)  # L1/L2 and R1/R2 overlap pairwise
    assert [(c["canonical_strand_probe"], c["probe_rank_within_te_record"]) for c in chosen] == [("L2", 1), ("R1", 2)]
    wider = cands + [_cand("5p", 50, "L3"), _cand("3p", 200, "R3")]
    assert len(_select_best_nonoverlapping(wider)) == 4


def test_finalizer_repairs_duplicate_sequences_within_a_record():
    fields = ("family\tte_id\tprobe_rank_within_te_record\tterminal_side\tk\tcanonical_strand_probe\t"
              "query_sequences_both_strands\tsource_start\tsource_end\tterminal_distance_bp\tpositionally_interpretable\t"
              "source_te_distinct_occurrence_count\treference_exact_occurrence_count\treference_class\tlibrary_role\n")
    hybrid = (fields + "F\tT1\t1\t5p\t21\tAACAACAACAACAACAACAAC\tx\t0\t21\t0\t1\t1\t0\tREFERENCE_ABSENT\tr\n"
              "F\tT1\t2\t3p\t21\tAACAACAACAACAACAACAAC\tx\t100\t121\t0\t1\t1\t0\tREFERENCE_ABSENT\tr\n").encode()
    audit = ("family\tte_id\tcanonical_strand_probe\tterminal_side\tsource_start\tsource_end\tk\tacquisition_safe\t"
             "source_exact_position_match\tpositionally_interpretable\treference_exact_occurrence_count\tterminal_distance_bp\n"
             "F\tT1\tAACAACAACAACAACAACAAC\t5p\t0\t21\t21\t1\t1\t1\t0\t0\nF\tT1\tAACAACAACAACAACAACAAC\t3p\t100\t121\t21\t1\t1\t1\t0\t0\n"
             "F\tT1\tGGTGGTGGTGGTGGTGGTGGT\t3p\t60\t81\t21\t1\t1\t1\t0\t5\n").encode()
    summary = b"te_id\tfamily\nT1\tF\nT2\tF\n"
    result = finalize_terminal(hybrid_codebook=hybrid, rebuilt_source_audit=audit,
                               gapfill_source_audit=b"te_id\tk\n", rebuilt_te_summary=summary)
    assert result.counts["records_repaired_for_duplicate_sequence"] == ["T1"]
    from evitea.languages.common import read_tsv_bytes

    probes = [r["canonical_strand_probe"] for r in read_tsv_bytes(result.codebook)]
    assert sorted(probes) == ["AACAACAACAACAACAACAAC", "GGTGGTGGTGGTGGTGGTGGT"]
    assert b"T2\tF\t0\tUNSUPPORTED_TERMINAL_SAFE" in result.record_status


def test_prepare_te_can_bind_built_languages(tmp_path, monkeypatch):
    import evitea.languages as languages
    from evitea.assets.te import load_prepared_te_assets, prepare_te_assets, verify_te_assets
    from evitea.inputs import load_family_map, validate_te_fasta
    from test_default_engine import INTERNAL, TERMINAL
    from test_phase3_asset_manifests import FIXTURES, _prepare_reference

    class Built:
        terminal_codebook = TERMINAL.encode()
        internal_only_codebook = INTERNAL.encode()

        def audit_tables(self):
            return {"FINAL_SAFE_TERMINAL_CODEBOOK.tsv": TERMINAL.encode()}

    seen = {}

    def fake_build(**kwargs):
        seen.update(kwargs)
        return Built()

    monkeypatch.setattr(languages, "build_languages", fake_build)
    reference = _prepare_reference(tmp_path)
    te = prepare_te_assets(validate_te_fasta(FIXTURES / "te.fa"), load_family_map(FIXTURES / "te_families.tsv"),
                           FIXTURES / "te_families.tsv", reference, tmp_path / "t", build_languages=True)
    assert te.languages["terminal_final"]["source"] == "built:evitea.languages.v1"
    assert seen["reference_fasta"] == reference.reference_fasta_path
    assert seen["repeatmasker_bed"] == reference.repeat_annotation_path
    assert (tmp_path / "t" / "languages" / "audit" / "FINAL_SAFE_TERMINAL_CODEBOOK.tsv").is_file()
    verify_te_assets(load_prepared_te_assets(tmp_path / "t"), reference)
