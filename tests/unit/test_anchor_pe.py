"""Paired-end anchor support behavior on small synthetic inputs."""

from __future__ import annotations

import pytest

from evitea.anchor_pe.alignments import (
    host_query_name,
    parse_host_sam,
    parse_reference_sam,
    terminal_softclips,
)
from evitea.anchor_pe.anchors import channel_a, channel_b_universe, score_candidates
from evitea.anchor_pe.clusters import annotate_corrected_clusters, cluster_mappings, interval_gap
from evitea.anchor_pe.codebook import CodeUnit, SafeCodebook, classify_pairs, one_mate_universe
from evitea.errors import InputContractError

PROBE = "ACGTACGTTGCAACGTTGCAA"  # 21-mer


def unit(uid, family="F", te="FBte1", rank=1):
    rc = PROBE.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return CodeUnit(uid, family, te, rank, "3p", 21, PROBE, (PROBE, rc), 10, 31, True)


def test_host_query_names_are_an_explicit_profile():
    assert host_query_name(7, 2, "pid") == "P7"
    assert host_query_name(7, 2, "aka-historical") == "p7|h2"
    with pytest.raises(InputContractError):
        host_query_name(7, 2, "other")


def test_production_host_read_names_are_pid_and_aka_is_reproduction_only():
    from evitea.anchor_pe.alignments import (
        HOST_READ_NAME_PROFILES,
        PRODUCTION_HOST_READ_NAMES,
        REPRODUCTION_ONLY_HOST_READ_NAMES,
    )
    from evitea.anchor_pe.engine import AnchorPEInputs

    assert PRODUCTION_HOST_READ_NAMES == "pid"
    assert REPRODUCTION_ONLY_HOST_READ_NAMES == ("aka-historical",)
    assert HOST_READ_NAME_PROFILES == ("pid", "aka-historical")
    assert AnchorPEInputs.__dataclass_fields__["host_read_names"].default == "pid"
    with pytest.raises(InputContractError):
        host_query_name(7, 2, "iso")  # alias of the default profile: refused


def test_classify_pairs_states_and_primary_multiblock():
    cb = SafeCodebook((unit(0), unit(1, rank=2), unit(2, family="G", te="FBte2")))
    hits = {5: {1: {0, 1}, 2: set()}, 6: {1: {0}, 2: {1}}, 7: {1: {0}, 2: {2}}, 8: {1: set(), 2: {2}}}
    rows, eligible = classify_pairs(hits, cb)
    states = {r["pair_id"]: r["eligibility_state"] for r in rows}
    assert states == {5: "PRIMARY_MULTI_BLOCK_COHERENT", 6: "BOTH_MATES_TE_POSITIVE",
                      7: "FAMILY_CONFLICT", 8: "NO_MULTI_RECORD_GE2_BLOCKS"}
    assert set(eligible) == {5, 8} and eligible[8]["host_mate"] == 1
    assert set(one_mate_universe(hits, eligible)) == {5, 8}


def sam(tmp_path, lines):
    path = tmp_path / "x.sam"
    path.write_text("".join("\t".join(map(str, x)) + "\n" for x in lines))
    return path


def test_host_sam_unique_high_needs_one_placement_and_mapq20(tmp_path):
    rec = lambda q, flag, pos, mapq: (q, flag, "chr2L", pos, mapq, "10M", "*", 0, 0, "A" * 10, "*")  # noqa: E731
    host = parse_host_sam(sam(tmp_path, [rec("P1", 0, 100, 42), rec("p2|h1", 16, 50, 19),
                                         rec("P3", 0, 10, 42), rec("P3", 256, 900, 1)]))
    assert host[1]["host_class"] == "UNIQUE_HIGH" and (host[1]["start"], host[1]["end"]) == (99, 109)
    assert host[2]["host_class"] == "MULTI_OR_LOW" and host[2]["strand"] == "-"
    assert host[3]["host_class"] == "MULTI_OR_LOW" and host[3]["reported_alignments"] == 2


def test_reference_ledger_evidence_complete_vs_te_unexplained(tmp_path):
    cb = SafeCodebook((unit(0),))
    te_read = "G" * 30 + PROBE + "T" * 49  # 100 bp, probe at 30..51
    host_read = "C" * 100
    reads = {1: (te_read, host_read), 2: (te_read, host_read), 3: (te_read, host_read)}
    eligible = {p: {"pair_id": p, "family": "F", "te_mate": 1, "host_mate": 2, "primary": False} for p in reads}
    def mate(q, flag, pos, pnext, cigar, seq):
        return (q, flag, "chr2L", pos, 42, cigar, "=", pnext, 400, seq, "*")

    lines = [  # a proper placement signature requires both mates: (chrom, min pos, max pos, |tlen|)
        mate("P1", 0x2 | 0x40, 100, 500, "100M", te_read), mate("P1", 0x2 | 0x80 | 0x10, 500, 100, "100M", host_read),
        mate("P2", 0x2 | 0x40, 100, 500, "40S60M", te_read), mate("P2", 0x2 | 0x80 | 0x10, 500, 100, "100M", host_read),
        mate("P3", 0x2 | 0x40, 100, 500, "100M", te_read),  # lone mate: not a proper signature
    ]
    ledger, records = parse_reference_sam(sam(tmp_path, lines), eligible, reads, cb.by_family())
    assert ledger[1]["reference_state"] == "EVIDENCE_COMPLETE_CLEAN_REFERENCE"
    assert ledger[2]["reference_state"] == "LOCAL_PROPER_TE_UNEXPLAINED"
    assert ledger[3]["reference_state"] == "NO_LOCAL_PROPER"
    assert records[2][0]["cigar"] == "40S60M"


def test_terminal_softclips_require_terminal_s():
    assert terminal_softclips("40S60M", 100) == ((0, 40), None)
    assert terminal_softclips("60M40S", 100) == (None, (60, 100))
    assert terminal_softclips("5H95M", 95) == (None, None)


def test_cluster_mappings_gap_ids_and_keep_all():
    rows = [
        {"pair_id": 1, "family": "F", "chrom": "c", "start": 0, "end": 100},      # mid 50
        {"pair_id": 2, "family": "F", "chrom": "c", "start": 500, "end": 600},    # mid 550: gap 500 -> same
        {"pair_id": 3, "family": "F", "chrom": "c", "start": 1050, "end": 1152},  # mid 1101: gap 551 -> new
    ]
    assert [c["pair_ids"] for c in cluster_mappings(rows, require_folds=True, prefix="X")] == ["1;2"]
    kept = cluster_mappings(rows, require_folds=True, prefix="X", keep_all=True)
    assert [(c["cluster_id"], c["pair_ids"], c["event_2p2fold"]) for c in kept] == [(1, "1;2", 1), (2, "3", 0)]
    assert interval_gap(0, 10, 11, 20) == 1 and interval_gap(0, 10, 5, 6) == 0


def test_novel_nonresident_annotation():
    cluster = {"cluster_id": 1, "family": "F", "chrom": "c", "start": 1000, "end": 1100}
    base = [{"family": "F", "chrom": "c", "start": 1601, "end": 1700}]
    reps = [{"family": "F", "chrom": "c", "start": 1600, "end": 1700}]
    (c,) = annotate_corrected_clusters([cluster], base, reps)
    assert c["novel"] and c["resident_compatible"] is True and not c["novel_nonresident"]


def test_channel_a_anchor_and_channel_b_completion():
    cb = SafeCodebook((unit(0),))
    te_seq = "G" * 40 + PROBE + "T" * 39  # probe at 40..61 overlaps the 50-bp leading clip
    targets = [{"cluster_id": 3, "family": "F", "chrom": "c", "start": 10_000, "end": 10_200}]
    ledger, records, host = {}, {}, {}
    for pid in (0, 1, 7):  # folds 0, 1, 2
        ledger[pid] = {"family": "F", "reference_state": "LOCAL_PROPER_TE_UNEXPLAINED"}
        records[pid] = [{"chrom": "c", "seq": te_seq, "cigar": "50S50M", "start": 10_100, "end": 10_150, "strand": "+"}]
    for pid, strand, start in ((10, "+", 9_500), (21, "+", 9_600), (12, "-", 10_300), (33, "-", 10_400)):
        ledger[pid] = {"family": "F", "reference_state": "NO_LOCAL_PROPER"}
        host[pid] = {"host_class": "UNIQUE_HIGH", "chrom": "c", "start": start, "end": start + 100, "strand": strand, "mapq": 42}
    witnesses, anchors = channel_a(targets, {p: 3 for p in (0, 1, 7)}, ledger, records, cb.by_family())
    assert len(witnesses) == 3
    assert [(a["breakpoint0"], a["anchor_host_side"], a["channel_A_pair_ids"]) for a in anchors] == [(10_100, "RIGHT", "0;1;7")]
    rows, (cand,) = score_candidates(anchors, channel_b_universe(ledger, host))
    assert cand.b_state == "BILATERAL" and cand.left_b_pair_ids == (10, 21) and cand.right_b_pair_ids == (12, 33)
    assert rows[0]["bilateral_B_2p2fold"] == 1


def test_decoder_strictly_decodes_a_terminal_record():
    pytest.importorskip("edlib")
    from evitea.anchor_pe.decoder import decode_fragments

    te = "A" * 5 + "CCGGTTAACC" * 3 + PROBE + "GATTACAGGC" * 2  # PROBE at 35..56, 3' terminal
    cb = SafeCodebook((CodeUnit(0, "F", "FBte1", 1, "3p", 21, PROBE, (PROBE,), 35, 56, True),))
    read = te[10:]
    (row,) = decode_fragments([{"pair_id": 3, "te_mate": 1, "host_mate": 2, "unit_ids": [0],
                                "old_primary": False, "te_read": read}], cb, {"FBte1": te}, {"FBte1": "F"})
    assert row["decoder_orientation"] == "+" and row["decoder_best_edit"] == 0
    assert row["decoder_terminal_side"] == "3p" and row["strict_fragment_decodable"] == 1
    assert row["masked_source_class"] == "MASKED_UNIQUE_RECORD_CONFIRMS_DECODER"
