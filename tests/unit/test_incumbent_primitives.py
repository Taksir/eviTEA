"""Clean incumbent primitives: alignment semantics, index access, parsing, geometry."""

from __future__ import annotations

import random

import pytest

from evitea.errors import ResumeStateError
from evitea.incumbent.alignment import hw_distance, hw_locations
from evitea.incumbent.coverage_panel import coverage_interval, merge_intervals, terminal_flags
from evitea.incumbent.events import IncumbentCall, IncumbentCallClass, legacy_tsv_bytes
from evitea.incumbent.host_mate_gate import parse_primary_host_sam, supported_host_cluster_pairs
from evitea.incumbent.recall_repair import (
    RescuePairEvidence,
    parse_rescue_host_sam,
    rescue_host_requests,
    supported_rescue_pairs,
)

from incumbent_fakes import paired_sample


def _edit(a: str, b: str) -> int:
    row = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        prev, row[0] = row[0], i
        for j, y in enumerate(b, 1):
            prev, row[j] = row[j], min(row[j] + 1, row[j - 1] + 1, prev + (x != y))
    return row[-1]


def _brute(query: str, target: str):
    best = len(query)
    per_end = {}
    for end in range(len(target)):
        costs = [(_edit(query, target[start : end + 1]), start) for start in range(end + 1)]
        value = min(c for c, _ in costs)
        per_end[end] = (value, min(s for c, s in costs if c == value))
        best = min(best, value)
    return best, tuple((per_end[e][1], e) for e in sorted(per_end) if per_end[e][0] == best)


def test_infix_alignment_matches_exhaustive_semantics():
    rng = random.Random(7)
    for _ in range(300):
        query = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 9)))
        target = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 16)))
        best, locations = _brute(query, target)
        for k in range(4):
            expected = best if best <= k else None
            assert hw_distance(query, target, k) == expected
            distance, found = hw_locations(query, target, k)
            assert distance == expected
            if expected is not None and best < len(query):
                assert found == locations


def test_msbwt_sample_provenance_and_recovery():
    pairs = [("ACGTACGTAA", "TTGCATTGCA"), ("GGGACGTCCC", "ACGTTTTAAA"), ("CATCATCATC", "GAGAGAGAGA")]
    index = paired_sample(pairs)
    index.verify_paired_layout()
    assert index.pair_count == 3 and index.record_count == 6
    assert index.pair_reads([2, 0]) == {0: pairs[0], 2: pairs[2]}
    dollars = index.dollar_ids("ACGT")
    files, pair_ids = index.file_and_pair(dollars)
    assert sorted(zip(pair_ids.tolist(), files.tolist())) == [(0, 0), (0, 0), (1, 0), (1, 1)]
    total, complete, located = index.seed_pairs("ACGT", 100)
    assert complete and located == {0, 1} and total == 4
    assert index.seed_pairs("ACGT", 3)[1] is False


def test_legacy_serialization_is_frozen_tab_lf():
    call = IncumbentCall("S", "V2R4_MATE_000001", IncumbentCallClass.MATE_RESCUE, "chr2L", 10, 12, "roo",
                         None, None, None, 2, 2, "two_sided_mate_cluster_rescue", 4.0, (1, 2, 3, 4))
    text = legacy_tsv_bytes([call]).decode()
    assert text.splitlines()[0].startswith("event_id\tchr\tstart\tend\tfamily")
    assert "\r" not in text and text.endswith("4.0\t1,2,3,4\n")


def test_window_geometry():
    assert coverage_interval(200, 10, 21) == (0, 10)
    assert coverage_interval(100, 10, 21) is None
    assert merge_intervals([(5, 9), (0, 4), None, (20, 21)]) == ((0, 9), (20, 21))
    assert terminal_flags(300, [(0, 21)]) == (True, False)
    assert terminal_flags(300, [(279, 21)]) == (False, True)


def test_primary_and_rescue_host_sam_semantics(tmp_path):
    primary = tmp_path / "p.sam"
    primary.write_text(
        "p1|h2\t0\tchr2L\t101\t42\t10M\t*\t0\t0\tA\tI\n"
        "p2|h1\t16\tchr2L\t400\t1\t5M2D5M\t*\t0\t0\tA\tI\n"
        "p2|h1\t256\tchrX\t9\t255\t10M\t*\t0\t0\tA\tI\n"
        "p3|h2\t0\tchr2L\t600\t3\t10M\t*\t0\t0\tA\tI\n"
        "p4|h2\t0\tchr2L\t1200\t40\t10M\t*\t0\t0\tA\tI\n",
        encoding="ascii",
    )
    rows = parse_primary_host_sam(primary, [1, 2, 3, 4, 5])
    assert [rows[p].status for p in (1, 2, 3, 4, 5)] == ["UNIQUE_HIGH", "MULTI", "MAPPED_LOWCONF", "UNIQUE_HIGH", "UNMAPPED"]
    assert rows[1].midpoint == 105 and rows[2].midpoint == 405 and rows[2].strand == "-"
    assert supported_host_cluster_pairs(rows) == set()
    rescue = tmp_path / "r.sam"
    rescue.write_text(
        "@HD\tVN:1.0\n"
        "pair_7\t0\tchr3R\t50\t30\t10M\t*\t0\t0\tA\tI\n"
        "pair_8\t0\tchr3R\t450\t30\t10M\t*\t0\t0\tA\tI\n"
        "pair_9\t0\tchr3R\t60\t10\t10M\t*\t0\t0\tA\tI\n",
        encoding="ascii",
    )
    mapped = parse_rescue_host_sam(rescue, [7, 8, 9, 10])
    assert [mapped[p].status for p in (7, 8, 9, 10)] == ["UNIQUE_HIGH", "UNIQUE_HIGH", "UNIQUE_LOW", "UNMAPPED"]
    assert supported_rescue_pairs(mapped) == {7, 8}
    evidence = {1: RescuePairEvidence(1, True, False, ("roo",)), 2: RescuePairEvidence(2, True, True, ("roo",)),
                3: RescuePairEvidence(3, False, True, ("412",))}
    assert rescue_host_requests(evidence) == {1: 1, 3: 0}


def test_stage_store_reuses_only_matching_inputs(tmp_path):
    from evitea.incumbent.engine import StageStore

    store = StageStore(tmp_path, {"sample": "A"})
    assert store.get("stage", lambda: [1, 2]) == [1, 2]
    assert StageStore(tmp_path, {"sample": "A"}).get("stage", lambda: pytest.fail("recomputed")) == [1, 2]
    with pytest.raises(ResumeStateError):
        StageStore(tmp_path, {"sample": "B"}).get("stage", lambda: [0])
