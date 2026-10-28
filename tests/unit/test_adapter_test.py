"""The library-free adapter test: places whose segments hold a sequencing-adapter core are set aside before the grouping."""

from __future__ import annotations

import random

from evitea.library_free import adapter as A
from evitea.library_free import priority as PR
from evitea.library_free import runtime as R

TRUSEQ = "AGATCGGAAGAGC"
NEXTERA = "CTGTCTCTTATACACATCT"
# A stretch of a TE-like sequence (no adapter core)
PLAIN = ("TGTTCACACATGAACACGAATATATTTAAAGACTTACAATTTTGGGCTCCGTTCATATCTTATGTAAATGAATCGAGAGCGATAAATTATATTTAGGATTTTG"
         "TTATCTAAGGCGACATGGGTGCATTGCTCAAAAACATGTAATTTAAGTGCACACTACATGAGTCAGTCACTTGAGATCGTTCCCCGCCTCCTAAAATAGTCC")


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_full_cores_anywhere_on_either_side():
    for core in (TRUSEQ, NEXTERA):
        assert A.adapter_segment("L", PLAIN[:50] + core + PLAIN[50:])
        assert A.adapter_segment("R", PLAIN[:50] + rc(core) + PLAIN[50:])
        assert A.adapter_segment("R", core + PLAIN)


def test_plain_segments_pass():
    assert not A.adapter_segment("L", PLAIN)
    assert not A.adapter_segment("R", PLAIN)
    assert not A.adapter_segment("L", "")
    rng = random.Random(7)
    hits = sum(A.adapter_segment(side, "".join(rng.choice("ACGT") for _ in range(300)))
               for side in "LR" for _ in range(500))
    assert hits <= 1


def test_shifted_core_at_the_junction_end():
    for k in range(1, A.SHIFT + 1):
        assert A.adapter_segment("L", TRUSEQ[k:] + PLAIN)
        assert A.adapter_segment("R", PLAIN + rc(TRUSEQ[k:]))
        assert A.adapter_segment("L", NEXTERA[k:] + PLAIN)
    # more than SHIFT letters inside the adapter, or the shifted core away from the junction end: not caught by the shift
    assert not A.adapter_segment("L", TRUSEQ[A.SHIFT + 1:] + PLAIN)
    assert not A.adapter_segment("L", PLAIN + TRUSEQ[1:])
    assert not A.adapter_segment("R", TRUSEQ[1:] + PLAIN)


def test_lower_case_segments():
    assert A.adapter_segment("L", (TRUSEQ[2:] + PLAIN).lower())


def test_mark_sets_aside_tested_stages_only():
    places = [dict(stage="S4_pass", lseg=TRUSEQ[1:] + PLAIN, rseg=PLAIN),
              dict(stage="S4_reject", lseg=PLAIN, rseg=PLAIN + rc(NEXTERA)),
              dict(stage="S4_pass", lseg=PLAIN, rseg=PLAIN),
              dict(stage="S4_abstain", lseg=TRUSEQ + PLAIN),
              dict(stage="S2_dead")]
    assert A.mark_adapter_places(places) == 2
    assert [p["stage"] for p in places] == ["S4_adapter", "S4_adapter", "S4_pass", "S4_abstain", "S2_dead"]


def test_set_aside_places_enter_no_duplicate_group_or_family():
    # dedupe ranks only these stages, and the grouping takes S4_pass places only
    assert A.STAGE not in PR.RANK


def test_saved_scan_pieces_do_not_depend_on_the_test():
    # the test runs after the pieces are loaded, so its module is not part of a piece's identity
    assert "adapter.py" not in R.CODE_FILES
