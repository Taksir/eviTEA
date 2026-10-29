"""The local heterozygous guard: a call with e above HOM_MAX is not added among confirmed places that all read as homozygous."""

from __future__ import annotations

import dataclasses
import math
import random

import pytest

from evitea import cli
from evitea.library_free import _p4_lib as P
from evitea.library_free import assemble as AS
from evitea.library_free import calling as CL
from evitea.library_free import het_guard as HG
from evitea.library_free import runtime as R
from evitea.library_free import single_copy as SC
from evitea.library_free.options import Options
from evitea.run_command import add_arguments, library_steps, no_library_steps

rng = random.Random(23)
NAN = float("nan")


def rnd(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def place(bL="100", bR="95", tsd="5"):
    return dict(bL=bL, bR=bR, tsd=tsd)


def test_confirmed_place():
    assert HG.confirmed("first", "ADDED", place(), 0.0) and HG.confirmed("first", "NEAR_EVITEA", place(), 0.4)
    assert not HG.confirmed("second", "ADDED", place(), 0.0)                       # first pass only
    for status in ("DROPPED_V2", "DROPPED_RULE_C", "DROPPED_CHECK", "DROPPED_RESIDENT_COPY", "DROPPED_SINGLE_COPY"):
        assert not HG.confirmed("first", status, place(), 0.0)
    assert not HG.confirmed("first", "ADDED", place(), NAN)                        # no measured e
    assert not HG.confirmed("first", "ADDED", place(bR=""), 0.0)                   # one junction
    for tsd, ok in (("1", False), ("2", True), ("20", True), ("21", False), ("", False), ("x", False)):
        assert HG.confirmed("first", "ADDED", place(tsd=tsd), 0.0) is ok


def test_only_a_heterozygous_looking_call_can_be_dropped():
    assert HG.dropped_candidate(P.HOM_MAX + 1e-6) and HG.dropped_candidate(0.9)
    assert not HG.dropped_candidate(P.HOM_MAX) and not HG.dropped_candidate(0.0) and not HG.dropped_candidate(NAN)


def _brute(places, pos):
    """The K places nearest to pos, ties to the left."""
    near = sorted(places, key=lambda x: (abs(x[0] - pos), x[0] > pos, x[1]))[:HG.K]
    return sum(h for _, h in near) / len(near)


def test_share_is_that_of_the_ten_nearest_places_ties_to_the_left():
    for trial in range(200):
        n = rng.choice([5, 6, 9, 10, 11, 15, 40])
        spots = [(x * 10, rng.random() < 0.4) for x in rng.sample(range(60), n)]         # many ties in distance
        profile = HG.Profile.build(("c1", p, h) for p, h in spots)
        profile = dataclasses.replace(profile, n=max(profile.n, HG.NMIN))             # past the sample minimum
        for pos in (0, 5, 125, 300, 595, 700):
            assert profile.share("c1", pos) == pytest.approx(_brute(spots, pos)), (trial, pos)


def test_fewer_than_five_places_on_the_contig_use_the_sample_share():
    items = [("c1", 100 * i, i % 2 == 0) for i in range(8)] + [("c2", 100 * i, False) for i in range(4)] + [("c3", 5, True)]
    p = HG.Profile.build(items)
    assert p.n == 13 and p.het == 5 and p.sample_share == pytest.approx(5 / 13)
    assert p.share("c2", 50) == pytest.approx(5 / 13) and p.share("c3", 5) == pytest.approx(5 / 13)     # fewer than 5 here
    assert p.share("absent", 5) == pytest.approx(5 / 13)
    assert p.share("c1", 50) == pytest.approx(4 / 8)                                                  # 8 places: all of them


def test_a_sample_with_fewer_than_ten_confirmed_places_has_no_guard():
    p = HG.Profile.build(("c1", 100 * i, False) for i in range(HG.NMIN - 1))
    assert not p.active and p.share("c1", 0) == 1.0 and not p.drop("c1", 0, 0.5)
    q = HG.Profile.build(("c1", 100 * i, False) for i in range(HG.NMIN))
    assert q.active and q.share("c1", 0) == 0.0 and q.drop("c1", 0, 0.5) and not q.drop("c1", 0, P.HOM_MAX)


def test_the_call_itself_counts_among_its_ten_nearest():
    # the call at 1000 is itself a confirmed, heterozygous-looking place: 3 of its 10 nearest (itself, 990, 1010) are heterozygous-looking
    others = [(990, True), (1010, True)] + [(p, False) for p in (980, 1020, 970, 1030, 960, 1040, 950, 1050)]
    with_self = HG.Profile.build(("c1", p, h) for p, h in others + [(1000, True)])
    without = HG.Profile.build(("c1", p, h) for p, h in others)
    assert with_self.share("c1", 1000) == pytest.approx(0.3) and not with_self.drop("c1", 1000, 0.4)
    assert without.share("c1", 1000) == pytest.approx(0.2) and without.drop("c1", 1000, 0.4)


def test_guard_parse_and_merge_view():
    assert HG.Guard.parse("on").on and not HG.Guard.parse("off").on
    with pytest.raises(ValueError):
        HG.Guard.parse("auto")
    g = HG.Guard(True)
    assert g.for_merge_view().fixed is None and g.for_merge_view().source == "main view"
    g.built = HG.Profile.build([("c1", 1, False)])
    m = g.for_merge_view()
    assert m.on and m.fixed is g.built and m.source == "main view"
    assert not HG.Guard(False).for_merge_view().on


GENOME = {"chrS": rnd(40_000), "chrA": rnd(40_000)}
D0 = 20.0


class Table:
    """Places rows of one pass and their segments (as tests/unit/test_single_copy.py builds them)."""

    def __init__(self):
        self.rows, self.segs, self.cluster = [], {}, 0

    def family(self, chrom, sites, e=0.0, ratio=1.0, tsd="5"):
        cl = str(self.cluster)
        self.cluster += 1
        body = rnd(200)
        for bL in sites:
            rs, T = rnd(150), str(int(round(ratio * D0)))
            self.rows.append(dict(
                chrom=chrom, first=str(bL - 10), last=str(bL), n_low="1", D=str(D0), stage="S4_pass", why="", bL=str(bL),
                bR=str(bL - 5), delta="5", tsd=tsd, two_sided="1", geno="het", fracL=str(e), fracR=str(e), rep="0.1",
                lstop="support", rstop="support", rstat="break", seg_keep="L,R", seg_reject="", dup="0", cluster=cl,
                track="family", fam_pair="1", fam_vtsd="1", target="", lseg_len="200", rseg_len="150", chkL_T=T,
                chkL_junction=rnd(60), chkR_T=T, chkR_junction=rnd(60), origin="", walked="", bL0=""))
            self.segs[(chrom, str(bL), cl)] = dict(lseg=body, rseg=rs)

    def get(self):
        return self.rows, self.segs


LIBRARY = (1000, 1010, 1020)
HOM = [3000 + 500 * i for i in range(9)]            # chrS: nine homozygous-looking confirmed places
HET_A = [2000 + 300 * i for i in range(6)]          # chrA: six heterozygous-looking confirmed places


def _tables():
    first, second = Table(), Table()
    first.family("chrS", list(LIBRARY))                          # library-confirmed, e 0: confirmed places
    first.family("chrS", HOM)                                    # ADDED, e 0: confirmed places
    first.family("chrS", [3250], e=0.3)                          # het-looking among homozygous ones: dropped
    first.family("chrS", [3750], e=P.HOM_MAX)                    # at the bound: stays
    first.family("chrS", [20_000], e=0.3, tsd="1")               # no TSD of 2 to 20: a candidate, not a confirmed place; dropped
    first.family("chrS", [30_000], e=0.3, ratio=2.0)             # fails the one-site check: stays DROPPED_CHECK
    first.family("chrA", HET_A, e=0.3)                           # a region of heterozygous-looking places: stays
    first.family("chrA", [2150], e=0.4)                          # among them: stays
    second.family("chrS", [15_000], e=0.3)                       # haplotype scan, among homozygous places: dropped
    second.family("chrA", [6000], e=0.3)                         # among heterozygous-looking places: stays
    return first.get(), second.get()


def _assemble(single_copy=None, guard=None):
    first, second = _tables()
    index = CL.CallIndex({"chrS": [(p, p + 8) for p in LIBRARY]})
    return AS.assemble(first, second, index, Options.previous().with_(check_scope="unconfirmed"), identity_filter=False,
                       reference=GENOME, single_copy=single_copy, het_guard=guard)


def _fields(rows, status=True):
    """The rows as text (nan equals nan), with or without the status."""
    return [repr(dataclasses.replace(r, status="") if not status else r) for r in rows]


def _status(rows):
    return {(r.scan, r.chrom, r.bL): r.status for r in rows}


def test_guard_changes_only_added_calls_and_keeps_regions_of_het_looking_places():
    base, base_summary = _assemble()
    b = _status(base)
    assert b[("first", "chrS", 3250)] == "ADDED" and b[("first", "chrS", 30_000)] == "DROPPED_CHECK"
    guard = HG.Guard(True)
    rows, summary = _assemble(guard=guard)
    s = _status(rows)
    changed = {k for k in b if b[k] != s[k]}
    assert changed == {("first", "chrS", 3250), ("first", "chrS", 20_000), ("second", "chrS", 15_000)}
    assert all(b[k] == "ADDED" and s[k] == HG.DROPPED for k in changed)
    assert s[("first", "chrS", 3750)] == "ADDED"                       # e at the bound is homozygous
    assert s[("first", "chrA", 2150)] == "ADDED" and s[("second", "chrA", 6000)] == "ADDED" and s[("first", "chrA", 2000)] == "ADDED"
    assert s[("first", "chrS", 30_000)] == "DROPPED_CHECK"
    info = summary["het_guard"]
    # confirmed: 3 library + 9 + 1 at the bound + the het-looking one at 3250 + 6 + 1 on chrA; the call without a TSD is not one
    assert info["mode"] == "on" and info["active"] and info["confirmed"] == 3 + 9 + 1 + 1 + 6 + 1
    assert info["heterozygous_looking"] == 1 + 6 + 1 and info["dropped"] == {"first": 2, "second": 1}
    assert info["dropped_per_contig"] == {"chrS": 3} and info["source"] == "own calls"
    assert info["near_bar"]["dropped"] == 0 and guard.built is not None and guard.built.n == info["confirmed"]
    # every other field of every row is unchanged, and the added-calls table loses exactly the dropped calls
    assert _fields(rows, status=False) == _fields(base, status=False)
    assert len(AS.added_calls_tsv_bytes(base).decode().splitlines()) - len(AS.added_calls_tsv_bytes(rows).decode().splitlines()) == 3
    assert summary["tolerance"] == base_summary["tolerance"]


def test_off_and_no_guard_give_the_rows_without_the_rule():
    base, _ = _assemble()
    rows, summary = _assemble(guard=HG.Guard(False))
    assert _fields(rows) == _fields(base) and summary["het_guard"] == {"mode": "off"}
    assert "het_guard" not in _assemble()[1]


def test_a_sample_with_few_confirmed_places_is_left_alone():
    old = HG.NMIN
    HG.NMIN = 30
    try:
        rows, summary = _assemble(guard=HG.Guard(True))
    finally:
        HG.NMIN = old
    assert _fields(rows) == _fields(_assemble()[0]) and summary["het_guard"]["active"] is False and summary["het_guard"]["dropped"] == {"first": 0, "second": 0}


def test_the_single_copy_rule_acts_first_and_changes_nothing_the_guard_counts():
    chrs = SC.Decision(SC.Choice.parse("chrS"))
    only_sc, _ = _assemble(single_copy=chrs)
    both, summary = _assemble(single_copy=chrs, guard=HG.Guard(True))
    plain, plain_summary = _assemble(guard=HG.Guard(True))
    b, o, p = _status(both), _status(only_sc), _status(plain)
    # the single-copy rule drops the het-looking added calls of chrS first; the guard never relabels them
    assert o[("first", "chrS", 3250)] == SC.DROPPED and b[("first", "chrS", 3250)] == SC.DROPPED
    assert b == o and not any(v == HG.DROPPED for v in b.values())     # nothing is left on chrS for the guard to drop
    assert any(v == HG.DROPPED for v in p.values())                    # without the single-copy rule it drops them itself
    # the confirmed places were taken before that drop: the same set, so the same share, with or without the single-copy rule
    assert summary["het_guard"]["confirmed"] == plain_summary["het_guard"]["confirmed"]
    assert summary["het_guard"]["heterozygous_looking"] == plain_summary["het_guard"]["heterozygous_looking"]
    assert p[("first", "chrA", 2150)] == b[("first", "chrA", 2150)] == "ADDED"


def test_the_merge_view_judges_its_calls_against_the_main_views_places():
    main = HG.Guard(True)
    rows, summary = _assemble(guard=main)
    merge_guard = main.for_merge_view()
    first, second = _tables()
    # the merge view's own places differ (here: only the chrS het-looking calls and the homozygous ones are present)
    first = ([r for r in first[0] if not (r["chrom"] == "chrA")], first[1])
    second = ([r for r in second[0] if r["chrom"] == "chrS"], second[1])
    index = CL.CallIndex({"chrS": [(p, p + 8) for p in LIBRARY]})
    m_rows, m_summary = AS.assemble(first, second, index, Options.previous().with_(check_scope="unconfirmed"), identity_filter=False,
                                    reference=GENOME, het_guard=merge_guard)
    info = m_summary["het_guard"]
    assert info["source"] == "main view" and info["confirmed"] == summary["het_guard"]["confirmed"]
    assert merge_guard.built is None                      # it built nothing of its own
    m = _status(m_rows)
    for key, status in _status(rows).items():
        if key in m and status == HG.DROPPED:
            assert m[key] == HG.DROPPED, key
    # without the main view's places the merge view alone would see a smaller sample
    alone, alone_summary = AS.assemble(first, second, index, Options.previous().with_(check_scope="unconfirmed"), identity_filter=False,
                                       reference=GENOME, het_guard=HG.Guard(True))
    assert alone_summary["het_guard"]["confirmed"] < info["confirmed"]


def test_dropped_calls_count_as_kept_in_the_family_records():
    from evitea.family_records import select

    assert HG.DROPPED in select.AWAY_KEPT


def test_saved_scan_pieces_do_not_depend_on_the_guard():
    assert "het_guard.py" not in R.CODE_FILES and "assemble.py" not in R.CODE_FILES and "runtime.py" not in R.CODE_FILES


def test_command_lines():
    parse = lambda *a: cli.build_parser().parse_args(["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o", *a])
    assert parse().het_guard == "on" and parse("--het-guard", "off").het_guard == "off"
    with pytest.raises(SystemExit):
        parse("--het-guard", "maybe")

    import argparse
    from pathlib import Path

    p = argparse.ArgumentParser()
    add_arguments(p)
    base = ["--reference", "ref.fa", "--read1", "a", "--read2", "b", "--output", "out", "--repeat-annotation", "rm.bed"]
    for flags, steps in ((["--no-library"], no_library_steps), (["--te-library", "te.fa", "--family-map", "m", "--library-free"], library_steps)):
        on, off = p.parse_args(base + flags), p.parse_args(base + flags + ["--het-guard", "off"])
        pick = lambda args: [w for name, w in steps(args, Path("out"), "s") if name == "library_free"][0]
        assert "--het-guard" not in pick(on) and pick(off)[-2:] == ["--het-guard", "off"]
