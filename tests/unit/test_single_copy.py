"""The single-copy rule: on a contig the sample carries in one copy, haplotype-scan calls and first-pass calls with the
heterozygous genotype are not added."""

from __future__ import annotations

import collections
import dataclasses
import random

import pytest

from evitea.library_free import _p4_lib as P
from evitea.library_free import assemble as AS
from evitea.library_free import calling as CL
from evitea.library_free import runtime as R
from evitea.library_free import single_copy as SC
from evitea.library_free.options import Options

rng = random.Random(11)


def rnd(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_choice_parse():
    assert SC.Choice.parse("auto") == SC.Choice("auto")
    assert SC.Choice.parse("none").mode == "none" and SC.Choice.parse("all").mode == "all"
    c = SC.Choice.parse("chrX, chrY")
    assert c.mode == "list" and c.names == ("chrX", "chrY") and c.as_text() == "chrX,chrY"
    for bad in ("", ",", "chrX,chrX"):
        with pytest.raises(ValueError):
            SC.Choice.parse(bad)


def test_dropped_calls():
    assert SC.dropped("second", 0.0) and SC.dropped("second", float("nan"))
    assert SC.dropped("first", P.HOM_MAX + 1e-6)
    assert not SC.dropped("first", P.HOM_MAX) and not SC.dropped("first", 0.0) and not SC.dropped("first", float("nan"))


def test_reference_copies_count_both_strands():
    rep = rnd(300)
    genome = {"c1": rnd(20_000) + rep + rnd(9_000) + P.rc(rep) + rnd(3_000).lower() + "N" * 40 + rnd(500),
              "c2": rnd(15_000) + rep[:80], "c3": rnd(20)}
    owner, words = SC._sample(genome, ["c1", "c2"], 2_000)
    old = SC.CHUNK
    SC.CHUNK = 5_000          # several steps per contig, words across step ends
    try:
        got = SC._reference_copies(genome, words)
    finally:
        SC.CHUNK = old
    seen = collections.Counter()
    for seq in genome.values():
        u = seq.upper()
        for i in range(len(u) - SC.K + 1):
            seen[u[i:i + SC.K]] += 1
    want = [seen[w] + seen[P.rc(w)] for w in words]
    assert list(got) == want
    assert set(want) >= {1, 2} and len(words) > 1_500 and set(owner) == {"c1", "c2"}


def test_measure_half_depth_contig():
    genome = {"chrA": rnd(150_000), "chrB": rnd(150_000), "chrX": rnd(150_000), "small": rnd(20_000)}
    rep = rnd(2_000)
    genome["chrA"] = genome["chrA"][:50_000] + rep + genome["chrA"][52_000:]
    genome["chrX"] = genome["chrX"][:50_000] + rep + genome["chrX"][52_000:]
    depth = {"chrA": 30, "chrB": 30, "chrX": 15, "small": 15}
    owner = {}
    for chrom, seq in genome.items():
        for i in range(len(seq) - SC.K + 1):
            owner.setdefault(seq[i:i + SC.K], set()).add(chrom)

    def count(w):
        where = owner.get(w, set()) | owner.get(P.rc(w), set())
        return sum(depth[c] for c in where) + (rng.random() < 0.2)

    profile = SC.measure(genome, count, list(genome), n_words=6_000)
    m = profile["contigs"]
    assert profile["genome_depth"] in (30, 30.5, 31)
    assert m["chrX"]["ratio"] == pytest.approx(0.5, abs=0.04) and m["chrA"]["ratio"] == pytest.approx(1.0, abs=0.04)
    assert m["small"]["words"] < SC.MIN_WORDS <= m["chrX"]["words"]
    # the shared stretch occurs twice in the reference: none of its words is counted
    assert profile["words"]["unique"] < profile["words"]["sampled"]


def _profile(ratios):
    return {"genome_depth": 30.0, "words": {}, "contigs": {c: {"depth": 30 * r, "words": 5_000, "ratio": r} for c, r in ratios.items()}}


def test_auto_needs_low_depth_and_mostly_homozygous_library_calls():
    d = SC.Decision(SC.Choice("auto"), profile=_profile({"chrX": 0.5, "chr1": 1.0, "chrY": 0.5, "chr2": 0.6}))
    near = ([("chrX", 0.0)] * 15 + [("chrX", 0.3)] * 1 + [("chr1", 0.0)] * 50 + [("chr1", 0.4)] * 20
            + [("chrY", 0.0)] * 9 + [("chr2", 0.0)] * 6 + [("chr2", 0.5)] * 6 + [("chrX", float("nan"))] * 30)
    flagged, info = d.decide(near, ["chrX", "chr1", "chrY", "chr2"])
    assert flagged == {"chrX"} and info["flagged"] == ["chrX"]
    assert info["contigs"]["chrX"]["near_evitea"] == 16 and info["contigs"]["chrX"]["near_evitea_het"] == 1
    assert info["contigs"]["chrY"]["low_depth"] and not info["contigs"]["chrY"]["guard"]   # fewer than 10 library calls
    assert info["contigs"]["chr2"]["low_depth"] and not info["contigs"]["chr2"]["guard"]   # half heterozygous: a pool
    assert not info["contigs"]["chr1"]["low_depth"]
    # the bound: 0.625 is in, a contig with too few counted words is not judged
    edge = SC.Decision(SC.Choice("auto"), profile=_profile({"chrX": SC.DEPTH_RATIO_MAX}))
    assert edge.decide([("chrX", 0.0)] * 10, ["chrX"])[0] == {"chrX"}
    few = _profile({"chrX": 0.5})
    few["contigs"]["chrX"]["words"] = SC.MIN_WORDS - 1
    assert SC.Decision(SC.Choice("auto"), profile=few).decide([("chrX", 0.0)] * 10, ["chrX"])[0] == frozenset()


def test_other_modes():
    near, contigs = [("chrX", 0.0)] * 20, ["chrX", "chr1"]
    assert SC.Decision(SC.Choice("none")).decide(near, contigs) == (frozenset(), {"mode": "none"})
    assert SC.Decision(SC.Choice("all")).decide([], contigs)[0] == {"chrX", "chr1"}
    assert SC.Decision(SC.Choice.parse("chr1")).decide(near, contigs)[0] == {"chr1"}
    flagged, info = SC.Decision(SC.Choice("auto"), reason="no msBWT index").decide(near, contigs)
    assert flagged == frozenset() and info["reason"] == "no msBWT index"
    flagged, info = SC.Decision.from_view(SC.Choice("auto"), ["chrX"], "main view").decide([], contigs)
    assert flagged == {"chrX"} and info["source"] == "main view"


GENOME = {"chrS": rnd(40_000), "chrA": rnd(40_000)}
D0 = 20.0


class Table:
    """Places rows of one pass and their segments (as tools/library_free_tests/test_assemble.py builds them)."""

    def __init__(self):
        self.rows, self.segs, self.cluster = [], {}, 0

    def family(self, chrom, sites, e=0.0, ratio=1.0):
        cl = str(self.cluster)
        self.cluster += 1
        body = rnd(200)
        for bL in sites:
            rs, T = rnd(150), str(int(round(ratio * D0)))
            self.rows.append(dict(
                chrom=chrom, first=str(bL - 10), last=str(bL), n_low="1", D=str(D0), stage="S4_pass", why="", bL=str(bL),
                bR=str(bL - 5), delta="5", tsd="5", two_sided="1", geno="het", fracL=str(e), fracR=str(e), rep="0.1",
                lstop="support", rstop="support", rstat="break", seg_keep="L,R", seg_reject="", dup="0", cluster=cl,
                track="family", fam_pair="1", fam_vtsd="1", target="", lseg_len="200", rseg_len="150", chkL_T=T,
                chkL_junction=rnd(60), chkR_T=T, chkR_junction=rnd(60), origin="", walked="", bL0=""))
            self.segs[(chrom, str(bL), cl)] = dict(lseg=body, rseg=rs)

    def get(self):
        return self.rows, self.segs


def _tables():
    first, second = Table(), Table()
    first.family("chrS", [1000, 1010, 1020])                     # library-confirmed: NEAR_EVITEA, set the tolerance
    first.family("chrS", [3000, 4000])                           # homozygous: stays ADDED
    first.family("chrS", [6000, 7000], e=0.12)                   # within HOM_MAX: stays ADDED
    first.family("chrS", [9000, 10000], e=0.3)                   # heterozygous: dropped on a single-copy contig
    first.family("chrS", [12000, 13000], e=0.3, ratio=2.0)       # fails the one-site check: stays DROPPED_CHECK
    first.family("chrA", [3000, 4000], e=0.3)                    # other contig: stays ADDED
    second.family("chrS", [20000, 21000], e=0.3)                 # haplotype scan: dropped on a single-copy contig
    second.family("chrA", [20000, 21000], e=0.3)
    return first.get(), second.get()


def _assemble(decision):
    first, second = _tables()
    index = CL.CallIndex({"chrS": [(p, p + 8) for p in (1000, 1010, 1020)]})
    return AS.assemble(first, second, index, Options.previous().with_(check_scope="unconfirmed"), identity_filter=False,
                       reference=GENOME, single_copy=decision)


def _status(rows):
    return {(r.scan, r.chrom, r.bL): r.status for r in rows}


def test_rule_changes_only_added_calls_on_flagged_contigs():
    base, base_summary = _assemble(None)
    b = _status(base)
    assert b[("first", "chrS", 9000)] == "ADDED" and b[("first", "chrS", 12000)] == "DROPPED_CHECK"
    assert b[("first", "chrS", 1000)] == "NEAR_EVITEA" and b[("second", "chrS", 20000)] == "ADDED"
    rows, summary = _assemble(SC.Decision(SC.Choice.parse("chrS")))
    s = _status(rows)
    changed = {k for k in b if b[k] != s[k]}
    assert changed == {("first", "chrS", 9000), ("first", "chrS", 10000), ("second", "chrS", 20000), ("second", "chrS", 21000)}
    assert all(b[k] == "ADDED" and s[k] == SC.DROPPED for k in changed)
    info = summary["single_copy"]
    assert info["dropped"] == {"first": 2, "second": 2} and info["near_bar"]["added"] == 2
    assert info["families_mostly_dropped"] == 2 and info["added_elsewhere_in_them"] == 0
    # every other field of every row is unchanged, and the added-calls table loses exactly the dropped calls
    assert [dataclasses.replace(r, status="") for r in rows] == [dataclasses.replace(r, status="") for r in base]
    added = AS.added_calls_tsv_bytes(rows).decode().splitlines()
    assert len(AS.added_calls_tsv_bytes(base).decode().splitlines()) - len(added) == 4
    assert summary["tolerance"] == base_summary["tolerance"]


def test_none_and_unflagged_auto_give_the_rows_without_the_rule():
    base, _ = _assemble(None)
    for d in (SC.Decision(SC.Choice("none")), SC.Decision(SC.Choice("auto"), profile=_profile({"chrS": 0.5, "chrA": 1.0}))):
        rows, summary = _assemble(d)       # auto: chrS has 3 library calls, fewer than the guard needs
        assert rows == base and "dropped" not in summary["single_copy"]


def test_auto_flags_through_the_guard():
    rows, summary = _assemble(SC.Decision(SC.Choice("auto"), profile=_profile({"chrS": 0.5, "chrA": 1.0})))
    assert summary["single_copy"]["contigs"]["chrS"]["near_evitea"] == 3
    old = SC.GUARD_MIN_CALLS
    SC.GUARD_MIN_CALLS = 3
    try:
        rows, summary = _assemble(SC.Decision(SC.Choice("auto"), profile=_profile({"chrS": 0.5, "chrA": 1.0})))
    finally:
        SC.GUARD_MIN_CALLS = old
    assert summary["single_copy"]["flagged"] == ["chrS"]
    assert _status(rows)[("second", "chrS", 20000)] == SC.DROPPED and _status(rows)[("second", "chrA", 20000)] == "ADDED"


def test_saved_scan_pieces_do_not_depend_on_the_rule():
    assert "single_copy.py" not in R.CODE_FILES and "assemble.py" not in R.CODE_FILES and "runtime.py" not in R.CODE_FILES
