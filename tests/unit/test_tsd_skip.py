"""The TSD skip of the family-identity filter: a family the filter would drop is kept when enough of its trusted members have a TSD."""

from __future__ import annotations

import argparse
import dataclasses
import functools
import random
from pathlib import Path

import pytest

from evitea import cli
from evitea.library_free import assemble as AS
from evitea.library_free import calling as CL
from evitea.library_free import runtime as R
from evitea.library_free import tsd_skip as TS
from evitea.library_free.options import Options
from evitea.run_command import add_arguments, library_steps, no_library_steps

rng = random.Random(41)


def rnd(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s: str, rate: float) -> str:
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in s)


GENOME = {"chrS": rnd(80_000)}
D0 = 20.0


def test_tsd_lengths_and_settings():
    for text, ok in (("2", True), ("20", True), ("5", True), ("1", False), ("21", False), ("", False), (None, False), ("x", False)):
        assert TS.has_tsd(text) is ok
    assert TS.Skip().on and TS.Skip().share == 0.3 and TS.Skip().minimum == 5
    for bad in (dict(share=0.0), dict(share=1.5), dict(minimum=0), dict(minimum=2.5)):
        with pytest.raises(ValueError):
            TS.Skip(**bad)


def test_judge_keeps_a_family_with_enough_trusted_members_and_enough_tsd():
    fams = {"a": (0.9, 10, ["5"] * 2 + ["1"] * 3),        # 5 trusted, 0.4
            "b": (0.9, 10, ["5"] * 3 + [""] * 7),         # 10 trusted, exactly 0.3
            "c": (0.9, 10, ["5"] * 4),                    # 4 trusted: too few
            "d": (0.9, 30, ["1"] * 20 + ["5"] * 2),       # 22 trusted, 0.09
            "e": (float("nan"), 3, [])}                   # none
    kept, report = TS.judge("first", fams, TS.Skip())
    assert kept == {"a", "b"}
    by = {r["cluster"]: r for r in report}
    assert by["b"]["tsd_share"] == 0.3 and by["b"]["kept"] and by["d"]["trusted"] == 22 and not by["d"]["kept"]
    assert by["e"]["tsd_share"] is None and by["e"]["c"] is None and by["a"]["calls"] == 10
    assert [r["cluster"] for r in report][:2] == ["d", "b"]                  # most trusted members first
    assert TS.judge("first", fams, TS.Skip(share=0.4))[0] == {"a"}
    assert TS.judge("first", fams, TS.Skip(minimum=4))[0] == {"a", "b", "c"}
    s = TS.summary(TS.Skip(), [report])
    assert s["judged"] == 5 and s["kept"] == 2 and s["without_trusted_member"] == 1 and len(s["families"]) == 4
    assert all("TSD share" in x for x in TS.lines(report, 5)) and len(TS.lines(report, 5)) == 3


class Table:
    """Places rows of one pass and their segments; ``diverse`` gives each member a copy of the family's segment with 10% changes (a family the filter drops)."""

    def __init__(self):
        self.rows, self.segs, self.cluster = [], {}, 0

    def family(self, sites, e=0.0, ratio=1.0, tsds=None, diverse=True):
        cl = str(self.cluster)
        self.cluster += 1
        body = rnd(200)
        for k, bL in enumerate(sites):
            tsd = (tsds[k] if tsds else "5")
            ls = mutate(body, 0.1) if diverse else body
            rs = rnd(150)
            T = str(int(round(ratio * D0)))
            self.rows.append(dict(
                chrom="chrS", first=str(bL - 10), last=str(bL), n_low="1", D=str(D0), stage="S4_pass", why="", bL=str(bL),
                bR=str(bL - 5), delta="5", tsd=str(tsd), two_sided="1", geno="het", fracL=str(e), fracR=str(e), rep="0.1",
                lstop="support", rstop="support", rstat="break", seg_keep="L,R", seg_reject="", dup="0", cluster=cl,
                track="family", fam_pair="1", fam_vtsd="1", target="", lseg_len=str(len(ls)), rseg_len=str(len(rs)),
                chkL_T=T, chkL_junction=rnd(60), chkR_T=T, chkR_junction=rnd(60), origin="", walked="", bL0=""))
            self.segs[("chrS", str(bL), cl)] = dict(lseg=ls, rseg=rs)
        return cl

    def get(self):
        return self.rows, self.segs


LIBRARY = (50_000, 50_010, 50_020, 50_030, 50_040)        # library calls near family N's places


@functools.cache
def _tables(second=True):
    first, scan = Table(), Table()
    first.family([1000 + 700 * i for i in range(6)])                                            # K: 6 trusted, all TSD
    first.family([10_000 + 700 * i for i in range(6)], tsds=["1"] * 6)                          # J: 6 trusted, no TSD
    first.family([20_000 + 700 * i for i in range(4)])                                          # F: 4 trusted: too few
    first.family([30_000 + 700 * i for i in range(6)], tsds=["5", "5", "1", "1", "1", "1"])     # M: 6 trusted, 0.33
    first.family([40_000 + 700 * i for i in range(6)], ratio=2.0)                               # U: fails the check: 0 trusted
    first.family(list(LIBRARY))                                                                 # N: near library calls
    first.family([60_000 + 700 * i for i in range(6)], diverse=False)                           # H: high identity: no filter
    if second:
        scan.family([70_000 + 900 * i for i in range(6)], e=0.3)                                # S: second-scan family, trusted
    return first.get(), (scan.get() if second else None)


def _assemble(skip, filter_on=True, scope="unconfirmed", second=True):
    first, scan = _tables(second)
    index = CL.CallIndex({"chrS": [(p, p + 8) for p in LIBRARY]})
    opts = Options.previous().with_(check_scope=scope)
    return AS.assemble(first, scan, index, opts, identity_filter=filter_on, reference=GENOME, tsd_skip=skip)


def _by_cluster(rows, scan="first"):
    out = {}
    for r in rows:
        if r.scan == scan:
            out.setdefault(r.cluster, set()).add(r.status)
    return out


def test_filter_alone_drops_the_diverged_families():
    rows, summary = _assemble(None)
    st = _by_cluster(rows)
    assert st["0"] == st["1"] == st["2"] == st["3"] == st["4"] == {"DROPPED_RULE_C"}
    assert st["6"] == {"ADDED"}                                                                    # high identity: not judged
    assert "rule_c_tsd_skip" not in summary


def test_skip_keeps_families_with_a_tsd_share_of_their_trusted_members():
    rows, summary = _assemble(TS.Skip())
    st = _by_cluster(rows)
    assert st["0"] == {"ADDED"}                                  # K: 6 trusted, all TSD
    assert st["1"] == {"DROPPED_RULE_C"}                         # J: no TSD
    assert st["2"] == {"DROPPED_RULE_C"}                         # F: too few trusted members
    assert st["3"] == {"ADDED"}                                  # M: 2 of 6 = 0.33, above the bar of 0.3
    assert st["4"] == {"DROPPED_RULE_C"}                         # U: no trusted member
    assert st["5"] == {"NEAR_EVITEA"}                            # N: trusted through the library calls; a kept family goes on to its other statuses
    assert st["6"] == {"ADDED"}
    info = summary["rule_c_tsd_skip"]
    assert info["on"] and info["share"] == 0.3 and info["minimum"] == 5
    fams = {(f["scan"], f["cluster"]): f for f in info["families"]}
    assert fams[("first", "0")]["trusted"] == 6 and fams[("first", "0")]["tsd_share"] == 1.0 and fams[("first", "0")]["kept"]
    assert fams[("first", "3")]["tsd"] == 2 and round(fams[("first", "3")]["tsd_share"], 2) == 0.33
    assert fams[("first", "2")]["trusted"] == 4 and not fams[("first", "2")]["kept"]
    assert ("first", "4") not in fams and info["without_trusted_member"] >= 1 and ("first", "6") not in fams
    # the bar is an option
    st4 = _by_cluster(_assemble(TS.Skip(share=0.4))[0])
    assert st4["0"] == {"ADDED"} and st4["3"] == {"DROPPED_RULE_C"}


def test_skip_on_the_haplotype_scan():
    off = _by_cluster(_assemble(None)[0], "second")
    on, summary = _assemble(TS.Skip())
    assert _by_cluster(on, "second")["0"] == {"ADDED"} and off["0"] == {"DROPPED_RULE_C"}
    assert any(f["scan"] == "second" and f["cluster"] == "0" and f["kept"] for f in summary["rule_c_tsd_skip"]["families"])


def test_skip_off_gives_the_rows_of_the_filter_alone():
    base, base_summary = _assemble(None)
    rows, summary = _assemble(TS.Skip(on=False))
    fields = lambda rs: [repr(r) for r in rs]
    assert fields(rows) == fields(base)                      # includes the check column: filter candidates are not sent to the check
    assert summary["rule_c_tsd_skip"]["note"] == "off" and summary["rule_c_tsd_skip"]["families"] == []
    assert summary["tolerance"] == base_summary["tolerance"]
    on, on_summary = _assemble(TS.Skip())
    assert on_summary["tolerance"] == base_summary["tolerance"]            # the tolerance still leaves the filter candidates out
    # filter candidates get a check value only with the skip on
    chk = lambda rs: {(r.scan, r.cluster): r.check for r in rs if r.rule_c_dropped and r.scan == "first"}
    assert all(v == "" for v in chk(base).values()) and any(v != "" for v in chk(on).values())


def test_without_the_filter_or_the_check_the_skip_does_nothing():
    base, _ = _assemble(None, filter_on=False)
    rows, summary = _assemble(TS.Skip(), filter_on=False)
    assert [repr(r) for r in rows] == [repr(r) for r in base]
    assert summary["rule_c_tsd_skip"]["note"] == "no family-identity filter" and summary["rule_c_tsd_skip"]["families"] == []
    rows, summary = _assemble(TS.Skip(), scope="none", second=False)
    assert summary["rule_c_tsd_skip"]["note"] == "no one-site check"
    assert [repr(r) for r in rows] == [repr(r) for r in _assemble(None, scope="none", second=False)[0]]


def test_only_the_filter_status_changes_for_a_kept_family():
    base, _ = _assemble(None)
    rows, _ = _assemble(TS.Skip())
    changed = [(a, b) for a, b in zip(base, rows) if a.status != b.status]
    assert changed and all(a.status == "DROPPED_RULE_C" for a, _ in changed)
    # every field other than the status and the check value is unchanged
    assert [dataclasses.replace(a, status="", check="") for a in base].__repr__() == [dataclasses.replace(b, status="", check="") for b in rows].__repr__()


def test_saved_scan_pieces_do_not_depend_on_the_skip():
    assert "tsd_skip.py" not in R.CODE_FILES


def test_command_lines_and_the_default_of_the_filter():
    parse = lambda *a: cli.build_parser().parse_args(["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o", *a])
    d = parse()
    assert d.rule_c == "on" and d.rule_c_tsd_skip == "on" and d.rule_c_tsd_share == 0.3 and d.rule_c_tsd_min == 5
    assert parse("--rule-c").rule_c == "on" and parse("--rule-c", "--adapter-test", "off").rule_c == "on"
    assert parse("--rule-c", "off").rule_c == "off" and parse("--rule-c-tsd-skip", "off").rule_c_tsd_skip == "off"
    o = parse("--rule-c-tsd-share", "0.4", "--rule-c-tsd-min", "8")
    assert o.rule_c_tsd_share == 0.4 and o.rule_c_tsd_min == 8
    with pytest.raises(SystemExit):
        parse("--rule-c", "maybe")


def test_every_route_runs_with_the_filter_on():
    p = argparse.ArgumentParser()
    add_arguments(p)
    base = ["--reference", "ref.fa", "--read1", "a", "--read2", "b", "--output", "out", "--repeat-annotation", "rm.bed"]
    lib = p.parse_args(base + ["--te-library", "te.fa", "--family-map", "m", "--library-free"])
    nolib = p.parse_args(base + ["--no-library"])
    step = lambda steps: [w for name, w in steps if name == "library_free"][0]
    parse = lambda words: cli.build_parser().parse_args(words)
    assert parse(step(library_steps(lib, Path("out"), "s"))).rule_c == "on"                   # the default
    assert parse(step(no_library_steps(nolib, Path("out"), "s"))).rule_c == "on"               # passed explicitly
