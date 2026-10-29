"""The library-free run flow on a synthetic sample (rt_driver.py).   Run: python -B tools/library_free_tests/test_runtime.py   (venv python)

U1  previous behaviour (options None) writes the previous version's files byte for byte (the previous package is run on the same sample).
U2  a regroup arm in the same work directory (only grouping / assembly switches differ, here also the check scope) reuses the saved pieces and gives the
    same places tables as a fresh run with those switches.
U3  EVITEA_PIECES_READONLY refuses a run whose detection switches differ (it would replace the saved pieces); the pieces stay.
U4  --detection-dir (saved tables) gives the same family calls as the run that wrote them.
U6  the command line: no switch gives the default switches, --previous-behaviour gives every rule off, a switch given with it still applies,
    polya-cross without anchored extension is refused.
U7  the repeat map (--repeat-map): with an empty library layer it feeds the resident-copy test (the summary's TE share is the map's); with
    library calls it is not used and the files equal the run without it; an explicit annotation wins over the map; with neither, the test is off.
U8  identity evidence and relative depth floor: with --identity-min-comparisons 0 and the floor off, the files equal those of a version
    without these switches (BASELINE_SRC = its src folder; skipped without it); the command line reads both switches.
U5  default switches: the extra outputs exist (second scan tables) and the summary carries the options and the assembly.
U9  the merge: Options.merge(), and the merge_view/ files of a default run, equal the merge's default files byte for byte on the
    elements sample (one and two haplotypes) and the cross sample (MERGE_SRC = the merge's src folder; skipped without it).
U10 the cross-species rules on the cross sample (synth_cross.py): each rule turns its scenario's place into a call (shared-word
    depth, repeat cutoff, poly-A tails, recurring sequence), the low-coverage block and the local copy stay out, the anchor's own
    share keeps the cutoff at REP_MAX, the places table names the rule, and a run from the saved tables gives both views' calls.
U11 the repeat cutoff from window counts, and the command line of the cross-species switches.
U12 single-copy contigs: the default (auto) measures the contig depth once (copy_profile.json in the work directory, reused) and
    flags nothing without library calls; none gives the same call tables; naming the contig changes only added calls (haplotype-scan
    calls and first-pass calls with e above HOM_MAX) to DROPPED_SINGLE_COPY, and the merge view takes the main view's contigs.
U13 local heterozygous guard: on by default and recorded in the summary of both views (the merge view from the main view's places);
    the sample has fewer than 10 confirmed places, so the guard is idle and --het-guard on and off give the call tables of the default run.
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
# PREVIOUS_SRC: the src folder of the previous version (a checkout of the version to compare against); U1 is skipped without it.
PREVIOUS = os.environ.get("PREVIOUS_SRC")
# BASELINE_SRC: the src folder of a version without the identity-evidence and relative-floor switches; U8 is skipped without it.
BASELINE = os.environ.get("BASELINE_SRC")
# MERGE_SRC: the src folder of the merge (the version before the cross-species rules); U9 is skipped without it.
MERGE = os.environ.get("MERGE_SRC")
FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def run(out, src, opts="-", extra=(), env_extra=None, seed=7):
    env = dict(os.environ, PYTHONPATH=str(src))
    env.update(env_extra or {})
    return subprocess.run([sys.executable, "-B", str(HERE / "rt_driver.py"), str(out), str(seed), opts, *extra], env=env,
                          capture_output=True, text=True)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:12]


def main():
    t = Path(tempfile.mkdtemp(prefix="rt_"))
    if PREVIOUS:
        fz, mf = t / "fz", t / "mf"
        a = run(fz, PREVIOUS)
        b = run(mf, ROOT / "src")
        check("U1 both runs finish", a.returncode == 0 and b.returncode == 0, (a.stderr[-300:], b.stderr[-300:]))
        names = ["library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv", "calls.library_free.tsv", "calls.combined.tsv"]
        same = {n: sha(fz / n) == sha(mf / n) for n in names if (fz / n).exists()}
        check("U1 previous behaviour is byte-identical to the previous version", all(same.values()) and len(same) >= 4, same)
    else:
        print("U1 skipped: PREVIOUS_SRC is not set")

    base = json.dumps({})        # default switches
    d0 = t / "d0"
    r = run(d0, ROOT / "src", base, ["--het"])
    check("U5 default switches run", r.returncode == 0, r.stderr[-400:])
    check("U5 second-scan tables written", (d0 / "library_free_second_places.tsv").exists() and (d0 / "library_free_second_segments.fa").exists())
    summ = json.loads((d0 / "library_free_summary.json").read_text())
    check("U5 summary carries options and assembly", "options" in summ and "assembly" in summ, list(summ))

    # U2: same work dir (rt_driver names it <out>_work), other grouping switches
    keep = json.dumps({"polya_tsd": True, "het_yield": False, "check_scope": "path-dropped"})
    # the output folder of the regroup arm differs but the work dir must be the same: patch the driver's work path
    drv = (HERE / "rt_driver.py").read_text().replace('work_dir=out.parent / (out.name + "_work")',
                                                    'work_dir=out.parent / (os.environ.get("WORKNAME", out.name) + "_work")')
    drv = drv.replace('work = out.parent / (out.name + "_in")', 'work = out.parent / (os.environ.get("WORKNAME", out.name) + "_in")')
    drv_path = t / "rt_driver_w.py"
    drv_path.write_text(drv.replace('HERE = os.path.dirname(os.path.abspath(__file__))', f'HERE = "{HERE}"'))
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), WORKNAME="d0")
    regroup = t / "d2"
    r2 = subprocess.run([sys.executable, "-B", str(drv_path), str(regroup), "7", keep, "--het"], env=env, capture_output=True, text=True)
    check("U2 regroup arm runs in the saved work directory", r2.returncode == 0, r2.stderr[-400:])
    fresh = t / "d3"
    env3 = dict(env, WORKNAME="d3")
    r3 = subprocess.run([sys.executable, "-B", str(drv_path), str(fresh), "7", keep, "--het"], env=env3, capture_output=True, text=True)
    check("U2 fresh run with the same switches runs", r3.returncode == 0, r3.stderr[-400:])
    same = {n: sha(regroup / n) == sha(fresh / n) for n in ("library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv",
                                                              "library_free_second_places.tsv")}
    check("U2 regroup from saved pieces equals a fresh run", all(same.values()), same)
    pk = sorted((t / "d0_work" / "pieces").glob("*.pkl"))
    mt = [p.stat().st_mtime_ns for p in pk]

    # U3
    ro = subprocess.run([sys.executable, "-B", str(drv_path), str(t / "d4"), "7", json.dumps({"walk_out": False}), "--het"],
                        env=dict(env, EVITEA_PIECES_READONLY="1"), capture_output=True, text=True)
    check("U3 read-only refuses a run that would scan again", ro.returncode != 0 and "EVITEA_PIECES_READONLY" in (ro.stderr + ro.stdout), ro.stderr[-200:])
    check("U3 the saved pieces are untouched", mt == [p.stat().st_mtime_ns for p in sorted((t / "d0_work" / "pieces").glob("*.pkl"))] and len(pk) > 0)
    ok = subprocess.run([sys.executable, "-B", str(drv_path), str(t / "d5"), "7", keep, "--het"],
                        env=dict(env, EVITEA_PIECES_READONLY="1"), capture_output=True, text=True)
    check("U3 read-only allows a regroup arm", ok.returncode == 0, ok.stderr[-200:])

    # U4
    dd = t / "d6"
    r4 = run(dd, ROOT / "src", base, ["--het", "--detection-dir", str(d0)])
    check("U4 run from saved tables", r4.returncode == 0, r4.stderr[-400:])
    check("U4 same family calls as the run that wrote the tables", sha(dd / "library_free_family_calls.tsv") == sha(d0 / "library_free_family_calls.tsv"))

    # U12
    def fam_rows(d):
        lines = (d / "library_free_family_calls.tsv").read_text().splitlines()
        head = lines[0].split("\t")
        return [dict(zip(head, x.split("\t"))) for x in lines[1:]]
    tables = ("library_free_family_calls.tsv", "calls.library_free.tsv", "calls.combined.tsv")
    sc0 = summ["assembly"].get("single_copy", {})
    prof = t / "d0_work" / "copy_profile.json"
    check("U12 auto by default: profile measured and kept, nothing flagged without library calls",
          sc0.get("mode") == "auto" and sc0.get("flagged") == [] and prof.is_file()
          and abs((sc0.get("contigs", {}).get("chrS", {}).get("ratio") or 0) - 1) < 0.01, sc0)
    pm = prof.stat().st_mtime_ns
    u_none, u_auto, u_list = t / "u_none", t / "u_auto", t / "u_list"
    runs = {d: subprocess.run([sys.executable, "-B", str(drv_path), str(d), "7", base, "--het", "--single-copy", mode],
                              env=dict(env, EVITEA_PIECES_READONLY="1"), capture_output=True, text=True)
            for d, mode in ((u_none, "none"), (u_auto, "auto"), (u_list, "chrS"))}
    check("U12 none / auto / named contig run from the saved pieces", all(x.returncode == 0 for x in runs.values()),
          [x.stderr[-300:] for x in runs.values()])
    check("U12 none and a second auto run give the call tables of the default run",
          all(sha(u / n) == sha(d0 / n) for u in (u_none, u_auto) for n in tables))
    check("U12 the profile is reused", prof.stat().st_mtime_ns == pm)
    a, b = fam_rows(d0), fam_rows(u_list)
    changed = [(x, y) for x, y in zip(a, b) if x != y]
    only = all({k for k in x if x[k] != y[k]} == {"status"} and x["status"] == "ADDED" and y["status"] == "DROPPED_SINGLE_COPY"
               and (x["scan"] == "second" or float(x["e"] or "nan") > 0.15) for x, y in changed)
    missed = [x for x in a if x["status"] == "ADDED" and (x["scan"] == "second" or float(x["e"] or "nan") > 0.15)]
    check("U12 a named contig turns exactly its heterozygous-looking added calls into DROPPED_SINGLE_COPY",
          len(a) == len(b) and only and len(changed) == len(missed) > 0, (len(changed), len(missed)))
    ms = json.loads((u_list / "merge_view" / "library_free_summary.json").read_text())["assembly"].get("single_copy", {}) \
        if (u_list / "merge_view").is_dir() else {"flagged": ["chrS"], "source": "main view"}
    check("U12 the merge view takes the main view's contigs", ms.get("flagged") == ["chrS"] and ms.get("source") == "main view", ms)

    # U13
    hg0 = summ["assembly"].get("het_guard", {})
    check("U13 on by default, recorded in the summary (idle: fewer than 10 confirmed places)",
          hg0.get("mode") == "on" and hg0.get("active") is False and hg0.get("confirmed", 99) < 10, hg0)
    mh = json.loads((d0 / "merge_view" / "library_free_summary.json").read_text())["assembly"].get("het_guard", {}) \
        if (d0 / "merge_view").is_dir() else {"source": "main view"}
    check("U13 the merge view takes the main view's places", mh.get("source") == "main view", mh)
    h_on, h_off = t / "h_on", t / "h_off"
    hruns = {d: subprocess.run([sys.executable, "-B", str(drv_path), str(d), "7", base, "--het", "--het-guard", mode],
                               env=dict(env, EVITEA_PIECES_READONLY="1"), capture_output=True, text=True)
             for d, mode in ((h_on, "on"), (h_off, "off"))}
    check("U13 on / off run from the saved pieces", all(x.returncode == 0 for x in hruns.values()), [x.stderr[-300:] for x in hruns.values()])
    check("U13 on and off give the call tables of the default run on this sample", all(sha(u / n) == sha(d0 / n) for u in (h_on, h_off) for n in tables))
    off = json.loads((h_off / "library_free_summary.json").read_text())["assembly"].get("het_guard")
    check("U13 off is recorded", off == {"mode": "off"}, off)

    sys.path.insert(0, str(ROOT / "src"))
    from evitea import cli
    from evitea.library_free.options import Options, from_args
    parse = lambda *a: cli.build_parser().parse_args(["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o",
                                                      "--msbwt-index", "m", *a])
    check("U6 no switch = default switches", from_args(parse()) == Options())
    check("U6 --previous-behaviour = every rule off", from_args(parse("--previous-behaviour")).is_previous)
    o = from_args(parse("--previous-behaviour", "--walk-out", "on"))
    check("U6 a switch given with it still applies", o.walk_out and not o.anchored_extension and not o.is_previous)
    check("U6 het scan and check options", from_args(parse("--het-scan", "second-walk", "--check-scope", "path-dropped", "--tol", "p95")).het_scan == "second-walk")
    try:
        from_args(parse("--anchored-extension", "off", "--polya-cross", "on"))
        check("U6 polya-cross without anchored extension refused", False)
    except ValueError:
        check("U6 polya-cross without anchored extension refused", True)

    # U7
    (t / "map.bed").write_text("chrS\t0\t5000\n")
    (t / "annot.bed").write_text("chrS\t0\t10000\n")
    guard = lambda d: json.loads((d / "library_free_summary.json").read_text())["assembly"].get("resident_guard")
    rep = lambda d: json.loads((d / "library_free_summary.json").read_text()).get("repeat_map")
    m0, m1, m2, m3 = t / "m_none", t / "m_map", t / "m_annot", t / "m_both"
    r0 = run(m0, ROOT / "src", base)
    r1 = run(m1, ROOT / "src", base, ["--repeat-map", str(t / "map.bed")])
    r2 = run(m2, ROOT / "src", base, ["--annotation", str(t / "annot.bed")])
    r3 = run(m3, ROOT / "src", base, ["--annotation", str(t / "annot.bed"), "--repeat-map", str(t / "map.bed")])
    check("U7 runs with an empty library layer finish", all(x.returncode == 0 for x in (r0, r1, r2, r3)), [x.stderr[-300:] for x in (r0, r1, r2, r3)])
    check("U7 no library call, no annotation, no map: the test is off", isinstance(guard(m0), str) and "off" in guard(m0), guard(m0))
    check("U7 no library call: the map feeds the test (TE share 5000 / 80000)", isinstance(guard(m1), dict) and guard(m1)["te_share"] == 0.0625
          and rep(m1) == {"path": str(t / "map.bed"), "used": True, "library_calls": 0}, (guard(m1), rep(m1)))
    check("U7 an annotation alone gives its own share", isinstance(guard(m2), dict) and guard(m2)["te_share"] == 0.125 and rep(m2) is None, guard(m2))
    check("U7 an annotation wins over the map", isinstance(guard(m3), dict) and guard(m3)["te_share"] == 0.125 and rep(m3)["used"] is False, (guard(m3), rep(m3)))
    p0, p1 = t / "p_none", t / "p_map"
    q0 = run(p0, ROOT / "src", base, ["--native", "px"])
    q1 = run(p1, ROOT / "src", base, ["--native", "px", "--repeat-map", str(t / "map.bed")])
    check("U7 runs with library calls finish", q0.returncode == 0 and q1.returncode == 0, (q0.stderr[-300:], q1.stderr[-300:]))
    names = ["library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv", "calls.library_free.tsv", "calls.combined.tsv"]
    same = {n: sha(p0 / n) == sha(p1 / n) for n in names + ["library_free_second_places.tsv"] if (p0 / n).exists()}
    check("U7 with library calls the map changes nothing", all(same.values()) and len(same) >= 5 and guard(p0) == guard(p1)
          and rep(p1)["used"] is False and rep(p1)["library_calls"] > 0, (same, rep(p1)))

    # U8
    names8 = ["library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv", "calls.library_free.tsv", "calls.combined.tsv",
              "library_free_second_places.tsv", "library_free_second_segments.fa"]
    if BASELINE:
        ob, on_ = t / "base8", t / "new8"
        rb = run(ob, BASELINE, base, ["--het"])
        rn = run(on_, ROOT / "src", json.dumps({"identity_min_comparisons": 0, "relative_depth_floor": False}), ["--het"])
        check("U8 both runs finish", rb.returncode == 0 and rn.returncode == 0, (rb.stderr[-300:], rn.stderr[-300:]))
        same = {n: sha(ob / n) == sha(on_ / n) for n in names8 if (ob / n).exists()}
        check("U8 limit 0 and the floor off: the files equal the version before the switches", all(same.values()) and len(same) == 7, same)
    else:
        print("U8 skipped: BASELINE_SRC is not set")
    from evitea import cli as _cli
    from evitea.library_free.options import Options as _O, from_args as _fa
    _p = lambda *a: _cli.build_parser().parse_args(["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o", "--msbwt-index", "m", *a])
    o8 = _fa(_p())
    check("U8 defaults: limit 3, floor on", o8.identity_min_comparisons == 3 and o8.relative_depth_floor is True)
    o9 = _fa(_p("--identity-min-comparisons", "0", "--relative-depth-floor", "on"))
    check("U8 the command line reads both switches", o9.identity_min_comparisons == 0 and o9.relative_depth_floor is True)
    check("U8 --previous-behaviour: limit 0, floor off", _fa(_p("--previous-behaviour")).identity_min_comparisons == 0 and not _fa(_p("--previous-behaviour")).relative_depth_floor)
    try:
        _O(identity_min_comparisons=-1)
        check("U8 a negative limit is refused", False)
    except ValueError:
        check("U8 a negative limit is refused", True)
    check("U8 the floor is part of the detection identity, the limit is not", _O(relative_depth_floor=False).detection_identity != _O().detection_identity
          and _O(identity_min_comparisons=0).detection_identity == _O().detection_identity)

    fl = summ["detection"].get("depth_floors", {})
    check("U8 floor on (default): the summary lists each contig's floor (at most 10)", list(fl) == ["chrS"] and 0 < fl["chrS"] <= 10, fl)
    f0 = t / "floor_off"
    rf = run(f0, ROOT / "src", json.dumps({"relative_depth_floor": False}))
    sf = json.loads((f0 / "library_free_summary.json").read_text()) if rf.returncode == 0 else {}
    check("U8 floor off: no floor entry in the summary", rf.returncode == 0 and "depth_floors" not in sf.get("detection", {"depth_floors": 0}),
          (list(sf.get("detection", {})), rf.stderr[-300:]))

    # U9
    outs = ["library_free_places.tsv", "library_free_segments.fa", "library_free_family_calls.tsv", "calls.library_free.tsv",
            "calls.combined.tsv", "calls.combined.vcf", "library_free_second_places.tsv", "library_free_second_segments.fa"]
    samples = (("elements", ["--het"], 7), ("elements", [], 7), ("cross", ["--sample", "cross"], 11))
    if MERGE:
        for k, (kind, extra, seed) in enumerate(samples):
            mg, mo, cr = t / f"u9_merge{k}", t / f"u9_mopt{k}", t / f"u9_cross{k}"
            ra = run(mg, MERGE, base, extra, seed=seed)
            rb = run(mo, ROOT / "src", "merge", extra, seed=seed)
            rc = run(cr, ROOT / "src", base, extra, seed=seed)
            check(f"U9 {kind}{' ' + extra[0] if extra else ''}: the three runs finish", ra.returncode == rb.returncode == rc.returncode == 0,
                  (ra.stderr[-300:], rb.stderr[-300:], rc.stderr[-300:]))
            same = {n: sha(mg / n) == sha(mo / n) for n in outs}
            check(f"U9 {kind}{' ' + extra[0] if extra else ''}: Options.merge() writes the merge's files", all(same.values()), same)
            same = {n: sha(mg / n) == sha(cr / "merge_view" / n) for n in outs}
            check(f"U9 {kind}{' ' + extra[0] if extra else ''}: merge_view/ holds the merge's files", all(same.values()), same)
    else:
        print("U9 skipped: MERGE_SRC is not set")

    # U10
    import csv
    sys.path.insert(0, str(HERE))
    import synth_cross as SX
    xs = SX.build(seed=11)["sites"]
    xa, xd = t / "x_anchor", t / "x_default"
    ra = run(xa, ROOT / "src", base, ["--sample", "cross", "--profile-min", "1000", "--repeat-anchor", "150,1"], seed=11)
    rd = run(xd, ROOT / "src", base, ["--sample", "cross"], seed=11)
    check("U10 cross sample runs finish", ra.returncode == 0 and rd.returncode == 0, (ra.stderr[-400:], rd.stderr[-400:]))

    def at(d, kind, k=0, view=""):
        """The non-duplicate place with a left break nearest the k-th site of ``kind`` (within 300 letters) and its family call."""
        pos = [x["pos"] for x in xs if x["kind"] == kind][k]
        rows = list(csv.DictReader(open(d / view / "library_free_places.tsv"), delimiter="\t"))
        near = [p for p in rows if p["bL"] and p["dup"] != "1" and abs(int(p["bL"]) - pos) <= 300]
        p = min(near, key=lambda p: abs(int(p["bL"]) - pos)) if near else {}
        calls = [c for c in csv.DictReader(open(d / view / "library_free_family_calls.tsv"), delimiter="\t")
                 if c["bL"] and abs(int(c["bL"]) - pos) <= 300 and c["status"] != "NOT_HET"]
        return p, calls

    def called(calls):
        return any(c["status"] == "ADDED" for c in calls)

    m, x = at(xa, "DIV", view="merge_view")[0], at(xa, "DIV")[0]
    check("U10 shared-word depth: the merge abstains on depth in the diverged block, the rule passes the place",
          "depth" in m.get("why", "") and x.get("stage") == "S4_pass" and x.get("cross_rule") == "shared_depth" and float(x["D"]) >= 10,
          (m.get("why"), x.get("stage"), x.get("D")))
    m, x = at(xa, "LOWCOV", view="merge_view")[0], at(xa, "LOWCOV")[0]
    check("U10 shared-word depth: a block read at low depth stays out", "depth" in m.get("why", "") and "depth" in x.get("why", "")
          and not x.get("cross_rule"), (m.get("why"), x.get("why"), x.get("cross_rule")))
    (m, mc), (x, xc) = at(xa, "REP", view="merge_view"), at(xa, "REP")
    check("U10 repeat cutoff: the merge abstains on repeat content, the rule calls the place",
          m.get("why") == "rep" and x.get("stage") == "S4_pass" and x.get("cross_rule") == "repeat_cutoff" and called(xc) and not called(mc),
          (m.get("why"), x.get("stage"), xc))
    d0x = at(xd, "REP")[0]
    sx = lambda d: json.loads((d / "library_free_summary.json").read_text())["detection"]
    check("U10 repeat cutoff: with the anchor's own share the cutoff is REP_MAX and the place stays out",
          sx(xd)["repeat_cutoff"]["cutoff"] == 0.5 and d0x.get("why") == "rep" and not d0x.get("cross_rule"), (sx(xd)["repeat_cutoff"], d0x.get("why")))
    check("U10 repeat cutoff: the summary gives the cutoff from the windows", 0.5 < sx(xa)["repeat_cutoff"]["cutoff"] < 1
          and sx(xa)["repeat_cutoff"]["windows"] == 150, sx(xa)["repeat_cutoff"])
    for k in (0, 1):
        m, x = at(xa, "POLYA", k, view="merge_view")[0], at(xa, "POLYA", k)[0]
        check(f"U10 poly-A tails ({'+-'[k]} strand): the merge rejects the tail side as short, the rule crosses the run and keeps both sides",
              "short" in m.get("seg_reject", "") and x.get("seg_keep") == "L,R" and x.get("cross_rule") == "polya_tails"
              and "a_x_" in (x.get("lstop", "") + x.get("rstop", "")), (m.get("seg_reject"), x.get("seg_keep"), x.get("lstop"), x.get("rstop")))
    (m, mc), (x, xc) = at(xa, "RECUR", view="merge_view"), at(xa, "RECUR")
    check("U10 recurring sequence: the merge rejects the copy beside its old copy as local, the rule calls it with its family",
          m.get("seg_reject") == "L:local,R:local" and x.get("stage") == "S4_pass" and x.get("cross_rule") == "recurring_sequence"
          and called(xc) and not called(mc), (m.get("seg_reject"), x.get("stage"), xc))
    x = at(xa, "LOCALDUP")[0]
    check("U10 recurring sequence: a local copy of the flank stays rejected", x.get("stage") == "S4_reject" and not x.get("cross_rule"),
          (x.get("stage"), x.get("cross_rule")))
    hdr = (xa / "library_free_places.tsv").read_text().split("\n", 1)[0].split("\t")
    hdr_m = (xa / "merge_view" / "library_free_places.tsv").read_text().split("\n", 1)[0].split("\t")
    check("U10 the places table names the rule; the merge view's has no such column", hdr[-1] == "cross_rule" and "cross_rule" not in hdr_m)
    by = sx(xa)["cross_rules"]["places_by_rule"]
    check("U10 the summary counts the places of each rule", all(by.get(f"{r}:S4_pass", 0) >= 1 for r in
                                                                ("shared_depth", "repeat_cutoff", "polya_tails", "recurring_sequence")), by)
    xr = t / "x_reuse"
    rr = run(xr, ROOT / "src", base, ["--sample", "cross", "--detection-dir", str(xa)], seed=11)
    check("U10 a run from the saved tables gives both views' family calls", rr.returncode == 0
          and sha(xr / "library_free_family_calls.tsv") == sha(xa / "library_free_family_calls.tsv")
          and sha(xr / "merge_view" / "library_free_family_calls.tsv") == sha(xa / "merge_view" / "library_free_family_calls.tsv"), rr.stderr[-300:])

    # U11
    import numpy as np
    from evitea.library_free import cross as X
    reps = np.array([0.1] * 90 + [0.6] * 6 + [0.7, 0.8, 0.9, 0.95])
    cut, sm = X.repeat_cutoff(reps, {"windows": 1000, "at_cutoff": 100})
    check("U11 a share at REP_MAX no larger than the anchor's keeps REP_MAX", cut == 0.5 and sm["at_rep_max"] == 10, sm)
    cut, sm = X.repeat_cutoff(reps, {"windows": 1000, "at_cutoff": 30})
    check("U11 a larger share: the cutoff is the window at the anchor's share from the top (3rd of 100)", cut == 0.8 and sm["windows_at_cutoff"] == 3, sm)
    cut, _ = X.repeat_cutoff(np.array([0.2] * 99 + [0.51]), {"windows": 1000, "at_cutoff": 5})
    check("U11 the cutoff is never below REP_MAX", cut == 0.51, cut)
    saved, X.REPEAT_ANCHOR = X.REPEAT_ANCHOR, None
    try:
        X.repeat_cutoff(reps)
        check("U11 genome-share without the anchor's counts is refused", False)
    except ValueError:
        check("U11 genome-share without the anchor's counts is refused", True)
    finally:
        X.REPEAT_ANCHOR = saved
    off = from_args(parse("--shared-word-depth", "off", "--repeat-cutoff", "fixed", "--polya-tails", "off", "--recurring-sequence", "off",
                          "--merge-view", "off"))
    check("U11 the four rules and the merge view off on the command line = Options.merge()", off == Options.merge(), off)
    check("U11 --previous-behaviour turns the cross-species rules off", not from_args(parse("--previous-behaviour")).any_cross)
    check("U11 the cross-species rules are part of the detection identity", Options().detection_identity != Options.merge().detection_identity
          and Options(merge_view=False).detection_identity == Options().detection_identity)
    for bad, msg in ((dict(anchored_extension=False, polya_cross=False), "polya-tails"), (dict(check_scope="none"), "recurring-sequence")):
        try:
            Options(**bad)
            check(f"U11 --{msg} without what it needs is refused", False)
        except ValueError:
            check(f"U11 --{msg} without what it needs is refused", True)

    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


main()
