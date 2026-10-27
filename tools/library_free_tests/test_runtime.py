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

    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


main()
