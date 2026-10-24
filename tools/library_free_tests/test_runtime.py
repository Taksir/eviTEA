"""The library-free run flow on a synthetic sample (rt_driver.py).   Run: python -B tests/test_runtime.py   (venv python)

U1  previous behaviour (options None) writes the previous version's files byte for byte (the previous package is run on the same sample).
U2  a regroup arm in the same work directory (only grouping / assembly switches differ, here also the check scope) reuses the saved pieces and gives the
    same places tables as a fresh run with those switches.
U3  EVITEA_PIECES_READONLY refuses a run whose detection switches differ (it would replace the saved pieces); the pieces stay.
U4  --detection-dir (saved tables) gives the same family calls as the run that wrote them.
U6  the command line: no switch gives the default switches, --previous-behaviour gives every rule off, a switch given with it still applies,
    polya-cross without anchored extension is refused.
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
# PREVIOUS_SRC: the src folder of the previous version (a checkout of the checkpoint branch); U1 is skipped without it.
PREVIOUS = os.environ.get("PREVIOUS_SRC")
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

    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


main()
