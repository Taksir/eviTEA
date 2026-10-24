"""Tests of the depth window (upper bound of the local-depth window).  Run: PYTHONPATH=src python tests/test_depth_window.py
Set PREVIOUS_SRC to the previous version's src folder to also compare against the previous package byte for byte."""
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(os.path.dirname(HERE)), "src")
sys.path.insert(0, HERE)

import synth  # noqa: E402
from evitea import depth_window  # noqa: E402
from evitea.breakpoints import walker  # noqa: E402
from evitea.library_free import _p4_lib as P  # noqa: E402
from evitea.library_free import detector as D  # noqa: E402

RESULTS = {}


def probe(src, depth, cap=None):
    env = dict(os.environ, PYTHONPATH=src, PYTHONDONTWRITEBYTECODE="1")
    env.pop(depth_window.ENV_VAR, None)
    args = [sys.executable, "-B", os.path.join(HERE, "probe_scan.py"), str(depth)] + ([cap] if cap is not None else [])
    return json.loads(subprocess.run(args, env=env, check=True, capture_output=True, text=True).stdout)


def t_parse():
    assert depth_window.parse_depth_cap(None) == math.inf
    assert depth_window.parse_depth_cap("none") == math.inf
    assert depth_window.parse_depth_cap("100") == 100.0
    assert depth_window.parse_depth_cap(250) == 250.0
    for bad in ("5", 9, "abc", "-3"):
        try:
            depth_window.parse_depth_cap(bad)
        except ValueError:
            continue
        raise AssertionError(bad)
    RESULTS["D1"] = "cap parsing: none/100/250 accepted; below the lower bound 10 and non-numbers refused"


def t_set_and_modules():
    assert walker.D_MAX == math.inf and P.D_MAX == math.inf and walker.D_MIN == P.D_MIN == 10
    depth_window.set_depth_cap(100)
    try:
        assert walker.D_MAX == 100.0 and P.D_MAX == 100.0 and os.environ[depth_window.ENV_VAR] == "100"
    finally:
        depth_window.set_depth_cap(None)
    assert walker.D_MAX == math.inf and P.D_MAX == math.inf and depth_window.ENV_VAR not in os.environ
    RESULTS["D2"] = "default has no upper bound, lower bound 10 in both modules; set_depth_cap(100) sets both and the environment, None undoes it"


def t_env_inherited():
    code = "from evitea.breakpoints import walker; from evitea.library_free import _p4_lib as P; print(walker.D_MAX, P.D_MAX)"
    env = dict(os.environ, PYTHONPATH=SRC, EVITEA_DEPTH_CAP="100", PYTHONDONTWRITEBYTECODE="1")
    out = subprocess.run([sys.executable, "-B", "-c", code], env=env, check=True, capture_output=True, text=True).stdout.split()
    assert out == ["100.0", "100.0"], out
    RESULTS["D3"] = "a worker process started with EVITEA_DEPTH_CAP=100 sees the bound in both modules"


def t_scan_in_window():
    """Depth inside the old window: places and segments byte-identical to the previous package, with and without a cap."""
    new_none, new_100 = probe(SRC, 30), probe(SRC, 30, "100")
    assert new_none["places_sha"] == new_100["places_sha"] and new_none["segments_sha"] == new_100["segments_sha"]
    assert new_none["stages"]["S4_pass"] == 5
    previous = os.environ.get("PREVIOUS_SRC")
    if previous:
        old = probe(previous, 30)
        assert (old["places_sha"], old["segments_sha"]) == (new_none["places_sha"], new_none["segments_sha"]), (old, new_none)
        RESULTS["D4"] = "depth 30 (D median 23): places and segments identical to the previous version with no cap and with --depth-cap 100 (5 S4_pass places)"
    else:
        RESULTS["D4"] = "depth 30: identical with no cap and with --depth-cap 100 (previous-version comparison skipped: PREVIOUS_SRC unset)"


def t_scan_above_window():
    """Depth 150 (D median 115): the old window abstains at every candidate, the new default calls them."""
    capped, free = probe(SRC, 150, "100"), probe(SRC, 150, "none")
    assert capped["stages"].get("S4_pass", 0) == 0 and capped["why"] == {"depth": 5}, capped
    assert free["stages"]["S4_pass"] == 5 and not free["why"], free
    low = probe(SRC, 30)
    assert free["segments_sha"] == low["segments_sha"], "the inserted DNA found at 150x differs from the 30x sequences"
    previous = os.environ.get("PREVIOUS_SRC")
    if previous:
        old = probe(previous, 150)
        assert (old["places_sha"], old["stages"]) == (capped["places_sha"], capped["stages"])
    RESULTS["D5"] = ("depth 150 (D median 115): --depth-cap 100 abstains on all 5 candidates (reason depth, identical to the previous version); "
                     "no cap gives 5 S4_pass places and the same inserted sequences as the 30x run")


def t_lower_bound():
    low = probe(SRC, 8, "none")
    assert low["stages"].get("S4_pass", 0) == 0, low
    RESULTS["D6"] = f"depth 8 (D median {low['d_median']}): no S4_pass place with no cap; the lower bound of 10 stays"


class _Half:
    """Single-strand view of a both-strand HashCounter (KmerCounter adds the reverse complement itself)."""

    def __init__(self, hc):
        self.hc = hc

    def __call__(self, seq):
        return self.hc.cnt(seq) / 2


def t_walker_window():
    from evitea.breakpoints.walker import KmerCounter, classify
    got = {}
    for depth in (30, 150):
        ref, sites, C = synth.build(depth=depth, glen=24000, sites=((3000, 1.0), (6000, 1.0), (9000, 1.0), (12000, 1.0), (15000, 1.0)))
        for cap in (100, None):
            depth_window.set_depth_cap(cap)
            try:
                r = classify(KmerCounter(_Half(C)), ref, 11000, 11010)
            finally:
                depth_window.set_depth_cap(None)
            got[(depth, cap)] = (r["cls"], r.get("why"))
    assert got[(30, 100)] == got[(30, None)] and got[(30, None)][0] != "UNINFORMATIVE", got
    assert got[(150, 100)] == ("UNINFORMATIVE", "depth"), got
    assert got[(150, None)][0] != "UNINFORMATIVE" or got[(150, None)][1] != "depth", got
    RESULTS["D7"] = f"breakpoint walk classify(): {got}"


def t_cli():
    from evitea import cli
    parser = cli.build_parser()
    for cmd, extra in (("call", []), ("refine-breakpoints", ["--run-dir", "r", "--msbwt-index", "m", "--reference-assets", "a",
                                                           "--te-assets", "b", "--output", "o"]),
                       ("library-free", ["--run-dir", "r", "--reference-assets", "a", "--output", "o"])):
        assert parser.parse_args([cmd] + extra).depth_cap == "none"
        assert parser.parse_args([cmd] + extra + ["--depth-cap", "100"]).depth_cap == "100"
    rc = cli.main(["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o", "--depth-cap", "5"])
    assert rc == 2, rc
    RESULTS["D8"] = "call, refine-breakpoints and library-free take --depth-cap (default none); a value below 10 exits with code 2"


def main():
    for t in (t_parse, t_set_and_modules, t_env_inherited, t_scan_in_window, t_scan_above_window, t_lower_bound, t_walker_window, t_cli):
        t()
    for k, v in RESULTS.items():
        print(k, v)
    print("ALL PASS")


if __name__ == "__main__":
    main()
