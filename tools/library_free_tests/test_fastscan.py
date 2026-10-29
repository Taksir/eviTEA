"""The scan-speed options of the library-free layer (--workers, --counter, --reference-cache).   Run: python -B tools/library_free_tests/test_fastscan.py   (venv python)

They change how fast a contig is scanned and nothing it yields. EXPECTED_CODE_HASH (a prefix of the detection code hash, for example the hash of the
version to keep calls identical with) makes H1 check it; the other checks need nothing but this tree. Checks that need pymsbwt are skipped without it.

C1  the fast counter equals the plain counts on every word asked of an exact index with interval searches (a suffix-array index of a read list):
    reference words, words with one letter changed, absent words, words of every length 1-60 and 100, walks, words with a letter other than A, C, G, T;
    first look (use False) and walks (use True); with caps small enough to empty both caches many times.
C2  an index without interval searches (the test fakes of rt_driver.py) gives the plain counts through the same class.
C3  the same on a real pymsbwt index.
C4  the self-check passes, counts the words it checked, and stops on a counter that is wrong by one on one word.
R1  reference tables: a cached folder holds exactly the freshly built arrays; RefMulti and RefWords answer the same on 10,000 random words, read through
    a memory map or from memory; a second call loads (no build); a folder written without the walk-out words gets them when a run needs them;
    an empty table (no multi-copy word) loads.
P1  scan_contig_tiled equals scan_contig_passes by value (every place, field, the info block and its key order; NaN equals NaN) on the synthetic
    samples, with 2 and 3 workers and tiles cut through candidate runs, for the default switches, the merge's, the previous behaviour's, the
    haplotype scan off and a sample with the cross-species rules on, with the plain counter, the fast counter (fallback and interval searches).
P2  the counter, tile and chunk sizes and the worker count change nothing else: workers 1 is scan_contig_passes itself; a worker error is raised.
W1  whole runs (rt_driver.py): workers 1 and 3 with both counters, on the exact-table index and on a real index, write every output file of both
    views byte for byte the same; the run summary is equal after taking out the timings and its speed entry.
W4  against the version without these files (FINAL_SRC = its src folder, for example a checkout of the build to keep calls identical with; skipped
    without it): its run and the new code's runs (fresh scan with 3 workers and the fast counter; reuse of its pieces; its reuse of the new
    pieces) write the same output files, and the saved pieces are equal by value.
W2  saved pieces are shared: pieces written by one setting are reused (no scan, EVITEA_PIECES_READONLY set) by every other setting, both ways.
S1  command line: library-free reads --workers, --counter and --reference-cache; `evitea run` gives the library-free step --workers from --threads;
    detect.sbatch gives the options (not to the previous version).
H1  the detection code hash (the identity of a saved piece) is unchanged by these files.
"""
import hashlib
import json
import math
import os
import pickle
import random
import subprocess
import sys
import tempfile
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))
import synth_cross as SX  # noqa: E402
import synth_elements as SE  # noqa: E402
import numpy as np  # noqa: E402
from evitea.library_free import _p4_lib as P  # noqa: E402
from evitea.library_free import detector as D  # noqa: E402
from evitea.library_free import fastcount as FC  # noqa: E402
from evitea.library_free import parallel as PA  # noqa: E402
from evitea.library_free import refcache as RC  # noqa: E402
from evitea.library_free import runtime as R  # noqa: E402
from evitea.library_free import scan as S  # noqa: E402
from evitea.library_free.options import Options  # noqa: E402
from evitea.library_free.walks import RefWords  # noqa: E402

try:
    from MUSCython import MultiStringBWTCython as MB
except ImportError:
    MB = None

FAIL, NCHK = [], [0]
TMP = Path(tempfile.mkdtemp(prefix="fastscan_"))


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail)[:600] if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def equal(a, b, path="x", out=None):
    """Equal by value: types, dict key order, list lengths, NaN equals NaN. Returns the first difference or None."""
    if type(a) is not type(b):
        return f"type {path}: {type(a).__name__} {type(b).__name__}"
    if isinstance(a, dict):
        if list(a) != list(b):
            return f"keys {path}: {list(a)[:12]} {list(b)[:12]}"
        for k in a:
            d = equal(a[k], b[k], f"{path}.{k}")
            if d:
                return d
        return None
    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return f"length {path}: {len(a)} {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            d = equal(x, y, f"{path}[{i}]")
            if d:
                return d
        return None
    if isinstance(a, np.ndarray):
        return None if np.array_equal(a, b) else f"array {path}"
    if isinstance(a, float) and math.isnan(a) and math.isnan(b):
        return None
    return None if a == b else f"value {path}: {a!r} {b!r}"


# ---------------------------------------------------------------- an exact index with interval searches
class SuffixIndex:
    """Backward search over the sorted suffixes of '$'-terminated reads, with pymsbwt's calls: countOccurrencesOfSeq and findIndicesOfStr
    (seq, givenRange as a tuple). Small, for tests."""

    def __init__(self, reads):
        sufs = []
        for r in reads:
            t = r + "$"
            sufs += [(t[i:], t[i - 1] if i else "$") for i in range(len(t))]
        sufs.sort(key=lambda x: x[0])
        self.total = len(sufs)
        bwt = [c for _, c in sufs]
        order = "$ACGT"
        self.fm, acc = {}, 0
        for c in order:
            self.fm[c] = acc
            acc += sum(1 for s, _ in sufs if s[0] == c)
        self.occ = {c: np.concatenate([[0], np.cumsum(np.array([x == c for x in bwt]))]) for c in order}
        self.steps = 0      # letters searched

    def _step(self, c, i):
        return self.fm[c] + int(self.occ[c][i])

    def findIndicesOfStr(self, seq, givenRange=None):
        if givenRange is not None and not isinstance(givenRange, tuple):
            raise TypeError("givenRange must be a tuple")
        l, h = (0, self.total) if givenRange is None else givenRange
        for c in reversed(seq):
            l, h = self._step(c, l), self._step(c, h)
            self.steps += 1
        return (l, h)

    def countOccurrencesOfSeq(self, seq, givenRange=None):
        l, h = (0, self.total) if givenRange is None else givenRange
        for c in reversed(seq):
            l, h = self._step(c, l), self._step(c, h)
            self.steps += 1
            if l == h:
                break
        return h - l


def random_reads(seed=3, glen=3000, n=900, rl=70):
    r = random.Random(seed)
    g = [r.choice("ACGT") for _ in range(glen)]
    block = g[200:330]
    for at in (900, 1500, 2200):                   # a repeat, so that some words have several occurrences
        g[at:at + len(block)] = block
    g = "".join(g)
    reads = []
    for _ in range(n):
        a = r.randrange(0, glen - rl)
        w = g[a:a + rl]
        if r.random() < 0.3:
            i = r.randrange(rl)
            w = w[:i] + r.choice([x for x in "ACGT" if x != w[i]]) + w[i + 1:]
        reads.append(w if r.random() < 0.5 else P.rc(w))
    return g, reads


def words_to_ask(g, seed=5):
    r = random.Random(seed)
    out = []
    for _ in range(400):
        a = r.randrange(0, len(g) - 101)
        for n in (1, 2, 5, 14, 18, 19, 25, 30, 31, 41, 60, 100):
            out.append(g[a:a + n])
    for w in list(out[::7]):
        i = r.randrange(len(w))
        out.append(w[:i] + r.choice("ACGT") + w[i + 1:])
    out += ["".join(r.choice("ACGT") for _ in range(31)) for _ in range(150)]            # absent
    out += ["ACGTN" + g[10:36], "N", "acgt", "A" * 31, "T" * 30, "AC" * 15, "AT" * 15]    # other letters, palindrome-like words
    return out + walk_words(g)


def walk_words(g, seed=6):
    """The words of walks: the four letters after a context, then the next context (the order of the walk)."""
    r = random.Random(seed)
    out = []
    for _ in range(60):
        a = r.randrange(0, len(g) - 80)
        ctx = g[a:a + 30]
        for i in range(a + 30, a + 70):
            out += [ctx + g[i], ctx] + [ctx + b for b in "ACGT"]
            ctx = ctx[1:] + g[i]
    return out


def counter_tests(label, bwt, g, expect_intervals, **limits):
    plain = FC.plain_count
    asked = words_to_ask(g)
    for use in (False, True):
        fast = FC.FastCounter(bwt, use, **limits)
        bad = next(((w, fast.cnt(w), plain(bwt, w)) for w in asked if fast.cnt(w) != plain(bwt, w)), None)
        check(f"{label}: {len(asked)} words, use={use}{' ' + str(limits) if limits else ''}", bad is None, bad)
        if use is True and fast.intervals != expect_intervals:
            check(f"{label}: interval searches {'used' if expect_intervals else 'not used'}", False, fast.intervals)
    again = FC.FastCounter(bwt, True, **limits)
    first = [again.cnt(w) for w in asked]
    check(f"{label}: a second ask of every word (caches) gives the same", first == [again.cnt(w) for w in asked])
    return again


def counter_checks():
    g, reads = random_reads()
    idx = SuffixIndex(reads)
    counter_tests("C1 suffix index", idx, g, True)
    counter_tests("C1 suffix index, tiny caps", idx, g, True, count_cap=40, interval_cap=90)
    walk = walk_words(g)
    idx.steps = 0
    fast = FC.FastCounter(idx, True)
    for w in walk:
        fast.cnt(w)
    fast_steps, idx.steps = idx.steps, 0
    plain = FC.FastCounter(idx, True)
    plain.intervals = False
    for w in walk:
        plain.cnt(w)
    check("C1 on the words of walks the fast counter searches fewer letters than the plain counter (strand-free cache kept in both)",
          fast_steps < 0.75 * idx.steps, (fast_steps, idx.steps))
    no_find = types.SimpleNamespace(countOccurrencesOfSeq=idx.countOccurrencesOfSeq)
    counter_tests("C2 index without interval searches", no_find, g, False)

    class Untyped(SuffixIndex):                                   # an older build: givenRange not accepted
        def findIndicesOfStr(self, seq, givenRange=None):
            if givenRange is not None:
                raise TypeError("untyped")
            return super().findIndicesOfStr(seq)

    old = Untyped.__new__(Untyped)
    old.__dict__.update(idx.__dict__)
    counter_tests("C2 index whose interval search raises TypeError", old, g, False)

    if MB is None:
        print("SKIP C3 (pymsbwt not installed)")
    else:
        import logging
        d = TMP / "real_idx"
        d.mkdir()
        MB.createMSBWTFromSeqs([r + "$" for r in reads], str(d), 1, True, logging.getLogger("t"))
        real = MB.loadBWT(str(d), useMemmap=True)
        counter_tests("C3 real pymsbwt index", real, g, True)
        counter_tests("C3 real pymsbwt index, tiny caps", real, g, True, count_cap=40, interval_cap=90)
        check("C3 the real index agrees with the suffix index on 300 words",
              all(FC.plain_count(real, w) == FC.plain_count(idx, w) for w in words_to_ask(g)[::40]))
    genome = {"c1": g, "c2": g[::-1]}
    info = FC.self_check(idx, genome, ["c1", "c2"], n=300)
    check("C4 the self-check passes and counts its words", info["words"] >= 300 and info["intervals"] is True, info)
    info = FC.self_check(no_find, genome, ["c1", "c2"], n=200)
    check("C4 the self-check on an index without interval searches", info["words"] >= 200 and info["intervals"] is False, info)

    class Wrong(SuffixIndex):                                     # the extension by the letter G is one too many (the interval route only)
        def findIndicesOfStr(self, seq, givenRange=None):
            l, h = super().findIndicesOfStr(seq, givenRange)
            return (l, h + 1) if seq == "G" and givenRange is not None else (l, h)

    wrong = Wrong.__new__(Wrong)
    wrong.__dict__.update(idx.__dict__)
    try:
        FC.self_check(wrong, {"c1": g}, ["c1"], n=4000, seed=1)
        check("C4 the self-check stops on a wrong count", False)
    except RuntimeError as exc:
        check("C4 the self-check stops on a wrong count and names the word", "self-check failed" in str(exc), str(exc)[:200])


# ---------------------------------------------------------------- reference tables
def refcache_checks():
    r = random.Random(9)
    unit = "".join(r.choice("ACGT") for _ in range(400))
    genome = {}
    for name, n in (("a", 6000), ("b", 4000), ("c", 1500)):
        s = list("".join(r.choice("ACGT") for _ in range(n)))
        for at in range(100, n - 500, 700):                  # copies of one unit, some on the other strand
            u = unit if r.random() < 0.5 else P.rc(unit)
            s[at:at + 400] = list(u)
        s[50:60] = list("NNNNNNNNNN")
        genome[name] = "".join(s)
    cache = TMP / "refcache"
    logs = []
    rm0 = D.reference_multi(genome)
    words0 = RefWords(genome)
    rm1, w1 = RC.load_or_build(genome, reference_sha256="r" * 64, code_sha256="c" * 64, cache_dir=cache, need_words=False, progress=logs.append)
    key = RC.cache_key("r" * 64, "c" * 64)
    names = sorted(p.name for p in (cache / key).iterdir())
    check("R1 first call builds and saves m21, m31 and about.json (no walk-out words asked)", names == ["about.json", "m21.npy", "m31.npy"] and w1 is None
          and not any(p.name.startswith(key + ".tmp") for p in cache.iterdir()), (names, logs))
    check("R1 the saved tables equal the built ones", all(np.array_equal(np.load(cache / key / f), arr) for f, arr in
                                                          (("m21.npy", rm0.m[P.REP_K]), ("m31.npy", rm0.m[P.K]))))
    check("R1 there are multi-copy words (the test is not empty)", len(rm0.m[P.REP_K]) > 0 and len(rm0.m[P.K]) > 0)
    logs.clear()
    rm2, w2 = RC.load_or_build(genome, reference_sha256="r" * 64, code_sha256="c" * 64, cache_dir=cache, need_words=True, progress=logs.append)
    check("R1 a later call loads the multi-copy words and adds the walk-out words", w2 is not None and (cache / key / "words31.npy").is_file()
          and any("loaded" in m for m in logs) and np.array_equal(np.asarray(w2.codes), words0.codes), logs)
    logs.clear()
    rm3, w3 = RC.load_or_build({}, reference_sha256="r" * 64, code_sha256="c" * 64, cache_dir=cache, need_words=True, progress=logs.append)
    check("R1 a third call loads everything (the genome is not read: an empty one is given)", np.array_equal(np.asarray(w3.codes), words0.codes)
          and isinstance(rm3.m[P.K], np.memmap) and isinstance(w3.codes, np.memmap) and not any("saved" in m for m in logs), (logs, type(w3.codes)))
    words = [g[i:i + n] for g in genome.values() for i in range(0, len(g) - 40, 13) for n in (21, 31)]
    words += ["".join(r.choice("ACGT") for _ in range(n)) for n in (21, 31) for _ in range(3000)]
    words += [unit[i:i + 31] for i in range(0, 300, 7)] + [P.rc(unit)[i:i + 21] for i in range(0, 300, 7)]
    words = [w for w in words if len(w) in (21, 31)]
    r.shuffle(words)
    words = words[:10000]
    same = all(rm0.word_multi(w) == rm3.word_multi(w) for w in words)
    check(f"R1 RefMulti.word_multi gives the same answer on {len(words)} words, mapped or built", same and any(rm0.word_multi(w) for w in words))
    check("R1 RefMulti.rep_frac gives the same answer", all(rm0.rep_frac(genome[c], x) == rm3.rep_frac(genome[c], x)
                                                           for c in genome for x in range(0, len(genome[c]), 211)))
    check("R1 RefWords.has gives the same answer", all(words0.has(w) == w3.has(w) for w in words if len(w) == 31) and
          any(w3.has(w) for w in words if len(w) == 31))
    codes = np.array([int(x) for x in r.sample(range(1 << 40), 500)] + [int(c) for c in rm0.m[P.K][:300]], dtype=np.uint64)
    check("R1 RefMulti._isin gives the same answer on codes", np.array_equal(rm0._isin(codes, P.K), rm3._isin(codes, P.K))
          and np.array_equal(rm0._isin(codes, P.REP_K), rm3._isin(codes, P.REP_K)))
    flat = {"x": "".join(r.choice("ACGT") for _ in range(2000))}
    rm4, w4 = RC.load_or_build(flat, reference_sha256="s" * 64, code_sha256="c" * 64, cache_dir=cache, need_words=True, progress=lambda m: None)
    rm5, w5 = RC.load_or_build({}, reference_sha256="s" * 64, code_sha256="c" * 64, cache_dir=cache, need_words=True, progress=lambda m: None)
    check("R1 a reference with no multi-copy word: the empty table loads and answers False", len(rm5.m[P.K]) == 0 and not rm5.word_multi(flat["x"][:31])
          and np.array_equal(np.asarray(w5.codes), RefWords(flat).codes))
    check("R1 another reference or code is another folder", RC.cache_key("s" * 64, "c" * 64) != key != RC.cache_key("r" * 64, "d" * 64))


# ---------------------------------------------------------------- tiles
def make_fixture(kind):
    if kind == "elements":
        S_ = SE.build(seed=7, ref_len=36_000, n_px=3, n_lx=2, n_uq=1, n_st=0)
        reads = S_["reads"]
    elif kind == "het":
        S_ = SE.build(seed=7, ref_len=36_000, n_px=3, n_lx=2, n_uq=1, n_st=0)
        r = random.Random(8)
        reads = []
        for _ in range(len(S_["reads"])):
            src = S_["sample"] if r.random() < 0.5 else S_["ref"]
            a = r.randrange(0, len(src) - 125)
            w = src[a:a + 125]
            reads.append(w if r.random() < 0.5 else P.rc(w))
    else:
        S_ = SX.build(seed=11, n_var=6)
        reads = S_["reads"]
    ref = S_["ref"]
    fa = TMP / f"{kind}.fa"
    fa.write_text(">chrS\n" + ref + "\n")
    return dict(ref=ref, reads=reads, fa=str(fa), genome={"chrS": ref}, rm=None, words=None)


def tile_parity(kind, optname, opts, fixture, configs, indexes):
    ref = fixture["ref"]
    if fixture["rm"] is None:
        fixture["rm"] = D.reference_multi(fixture["genome"])
        fixture["words"] = RefWords(fixture["genome"])
        fixture["fw"] = P.MappyFirewall(fixture["fa"], None)
    rm, words, fw = fixture["rm"], fixture["words"], fixture["fw"]
    words = words if opts.walk_out else None
    key = (kind, optname)
    serial = S.scan_contig_passes(P.Counter(indexes["plain"], cache=False), ref, "chrS", rm, fw, opts, words, P.REP_MAX)
    want = (serial[0], serial[1], {k: v for k, v in serial[2].items() if k not in ("walk_log", "cross_log")})
    check(f"P1 {kind}/{optname}: the scan makes places ({len(serial[0])} first pass, {len(serial[1] or [])} second)", len(serial[0]) > 0, key)
    cands = [p for p in serial[0] if "first" in p]
    for index, counter, workers, tile, chunk in configs:
        bwt = indexes[index]
        make = ((lambda use: FC.FastCounter(bwt, use)) if counter == "fast" else (lambda use: P.Counter(bwt, cache=use)))
        t = tile if tile else (cands[0]["first"] + 1 if cands else 5000)
        got = PA.scan_contig_tiled(make, ref, "chrS", rm, fw, opts, words, P.REP_MAX, workers=workers, tile=t, s1_chunk=chunk)
        edges = sum(1 for p in cands if p["first"] // t != p["last"] // t)
        check(f"P1 {kind}/{optname}: {index} index, {counter} counter, {workers} workers, tile {t}, {edges} places across a tile edge",
              equal(want, got) is None, equal(want, got))
    return serial


def tile_checks():
    fixtures = {k: make_fixture(k) for k in ("elements", "het", "cross")}
    for kind, fx in fixtures.items():
        C = SE.SubCounter(fx["reads"])
        fx["bwts"] = {"plain": types.SimpleNamespace(countOccurrencesOfSeq=C._occ)}
        if MB is not None:
            import logging
            d = TMP / f"real_{kind}"
            d.mkdir()
            MB.createMSBWTFromSeqs([w.upper() + "$" for w in fx["reads"]], str(d), 1, True, logging.getLogger("t"))
            fx["bwts"]["real"] = MB.loadBWT(str(d), useMemmap=True)
    plan = [
        ("elements", "default", Options(), [("plain", "plain", 2, 13000, 1500), ("plain", "fast", 3, 5000, 700)] +
         ([("real", "fast", 3, 5000, 700), ("real", "plain", 2, 0, 5000)] if MB else [])),
        ("het", "default", Options(), [("plain", "fast", 3, 5000, 900), ("plain", "plain", 2, 0, 3000)] + ([("real", "fast", 2, 13000, 900)] if MB else [])),
        ("elements", "merge", Options.merge(), [("plain", "plain", 2, 13000, 1500), ("plain", "fast", 3, 5000, 700)]),
        ("het", "previous", Options.previous(), [("plain", "plain", 3, 5000, 900), ("plain", "fast", 2, 13000, 900)]),
        ("het", "het scan off", Options(het_scan="off"), [("plain", "plain", 3, 5000, 900), ("plain", "fast", 2, 0, 900)]),
        ("cross", "cross rules on", Options(), [("plain", "plain", 3, 13000, 3000), ("plain", "fast", 2, 5000, 2000)] + ([("real", "fast", 3, 0, 3000)] if MB else [])),
    ]
    for kind, name, opts, configs in plan:
        fx = fixtures[kind]
        tile_parity(kind, name, opts, fx, configs, fx["bwts"])
    # P2
    fx = fixtures["elements"]
    bwt = fx["bwts"]["plain"]
    fx_rm, fx_fw = fx["rm"], fx["fw"]
    one = PA.scan_contig_tiled(lambda use: P.Counter(bwt, cache=use), fx["ref"], "chrS", fx_rm, fx_fw, Options(), fx["words"], P.REP_MAX, workers=1)
    direct = S.scan_contig_passes(P.Counter(bwt, cache=False), fx["ref"], "chrS", fx_rm, fx_fw, Options(), fx["words"], P.REP_MAX)
    check("P2 workers 1 returns what scan_contig_passes returns (walk_log and cross_log included)", "walk_log" in one[2] and equal(one, direct) is None)

    def bad_counter(use):
        if use:
            raise RuntimeError("worker failed on purpose")
        return P.Counter(bwt, cache=False)

    try:
        PA.scan_contig_tiled(bad_counter, fx["ref"], "chrS", fx_rm, fx_fw, Options(), fx["words"], P.REP_MAX, workers=2, tile=5000)
        check("P2 a worker error is raised", False)
    except Exception as exc:  # noqa: BLE001
        check("P2 a worker error is raised (no pieces from a failed pool)", "on purpose" in str(exc), repr(exc)[:200])
    check("P2 the shared state is emptied after a run, also after an error", PA._SHARED == {})
    timing = {}
    PA.scan_contig_tiled(lambda use: P.Counter(bwt, cache=use), fx["ref"], "chrS", fx_rm, fx_fw, Options(), fx["words"], P.REP_MAX, workers=2, tile=13000,
                         timing=timing)
    check("P2 timing gives the seconds of the first look and the walks", set(timing) == {"first_look", "walks"} and all(v > 0 for v in timing.values()), timing)
    a = {"x": {"n": 1, "mode": "m"}, "k": 2}
    try:
        PA._join_stats(a, {"x": {"n": 1, "mode": "other"}, "k": 3})
        check("P2 tile stats that disagree on a copied key are refused", False)
    except RuntimeError:
        check("P2 tile stats that disagree on a copied key are refused", True)
    check("P2 stats are summed with key order kept", list(PA._join_stats({"b": 1, "a": {"z": 1, "y": 2}}, {"b": 2, "a": {"z": 3, "y": 4}})) == ["b", "a"]
          and PA._join_stats({"b": 1, "a": {"z": 1, "y": 2}}, {"b": 2, "a": {"z": 3, "y": 4}}) == {"b": 3, "a": {"z": 4, "y": 6}})


# ---------------------------------------------------------------- whole runs
def driver(out, extra=(), sample_args=("--het",), env=None, seed=7):
    return driver_at(str(ROOT / "src"), out, extra, sample_args, env, seed)


def driver_at(src, out, extra=(), sample_args=("--het",), env=None, seed=7):
    """rt_driver.py with the package of the src folder ``src`` on the path."""
    e = dict(os.environ, PYTHONPATH=str(src))
    e.update(env or {})
    return subprocess.run([sys.executable, "-B", str(HERE / "rt_driver.py"), str(out), str(seed), "{}", "--small", *sample_args, *extra], env=e,
                          capture_output=True, text=True)


def scrub(x):
    """A run summary without the timings, the memory, the speed entry and the paths of the run folders."""
    drop = {"seconds", "rule_c_seconds", "peak_rss_gb", "speed", "scan_seconds", "source_run", "output"}
    if isinstance(x, dict):
        return {k: scrub(v) for k, v in x.items() if k not in drop}
    if isinstance(x, list):
        return [scrub(v) for v in x]
    return x


def folder_state(out):
    files, summaries = {}, {}
    for p in sorted(Path(out).rglob("*")):
        if p.is_file():
            rel = str(p.relative_to(out))
            if p.name == R.SUMMARY:
                summaries[rel] = scrub(json.loads(p.read_text()))
            else:
                files[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return files, summaries


def run_checks():
    idx_fake = TMP / "idx_fake"
    idx_real = TMP / "idx_real"
    settings = [("workers1_plain", ["--workers", "1", "--counter", "plain"]), ("workers3_fast", ["--workers", "3", "--counter", "fast"]),
                ("workers3_plain", ["--workers", "3", "--counter", "plain"]), ("workers2_fast", ["--workers", "2", "--counter", "fast"])]
    kinds = [("fake", ["--index-dir", str(idx_fake)])] + ([("real", ["--index", "real", "--index-dir", str(idx_real)])] if MB else [])
    for kind, kextra in kinds:
        results = {}
        for name, extra in settings:
            out = TMP / f"w1_{kind}_{name}"
            r = driver(out, extra + kextra + ["--reference-cache", str(TMP / f"cache_{kind}")])
            check(f"W1 {kind} index, {name}: the run finishes", r.returncode == 0, r.stderr[-500:])
            if r.returncode == 0:
                results[name] = folder_state(out)
        base = results.get("workers1_plain")
        for name, state in results.items():
            if base is None or name == "workers1_plain":
                continue
            check(f"W1 {kind} index, {name}: all {len(state[0])} output files of both views equal those of workers 1 / plain",
                  state[0] == base[0] and len(base[0]) > 8 and "merge_view/calls.combined.tsv" in base[0], [k for k in base[0] if state[0].get(k) != base[0][k]])
            check(f"W1 {kind} index, {name}: the run summaries (without timings and speed) are equal", state[1] == base[1] and len(base[1]) == 2,
                  [k for k in base[1] if state[1].get(k) != base[1][k]])
        out = TMP / f"w1_{kind}_workers3_fast"
        summ = json.loads((out / R.SUMMARY).read_text()) if (out / R.SUMMARY).is_file() else {}
        sp = summ.get("speed", {})
        check(f"W1 {kind} index: the summary's speed entry records workers, counter, tile, cache and seconds",
              sp.get("workers") == 3 and sp.get("counter") == "fast" and sp.get("tile") == PA.TILE and sp.get("self_check_words", 0) >= 2000
              and set(sp.get("seconds", {})) >= {"setup", "scan", "first_look", "walks"} and "speed" not in json.loads(
                  (out / "merge_view" / R.SUMMARY).read_text()), sp)
    # W2: pieces written by each setting are reused by the others
    work = TMP / "w2_work"
    idx = TMP / "idx_w2"
    first = None
    for step, (name, extra) in enumerate([settings[0], settings[1], settings[2], settings[0]]):
        out = TMP / f"w2_{step}_{name}"
        env = {"EVITEA_PIECES_READONLY": "1"} if step else None
        r = driver(out, extra + ["--work-dir", str(work), "--index-dir", str(idx), "--reference-cache", str(TMP / "cache_w2")], env=env)
        check(f"W2 step {step} {name}: " + ("scans and writes the pieces" if not step else "reuses the saved pieces (read-only: a scan would stop it)"),
              r.returncode == 0, r.stderr[-500:])
        if r.returncode == 0:
            state = folder_state(out)
            if first is None:
                first = state
            else:
                check(f"W2 step {step} {name}: the output files equal those of the run that wrote the pieces", state[0] == first[0],
                      [k for k in first[0] if state[0].get(k) != first[0][k]])
    # the reverse: pieces written with 3 workers and the fast counter, then read by one plain process
    work2, idx2 = TMP / "w2b_work", TMP / "idx_w2b"
    ra = driver(TMP / "w2b_a", settings[1][1] + ["--work-dir", str(work2), "--index-dir", str(idx2)])
    rb = driver(TMP / "w2b_b", settings[0][1] + ["--work-dir", str(work2), "--index-dir", str(idx2)], env={"EVITEA_PIECES_READONLY": "1"})
    check("W2 pieces written with 3 workers and the fast counter are read by one plain process", ra.returncode == 0 and rb.returncode == 0
          and folder_state(TMP / "w2b_a")[0] == folder_state(TMP / "w2b_b")[0], (ra.stderr[-300:], rb.stderr[-300:]))
    pieces = sorted((work2 / "pieces").glob("*.pkl"))
    check("W2 the saved piece is one file for the contig", len(pieces) == 1, pieces)
    # tools/compare_scan.py: the comparison of the Sycamore checks, on pieces of two settings with one identity
    import compare_scan as CS
    work3 = TMP / "w2c_work"
    rc = driver(TMP / "w2c", settings[0][1] + ["--work-dir", str(work3), "--index-dir", str(idx2)])
    got = CS.compare_pieces(work2, work3)
    check("W3 compare_scan: the piece of 3 workers / fast and the piece of 1 worker / plain are EQUAL by value, identity included",
          rc.returncode == 0 and got and all(state == "EQUAL" for _, state, _ in got), got)
    got = CS.compare_outputs(TMP / "w2b_a", TMP / "w2c")
    check("W3 compare_scan: the output folders of the two runs are EQUAL file by file", got and all(state == "EQUAL" for _, state, _ in got) and len(got) > 16,
          [g for g in got if g[1] != "EQUAL"])
    with (work3 / "pieces" / "chrS.pkl").open("rb") as handle:
        piece = pickle.load(handle)
    piece["first"][0]["D"] = 12345.0
    other = TMP / "w2c_other"
    (other / "pieces").mkdir(parents=True)
    with (other / "pieces" / "chrS.pkl").open("wb") as handle:
        pickle.dump(piece, handle)
    (other / "pieces" / "chrS.complete.json").write_text((work3 / "pieces" / "chrS.complete.json").read_text())
    got = CS.compare_pieces(work3, other)
    check("W3 compare_scan: a changed field is DIFFERENT and named", got[0][1] == "DIFFERENT" and "first[0].D" in (got[0][2] or ""), got)
    mutated = TMP / "w2c_mutated"
    import shutil
    shutil.copytree(TMP / "w2c", mutated)
    target = mutated / "calls.combined.tsv"
    target.write_text(target.read_text() + "x\n")
    got = CS.compare_outputs(TMP / "w2c", mutated)
    check("W3 compare_scan: a changed output file is DIFFERENT", [g[0] for g in got if g[1] == "DIFFERENT"] == ["calls.combined.tsv"], got)


def final_src_checks():
    final = os.environ.get("FINAL_SRC")
    if not final:
        print("SKIP W4 (FINAL_SRC not set)")
        return
    import compare_scan as CS
    idx = TMP / "idx_w4"
    w1, w2 = TMP / "w4_work_old", TMP / "w4_work_new"
    common = ["--index-dir", str(idx)]
    old = driver_at(final, TMP / "w4_old", common + ["--work-dir", str(w1)])
    reuse = driver_at(str(ROOT / "src"), TMP / "w4_new_reuse", common + ["--work-dir", str(w1), "--workers", "3", "--counter", "fast"],
                      env={"EVITEA_PIECES_READONLY": "1"})
    fresh = driver_at(str(ROOT / "src"), TMP / "w4_new_fresh", common + ["--work-dir", str(w2), "--workers", "3", "--counter", "fast",
                                                                    "--reference-cache", str(TMP / "w4_cache")])
    back = driver_at(final, TMP / "w4_old_reuse", common + ["--work-dir", str(w2)], env={"EVITEA_PIECES_READONLY": "1"})
    check("W4 the four runs finish (the new code reuses the old pieces without a scan, the old code reuses the new ones)",
          all(r.returncode == 0 for r in (old, reuse, fresh, back)), [r.stderr[-300:] for r in (old, reuse, fresh, back)])
    for name in ("w4_new_reuse", "w4_new_fresh", "w4_old_reuse"):
        got = CS.compare_outputs(TMP / "w4_old", TMP / name)
        check(f"W4 {name}: all {len(got)} output files equal those of the final build's run", got and all(s == "EQUAL" for _, s, _ in got) and len(got) > 16,
              [g for g in got if g[1] != "EQUAL"])
    got = CS.compare_pieces(w1, w2)
    check("W4 the piece the final build wrote and the piece of the new code are EQUAL by value", got and all(s == "EQUAL" for _, s, _ in got), got)


def command_line_checks():
    import argparse
    from evitea import cli, run_command

    parser = cli.build_parser()
    base = ["library-free", "--run-dir", "r", "--reference-assets", "a", "--output", "o", "--sample-assets", "s"]
    a = parser.parse_args(base)
    check("S1 library-free defaults: one worker, plain counter, no reference cache folder", (a.workers, a.counter, a.reference_cache) == (1, "plain", None))
    a = parser.parse_args(base + ["--workers", "8", "--counter", "fast", "--reference-cache", "/x/cache"])
    check("S1 library-free reads --workers, --counter and --reference-cache", (a.workers, a.counter, str(a.reference_cache)) == (8, "fast", "/x/cache"))
    try:
        with open(os.devnull, "w") as sink:
            old = sys.stderr
            sys.stderr = sink
            try:
                parser.parse_args(base + ["--counter", "other"])
            finally:
                sys.stderr = old
        check("S1 an unknown counter is refused", False)
    except SystemExit:
        check("S1 an unknown counter is refused", True)
    ns = argparse.Namespace(threads=0, het_guard="on", keep_resident_copies=False, contigs=None, te_library="l", family_map="f", reference="g",
                            library_free=True, repeat_annotation=None, workers=None, read1="1", read2="2", read_length=None,
                            fragment_q95=None, fragment_q99=None, primary_contigs=None)
    steps = dict(run_command.library_steps(ns, Path("/o"), "S"))
    lf = steps["library_free"]
    check("S1 evitea run: --threads 0 gives the library-free step --workers 0 (every CPU)", lf[lf.index("--workers") + 1] == "0", lf)
    ns.threads = 6
    lf = dict(run_command.library_steps(ns, Path("/o"), "S"))["library_free"]
    check("S1 evitea run: --threads 6 gives --workers 6", lf[lf.index("--workers") + 1] == "6")
    ns.repeat_annotation = None
    ns.library_free = False
    nl = dict(run_command.no_library_steps(ns, Path("/o"), "S"))["library_free"]
    check("S1 evitea run without a library gives the library-free step --workers too", nl[nl.index("--workers") + 1] == "6")
    text = (ROOT / "tools" / "slurm" / "detect.sbatch").read_text()
    check("S1 detect.sbatch gives --workers (WORKERS or the job's CPUs), --counter and --reference-cache, not to the previous version",
          "--workers ${WORKERS:-${SLURM_CPUS_PER_TASK:-1}}" in text and "--counter $COUNTER" in text and "--reference-cache $REFCACHE" in text
          and 'if [ -z "${PREVIOUS:-}" ]; then' in text)


def hash_check():
    h = R._code_sha256()
    print("detection code hash", h[:8])
    want = os.environ.get("EXPECTED_CODE_HASH")
    if want:
        check(f"H1 the detection code hash is {want}", h.startswith(want), h[:16])
    else:
        print("SKIP H1 (EXPECTED_CODE_HASH not set)")
    here = Path(R.__file__).parent
    check("H1 the new files are not among the files of the identity", not {"fastcount.py", "parallel.py", "refcache.py"} & set(R.CODE_FILES))
    check("H1 the identity of a piece does not mention the scan-speed options", "workers" not in json.dumps(Options().detection_identity) and
          not any(w in (here / "options.py").read_text() for w in ("workers", "fastcount", "refcache")))


if __name__ == "__main__":
    sections = sys.argv[1:] or ["counter", "ref", "tiles", "runs", "final", "cli", "hash"]
    for name, fn in (("counter", counter_checks), ("ref", refcache_checks), ("tiles", tile_checks), ("runs", run_checks), ("final", final_src_checks), ("cli", command_line_checks),
                     ("hash", hash_check)):
        if name in sections:
            print(f"== {name}", flush=True)
            fn()
    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)
