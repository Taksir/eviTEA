"""The repeat map (src/evitea/repeat_map.py).   Run: PYTHONPATH=src python -B tools/library_free_tests/test_repeat_map.py

M1  hit lines to intervals: the window's own locus, hits below 80% identity or 100 aligned bases are left out; a reversed subject
    is read in forward coordinates; subject sides are marked only on the chosen contigs; an unknown subject is refused.
M2  windows: step 250, a window of more than half N is left out, a contig shorter than a window has none.
M3  a reference with a planted repeat (BLAST+ on PATH; skipped without it): the map covers the planted copies (forward,
    diverged and reverse-complement) and nothing else, does not depend on the number of parts, refuses to replace a map, and a
    second build in the same folder with its saved BLAST work gives the same file.
M4  the command: ``evitea build-repeat-map`` writes the map and the summary.
"""
import json
import random
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from evitea import cli  # noqa: E402
from evitea import repeat_map as RM  # noqa: E402
from evitea.errors import InputContractError  # noqa: E402

FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def main():
    # M1
    hits = [
        "chrA|0\tchrA\t100.0\t500\t1\t500\t1\t500",            # own locus: left out
        "chrA|1000\tchrA\t99.0\t300\t11\t310\t2000\t2299",    # a second place on the same contig
        "chrA|2000\tchrB\t79.9\t300\t1\t300\t100\t399",       # identity below 80
        "chrA|2000\tchrB\t95.0\t99\t1\t99\t100\t198",         # shorter than 100
        "chrA|3000\tlcl|chrB\t90.0\t200\t101\t300\t900\t701", # reversed subject, lcl prefix
        "chrA|4000\tchrC\t92.0\t150\t1\t150\t50\t199",        # subject on a contig that is not chosen
    ]
    iv = RM.intervals_from_hits(hits, ["chrA", "chrB"])
    merged = RM.merge_intervals(iv)
    check("M1 intervals", merged == {"chrA": [(1010, 1310), (1999, 2299), (3100, 3300), (4000, 4150)], "chrB": [(700, 900)]}, merged)
    try:
        RM.intervals_from_hits(["chrA|0\tchrZ\t90.0\t200\t1\t200\t1\t200"], ["chrA"], known=["chrA", "chrB"])
        check("M1 an unknown subject is refused", False)
    except InputContractError:
        check("M1 an unknown subject is refused", True)

    # M2
    t = Path(tempfile.mkdtemp(prefix="rmap_"))
    r = random.Random(11)
    seq = lambda n: "".join(r.choice("ACGT") for _ in range(n))
    fa = t / "w.fa"
    body = seq(1200) + "N" * 600 + seq(500)
    fa.write_text(">long extra words\n" + body + "\n>short\n" + seq(400) + "\n")
    win = list(RM.window_records(fa, ["long", "short"]))
    starts = [s for _, s, _ in win]
    ok_n = (all(w.count("N") <= 250 for _, _, w in win) and not {1000, 1250, 1500} & set(starts) and {0, 750, 1750} <= set(starts)
            and starts == sorted(starts))
    check("M2 windows of 500 at step 250; N-rich windows and short contigs have none", ok_n and all(c == "long" for c, _, _ in win)
          and all(len(w) == 500 for _, _, w in win), starts)

    # M3
    elem = seq(700)
    div = list(elem)
    for k in range(0, 700, 20):                  # 5% substitutions
        div[k] = {"A": "C", "C": "G", "G": "T", "T": "A"}[div[k]]
    div = "".join(div)

    def contig(parts):
        return "".join(parts)

    a = contig([seq(1000), elem, seq(1500)])                 # element at 1000
    b = contig([seq(3000), div, seq(1000)])                  # diverged copy at 3000
    c = contig([seq(2000), rc(elem), seq(2000)])             # reverse complement at 2000
    d = seq(4000)                                            # nothing repeated, not chosen
    ref = t / "ref.fa"
    ref.write_text("".join(f">{n}\n{s}\n" for n, s in (("chrA", a), ("chrB", b), ("chrC", c), ("chrD", d))))
    planted = {"chrA": (1000, 1700), "chrB": (3000, 3700), "chrC": (2000, 2700)}
    if not all(shutil.which(x) for x in RM.TOOLS):
        print("M3 skipped: BLAST+ (makeblastdb, blastn) is not on PATH")
    else:
        beds = {}
        for parts in (1, 3):
            out = t / f"map{parts}"
            s = RM.build_repeat_map(reference_fasta=ref, output=out, contigs=["chrA", "chrB", "chrC"], parts=parts, progress=lambda m: None)
            beds[parts] = (out / "repeat_map.bed").read_text()
        check("M3 the map does not depend on the number of parts", beds[1] == beds[3], (beds[1], beds[3]))
        rows = [line.split("\t") for line in beds[1].splitlines()]
        got = {}
        for contig_, x, y in rows:
            got.setdefault(contig_, []).append((int(x), int(y)))
        cover = {n: sum(min(y, hi) - max(x, lo) for x, y in got.get(n, ()) if min(y, hi) > max(x, lo)) for n, (lo, hi) in planted.items()}
        check("M3 every planted copy is covered (at least 90%)", all(v >= 630 for v in cover.values()), cover)
        outside = {n: [(x, y) for x, y in v if x < planted[n][0] - 40 or y > planted[n][1] + 40] for n, v in got.items()}
        check("M3 nothing outside the planted copies, and chrD (not chosen, no repeat) has no interval",
              not any(outside.values()) and "chrD" not in got, outside)
        check("M3 summary counts", s["intervals"] == len(rows) and s["bases"] == sum(int(y) - int(x) for _, x, y in rows)
              and s["contigs"] == ["chrA", "chrB", "chrC"], s)
        try:
            RM.build_repeat_map(reference_fasta=ref, output=t / "map1", contigs=["chrA", "chrB", "chrC"], parts=1, progress=lambda m: None)
            check("M3 an existing map is not replaced", False)
        except InputContractError:
            check("M3 an existing map is not replaced", True)
        (t / "map3" / "repeat_map.bed").unlink()
        RM.build_repeat_map(reference_fasta=ref, output=t / "map3", contigs=["chrA", "chrB", "chrC"], parts=1, progress=lambda m: None)
        check("M3 a new build with saved BLAST work gives the same map", (t / "map3" / "repeat_map.bed").read_text() == beds[1])
        try:
            RM.build_repeat_map(reference_fasta=ref, output=t / "map3", contigs=["chrA", "chrB"], progress=lambda m: None)
            check("M3 saved BLAST work of other contigs is refused", False)
        except InputContractError:
            check("M3 saved BLAST work of other contigs is refused", True)

        # M4
        cli_out = t / "cli_map"
        code = cli.main(["build-repeat-map", "--reference-fasta", str(ref), "--output", str(cli_out), "--contigs", "chrA,chrB,chrC", "--parts", "2"])
        summ = json.loads((cli_out / "repeat_map_summary.json").read_text())
        check("M4 command writes the map and the summary", code == 0 and (cli_out / "repeat_map.bed").read_text() == beds[1]
              and summ["constants"]["window"] == 500 and summ["constants"]["min_identity"] == 80.0, summ)

    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


main()
