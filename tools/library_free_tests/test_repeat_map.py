"""The repeat map (src/evitea/repeat_map.py).   Run: PYTHONPATH=src python -B tools/library_free_tests/test_repeat_map.py

M1  hit lines to intervals: the window's own locus, hits below 80% identity or 100 aligned bases are left out; a reversed subject
    is read in forward coordinates; subject sides are marked only on the chosen contigs; blastn's wrapped ids of accession-like names
    (ref|NW_001845059.1|, gi|..|ref|..|) are read as the reference's names; an unknown subject is refused.
M2  windows: step 250, a window of more than half N is left out, a contig shorter than a window has none.
M3  a reference with a planted repeat (BLAST+ on PATH; skipped without it): the map covers the planted copies (forward,
    diverged and reverse-complement) and nothing else, does not depend on the number of parts, refuses to replace a map, and a
    second build in the same folder with its saved BLAST work gives the same file.
M4  the command: ``evitea build-repeat-map`` writes the map and the summary.
M5  default contigs: the primary chromosomes and every contig of at least 1 Mb, in reference order. On dm6 and mm10 lengths the
    list is that of the contigs of at least 1 Mb; on yeast lengths it is all 16 chromosomes (the mitochondrion stays out); the
    map's default reads the same rule from a FASTA.
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
    names = {"chrA", "chrB", "NW_001845059.1", "NW_001845060.1"}
    acc = RM.intervals_from_hits(["chrA|1000\tref|NW_001845059.1|\t95.0\t300\t1\t300\t50\t349",
                                  "chrA|3000\tgi|123|ref|NW_001845060.1|\t95.0\t300\t1\t300\t500\t201",
                                  "chrA|5000\tgb|chrB|\t95.0\t200\t1\t200\t10\t209"], ["chrA", "chrB"], known=names)
    check("M1 wrapped accession ids are read as the reference's names; hits to contigs not chosen mark the query side only",
          RM.merge_intervals(acc) == {"chrA": [(1000, 1300), (3000, 3300), (5000, 5200)], "chrB": [(9, 209)]}, RM.merge_intervals(acc))
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
    unplaced = seq(300) + elem[:400] + seq(300)                 # an unplaced contig with an NCBI-style name, not chosen
    ref.write_text("".join(f">{n}{' desc' if n.startswith('NW') else ''}\n{s}\n"
                           for n, s in (("chrA", a), ("chrB", b), ("NW_001845059.1", unplaced), ("chrC", c), ("chrD", d))))
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
        bed6 = (t / "map1" / RM.BED6).read_text()
        check("M3 the BED6 copy holds the same intervals, named repeat, score 0, no strand",
              bed6 == "".join(f"{line}\trepeat\t0\t.\n" for line in beds[1].splitlines()), bed6[:200])
        from evitea.inputs import load_repeat_annotation, validate_reference_fasta
        loaded = load_repeat_annotation(t / "map1" / RM.BED6, validate_reference_fasta(ref))
        check("M3 prepare-reference reads the BED6 copy as a repeat annotation", len(loaded.records) == len(beds[1].splitlines()))
        rows = [line.split("\t") for line in beds[1].splitlines()]
        got = {}
        for contig_, x, y in rows:
            got.setdefault(contig_, []).append((int(x), int(y)))
        cover = {n: sum(min(y, hi) - max(x, lo) for x, y in got.get(n, ()) if min(y, hi) > max(x, lo)) for n, (lo, hi) in planted.items()}
        check("M3 every planted copy is covered (at least 90%)", all(v >= 630 for v in cover.values()), cover)
        outside = {n: [(x, y) for x, y in v if x < planted[n][0] - 40 or y > planted[n][1] + 40] for n, v in got.items()}
        check("M3 nothing outside the planted copies, and chrD (not chosen, no repeat) has no interval",
              not any(outside.values()) and "chrD" not in got, outside)
        check("M3 the NCBI-style unplaced contig (not chosen) has no interval, and its copy still marks the other contigs",
              "NW_001845059.1" not in got and cover["chrA"] >= 630, got)
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

    # M5
    from evitea.library_free.runtime import MIN_DEFAULT_CONTIG, default_contigs

    def old_rule(lengths):
        return tuple(name for name, n in lengths.items() if n >= MIN_DEFAULT_CONTIG)

    dm6 = {"chr2L": 23_513_712, "chr2R": 25_286_936, "chr3L": 28_110_227, "chr3R": 32_079_331, "chr4": 1_348_131,
           "chrM": 19_524, "chrX": 23_542_271, "chrY": 3_667_352, "chrX_DS483648v1_random": 28_022, "NW_001845059.1": 7_540}
    mm10 = {f"chr{c}": n for c, n in zip([*range(1, 20), "X", "Y"],
            [195_471_971, 182_113_224, 160_039_680, 156_508_116, 151_834_684, 149_736_546, 145_441_459, 129_401_213,
             124_595_110, 130_694_993, 122_082_543, 120_129_022, 120_421_639, 124_902_244, 104_043_685, 98_207_768,
             94_987_271, 90_702_639, 61_431_566, 171_031_299, 91_744_698])}
    mm10.update({"chrM": 16_299, "chr1_GL456211_random": 241_735, "chrUn_JH584304": 114_452})
    yeast = dict(zip(["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI",
                      "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrM"],
                     [230_218, 813_184, 316_620, 1_531_933, 576_874, 270_161, 1_090_940, 562_643, 439_888, 745_751, 666_816,
                      1_078_177, 924_431, 784_333, 1_091_291, 948_066, 85_779]))
    check("M5 dm6: the default is the contigs of at least 1 Mb", default_contigs(dm6) == old_rule(dm6)
          == ("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"), default_contigs(dm6))
    check("M5 mm10: the default is the contigs of at least 1 Mb (chr1-19, X, Y)",
          default_contigs(mm10) == old_rule(mm10) == tuple(f"chr{c}" for c in [*range(1, 20), "X", "Y"]), default_contigs(mm10))
    check("M5 yeast: all 16 chromosomes in reference order, no mitochondrion (4 reach 1 Mb)",
          default_contigs(yeast) == tuple(list(yeast)[:16]) and len(old_rule(yeast)) == 4, default_contigs(yeast))
    big_mito = {"chr1": 30_000_000, "chrM": 1_200_000, "scaffold_9": 2_000_000}
    check("M5 every contig of at least 1 Mb stays in the default", default_contigs(big_mito) == ("chr1", "chrM", "scaffold_9"),
          default_contigs(big_mito))
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "small.fa"
        fa.write_text("".join(f">{n}\n{'ACGT' * (k // 4)}\n" for n, k in
                              [("chrI", 1000), ("chrII", 900), ("chrM", 400), ("chrIII", 152), ("tig_7", 48)]))
        check("M5 the map's default reads the rule from the FASTA", RM.default_contigs(fa) == ("chrI", "chrII", "chrIII"),
              RM.default_contigs(fa))

    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


main()
