"""Family records: flagging, the record guard, the record builder, the build outputs, and the second-pass calls.
Run: PYTHONPATH=src python -B tests/test_records.py

F1  flag rules kept / before-filters / both on crafted rows (away >= 2 and away > near; scans apart; statuses counted as stated).
F2  guard threshold = max(0.5, 2 x TE share); inside share from the annotation or the table's column.
F3  build_record: mixed-orientation members give one record with both ends; fewer than two usable members give none.
F4  join_ends joins on an exact overlap, otherwise separates by SPACER N.
F5  word_copies finds the copies of a record in the reference (diverged copies, both strands), none in unrelated sequence.
F6  build outputs: library, map (names), padded map >= 100 entries and parses, annotation, report columns, sha file; second-scan names.
F7  guard on: a family of old copies gets no record (listed), guard off builds it.
F8  add_record_calls: padding calls dropped, own-copy calls dropped, calls near the existing set not added, close calls counted once.
F9  the two commands parse and run through the CLI.
"""
import csv
import os
import random
import sys
import tempfile
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

try:
    import refs  # reference implementations: development only
    BANNED = refs.banned_words()
except Exception:
    BANNED = None
from evitea.family_records import build as B, run as RUN, select as S  # noqa: E402
from evitea.library_free.resident import CopyIndex  # noqa: E402

FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


rng = random.Random(11)
rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
rc = B.rc


def mutate(s, rate):
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in s)


def row(cluster, status, chrom="chrS", bL=1000, scan="first", near=0, rc_drop=0, inside=""):
    return dict(chrom=chrom, bL=str(bL), cluster=str(cluster), status=status, near_evitea=str(near), rule_c_dropped=str(rc_drop),
                scan=scan, family_inside=inside)


def members(elem, n, ln=90, rate=0.02):
    out = []
    for k in range(n):
        e = mutate(elem, rate)
        l, r = e[:ln], e[-ln:]
        out.append((l, r) if k % 2 == 0 else (rc(r), rc(l)))
    return out


def contains_either(seq, piece):
    return piece in seq or rc(piece) in seq


def write_fa(path, d):
    Path(path).write_text("".join(">%s\n%s\n" % (k, v) for k, v in d.items()))


def main():
    # F1
    rows = [row(1, "ADDED", bL=1000), row(1, "ADDED", bL=2000), row(1, "NEAR_EVITEA", bL=3000, near=1),
            row(2, "ADDED", bL=5000), row(2, "NEAR_EVITEA", bL=6000, near=1),
            row(3, "DROPPED_V2", bL=7000), row(3, "DROPPED_CHECK", bL=8000), row(3, "NEAR_EVITEA", bL=9000, near=1),
            row(4, "DROPPED_RESIDENT_COPY", bL=1100), row(4, "ADDED", bL=1200),
            row(5, "ADDED", bL=1300, scan="second"), row(5, "NEAR_FIRST_PASS", bL=1400, scan="second"),
            row(6, "DROPPED_RULE_C", bL=1500, rc_drop=1), row(6, "DROPPED_RULE_C", bL=1600, rc_drop=1),
            row(7, "ADDED", bL=1700, scan="second"), row(7, "ADDED", bL=1800, scan="second"), row(7, "NOT_HET", bL=1900, scan="second")]
    k = S.flagged_families(rows, "kept")
    b = S.flagged_families(rows, "before-filters")
    both = S.flagged_families(rows, "both")
    check("F1 kept: cluster 1 (2 away, 1 near) and 4 (dropped resident + added) and second-scan 7", list(k) == [("first", "1"), ("first", "4"), ("second", "7")], list(k))
    check("F1 kept: cluster 2 (1 away) and 5 (1 away, 1 near) not flagged", ("first", "2") not in k and ("second", "5") not in k)
    check("F1 before-filters: cluster 3 (two away before the filters) flagged, the family-identity filter dropped 6 not", ("first", "3") in b and ("first", "6") not in b, list(b))
    check("F1 both is the union", set(both) == set(k) | set(b))
    check("F1 flagged_by lists the rule", both[("first", "3")]["flagged_by"] == ["before-filters"] and both[("first", "1")]["flagged_by"] == ["kept", "before-filters"])
    check("F1 a second-scan NOT_HET call is not counted", b[("second", "7")]["before-filters"] == (2, 0))
    try:
        S.flagged_families(rows, "nope")
        check("F1 unknown rule refused", False)
    except ValueError:
        check("F1 unknown rule refused", True)

    # F2
    check("F2 threshold floor 0.5", S.guard_threshold(0.1) == 0.5 and S.guard_threshold(None) == 0.5)
    check("F2 threshold 2 x share above 0.25", abs(S.guard_threshold(0.4) - 0.8) < 1e-9)
    cp = CopyIndex({"chrS": [(900, 1100)]})
    rr = [row(1, "ADDED", bL=1000), row(1, "ADDED", bL=5000)]
    check("F2 inside share from the annotation", S.inside_share(rr, cp) == 0.5)
    check("F2 inside share from the column", S.inside_share([row(1, "ADDED", inside="0.7500")], None) == 0.75)
    check("F2 unknown without either", S.inside_share([row(1, "ADDED")], None) is None)

    # F3, F4
    elem = rnd(420)
    rec = B.build_record(members(elem, 6))
    check("F3 record built from 6 members", rec is not None and rec.placed >= 5, rec)
    check("F3 both ends in the record, spacer between", rec and contains_either(rec.seq, elem[10:70]) and rec.joined == "spacer" and "N" * 100 in rec.seq)
    ident = lambda x, y: sum(p == q for p, q in zip(x, y)) / max(1, len(y))
    a_, b_ = rec.seq[:rec.a_len], rec.seq[-rec.b_len:]
    fwd = ident(a_, elem[:90]) >= 0.9 and ident(b_, elem[-90:]) >= 0.9
    rev = ident(a_, rc(elem[-90:])) >= 0.9 and ident(b_, rc(elem[:90])) >= 0.9
    check("F3 end A and end B are the two ends of the element, in one orientation", fwd or rev, (fwd, rev))
    check("F3 one usable member gives none", B.build_record(members(elem, 1)) is None)
    a = rnd(80)
    b2 = a[40:] + rnd(50)
    j, how = B.join_ends(a, b2)
    check("F4 overlap joined", how == "overlap" and j == a[:40] + b2, how)
    j, how = B.join_ends(rnd(60), rnd(60))
    check("F4 spacer otherwise", how == "spacer" and "N" * B.SPACER in j)

    # reference with diverged copies of elem on both strands
    ref_parts = [rnd(30_000)]
    for k in range(5):
        c = mutate(elem, 0.04)
        ref_parts += [c if k % 2 == 0 else rc(c), rnd(8_000)]
    ref = "".join(ref_parts)
    other = rnd(20_000)
    tmp = Path(tempfile.mkdtemp(prefix="rec_"))
    ref_fa = tmp / "ref.fa"
    write_fa(ref_fa, {"chrS": ref, "chrU": other})
    cps = B.word_copies({"FBte9900001": elem, "FBte9900002": rnd(300)}, ref_fa)
    mine = [c for c in cps if c[3] == "FBte9900001"]
    check("F5 copies of the record found (>= 4 of 5 at 4% divergence)", len([c for c in mine if c[0] == "chrS"]) >= 4, len(mine))
    check("F5 none for an unrelated record, none on the unrelated contig", not [c for c in cps if c[3] == "FBte9900002"] and not [c for c in mine if c[0] == "chrU"])

    # F6 build outputs on a table of two first-pass families and one second-scan family
    el2, el3 = rnd(420), rnd(380)
    fam_rows, seg = [], {}
    def add_family(scan, cl, elem_, n, start):
        mem = members(elem_, n)
        for k in range(n):
            bL = start + 3000 * k
            fam_rows.append(row(cl, "ADDED", bL=bL, scan=scan))
            l, r = mem[k]
            seg[(scan, "chrS", str(bL), str(cl))] = (l, r)
    add_family("first", 1, elem, 5, 1000)
    add_family("first", 2, el2, 4, 200_000)
    add_family("second", 9, el3, 4, 400_000)
    fam_rows.append(row(3, "NEAR_EVITEA", bL=700_000, near=1))
    fam_tsv = tmp / "fc.tsv"
    cols = ["chrom", "bL", "cluster", "status", "near_evitea", "rule_c_dropped", "scan", "family_inside"]
    fam_tsv.write_text("\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n" for r in fam_rows))
    for scan, name in (("first", "seg1.fa"), ("second", "seg2.fa")):
        (tmp / name).write_text("".join(">%s:%s:lseg:cl%s\n%s\n>%s:%s:rseg:cl%s\n%s\n" % (c, b, cl, l, c, b, cl, r)
                                        for (sc, c, b, cl), (l, r) in seg.items() if sc == scan))
    # the real library builder needs a real reference; the end-word step is replaced by a stub here (the real one needs a
    # full reference)
    real_languages = B.languages
    B.languages = lambda te, fm, bed, ref: (b"te_id\tword\n" + b"".join(("%s\tAAAA\n" % ln.strip()[1:]).encode() for ln in open(te) if ln.startswith(">")),
                                            b"te_id\tword\n")
    try:
        out = tmp / "out"
        summ = RUN.build_family_records(family_calls=[fam_tsv], segments=tmp / "seg1.fa", second_segments=tmp / "seg2.fa",
                                        reference_fasta=ref_fa, output=out, progress=lambda m: None)
        names = [l.split("\t")[:2] for l in (out / "family_map.tsv").read_text().splitlines()[1:]]
        check("F6 three records, names by scan", [n[1] for n in names] == ["LF1", "LF2", "LFh9"] and [n[0] for n in names] == ["FBte9900001", "FBte9900002", "FBte9700009"], names)
        pm = (out / "family_map_padded.tsv").read_text().splitlines()
        check("F6 padded map has 100 entries (3 records + 97 inert)", len(pm) - 1 == 100 and pm[-1].startswith("FBte9800097\tPAD097"), len(pm))
        from evitea.anchor_pe.clusters import load_family_map
        mp, fams = load_family_map(RUN.read_tsv(out / "family_map_padded.tsv"))
        check("F6 padded map parses as a family map", len(mp) == 100)
        if BANNED is not None:
            check("F6 no competitor word in any output", not any(w in (out / f).read_text().lower() for f in os.listdir(out) for w in BANNED))
        rep = RUN.read_tsv(out / "records_report.tsv")
        check("F6 report lists every flagged family with its status", [r["status"] for r in rep] == ["built", "built", "built"] and rep[0]["flagged_by"] == "kept", rep[0])
        check("F6 annotation has copies of record 1", any(l.endswith("FBte9900001") for l in (out / "annotation.bed").read_text().splitlines()))
        sh = (out / "sha12.txt").read_text()
        check("F6 sha file lists records and files", sh.startswith("records 3") and "te_library.fa " in sh)
        check("F6 summary", summ["built"] == 3 and summ["skipped_guard"] == [])

        # F7 guard: annotate record-1 sites as old copies
        bed = tmp / "rm.bed"
        bed.write_text("".join("chrS\t%d\t%d\tX\n" % (int(r_["bL"]) - 50, int(r_["bL"]) + 50) for r_ in fam_rows if r_["cluster"] == "1"))
        out2 = tmp / "out_guard"
        s2 = RUN.build_family_records(family_calls=[fam_tsv], segments=tmp / "seg1.fa", second_segments=tmp / "seg2.fa",
                                      reference_fasta=ref_fa, output=out2, repeat_annotation=bed, progress=lambda m: None)
        check("F7 guard skips the family of old copies and lists it", s2["skipped_guard"] == ["LF1"] and s2["built"] == 2, s2)
        check("F7 skipped family is not in the library", "FBte9900001" not in (out2 / "te_library.fa").read_text())
        out3 = tmp / "out_noguard"
        s3 = RUN.build_family_records(family_calls=[fam_tsv], segments=tmp / "seg1.fa", second_segments=tmp / "seg2.fa",
                                      reference_fasta=ref_fa, output=out3, repeat_annotation=bed, record_guard=False, progress=lambda m: None)
        check("F7 guard off builds it", s3["built"] == 3 and s3["skipped_guard"] == [])
        # F7b the repeat map stands in for the annotation; an annotation given wins over it
        out5 = tmp / "out_map"
        s5 = RUN.build_family_records(family_calls=[fam_tsv], segments=tmp / "seg1.fa", second_segments=tmp / "seg2.fa",
                                      reference_fasta=ref_fa, output=out5, repeat_map=bed, progress=lambda m: None)
        check("F7b the repeat map gives the guard of the annotation", s5["skipped_guard"] == ["LF1"] and s5["built"] == 2
              and s5["record_guard"]["repeat_map"] == str(bed) and s5["record_guard"]["te_share"] == s2["record_guard"]["te_share"], s5)
        none_bed = tmp / "none.bed"
        none_bed.write_text("chrS\t0\t10\tX\n")
        s6 = RUN.build_family_records(family_calls=[fam_tsv], segments=tmp / "seg1.fa", second_segments=tmp / "seg2.fa",
                                      reference_fasta=ref_fa, output=tmp / "out_both", repeat_annotation=none_bed, repeat_map=bed,
                                      progress=lambda m: None)
        check("F7b an annotation given wins over the map", s6["built"] == 3 and "repeat_map" not in s6["record_guard"], s6)
        out4 = tmp / "out_copy"
        RUN.build_family_records(family_calls=[fam_tsv], segments=tmp / "seg1.fa", second_segments=tmp / "seg2.fa",
                                 reference_fasta=ref_fa, output=out4, copy_report=True, progress=lambda m: None)
        cr = RUN.read_tsv(out4 / "copy_report.tsv")
        check("F7 copy report: absent share and copies per record", len(cr) == 3 and float(cr[0]["reference_absent_share"]) < 0.9 and int(cr[0]["reference_copies"]) >= 4, cr[0])
        check("F7 copy report off by default", not (out / "copy_report.tsv").exists())
    finally:
        B.languages = real_languages

    # F8 add_record_calls
    calls = tmp / "calls.native.tsv"
    ch = ["chrom", "family", "start", "end", "total_support"]
    crow = [("chrS", "LF1", 100_000, 100_000, 9), ("chrS", "LF1", 100_200, 100_200, 12),     # close pair: one added (higher support)
            ("chrS", "LF1", 300_000, 300_000, 5),                                              # near an existing call
            ("chrS", "LF2", 30_100, 30_100, 7),                                                # inside an own copy of LF2
            ("chrS", "PAD003", 500_000, 500_000, 4),                                            # padding
            ("chrS", "LF2", 600_000, 600_000, 3)]
    calls.write_text("\t".join(ch) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in crow))
    base = tmp / "base.tsv"
    base.write_text("event_id\tchr\tstart\tend\tfamily\tsource\nE1\tchrS\t300100\t300100\troo\tEVITEA\n")
    rec_dir = tmp / "rec"
    rec_dir.mkdir()
    (rec_dir / "family_map.tsv").write_text(RUN.MAP_HEADER + "FBte9900001\tLF1\t\t\t\t\nFBte9900002\tLF2\t\t\t\t\nFBte9800003\tPAD003\t\t\t\t\n")
    (rec_dir / "annotation.bed").write_text("chrS\t30000\t30500\tFBte9900002\n")
    rp = RUN.add_record_calls(calls=calls, base_calls=base, annotation=rec_dir / "annotation.bed", family_map=rec_dir / "family_map.tsv", output=tmp / "add")
    check("F8 report counts", (rp["record_calls"], rp["inside_own_copy"], rp["near_existing"], rp["same_place"], rp["added"], rp["padding_calls"]) == (5, 1, 1, 1, 2, 1), rp)
    add = RUN.read_tsv(tmp / "add" / "record_calls_added.tsv")
    check("F8 of a close pair the call nearer the contig start is added", sorted((r["family"], r["start"]) for r in add) == [("LF1", "100000"), ("LF2", "600000")], add)
    comb = RUN.read_tsv(tmp / "add" / "calls.combined.tsv")
    check("F8 combined = existing + added with source RECORD", len(comb) == 3 and comb[0]["source"] == "EVITEA" and comb[-1]["source"] == "RECORD")

    # F9 CLI
    from evitea import cli
    rc_ = cli.main(["add-record-calls", "--calls", str(calls), "--base-calls", str(base), "--records", str(rec_dir), "--output", str(tmp / "add2")])
    check("F9 add-record-calls through the CLI", rc_ == 0 and (tmp / "add2" / "record_calls_report.json").exists())
    p = cli.build_parser()
    ns = p.parse_args(["build-family-records", "--family-calls", str(fam_tsv), "--reference-fasta", str(ref_fa), "--output", str(tmp / "x")])
    check("F9 build-family-records defaults", ns.flag_rule == "kept" and ns.record_guard == "on" and ns.copy_report == "off" and ns.repeat_map is None)
    ns2 = p.parse_args(["build-family-records", "--family-calls", str(fam_tsv), "--reference-fasta", str(ref_fa), "--output", str(tmp / "x"),
                        "--repeat-map", str(tmp / "m.bed")])
    check("F9 --repeat-map is read", ns2.repeat_map == tmp / "m.bed")

    print("\n%d checks, %d failed" % (NCHK[0], len(FAIL)))
    sys.exit(1 if FAIL else 0)


main()
