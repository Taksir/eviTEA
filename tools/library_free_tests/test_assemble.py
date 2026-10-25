"""Family-call assembly on crafted tables: the one-site check, check scopes, the resident-copy test, the join of the first pass
with the haplotype scan, and the family-identity filter on trusted members.
Run: PYTHONPATH=src python -B tools/library_free_tests/test_assemble.py

A1  scope none: assemble gives calling.family_calls' rows (status, e, c) on the same tables.
A2  unconfirmed: a junction at twice the single-copy depth is dropped (DROPPED_CHECK); one at single-copy depth is added.
A3  unconfirmed: a junction whose k-mers are in the reference is dropped at single-copy depth.
A4  scope path-dropped: a call with e below T_E is not checked; one above is, and e above REFERENCE_PATH_CEILING stays dropped by the reference-path filter.
A5  tolerance: the p99 of the library-confirmed sides; own-calls source and the fallback when none is confirmed.
A6  second scan: NOT_HET / NEAR_FIRST_PASS / DROPPED_CHECK / ADDED by the rules; a call within 500 bp of a first-pass family call is not added.
A7  family-identity filter on trusted members keeps a family that its own grouping dropped because of untrusted diverged members.
A8  resident-copy test: family and per-call modes; off without an annotation.
A9  the check refuses tables with no recorded junctions.
A10 identity evidence (--identity-min-comparisons): a family whose plain-extension members give one comparison below the cutoff but whose
    members as a whole give many above it is kept at 3 and dropped at 0. A family with one comparison in all is not judged at 3 (kept)
    and judged at 0 (dropped). The same holds for the family-identity filter on trusted members of the haplotype scan.
"""
import os
import random
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from evitea.errors import InputContractError  # noqa: E402
from evitea.library_free import _p4_lib as P, calling as CL, assemble as AS, check as CK  # noqa: E402
from evitea.library_free.options import Options  # noqa: E402
from evitea.library_free.resident import CopyIndex  # noqa: E402

FAIL, NCHK = [], [0]


def check(name, ok, detail=""):
    NCHK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


rng = random.Random(3)
rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
GENOME = {"chrS": rnd(60_000)}
D0 = 20.0


def mutate(s, rate):
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in s)


class Table:
    """Places rows of one pass and their segments."""

    def __init__(self):
        self.rows, self.segs, self.cluster = [], {}, 0

    def family(self, sites, e=0.0, ratio=1.0, junction=None, lseg=None, rseg=None, tsd=5, stage_extra=None, two=True):
        """One family: a list of (bL) sites sharing a cluster id; returns the cluster id."""
        cl = str(self.cluster)
        self.cluster += 1
        body = lseg or rnd(200)
        for k, bL in enumerate(sites):
            ls = body if k == 0 else mutate(body, 0.0)
            rs = rseg or rnd(150)
            T = int(round(ratio * D0))
            junc = (junction[k] if isinstance(junction, list) else junction) or rnd(60)
            row = dict(chrom="chrS", first=str(bL - 10), last=str(bL), n_low="1", D=str(D0), stage="S4_pass", why="", bL=str(bL),
                       bR=str(bL - tsd) if two else "", delta=str(tsd) if two else "", tsd=str(tsd) if two else "",
                       two_sided="1" if two else "0", geno="het", fracL=str(e), fracR=str(e) if two else "", rep="0.1",
                       lstop="support", rstop="support" if two else "", rstat="break" if two else "through", seg_keep="L,R" if two else "L",
                       seg_reject="", dup="0", cluster=cl, track="family", fam_pair="1", fam_vtsd="1", target="",
                       lseg_len=str(len(ls)), rseg_len=str(len(rs) if two else 0),
                       chkL_T=str(T), chkL_junction=junc, chkR_T=str(T) if two else "", chkR_junction=rnd(60) if two else "",
                       origin="", walked="", bL0="")
            if stage_extra:
                row.update(stage_extra)
            self.rows.append(row)
            self.segs[("chrS", str(bL), cl)] = dict(lseg=ls, **({"rseg": rs} if two else {}))
        return cl

    def get(self):
        return self.rows, self.segs


def index(*pos):
    return CL.CallIndex({"chrS": [(p, p + 8) for p in pos]})


def run(first, second=None, idx=None, **kw):
    opts = kw.pop("opts", Options.previous().with_(check_scope="unconfirmed"))
    rc = kw.pop("identity_filter", False)
    return AS.assemble(first.get(), second.get() if second else None, idx or index(), opts, identity_filter=rc, reference=GENOME,
                       aligner=kw.pop("aligner", None), **kw)


def status(rows, scan="first"):
    return {(r.bL): r.status for r in rows if r.scan == scan}


def main():
    # ---------------------------------------------------------------- A1
    t = Table()
    t.family([1000, 6000], e=0.0)
    t.family([2000, 7000], e=0.30)
    t.family([3000, 8000], e=0.0)
    idx = index(3000)
    rows, _ = AS.assemble(t.get(), None, idx, Options.previous(), identity_filter=True, reference=None)
    fz = CL.family_calls(t.rows, idx, CL.family_similarity(t.rows, t.segs), rule_c=True)
    check("A1 scope none = calling.family_calls (status, e, c, near)",
          [(r.status, r.e, r.c, r.near_evitea, r.bL) for r in rows] == [(r.status, r.e, r.c, r.near_evitea, r.bL) for r in fz]
          and {r.status for r in rows} >= {"ADDED", "DROPPED_V2", "NEAR_EVITEA"}, sorted({r.status for r in rows}))
    check("A1 scope none: base columns only in the table", AS.family_calls_tsv_bytes(rows, extended=False)
          == CL.family_calls_tsv_bytes(fz))

    # ---------------------------------------------------------------- A2, A3
    t = Table()
    t.family([1000, 1010], ratio=1.0)                 # library-confirmed pair (near native): sets the tolerance
    t.family([5000, 9000], ratio=1.0)                 # single-copy depth, novel junction
    t.family([12000, 15000], ratio=2.0)               # twice the single-copy depth
    ref_junc = GENOME["chrS"][30000:30060]
    t.family([20000, 22000], ratio=1.0, junction=ref_junc)   # junction k-mers are in the reference
    idx = index(1000, 1010)
    # tolerance sides: ratios 1.0 for the library-confirmed calls
    rows, summ = run(t, idx=idx)
    st = status(rows)
    check("A2 single-copy depth, novel junction: ADDED", st[5000] == "ADDED" and st[9000] == "ADDED", st)
    check("A2 twice the single-copy depth: DROPPED_CHECK", st[12000] == "DROPPED_CHECK" and st[15000] == "DROPPED_CHECK", st)
    check("A3 junction words in the reference: DROPPED_CHECK", st[20000] == "DROPPED_CHECK", st)
    check("A2 the library-confirmed pair is NEAR_EVITEA and is not checked", st[1000] == "NEAR_EVITEA")
    check("A2 summary: tolerance from the library-confirmed sides", summ["tolerance"]["source"] == "library" and summ["tolerance"]["sides"] == 4,
          summ["tolerance"])

    # ---------------------------------------------------------------- A4
    t = Table()
    t.family([1000, 1010], ratio=1.0)
    t.family([5000, 9000], e=0.02, ratio=2.0)           # e < T_E, bad depth
    t.family([12000, 15000], e=0.30, ratio=2.0)         # e in range, bad depth
    t.family([20000, 22000], e=0.30, ratio=1.0)         # e in range, good depth
    t.family([25000, 27000], e=0.95, ratio=1.0)         # e above REFERENCE_PATH_CEILING
    rows, _ = run(t, idx=index(1000, 1010), opts=Options.previous().with_(check_scope="path-dropped"))
    st = status(rows)
    check("A4 path-dropped: e < T_E is not checked", st[5000] == "ADDED", st)
    check("A4 path-dropped: e in [T_E, REFERENCE_PATH_CEILING] is checked", st[12000] == "DROPPED_CHECK" and st[20000] == "ADDED", st)
    check("A4 path-dropped: e above REFERENCE_PATH_CEILING stays dropped by the reference-path filter", st[25000] == "DROPPED_V2", st)
    rows, _ = run(t, idx=index(1000, 1010))
    st = status(rows)
    check("A4 unconfirmed: e < T_E is checked too", st[5000] == "DROPPED_CHECK" and st[20000] == "ADDED", st)

    # ---------------------------------------------------------------- A5
    t = Table()
    t.family([5000, 9000], ratio=1.0)
    t.family([12000, 15000], ratio=1.6)
    rows, summ = run(t, idx=index())
    check("A5 no library-confirmed side: the tolerance comes from the layer's own two-sided calls",
          summ["tolerance"]["source"] == "own-calls" and summ["tolerance"]["sides"] == 8, summ["tolerance"])
    st = status(rows)
    check("A5 own-calls tolerance (p99 of 1.0 and 1.6) keeps both families", st[5000] == "ADDED" and st[12000] == "ADDED", st)
    t = Table()
    t.family([1000, 1010], ratio=1.0)
    t.family([5000, 9000], ratio=1.0)
    rows, summ = run(t, idx=index(1000, 1010), opts=Options.previous().with_(check_scope="unconfirmed", tol_source="own-calls"))
    check("A5 --tol-source own-calls ignores the library-confirmed calls", summ["tolerance"]["source"] == "own-calls")

    # ---------------------------------------------------------------- A6
    f = Table()
    f.family([1000, 1010], ratio=1.0)
    f.family([30000, 30500], ratio=1.0)                 # first-pass family call at 30000
    h = Table()
    h.family([5000, 9000], e=0.30, ratio=1.0)          # het site absent from the first pass
    h.family([12000, 15000], e=0.30, ratio=2.0)        # fails the check
    h.family([30100, 41000], e=0.30, ratio=1.0)        # one place within 500 bp of a first-pass call: that call is NEAR_FIRST_PASS
    h.family([45000, 47000], e=0.0, ratio=1.0)         # no allele-share evidence
    h.family([50000, 52000], e=0.95, ratio=1.0)
    rows, _ = run(f, h, idx=index(1000, 1010))
    sh = status(rows, "second")
    check("A6 second scan: new het site ADDED", sh[5000] == "ADDED" and sh[9000] == "ADDED", sh)
    check("A6 second scan: fails the check -> DROPPED_CHECK", sh[12000] == "DROPPED_CHECK", sh)
    check("A6 second scan: within 500 bp of a first-pass family call -> NEAR_FIRST_PASS (the family's other call is judged apart)",
          sh[30100] == "NEAR_FIRST_PASS" and sh[41000] == "ADDED", sh)
    check("A6 second scan: e < T_E -> NOT_HET; e > REFERENCE_PATH_CEILING -> DROPPED_V2", sh[45000] == "NOT_HET" and sh[50000] == "DROPPED_V2", sh)
    added = AS.added_calls_tsv_bytes(rows).decode().splitlines()
    check("A6 added calls name second-scan families LF_cluster_h<N>", any("LF_cluster_h" in l for l in added) and
          any("\tLF_cluster_0\t" in l for l in added) or len([l for l in added if "LF_cluster_h" in l]) >= 2, added[:4])
    try:
        AS.assemble(f.get(), h.get(), index(), Options.previous(), identity_filter=False, reference=GENOME)
        ok = False
    except InputContractError:
        ok = True
    check("A6 a second scan with check scope none is refused", ok)

    # ---------------------------------------------------------------- A7
    base = rnd(220)
    f = Table()
    f.family([1000, 1010], ratio=1.0)
    h = Table()
    cl = h.family([5000, 9000], e=0.30, ratio=1.0, lseg=base)       # two trusted members, identical segments
    # two more places in the same family from untrusted calls (e below T_E) with diverged segments
    cl_row = [r for r in h.rows if r["cluster"] == cl]
    for bL, d in ((13000, 0.12), (17000, 0.12)):
        ls = mutate(base, d)
        row = dict(cl_row[0], bL=str(bL), bR=str(bL - 5), first=str(bL - 10), last=str(bL), fracL="0.0", fracR="0.0",
                   lseg_len=str(len(ls)))
        h.rows.append(row)
        h.segs[("chrS", str(bL), cl)] = dict(lseg=ls, rseg=h.segs[("chrS", "5000", cl)]["rseg"])
    # make the right segment length match (rseg of the first place reused)
    for r in h.rows:
        r["rseg_len"] = str(len(h.segs[("chrS", r["bL"], r["cluster"])]["rseg"]))
        r["lseg_len"] = str(len(h.segs[("chrS", r["bL"], r["cluster"])]["lseg"]))
    rows, _ = run(f, h, idx=index(1000, 1010), identity_filter=True)
    own = CL.family_similarity(h.rows, h.segs)
    own_c = own[cl][0]
    sh = status(rows, "second")
    check("A7 setup: the family's own the family-identity filter drops it (median identity below the limit)", own_c <= CL.RULE_C_MAX, round(own_c, 3))
    check("A7 family-identity filter on trusted members: the two trusted calls are kept", sh[5000] == "ADDED" and sh[9000] == "ADDED", sh)
    check("A7 the untrusted calls (e below T_E) are NOT_HET", sh[13000] == "NOT_HET" and sh[17000] == "NOT_HET", sh)

    # ---------------------------------------------------------------- A8
    t = Table()
    t.family([1000, 1010], ratio=1.0)
    t.family([5000, 9000, 12000], ratio=1.0)           # family inside copies
    t.family([20000, 22000, 24000], ratio=1.0)         # family outside copies
    t.family([40000, 42000, 44000, 46000, 48000], ratio=1.0)   # one of five inside
    bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
    for s in (4900, 8900, 11900, 39900):
        bed.write(f"chrS\t{s}\t{s + 400}\tTE\t0\t+\n")
    # TE share: add a big block elsewhere so that the share is about 10%
    bed.write("chrS\t50000\t55000\tTE\t0\t+\n")
    bed.close()
    cx = CopyIndex.from_file(bed.name, contigs={"chrS"})
    share = cx.share(60_000)
    rows, summ = run(t, idx=index(1000, 1010), copies=cx, genome_bp=60_000, opts=Options.previous().with_(check_scope="unconfirmed", resident_guard="family"))
    st = status(rows)
    check("A8 family mode: calls inside copies of an enriched family are dropped", st[5000] == "DROPPED_RESIDENT_COPY"
          and st[9000] == "DROPPED_RESIDENT_COPY" and st[12000] == "DROPPED_RESIDENT_COPY", (round(share, 3), st))
    check("A8 family mode: a family outside copies, and the one inside call of a family below the enrichment limit, stay",
          st[20000] == "ADDED" and st[40000] == "ADDED", st)
    rows, _ = run(t, idx=index(1000, 1010), copies=cx, genome_bp=60_000, opts=Options.previous().with_(check_scope="unconfirmed", resident_guard="per-call"))
    st = status(rows)
    check("A8 per-call mode: every call inside a copy is dropped", st[40000] == "DROPPED_RESIDENT_COPY" and st[42000] == "ADDED", st)
    rows, summ = run(t, idx=index(1000, 1010), opts=Options.previous().with_(check_scope="unconfirmed", resident_guard="family"))
    check("A8 no annotation: the test is off", status(rows)[5000] == "ADDED" and "no repeat annotation" in str(summ["resident_guard"]))

    # ---------------------------------------------------------------- A9
    t = Table()
    t.family([5000, 9000])
    for r in t.rows:
        for k in ("chkL_T", "chkL_junction", "chkR_T", "chkR_junction"):
            r.pop(k)
    try:
        run(t)
        ok = False
    except InputContractError:
        ok = True
    check("A9 tables without recorded junctions are refused when a check is asked for", ok)
    # ---------------------------------------------------------------- A10
    def crafted(sites, frozen, body, rates):
        """A family of two-sided calls: member k has the left segment mutate(body, rates[k]). Its extension stopped at a plain-extension end
        (support) for the indices in ``frozen`` and at an anchored end (a_x) for the others."""
        t = Table()
        cl = t.family(sites)
        for k, (row, bL) in enumerate(zip(t.rows, sites)):
            ls = body if rates[k] == 0 else mutate(body, rates[k])
            t.segs[("chrS", str(bL), cl)]["lseg"] = ls
            row["lseg_len"] = str(len(ls))
            row["lstop"] = row["rstop"] = "support" if k in frozen else "a_x"
        return t, cl

    body = rnd(200)
    sites = [1000, 3000, 5000, 7000, 9000]
    pri = Options.previous().with_(first_pass_priority=True, check_scope="none")
    for min_ in (0, 3):
        t, cl = crafted(sites, {0, 1}, body, [0, 0.06, 0, 0, 0])
        rows, summ = run(t, opts=pri.with_(identity_min_comparisons=min_), identity_filter=True)
        st = {r.status for r in rows}
        c, used = rows[0].c, rows[0].c_used
        if min_ == 0:
            check("A10 limit 0: two frozen members give one comparison below the cutoff; the family is dropped",
                  st == {"DROPPED_RULE_C"} and used == 1 and c <= CL.RULE_C_MAX, (st, c, used))
        else:
            check("A10 limit 3: the frozen members give one comparison, so all members (10 comparisons, median above the cutoff) judge; the family is kept",
                  st == {"ADDED"} and used == 10 and c > CL.RULE_C_MAX and summ["identity_min_comparisons"] == 3, (st, c, used))
    t, cl = crafted(sites, {0, 1, 2}, body, [0, 0.06, 0.06, 0, 0])
    rows, _ = run(t, opts=pri.with_(identity_min_comparisons=3), identity_filter=True)
    check("A10 limit 3: three frozen members give three comparisons and are used (median below the cutoff: dropped)",
          {r.status for r in rows} == {"DROPPED_RULE_C"} and rows[0].c_used == 3, (rows[0].c, rows[0].c_used))
    two = [1000, 3000]
    for min_ in (0, 3):
        t, cl = crafted(two, {0, 1}, body, [0, 0.06])
        rows, summ = run(t, opts=pri.with_(identity_min_comparisons=min_), identity_filter=True)
        if min_ == 0:
            check("A10 limit 0: a family of two members is judged on its one comparison (dropped)", {r.status for r in rows} == {"DROPPED_RULE_C"})
        else:
            check("A10 limit 3: a family with one comparison is not judged (c is NaN), the filter keeps it, and the summary counts it",
                  {r.status for r in rows} == {"ADDED"} and rows[0].c != rows[0].c and rows[0].c_used == 1 and summ["families_not_judged"] == 1,
                  (rows[0].c, rows[0].c_used, summ.get("families_not_judged")))
    # trusted members of the haplotype scan
    for min_ in (0, 3):
        t, cl = crafted(sites, {0, 1}, body, [0, 0.06, 0, 0, 0])
        out = AS._trusted_identity_filter(types.SimpleNamespace(places=t.rows, segments=t.segs), list(range(len(t.rows))), True, False, None, min_)
        dropped = {v["dropped2"] for v in out.values()}
        check(f"A10 trusted members, limit {min_}: " + ("one frozen comparison below the cutoff drops the family" if min_ == 0
                                                     else "all members judge and keep it"),
              dropped == ({True} if min_ == 0 else {False}), dropped)

    print("RESULT", ("FAIL " + ",".join(FAIL)) if FAIL else "ALL PASS (%d/%d)" % (NCHK[0], NCHK[0]))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
