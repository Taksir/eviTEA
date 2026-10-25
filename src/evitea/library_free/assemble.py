"""Family calls of the library-free layer: the reference-path filter, the family-identity filter, the one-site check, the resident-copy test,
and the join of the first pass with the haplotype scan.

Input: the grouped tables of each pass as text fields (``detector.parse_places`` / ``parse_segments``), so a saved detection can be
assembled again with other switches. The first pass is judged as ``calling.family_calls`` judges it (reference-path filter,
the family-identity filter, near an EviTEA call), with these changes by switch:

check scope   none         the previous version's rule: a call with reference-path fraction e >= T_E is dropped.
              path-dropped   a call that the reference-path filter drops (T_E <= e <= REFERENCE_PATH_CEILING) is kept when it passes the one-site check (``check.py``).
              unconfirmed  every library-free call the library layer does not make (no EviTEA call within 500 bp) must pass the
                           one-site check, and e <= REFERENCE_PATH_CEILING. e stays in the output as the allele-share column.
Calls are judged in this order: DROPPED_V2 (e above REFERENCE_PATH_CEILING, or e >= T_E with scope none), DROPPED_RULE_C, NEAR_EVITEA, DROPPED_CHECK,
DROPPED_RESIDENT_COPY, ADDED.

Haplotype scan. Its family calls are judged on their own grouping and join the first pass only where the first pass has no family call
(of any status) within 500 bp (MERGE_D), and only if they show reads of both alleles (T_E <= e <= REFERENCE_PATH_CEILING) and pass the one-site check.
The family-identity filter for them is decided on trusted members: the places of calls that lie within 500 bp of a first-pass family call or that pass
the check with T_E <= e <= REFERENCE_PATH_CEILING are grouped again among themselves, the family rule and the family-identity filter are applied to those clusters, and a
call that the family-identity filter dropped in its own grouping is kept when this second decision keeps it. Status of such a call: NOT_HET (no
allele-share evidence: left to the first pass), DROPPED_V2, NEAR_FIRST_PASS, DROPPED_CHECK, DROPPED_RULE_C, NEAR_EVITEA,
DROPPED_RESIDENT_COPY, ADDED.
"""

from __future__ import annotations

import bisect
import math
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import astuple, dataclass, fields

from ..errors import InputContractError
from . import RULE_C_MAX, T_E
from . import _p4_lib as P
from . import calling as CL
from . import check as CK
from . import priority as PR
from .extension import is_anchored
from .options import Options

__all__ = ["REFERENCE_PATH_CEILING", "FamilyCallX", "assemble", "added_calls_tsv_bytes", "family_calls_tsv_bytes", "COLUMNS_X"]

REFERENCE_PATH_CEILING = 0.9
NEAR_FIRST_BP = P.MERGE_D


@dataclass(frozen=True, slots=True)
class FamilyCallX(CL.FamilyCall):
    scan: str = "first"
    origin: str = ""
    check: str = ""
    ratio_L: str = ""
    ratio_R: str = ""
    inside_copy: str = ""
    family_inside: str = ""


BASE_COLUMNS = CL.FAMILY_CALL_COLUMNS
COLUMNS_X = tuple(f.name for f in fields(FamilyCallX))


@dataclass
class _Pass:
    scan: str
    places: list
    segments: dict
    calls: list            # indexes of family-call rows
    e: dict
    c: dict                # cluster -> (c, pairs, used)
    near: dict             # place index -> bool (EviTEA call within 500 bp)
    sides: dict            # place index -> list[SideInput]
    check: dict            # place index -> 1/0 (only for measured calls)


def _positions(rows):
    by = defaultdict(list)
    for r in rows:
        by[r["chrom"]].append(int(r["bL"]))
    return {c: sorted(v) for c, v in by.items()}


def _within(pos, chrom, x, d):
    v = pos.get(chrom)
    if not v:
        return False
    j = bisect.bisect_left(v, x - d)
    return j < len(v) and v[j] <= x + d


def _fnum(v):
    return CL._num(v)


def _judged(used, all_members, minimum, aligner):
    """(c, pairs, usable pairs) of a family from at least ``minimum`` usable member comparisons: the members ``used`` (the frozen-extension
    members when the first-pass priority picked them) when they give that many, otherwise all members; with fewer the family is not judged
    (c = NaN, which the family-identity filter never drops).  ``minimum`` 0 is the plain measure on ``used``."""
    r = CL.feat_c(used, aligner)
    if minimum == 0 or r[2] >= minimum:
        return r
    r = CL.feat_c(all_members, aligner) if all_members != used else r
    if r[2] >= minimum:
        return r
    return (float("nan"), r[1], r[2])


def _family_similarity(places, segments, priority, aligner, minimum=0):
    members = PR.cluster_members(places, segments, priority=priority)
    allm = PR.cluster_members(places, segments, priority=False) if priority and minimum else members
    wanted = sorted({p["cluster"] for p in places if CL.is_family_call(p)}, key=int)
    return {c: _judged(members.get(c, []), allm.get(c, []), minimum, aligner) for c in wanted}


def _prepare(scan, places, segments, index, priority, aligner, minimum=0):
    calls = [i for i, p in enumerate(places) if CL.is_family_call(p)]
    e = {i: CL.e_value(places[i]["fracL"], places[i]["fracR"], places[i]["rstat"]) for i in calls}
    c = _family_similarity(places, segments, priority, aligner, minimum)
    near = {i: index.near(places[i]["chrom"], int(places[i]["bL"])) for i in calls}
    sides = {i: CK.side_inputs(places[i]) for i in calls}
    return _Pass(scan, places, segments, calls, e, c, near, sides, {})


def _identity_filter_of(pa, i):
    c = pa.c.get(pa.places[i]["cluster"], (float("nan"), 0, 0))
    return c, CL.rule_c_dropped(c[0])


def _need_columns(pa):
    for i in pa.calls:
        if "chkL_T" not in pa.places[i] and "chkR_T" not in pa.places[i]:
            raise InputContractError("the places table has no one-site check columns (chkL_T, chkR_T); "
                                     "run the detection with --check-scope path-dropped or unconfirmed")
        return


def _tolerance(first: _Pass, opts: Options, identity_filter: bool):
    lib, own = [], []
    for i in first.calls:
        for s in first.sides[i]:
            if first.near[i]:
                lib.append(s.ratio)
        p = first.places[i]
        c, dc = _identity_filter_of(first, i)
        e = first.e[i]
        if p["two_sided"] == "1" and not CL.v2_dropped(e) and not (identity_filter and dc):
            own.extend(s.ratio for s in first.sides[i])
    return CK.tolerance(lib, own, opts.tol_source)


def _measure(passes, ref_index, wanted, tol, opts):
    """Novelty for the wanted (pass, place index) calls; sets ``check`` on each pass."""
    junctions = {}
    for pa, idxs in wanted:
        for i in idxs:
            for s in pa.sides[i]:
                junctions[(pa.scan, i, s.side)] = s.junction
    nov = CK.novelty(junctions, ref_index) if junctions else {}
    t = tol.value(opts.tol) if tol is not None else math.inf
    for pa, idxs in wanted:
        for i in idxs:
            novs = {s.side: nov[(pa.scan, i, s.side)] for s in pa.sides[i]}
            pa.check[i] = CK.call_pass(pa.sides[i], novs, t)
    return len(junctions)


def _row(pa, i, status, extra=None):
    p = pa.places[i]
    c, dc = _identity_filter_of(pa, i)
    e = pa.e[i]
    base = (i, p["chrom"], int(p["bL"]), p["bR"], p["tsd"], p["cluster"], p["fracL"], p["fracR"], p["rstat"], e, c[0], c[1], c[2],
            CL.v2_dropped(e), dc, pa.near[i], status)
    ratios = {s.side: s.ratio for s in pa.sides[i]}
    chk = pa.check.get(i)
    x = dict(scan=pa.scan, origin=p.get("origin", "") or "", check="" if chk is None else str(chk),
             ratio_L=("" if "L" not in ratios else f"{ratios['L']:.4f}"), ratio_R=("" if "R" not in ratios else f"{ratios['R']:.4f}"))
    if extra:
        x.update(extra)
    return FamilyCallX(*base, **x)


def _first_status(pa, i, opts, identity_filter, scope):
    e = pa.e[i]
    dv = CL.v2_dropped(e)
    _c, dc = _identity_filter_of(pa, i)
    if scope == "none":
        if dv:
            return "DROPPED_V2"
    elif not math.isnan(e) and e > REFERENCE_PATH_CEILING:
        return "DROPPED_V2"
    if identity_filter and dc:
        return "DROPPED_RULE_C"
    if pa.near[i]:
        return "NEAR_EVITEA"
    if (scope == "unconfirmed") or (scope == "path-dropped" and dv):
        if pa.check.get(i) != 1:
            return "DROPPED_CHECK"
    return "ADDED"


def _trusted_identity_filter(second: _Pass, trusted: list[int], priority: bool, polya_tsd: bool, aligner, minimum=0):
    """Group the trusted places among themselves; (cluster2, family vtsd2, c2 tuple, the family-identity filter dropped2) per place index."""
    G, back = [], []
    for i in trusted:
        p = second.places[i]
        d = second.segments.get((p["chrom"], p["bL"], p["cluster"]), {})
        keep = [x for x in p["seg_keep"].split(",") if x]
        segs = [d[k] for k, side in (("lseg", "L"), ("rseg", "R")) if side in keep and k in d]
        tsd = p.get("tsd", "")
        stop = {"lseg": p.get("lstop"), "rseg": p.get("rstop")}
        allsegs = [d[k] for k in ("lseg", "rseg") if len(d.get(k, "")) >= PR.SEG_MIN]
        fzsegs = [d[k] for k in ("lseg", "rseg") if len(d.get(k, "")) >= PR.SEG_MIN and not is_anchored(stop[k])]
        G.append(dict(segs=segs, two_sided=int(p["two_sided"] or 0), tsd=(int(tsd) if tsd not in ("", None) else None),
                      lseg=d.get("lseg", ""), rseg=d.get("rseg", ""), all_segs=allsegs, fz_segs=fzsegs))
        back.append(i)
    if not G:
        return {}
    cl, _ = P.group(G)
    fam = PR.families(G, cl, polya_tsd=polya_tsd)
    mem = defaultdict(list)
    for g, c in zip(G, cl):
        mem[c].append(g)
    c2 = {}
    for c, f in fam.items():
        if not f["vtsd"]:
            continue
        fz = [g["fz_segs"] for g in mem[c]]
        every = [g["all_segs"] for g in mem[c]]
        use = fz if (priority and sum(1 for m in fz if m) >= 2) else every
        c2[c] = _judged(use, every, minimum, aligner)
    out = {}
    for k, i in enumerate(back):
        c = cl[k]
        f = fam[c]
        cc = c2[c][0] if f["vtsd"] else float("nan")
        out[i] = dict(cluster2=c, vtsd2=f["vtsd"], dropped2=(not math.isnan(cc)) and CL.rule_c_dropped(cc) if f["vtsd"] else False)
    return out


def assemble(first, second, index: CL.CallIndex, opts: Options, *, identity_filter: bool, reference=None, copies=None, genome_bp: int = 0,
             aligner: Callable | None = None):
    """(rows, summary). ``first`` / ``second`` = (places rows, segments dict) of each pass (``second`` None without a haplotype
    scan); ``reference`` = dict contig -> sequence (needed for the check); ``copies`` = ``resident.CopyIndex`` (None: no
    resident-copy test); ``genome_bp`` = bases of the reference (the TE share's denominator)."""
    scope = opts.check_scope
    if second is not None and scope == "none":
        raise InputContractError("the haplotype scan joins calls through the one-site check; use --check-scope path-dropped or "
                                 "unconfirmed with --het-scan")
    minimum = opts.identity_min_comparisons
    F = _prepare("first", first[0], first[1], index, opts.first_pass_priority, aligner, minimum)
    H = _prepare("second", second[0], second[1], index, opts.first_pass_priority, aligner, minimum) if second is not None else None
    summary: dict = {"check_scope": scope, "first_pass_family_calls": len(F.calls),
                     "second_scan_family_calls": len(H.calls) if H else 0,
                     "families_with_similarity": len(F.c) + (len(H.c) if H else 0)}
    if minimum:
        summary["identity_min_comparisons"] = minimum
        summary["families_not_judged"] = sum(math.isnan(v[0]) for v in F.c.values()) + (sum(math.isnan(v[0]) for v in H.c.values()) if H else 0)
    tol = None
    if scope != "none":
        if F.calls:
            _need_columns(F)
        if H is not None and H.calls:
            _need_columns(H)
        if reference is None:
            raise InputContractError("the one-site check needs the reference sequence")
        tol = _tolerance(F, opts, identity_filter)
        summary["tolerance"] = None if tol is None else dict(p95=tol.p95, p99=tol.p99, source=tol.source, sides=tol.n_sides, used=opts.tol)
        wanted = []
        need = []
        for i in F.calls:
            e = F.e[i]
            _c, dc = _identity_filter_of(F, i)
            if (not math.isnan(e) and e > REFERENCE_PATH_CEILING) or (identity_filter and dc) or F.near[i]:
                continue
            if scope == "unconfirmed" or (scope == "path-dropped" and CL.v2_dropped(e)):
                need.append(i)
        wanted.append((F, need))
        if H is not None:
            hn = [i for i in H.calls if not math.isnan(H.e[i]) and T_E <= H.e[i] <= REFERENCE_PATH_CEILING]
            wanted.append((H, hn))
        summary["junctions_measured"] = _measure([F] + ([H] if H else []), CK.RefIndex(reference), wanted, tol, opts)
    # first pass statuses
    status = {i: _first_status(F, i, opts, identity_filter, scope) for i in F.calls}
    # second scan
    hstatus, hextra = {}, {}
    if H is not None:
        fpos = _positions([F.places[i] for i in F.calls])
        member = {i: _within(fpos, H.places[i]["chrom"], int(H.places[i]["bL"]), NEAR_FIRST_BP) for i in H.calls}
        trusted = [i for i in H.calls if member[i] or (T_E <= H.e[i] <= REFERENCE_PATH_CEILING and H.check.get(i) == 1)]
        rct = _trusted_identity_filter(H, trusted, opts.first_pass_priority, opts.polya_tsd, aligner, minimum) if identity_filter else {}
        for i in H.calls:
            e = H.e[i]
            _c, dc = _identity_filter_of(H, i)
            if math.isnan(e) or e < T_E:
                hstatus[i] = "NOT_HET"
            elif e > REFERENCE_PATH_CEILING:
                hstatus[i] = "DROPPED_V2"
            elif member[i]:
                hstatus[i] = "NEAR_FIRST_PASS"
            elif H.check.get(i) != 1:
                hstatus[i] = "DROPPED_CHECK"
            elif identity_filter and dc and not (i in rct and rct[i]["vtsd2"] and not rct[i]["dropped2"]):
                hstatus[i] = "DROPPED_RULE_C"
            elif H.near[i]:
                hstatus[i] = "NEAR_EVITEA"
            else:
                hstatus[i] = "ADDED"
    # resident-copy test on the calls that would be added
    inside_info: dict = {}
    if opts.resident_guard != "off" and copies is not None:
        from . import resident as RS

        share = copies.share(genome_bp)
        summary["resident_guard"] = dict(mode=opts.resident_guard, te_share=round(share, 4), threshold=round(RS.FAMILY_ENRICHMENT * share, 4))
        for pa, st in ((F, status), (H, hstatus)):
            if pa is None:
                continue
            fr = RS.resident_flags([pa.places[i] for i in pa.calls], copies, share, opts.resident_guard)
            for i, (ins, fam, drop) in zip(pa.calls, fr):
                inside_info[(pa.scan, i)] = (ins, fam)
                if drop and st[i] == "ADDED":
                    st[i] = "DROPPED_RESIDENT_COPY"
    elif opts.resident_guard != "off":
        summary["resident_guard"] = "off: no repeat annotation"
    rows = []
    for pa, st in ((F, status), (H, hstatus)):
        if pa is None:
            continue
        for i in pa.calls:
            ins = inside_info.get((pa.scan, i))
            extra = {} if ins is None else dict(inside_copy=str(int(ins[0])), family_inside=f"{ins[1]:.4f}")
            rows.append(_row(pa, i, st[i], extra))
    summary["statuses"] = {s: dict(sorted(_count([r for r in rows if r.scan == s]).items())) for s in ("first", "second")}
    return tuple(rows), summary


def _count(rows):
    out = defaultdict(int)
    for r in rows:
        out[r.status] += 1
    return out


def family_calls_tsv_bytes(rows, *, extended: bool) -> bytes:
    if not extended:
        return CL.family_calls_tsv_bytes([CL.FamilyCall(*astuple(r)[:len(BASE_COLUMNS)]) for r in rows])
    return CL._tsv(COLUMNS_X, ([CL._cell(v) for v in astuple(r)] for r in rows))


def added_calls_tsv_bytes(rows, *, second_prefix: str = "h") -> bytes:
    """The added calls of both passes; the call interval is [bL, bL + 1), as in the previous scored union. Second-scan families are
    named ``LF_cluster_h<N>``."""
    seen: dict[str, int] = defaultdict(int)
    body = []
    for r in rows:
        if r.status != "ADDED":
            continue
        base = f"LF_{r.chrom}_{r.bL}"
        seen[base] += 1
        eid = base if seen[base] == 1 else f"{base}_{seen[base]}"
        fam = f"LF_cluster_{second_prefix}{r.cluster}" if r.scan == "second" else f"LF_cluster_{r.cluster}"
        body.append([eid, r.chrom, r.bL, r.bL + 1, fam, r.cluster, r.bL, r.bR, r.tsd, CL._cell(r.e), CL._cell(r.c)])
    return CL._tsv(CL.ADDED_CALL_COLUMNS, body)
