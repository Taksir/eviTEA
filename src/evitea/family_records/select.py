"""Which families get a record, and the record guard.

A family is flagged when it has at least FAMILY_MIN calls away from the library layer's calls and more away than near. The calls
counted depend on the flag rule:

kept            away = calls that the filters kept or the resident-copy test, the single-copy rule or the local heterozygous guard
                dropped last (ADDED, DROPPED_RESIDENT_COPY, DROPPED_SINGLE_COPY, DROPPED_HOM_REGION);
                near = calls near a library call (NEAR_EVITEA). A haplotype-scan call near a first-pass call (NEAR_FIRST_PASS) is
                neither: the first-pass call may itself belong to a family the library lacks.
before-filters  away = calls that are not near a library call and that the family-identity filter does not drop, regardless of the
                reference-path filter, the one-site check and the resident-copy test; near as above.
both            a family flagged by either rule.

Cluster numbers belong to a scan, so a family is named by (scan, cluster).

Record guard    a family whose calls lie inside annotated reference TE copies more than twice as often as the genome's TE fraction
                (and in more than half of its calls) is a family of old copies, not of non-reference insertions: it gets no record.
"""

from __future__ import annotations

from collections import defaultdict

from ..library_free.resident import FAMILY_ENRICHMENT, CopyIndex

FAMILY_MIN = 2
GUARD_FLOOR = 0.5
FLAG_RULES = ("kept", "before-filters", "both")
FlagRule = str

AWAY_KEPT = ("ADDED", "DROPPED_RESIDENT_COPY", "DROPPED_SINGLE_COPY", "DROPPED_HOM_REGION")
NEAR_STATUS = ("NEAR_EVITEA",)
NOT_COUNTED = ("NOT_HET", "NEAR_FIRST_PASS")        # near a first-pass call, which may itself be a call of a family the library lacks


def truthy(v) -> bool:
    return str(v).strip().lower() in ("1", "true")


def scan_of(row) -> str:
    return row.get("scan") or "first"


def _counts(rows, rule):
    if rule == "kept":
        away = sum(r["status"] in AWAY_KEPT for r in rows)
        near = sum(r["status"] in NEAR_STATUS for r in rows)
    else:
        away = sum(r["status"] not in NOT_COUNTED and r["status"] not in NEAR_STATUS and not truthy(r["near_evitea"])
                   and not truthy(r["rule_c_dropped"]) for r in rows)
        near = sum(r["status"] in NEAR_STATUS or truthy(r["near_evitea"]) for r in rows)
    return away, near


def flagged_families(rows, rule: FlagRule = "kept"):
    """{(scan, cluster): {'rows': [...], 'kept': (away, near), 'before-filters': (away, near), 'flagged_by': [...]}} for the
    flagged families, in scan then cluster order."""
    if rule not in FLAG_RULES:
        raise ValueError(f"flag rule {rule!r}: use one of {FLAG_RULES}")
    byc = defaultdict(list)
    for r in rows:
        byc[(scan_of(r), r["cluster"])].append(r)
    out = {}
    for key in sorted(byc, key=lambda k: (k[0] != "first", int(k[1]))):
        rr = byc[key]
        counts = {m: _counts(rr, m) for m in ("kept", "before-filters")}
        by = [m for m in ("kept", "before-filters") if counts[m][0] >= FAMILY_MIN and counts[m][0] > counts[m][1]]
        use = by if rule == "both" else [m for m in by if m == rule]
        if use:
            out[key] = dict(rows=rr, flagged_by=use, **counts)
    return out


def guard_threshold(te_share: float | None) -> float:
    """The family-inside fraction above which a family gets no record. An unknown TE fraction counts as 0 (the floor alone)."""
    return max(GUARD_FLOOR, FAMILY_ENRICHMENT * (te_share or 0.0))


def inside_share(rows, copies: CopyIndex | None) -> float | None:
    """Fraction of the family's calls (any status) whose break lies inside an annotated copy: measured with the annotation when
    given, else read from the family-calls table's ``family_inside`` column; None when neither exists."""
    if copies is not None:
        k = sum(copies.inside(r["chrom"], int(r["bL"])) for r in rows)
        return k / len(rows)
    vals = [r.get("family_inside", "") for r in rows]
    if all(v not in ("", None) for v in vals):
        return float(vals[0])
    return None
