"""500-bp host-midpoint clustering, novelty masks and residency."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

from ..constants import MIN_FOLDS, MIN_PAIRS
from ..errors import AssetCompatibilityError
from ..family_table import family_pairs, is_two_column_table
from .codebook import fold_of

__all__ = [
    "CLUSTER_GAP",
    "cluster_mappings",
    "interval_gap",
    "family_map_rows",
    "load_family_map",
    "parse_repeatmasker",
    "unmapped_records",
    "relaxed_provisional_clusters",
    "annotate_corrected_clusters",
]

CLUSTER_GAP = 500
_FAR = 10**12


def cluster_mappings(rows: Iterable[Mapping], *, require_folds: bool, prefix: str, keep_all: bool = False) -> list[dict]:
    """Single-linkage on host midpoints (gap <= 500) per (family, chrom).

    Cluster ids are assigned 1.. in (family, chrom) order, then midpoint order,
    over emitted clusters only.
    """

    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        item = dict(row)
        item["midpoint"] = (int(row["start"]) + int(row["end"])) // 2
        grouped[(item["family"], item["chrom"])].append(item)
    out: list[dict] = []
    cid = 0
    for (family, chrom), values in sorted(grouped.items()):
        values.sort(key=lambda r: (int(r["midpoint"]), int(r["pair_id"])))

        def flush(members: list[dict]) -> None:
            nonlocal cid
            if not members:
                return
            pids = sorted({int(x["pair_id"]) for x in members})
            folds = sorted({fold_of(p) for p in pids})
            supported = len(pids) >= MIN_PAIRS and (not require_folds or len(folds) >= MIN_FOLDS)
            if not supported and not keep_all:
                return
            mids = [int(x["midpoint"]) for x in members]
            cid += 1
            out.append({
                "cluster_id": cid, "cluster_label": f"{prefix}_{cid:05d}", "family": family, "chrom": chrom,
                "start": min(mids), "end": max(mids), "pair_count": len(pids), "fold_count": len(folds),
                "event_2p2fold": int(len(pids) >= MIN_PAIRS and len(folds) >= MIN_FOLDS),
                "pair_ids": ";".join(map(str, pids)), "folds": ";".join(map(str, folds)),
            })

        current: list[dict] = []
        last = None
        for row in values:
            mid = int(row["midpoint"])
            if current and mid - int(last) > CLUSTER_GAP:
                flush(current)
                current = []
            current.append(row)
            last = mid
        flush(current)
    return out


TWO_COLUMN_RECORD = "te_record_id"


def interval_gap(a0: int, a1: int, b0: int, b1: int) -> int:
    if max(a0, b0) <= min(a1, b1):
        return 0
    return b0 - a1 if a1 < b0 else a0 - b1


def _fbte(text: str) -> str | None:
    match = re.search(r"\bFBte\d+\b", text or "", re.I)
    return None if not match else "FBte" + re.sub(r"\D", "", match.group(0))


def _canonical(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (text or "").strip().lower())


def load_family_map(rows: list[Mapping[str, str]]) -> tuple[dict[str, str], set[str]]:
    """``(record -> family, families)`` from the rows of ``family_map_rows``: FlyBase ids parsed as FlyBase ids
    (``FBte`` + digits), the record names of a two-column table (column ``te_record_id``) taken as written."""

    if not rows:
        raise AssetCompatibilityError("family map is empty")
    header = list(rows[0])
    plain = TWO_COLUMN_RECORD in header
    fb = TWO_COLUMN_RECORD if plain else next((h for h in header if "fbte" in h.lower()), None)
    fam = next((h for h in header if "family" in h.lower()), None)
    if not fb or not fam:
        raise AssetCompatibilityError(f"cannot identify family-map columns: {header}")
    mapping: dict[str, str] = {}
    for row in rows:
        te = (row.get(fb, "") or "").strip() if plain else _fbte(row.get(fb, ""))
        family = (row.get(fam, "") or "").strip()
        if te and family:
            mapping[te] = family
    if not mapping:
        raise AssetCompatibilityError("no family-map entries parsed")
    return mapping, set(mapping.values())


def family_map_rows(path: str | Path) -> list[dict[str, str]]:
    """Rows of the family map as dicts: the FlyBase table's own columns, or ``te_record_id`` / ``family`` for a two-column
    ``record<TAB>family`` table (``evitea.family_table``)."""

    if is_two_column_table(path):
        return [{TWO_COLUMN_RECORD: record, "family": family} for record, family in family_pairs(path)]
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unmapped_records(records, mapping: Mapping[str, str]) -> int:
    """How many library records have no family in the map."""

    return sum(1 for record in records if record not in mapping)


def parse_repeatmasker(path: str | Path, fbmap: Mapping[str, str], known: set[str]) -> list[dict]:
    norm: dict[str, str] = {}
    for family in known:
        norm.setdefault(_canonical(family), family)
    out: list[dict] = []
    with Path(path).open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            x = line.rstrip("\r\n").split("\t")
            if len(x) < 3:
                continue
            name = x[3] if len(x) >= 4 else ""
            te = _fbte(name) or next(
                (t for t in (name, name.split("#")[0], name.split("|")[0], name.split(":")[0]) if t in fbmap), None)
            family = fbmap.get(te, "") if te else ""
            if not family:
                for token in (name, name.split("#")[0], name.split("|")[0], name.split(":")[0]):
                    q = norm.get(_canonical(token))
                    if q:
                        family = q
                        break
            out.append({"chrom": x[0], "start": int(x[1]), "end": int(x[2]), "family": family})
    return out


def relaxed_provisional_clusters(
    primary_clusters: list[dict],
    old_primary_unique: list[dict],
    decoded_unique: list[dict],
    strict_new: Mapping[int, Mapping],
    decoded_host: Mapping[int, Mapping],
) -> list[dict]:
    """Delayed-aggregation novelty mask from strictly decoded fragments outside the primary set.

    Combine decoded UNIQUE_HIGH fragments with the primary set. Keep clusters
    with a decoded member that lie > 500 bp from every same-family primary
    cluster. Retain those with >= 2 masked-family-confirmed decoded
    UNIQUE_HIGH members.
    """

    combined = cluster_mappings(old_primary_unique + decoded_unique, require_folds=False, prefix="FRAGMENT_DECODED")
    by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for c in primary_clusters:
        by_key[(c["family"], c["chrom"])].append(c)
    strict_ids = set(strict_new)
    out: list[dict] = []
    for c in combined:
        members = sorted({int(x) for x in str(c["pair_ids"]).split(";") if x} & strict_ids)
        if not members:
            continue
        gap = min((interval_gap(int(c["start"]), int(c["end"]), int(b["start"]), int(b["end"]))
                   for b in by_key[(c["family"], c["chrom"])]), default=_FAR)
        if gap <= CLUSTER_GAP:
            continue
        confirmed = [
            pid for pid in members
            if str(strict_new[pid].get("masked_family_confirmed", "0")) == "1"
            and decoded_host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"
        ]
        if len(confirmed) < MIN_PAIRS:
            continue
        folds = sorted({fold_of(p) for p in confirmed})
        out.append({
            "cluster_id": c["cluster_id"], "cluster_label": f"RELAXED_{int(c['cluster_id']):05d}",
            "family": c["family"], "chrom": c["chrom"], "start": c["start"], "end": c["end"],
            "pair_count": len(confirmed), "fold_count": len(folds),
            "pair_ids": ";".join(map(str, confirmed)), "folds": ";".join(map(str, folds)),
        })
    return out


def annotate_corrected_clusters(clusters: list[dict], baseline: list[dict], repeats: list[dict]) -> list[dict]:
    """``novel`` = > 500 bp from every same-family baseline cluster;
    ``resident_compatible`` = <= 500 bp from a same-family RepeatMasker TE."""

    base: dict[tuple[str, str], list[dict]] = defaultdict(list)
    reps: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for b in baseline:
        base[(b["family"], b["chrom"])].append(b)
    for r in repeats:
        if r["family"]:
            reps[(r["family"], r["chrom"])].append(r)
    out = []
    for c in clusters:
        key, s, e = (c["family"], c["chrom"]), int(c["start"]), int(c["end"])
        nb = min((interval_gap(s, e, int(b["start"]), int(b["end"])) for b in base[key]), default=_FAR)
        nr = min((interval_gap(s, e, int(r["start"]), int(r["end"])) for r in reps[key]), default=_FAR)
        item = dict(c)
        item["nearest_baseline_same_family_bp"] = "" if nb == _FAR else nb
        item["nearest_repeatmasker_same_family_bp"] = "" if nr == _FAR else nr
        item["novel"] = nb > CLUSTER_GAP
        item["resident_compatible"] = nr <= CLUSTER_GAP
        item["novel_nonresident"] = item["novel"] and not item["resident_compatible"]
        out.append(item)
    return out
