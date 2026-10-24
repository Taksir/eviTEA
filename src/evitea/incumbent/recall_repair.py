"""Candidate-recall repair for under-covered panel groups (lineage stage G1H-R).

Groups whose pair-aware coverage stays below 95% are revisited using wider
exact G1A anchors of their family (family-unique, not low-complexity,
nonempty).  Anchors of width <= 1024 act as a geometric proxy for what the
compact screen already covers; for each group the smallest cap of
1536..16384 whose anchors can reach 98% window coverage (else 16384) bounds a
greedy choice by new windows per interval row.  Only pairs not already routed
by the pair-aware panel are enumerated; new one-mate pairs send their opposite
mate to Bowtie2 (``pair_{id}``, ``-k 2 --very-sensitive``) and UNIQUE_HIGH
mates are clustered (gap 500 bp, >= 2 pairs).  Supported pairs are added to the
working candidate pairs.

Groups are matched to anchors by family label, so every record of a family
contributes anchor coordinates to each same-family group, exactly as frozen.
"""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..errors import BackendCompatibilityError
from .anchor_census import AnchorCatalogue
from .host_mate_gate import HostMapping
from .pair_aware_panel import PairAwareGroup

__all__ = [
    "MAX_RESCUE_INTERVAL_MASS",
    "RescueAnchor",
    "RescuePairEvidence",
    "check_rescue_interval_mass",
    "plan_rescue_anchors",
    "enumerate_rescue_pairs",
    "rescue_host_requests",
    "RESCUE_HOST_OPTIONS",
    "RESCUE_READ_NAME_PROFILES",
    "rescue_query_name",
    "parse_rescue_host_sam",
    "supported_rescue_pairs",
]

UNDERCOVERAGE_THRESHOLD = 0.95
RESCUE_TARGET_COVERAGE = 0.98
BASE_INTERVAL_CAP = 1024
CAP_LADDER = (1536, 2048, 3072, 4096, 8192, 16384)
MAX_RESCUE_INTERVAL_MASS = 1_500_000  # frozen default; a resource limit, settable with --max-rescue-interval-mass
READ_WINDOW = 125
RESCUE_MAPQ = 20
RESCUE_CLUSTER_GAP = 500
RESCUE_CLUSTER_SUPPORT = 2
RESCUE_HOST_OPTIONS = ("--very-sensitive", "-k", "2", "--no-unal")

# Bowtie2 seeds its per-read tie-breaking from the query name, so the rescue
# mapping's names are part of the reproduction identity:
# * ``pair``    - ``pair_{id}``: production default; the AKA017 historical run.
# * ``numeric`` - ``{id}``: reproduction only; the ISO1 historical run (Linux /
#   WSL Bowtie2; the research notes report that native Windows Bowtie2 2.5.5
#   rejects all-numeric FASTA names).
RESCUE_READ_NAME_PROFILES = ("pair", "numeric")


def rescue_query_name(pair_id: int, profile: str = "pair") -> str:
    if profile == "pair":
        return f"pair_{pair_id}"
    if profile == "numeric":
        return str(int(pair_id))
    raise ValueError(f"unknown rescue read-name profile {profile!r}")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(slots=True)
class RescueAnchor:
    query_id: str
    query_seq: str
    k: int
    width: int
    family: str


@dataclass(frozen=True, slots=True)
class RescuePairEvidence:
    pair_id: int
    mate1_hit: bool
    mate2_hit: bool
    families: tuple[str, ...]


def _window_interval(start: int, k: int, n_windows: int) -> tuple[int, int] | None:
    lo = max(0, start + k - READ_WINDOW)
    hi = min(n_windows - 1, start)
    return None if lo > hi else (lo, hi)


def _mean(mask: np.ndarray) -> float:
    return float(mask.mean()) if len(mask) else 0.0


def check_rescue_interval_mass(anchors: Sequence[RescueAnchor], limit: int) -> None:
    """Fail closed when the selected anchors' msBWT interval rows exceed the limit (a resource limit only)."""

    mass = sum(anchor.width for anchor in anchors)
    if mass > limit:
        raise BackendCompatibilityError(
            f"rescue anchor interval mass {mass:,} exceeds the {limit:,}-row budget (--max-rescue-interval-mass)"
        )


def plan_rescue_anchors(
    catalogue: AnchorCatalogue,
    groups: Sequence[PairAwareGroup],
    *,
    max_interval_mass: int = MAX_RESCUE_INTERVAL_MASS,
) -> list[RescueAnchor]:
    queries = catalogue.by_id()
    eligible: dict[tuple[str, str], list[tuple[int, int, object]]] = defaultdict(list)
    for source in catalogue.sources:
        if not source.queryable:
            continue
        query = queries.get(source.query_id)
        if query is None or query.width <= 0 or query.low_complexity or query.source_family_count != 1:
            continue
        eligible[(source.family, source.orientation)].append((source.start, source.k, query))

    selected: dict[str, RescueAnchor] = {}
    low = [g for g in groups if g.selected_coverage < UNDERCOVERAGE_THRESHOLD]
    for group in sorted(low, key=lambda g: (g.selected_coverage, g.family, g.orientation)):
        n = group.window_count
        rows = eligible.get((group.family, group.orientation), [])
        covered = np.zeros(n, dtype=np.bool_)
        for start, k, query in rows:
            if query.width <= BASE_INTERVAL_CAP:
                iv = _window_interval(start, k, n)
                if iv is not None:
                    covered[iv[0] : iv[1] + 1] = True
        selected_cap = CAP_LADDER[-1]
        achievable = _mean(covered)
        for cap in CAP_LADDER:
            trial = covered.copy()
            for start, k, query in rows:
                if BASE_INTERVAL_CAP < query.width <= cap:
                    iv = _window_interval(start, k, n)
                    if iv is not None:
                        trial[iv[0] : iv[1] + 1] = True
            achievable = _mean(trial)
            if achievable + 1e-15 >= RESCUE_TARGET_COVERAGE:
                selected_cap = cap
                break
        target = min(RESCUE_TARGET_COVERAGE, achievable)
        windows: dict[str, list[np.ndarray]] = {}
        first: dict[str, object] = {}
        for start, k, query in rows:
            if not BASE_INTERVAL_CAP < query.width <= selected_cap:
                continue
            first.setdefault(query.query_id, query)
            iv = _window_interval(start, k, n)
            if iv is not None:
                windows.setdefault(query.query_id, []).append(np.arange(iv[0], iv[1] + 1, dtype=np.int32))
        options = [
            (qid, first[qid], np.unique(np.concatenate(parts))) for qid, parts in windows.items()
        ]
        while _mean(covered) + 1e-15 < target:
            best = best_key = None
            for position, (qid, query, idx) in enumerate(options):
                gain = int((~covered[idx]).sum())
                if gain <= 0:
                    continue
                width = max(1, query.width)
                rank = (gain / width, gain, -width, query.k, qid)
                if best_key is None or rank > best_key:
                    best_key, best = rank, position
            if best is None:
                break
            qid, query, idx = options.pop(best)
            covered[idx] = True
            anchor = selected.get(qid)
            if anchor is None:
                selected[qid] = RescueAnchor(qid, query.query_seq.strip().upper().replace("U", "T"), query.k, query.width, group.family)
            elif anchor.family != group.family:
                raise BackendCompatibilityError(f"family-unique rescue anchor {qid} selected for two families")
    anchors = sorted(selected.values(), key=lambda a: (a.family, a.width, a.query_id))
    check_rescue_interval_mass(anchors, max_interval_mass)
    return anchors


def enumerate_rescue_pairs(anchors: Sequence[RescueAnchor], routed_pairs: np.ndarray, index) -> dict[int, RescuePairEvidence]:
    routed = set(int(x) for x in routed_pairs.tolist())
    masks: dict[int, int] = {}
    families: dict[int, set[str]] = defaultdict(set)
    for anchor in anchors:
        dollars = index.dollar_ids(anchor.query_seq, cache=False)
        if len(dollars) != anchor.width:
            raise BackendCompatibilityError(f"interval width changed for rescue anchor {anchor.query_id}")
        files, pairs = index.file_and_pair(dollars)
        for file_id, pair_id in zip(files.tolist(), pairs.tolist()):
            if pair_id in routed:
                continue
            masks[pair_id] = masks.get(pair_id, 0) | (1 << file_id)
            families[pair_id].add(anchor.family)
    return {
        pid: RescuePairEvidence(pid, bool(masks[pid] & 1), bool(masks[pid] & 2), tuple(sorted(families[pid])))
        for pid in sorted(masks)
    }


def rescue_host_requests(evidence: Mapping[int, RescuePairEvidence]) -> dict[int, int]:
    out: dict[int, int] = {}
    for pid, row in evidence.items():
        if row.mate1_hit and not row.mate2_hit:
            out[pid] = 1
        elif row.mate2_hit and not row.mate1_hit:
            out[pid] = 0
    return out


def parse_rescue_host_sam(sam: Path, expected: Iterable[int]) -> dict[int, HostMapping]:
    records: dict[int, list[tuple]] = defaultdict(list)
    with Path(sam).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            pair_id = int(qname[len("pair_") :]) if qname.startswith("pair_") else int(qname)
            flag = int(fields[1])
            if flag & 0x4:
                continue
            pos = int(fields[3])
            span = sum(int(n) for n, op in _CIGAR.findall(fields[5]) if op in "MDN=X") if fields[5] != "*" else 0
            midpoint = pos + ((span - 1) // 2 if span > 0 else 0)
            records[pair_id].append(
                (fields[2], pos, midpoint, int(fields[4]), "-" if flag & 0x10 else "+", fields[5], bool(flag & 0x100))
            )
    out: dict[int, HostMapping] = {}
    for pair_id in sorted(set(int(x) for x in expected)):
        rows = records.get(pair_id, [])
        if not rows:
            out[pair_id] = HostMapping(pair_id, "UNMAPPED", 0, "", None, None, None, "", "")
            continue
        primary = next((row for row in rows if not row[6]), rows[0])
        if len(rows) == 1 and primary[3] >= RESCUE_MAPQ:
            status = "UNIQUE_HIGH"
        elif len(rows) == 1:
            status = "UNIQUE_LOW"
        else:
            status = "MULTI"
        out[pair_id] = HostMapping(pair_id, status, len(rows), primary[0], primary[1], primary[2], primary[3], primary[4], primary[5])
    return out


def supported_rescue_pairs(mapping: Mapping[int, HostMapping]) -> set[int]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in mapping.values():
        if row.status == "UNIQUE_HIGH":
            by_chrom[row.chrom].append((int(row.midpoint), row.pair_id))
    kept: set[int] = set()
    for chrom in sorted(by_chrom):
        current: list[tuple[int, int]] = []
        last = None
        for value in sorted(by_chrom[chrom]) + [None]:
            if value is not None and (last is None or value[0] - last <= RESCUE_CLUSTER_GAP):
                current.append(value)
                last = value[0]
                continue
            pairs = {pid for _, pid in current}
            if len(pairs) >= RESCUE_CLUSTER_SUPPORT:
                kept.update(pairs)
            if value is not None:
                current, last = [value], value[0]
    return kept
