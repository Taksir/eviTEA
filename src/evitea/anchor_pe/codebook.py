"""Final safe terminal codebook, its msBWT projection, and pair classification."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

from ..constants import N_FOLDS
from ..errors import AssetCompatibilityError, InputContractError

__all__ = [
    "CodeUnit",
    "SafeCodebook",
    "PairHits",
    "load_safe_codebook",
    "project_codebook",
    "classify_pairs",
    "one_mate_universe",
    "revcomp",
    "fold_of",
]

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def fold_of(pair_id: int) -> int:
    return int(pair_id) % N_FOLDS


def _as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def _as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "t"}


@dataclass(frozen=True, slots=True)
class CodeUnit:
    unit_id: int
    family: str
    te_id: str
    probe_rank: int
    terminal_side: str
    k: int
    canonical: str
    sequences: tuple[str, ...]
    source_start: int
    source_end: int
    positionally_interpretable: bool


@dataclass(frozen=True, slots=True)
class SafeCodebook:
    units: tuple[CodeUnit, ...]

    def by_family(self) -> dict[str, list[CodeUnit]]:
        out: dict[str, list[CodeUnit]] = defaultdict(list)
        for unit in self.units:
            out[unit.family].append(unit)
        return out

    def record_families(self) -> dict[str, str]:
        out: dict[str, str] = {}
        for unit in self.units:
            if unit.te_id and unit.family:
                out[unit.te_id] = unit.family
        return out

    def probe_index(self) -> dict[tuple[str, str], list[CodeUnit]]:
        """Units keyed by ``(te_id, oriented probe)`` over both strands."""

        index: dict[tuple[str, str], list[CodeUnit]] = defaultdict(list)
        for unit in self.units:
            if not unit.te_id or not unit.canonical:
                continue
            index[(unit.te_id, unit.canonical)].append(unit)
            index[(unit.te_id, revcomp(unit.canonical).upper())].append(unit)
        return index


def load_safe_codebook(path: str | Path) -> SafeCodebook:
    """Load the final safe terminal codebook (one unit per row)."""

    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise AssetCompatibilityError("safe terminal codebook is empty")
    units = []
    for index, row in enumerate(rows):
        seqs = tuple(dict.fromkeys(
            x.strip().upper() for x in str(row.get("query_sequences_both_strands", "")).split(";") if x.strip()
        ))
        if not seqs:
            raise AssetCompatibilityError(f"codebook row {index + 1} has no query sequences")
        canonical = str(row.get("canonical_strand_probe", "")).strip().upper()
        units.append(
            CodeUnit(
                unit_id=index,
                family=row["family"],
                te_id=str(row.get("te_id", "")).strip(),
                probe_rank=int(row["probe_rank_within_te_record"]),
                terminal_side=str(row.get("terminal_side", "")),
                k=_as_int(row.get("k"), 0) or len(canonical),
                canonical=canonical,
                sequences=seqs,
                source_start=_as_int(row.get("source_start"), -1),
                source_end=_as_int(row.get("source_end"), -1),
                positionally_interpretable=_as_bool(row.get("positionally_interpretable", "0")),
            )
        )
    return SafeCodebook(tuple(units))


# pair_id maps to {1: unit ids on mate 1, 2: unit ids on mate 2}
PairHits = dict[int, dict[int, set[int]]]


def project_codebook(
    sample,
    codebook: SafeCodebook,
    *,
    progress: Callable[[str], None] | None = None,
) -> tuple[PairHits, list[dict[str, object]]]:
    """Exact projection of every oriented probe; each read counted once per unit.

    ``sample`` is an ``incumbent.MsbwtSample`` (interval, dollar id, provenance).
    """

    pair_hits: PairHits = defaultdict(lambda: {1: set(), 2: set()})
    observation: list[dict[str, object]] = []
    for n, unit in enumerate(codebook.units, 1):
        raw_rows = 0
        seen_dollars: set[int] = set()
        pairs: set[int] = set()
        mates: Counter[int] = Counter()
        for query in unit.sequences:
            lo, hi = sample.interval(query, cache=False)
            raw_rows += hi - lo
            if hi == lo:
                continue
            dollars = [int(d) for d in sample.dollar_ids(query, cache=False).tolist()]
            fresh = []
            for dollar in dollars:
                if dollar not in seen_dollars:
                    seen_dollars.add(dollar)
                    fresh.append(dollar)
            if not fresh:
                continue
            files, pair_ids = sample.file_and_pair(fresh)
            for file_id, pair_id in zip(files.tolist(), pair_ids.tolist()):
                if not 0 <= int(pair_id) < sample.pair_count:
                    raise AssetCompatibilityError(f"pair id outside the sample: {pair_id}")
                mate = int(file_id) + 1
                pair_hits[int(pair_id)][mate].add(unit.unit_id)
                pairs.add(int(pair_id))
                mates[mate] += 1
        observation.append({
            "unit_id": unit.unit_id, "family": unit.family, "te_id": unit.te_id,
            "probe_rank": unit.probe_rank, "terminal_side": unit.terminal_side,
            "raw_msbwt_rows": raw_rows, "distinct_dollar_ids": len(seen_dollars),
            "distinct_pair_ids": len(pairs), "mate1_dollar_ids": mates[1],
            "mate2_dollar_ids": mates[2], "observed": int(bool(seen_dollars)),
        })
        if progress and (n % 50 == 0 or n == len(codebook.units)):
            progress(f"codebook projection {n}/{len(codebook.units)}; pairs={len(pair_hits):,}")
    return {pid: {1: set(v[1]), 2: set(v[2])} for pid, v in pair_hits.items()}, observation


def classify_pairs(pair_hits: Mapping[int, Mapping[int, set[int]]], codebook: SafeCodebook) -> tuple[list[dict], dict[int, dict]]:
    """Signature rows for every hit pair and the one-mate family-resolved set."""

    units = {u.unit_id: u for u in codebook.units}
    rows: list[dict] = []
    eligible: dict[int, dict] = {}
    for pid in sorted(pair_hits):
        m1, m2 = pair_hits[pid][1], pair_hits[pid][2]
        families = {units[i].family for i in m1 | m2}
        family = next(iter(families)) if len(families) == 1 else ""
        te_mate = 0
        primary = False
        if len(families) != 1:
            state = "FAMILY_CONFLICT"
        elif m1 and m2:
            state = "BOTH_MATES_TE_POSITIVE"
        else:
            te_mate = 1 if m1 else 2
            selected = m1 if m1 else m2
            by_te: dict[str, set[int]] = defaultdict(set)
            for uid in selected:
                by_te[units[uid].te_id].add(uid)
            primary = any(len(v) >= 2 for v in by_te.values())
            state = "PRIMARY_MULTI_BLOCK_COHERENT" if primary else "NO_MULTI_RECORD_GE2_BLOCKS"
            eligible[pid] = {
                "pair_id": pid, "family": family, "te_mate": te_mate,
                "host_mate": 2 if te_mate == 1 else 1, "unit_ids": sorted(selected),
                "primary": primary, "fold": fold_of(pid),
            }
        rows.append({
            "pair_id": pid, "fold": fold_of(pid), "mate1_units": ";".join(map(str, sorted(m1))),
            "mate2_units": ";".join(map(str, sorted(m2))), "families": ";".join(sorted(families)),
            "family": family, "te_mate": te_mate, "host_mate": (2 if te_mate == 1 else 1) if te_mate else 0,
            "primary_multiblock": int(primary), "eligibility_state": state,
        })
    return rows, eligible


def one_mate_universe(pair_hits: Mapping[int, Mapping[int, set[int]]], eligible: Mapping[int, dict]) -> dict[int, dict]:
    """Broader one-mate probe-positive universe used by the relaxed novelty mask.

    Includes one-mate family-conflict pairs (not only the family-resolved set).
    Insertion order follows ``pair_hits``.
    """

    out: dict[int, dict] = {}
    for pid, hits in pair_hits.items():
        nonempty = [m for m in (1, 2) if hits[m]]
        if len(nonempty) != 1:
            continue
        te_mate = nonempty[0]
        out[pid] = {
            "pair_id": pid, "te_mate": te_mate, "host_mate": 2 if te_mate == 1 else 1,
            "unit_ids": sorted(hits[te_mate]),
            "old_primary": bool(pid in eligible and eligible[pid]["primary"]),
        }
    if any(pid < 0 for pid in out):
        raise InputContractError("negative pair id")
    return out
