"""TE-library evidence primitives of the library layer.

Every library layer step uses these primitives:

* whole-segment infix alignment of a TE-facing sequence anywhere in each TE
  consensus record, either strand, at most three edits;
* collapse to the best distance per biological family (raw family-map labels,
  no alias normalization at this level);
* a unique best family is required; equal-best families are rejected;
* side clustering of qualified breakpoints (single link, radius 8 bp, modal
  coordinate with ties to the larger coordinate).

It also builds the candidate-family 21-mer panels.  A record joins a family
panel when its first header token, its mapped family, or a substring match of
the family label selects it.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..errors import InputContractError
from ..family_table import family_pairs, is_two_column_table
from .alignment import hw_distance
from .sequence import fasta_records_full_header, revcomp, shannon_entropy_acgt

MAX_TE_EDITS = 3
TE_INDEX_KMER = 8
SIDE_CLUSTER_RADIUS = 8
_FBTE = re.compile(r"FBte\d+", re.IGNORECASE)
FROZEN_ALIASES = {"rooA": "roo", "Rt1a": "Rt1b"}

__all__ = [
    "MAX_TE_EDITS",
    "SIDE_CLUSTER_RADIUS",
    "normalize_family",
    "load_record_family_map",
    "load_heuristic_family_map",
    "TERecord",
    "load_te_records",
    "LocalHit",
    "LocalTEIndex",
    "family_best_classify",
    "QualifiedObservation",
    "SideCluster",
    "cluster_observations",
    "mode_hi",
    "load_family_kmer_panels",
]


def normalize_family(label: str) -> str:
    """Apply only the two built-in aliases (``rooA`` to ``roo``, ``Rt1a`` to ``Rt1b``)."""

    return FROZEN_ALIASES.get(str(label), str(label))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_record_family_map(path: str | Path) -> dict[str, str]:
    """Read ``FBte_ID`` to ``FlyBase_Family`` (raw labels) with column fallbacks, or a two-column ``record<TAB>family`` table."""

    if is_two_column_table(path):
        return dict(family_pairs(path))
    out: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            record = (row.get("FBte_ID") or row.get("fbte_id") or row.get("FBte") or "").strip()
            family = (
                row.get("FlyBase_Family") or row.get("family") or row.get("Family") or ""
            ).strip()
            if record and family:
                out[record] = family
    return out


def load_heuristic_family_map(path: str | Path) -> dict[str, str]:
    """Read the family map: first ID-like column and first ``family`` column, or a two-column ``record<TAB>family`` table."""

    if is_two_column_table(path):
        return dict(family_pairs(path))
    rows = _read_tsv(Path(path))
    if not rows:
        return {}
    fields = list(rows[0])
    id_fields = [
        name for name in fields if "fbte" in name.lower() or name.lower() in {"te_id", "id", "element", "name"}
    ]
    family_fields = [name for name in fields if "family" in name.lower()]
    id_field = id_fields[0] if id_fields else fields[0]
    family_field = family_fields[0] if family_fields else fields[-1]
    out: dict[str, str] = {}
    for row in rows:
        key = (row.get(id_field) or "").strip()
        value = (row.get(family_field) or "").strip()
        if key and value:
            out[key] = value
    return out


def _header_family(header: str, record_id: str, family_map: Mapping[str, str]) -> str:
    if record_id in family_map:
        return family_map[record_id]
    match = re.search(r"(?:^|[;\s])name=([^;\s]+)", header)
    if match:
        name = match.group(1)
        return name.split("\\")[-1] if "\\" in name else name
    tokens = header.split()
    if len(tokens) > 1 and "=" not in tokens[1]:
        return tokens[1].split("\\")[-1]
    return record_id


@dataclass(frozen=True, slots=True)
class TERecord:
    record_id: str
    family: str
    sequence: str


def load_te_records(path: str | Path, family_map: Mapping[str, str]) -> tuple[TERecord, ...]:
    """Parse TE consensus records, keyed by the ``FBte`` ID or first header token."""

    records: list[TERecord] = []
    header: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if header is None:
            return
        sequence = "".join(chunks).upper()
        match = _FBTE.search(header)
        record_id = match.group(0) if match else header.split()[0]
        if sequence:
            records.append(TERecord(record_id, _header_family(header, record_id, family_map), sequence))

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:]
                chunks = []
            else:
                chunks.append(re.sub("[^ACGTNacgtn]", "", line))
    finish()
    if not records:
        raise InputContractError(f"no TE FASTA records parsed from {path}")
    return tuple(records)


@dataclass(frozen=True, slots=True)
class LocalHit:
    family: str
    record_id: str
    edit_distance: int
    orientation: str


class LocalTEIndex:
    """Exact infix TE-library matcher with an 8-mer candidate prefilter.

    With at most three edits every qualifying query of length >= 35 shares an
    exact 8-mer with its target, so the prefilter never changes a result.
    """

    def __init__(self, records: Sequence[TERecord]):
        self.records: list[tuple[str, str, str, str]] = []
        self._kmers: dict[str, set[int]] = defaultdict(set)
        for record_index, record in enumerate(records):
            forward = record.sequence.upper().replace("N", "A")
            reverse = revcomp(forward).replace("N", "A")
            self.records.append((record.record_id, record.family, forward, reverse))
            seen: set[str] = set()
            for text in (forward, reverse):
                for offset in range(0, len(text) - TE_INDEX_KMER + 1):
                    kmer = text[offset : offset + TE_INDEX_KMER]
                    if len(kmer) == TE_INDEX_KMER and set(kmer) <= set("ACGT"):
                        seen.add(kmer)
            for kmer in seen:
                self._kmers[kmer].add(record_index)

    def candidate_ids(self, query: str) -> set[int]:
        normalized = query.upper().replace("N", "A")
        out: set[int] = set()
        for offset in range(0, len(normalized) - TE_INDEX_KMER + 1):
            out.update(self._kmers.get(normalized[offset : offset + TE_INDEX_KMER], ()))
        return out

    def record_hits(self, query: str) -> list[tuple[int, str, str, int, str]]:
        """``(record_index, record_id, family, edit_distance, orientation)`` rows."""

        normalized = query.upper().replace("N", "A")
        hits: list[tuple[int, str, str, int, str]] = []
        for record_index in sorted(self.candidate_ids(normalized)):
            record_id, family, forward, reverse = self.records[record_index]
            best: int | None = None
            orientation = ""
            for label, target in (("fwd", forward), ("rc", reverse)):
                distance = hw_distance(normalized, target, MAX_TE_EDITS)
                if distance is not None and (best is None or distance < best):
                    best = distance
                    orientation = label
                    if distance == 0:
                        break
            if best is not None:
                hits.append((record_index, record_id, family, best, orientation))
        return hits


def family_best_classify(
    hits: Iterable[tuple[int, str, str, int, str]],
) -> tuple[LocalHit | None, dict[str, object]]:
    """Return unique best biological family among TE-record hits (raw labels)."""

    rows = list(hits)
    if not rows:
        return None, {"reason": "NO_TE_MATCH_LE3", "best_distance": "", "best_families": "", "best_fbte": ""}
    best_by_family: dict[str, dict[str, object]] = {}
    for record_index, record_id, family, distance, orientation in rows:
        old = best_by_family.get(family)
        if old is None or distance < old["distance"] or (
            distance == old["distance"] and record_index < old["record_index"]
        ):
            best_by_family[family] = {
                "distance": distance,
                "record_id": record_id,
                "orientation": orientation,
                "record_index": record_index,
            }
    best = min(int(value["distance"]) for value in best_by_family.values())
    winners = sorted(family for family, value in best_by_family.items() if value["distance"] == best)
    if len(winners) != 1:
        return None, {
            "reason": "EQUAL_BEST_FAMILY_TIE",
            "best_distance": best,
            "best_families": ",".join(winners),
            "best_fbte": ",".join(sorted(str(best_by_family[f]["record_id"]) for f in winners)),
        }
    family = winners[0]
    value = best_by_family[family]
    return (
        LocalHit(family, str(value["record_id"]), int(value["distance"]), str(value["orientation"])),
        {"reason": "QUALIFIED", "best_distance": best, "best_families": family, "best_fbte": value["record_id"]},
    )


@dataclass(frozen=True, slots=True)
class QualifiedObservation:
    read_id: str
    chrom: str
    side: str
    breakpoint: int
    clip_len: int
    family: str
    fbte_id: str
    edit_distance: int
    orientation: str


@dataclass(slots=True)
class SideCluster:
    chrom: str
    family: str
    side: str
    coord: int
    reads: set[str]
    positions: list[int]
    obs: list[QualifiedObservation] = field(default_factory=list)

    @property
    def support(self) -> int:
        return len(self.reads)


def mode_hi(values: Sequence[int]) -> int:
    counts = Counter(values)
    top = max(counts.values())
    return max(value for value, count in counts.items() if count == top)


def _finish_cluster(chrom: str, family: str, side: str, rows: list[QualifiedObservation]) -> SideCluster:
    positions = [row.breakpoint for row in rows]
    return SideCluster(chrom, family, side, mode_hi(positions), {row.read_id for row in rows}, positions, list(rows))


def cluster_observations(observations: Iterable[QualifiedObservation]) -> list[SideCluster]:
    buckets: dict[tuple[str, str, str], list[QualifiedObservation]] = defaultdict(list)
    seen: set[tuple[str, str, str, str, int]] = set()
    for row in observations:
        key = (row.read_id, row.chrom, row.family, row.side, row.breakpoint)
        if key in seen:
            continue
        seen.add(key)
        buckets[(row.chrom, row.family, row.side)].append(row)
    out: list[SideCluster] = []
    for (chrom, family, side), rows in sorted(buckets.items()):
        rows.sort(key=lambda item: (item.breakpoint, item.read_id))
        current: list[QualifiedObservation] = []
        for row in rows:
            if not current:
                current = [row]
                continue
            if row.breakpoint - current[-1].breakpoint <= SIDE_CLUSTER_RADIUS:
                current.append(row)
            else:
                out.append(_finish_cluster(chrom, family, side, current))
                current = [row]
        if current:
            out.append(_finish_cluster(chrom, family, side, current))
    return out


def load_family_kmer_panels(
    te_fasta: str | Path,
    family_map_path: str | Path,
    wanted_families: Iterable[str],
    k: int,
    min_entropy: float,
) -> tuple[dict[str, set[str]], dict[str, object]]:
    """Build candidate-family k-mer panels by the module's record-selection rule."""

    family_map = load_heuristic_family_map(family_map_path)
    wanted = set(wanted_families)
    panels: dict[str, set[str]] = {family: set() for family in wanted}
    records = 0
    mapped = 0
    for header, sequence in fasta_records_full_header(te_fasta):
        records += 1
        token = header.split()[0] if header.split() else header
        mapped_family = family_map.get(header, family_map.get(token))
        targets: set[str] = set()
        if token in wanted:
            targets.add(token)
        if mapped_family in wanted:
            targets.add(mapped_family)
        for family in wanted:
            if family and family in header:
                targets.add(family)
        if not targets:
            continue
        mapped += 1
        forward = sequence.upper()
        for text in (forward, revcomp(forward)):
            for offset in range(0, len(text) - k + 1):
                kmer = text[offset : offset + k]
                if "N" in kmer or shannon_entropy_acgt(kmer) < min_entropy:
                    continue
                for family in targets:
                    panels[family].add(kmer)
    missing = sorted(family for family, kmers in panels.items() if not kmers)
    return panels, {"te_records": records, "mapped_records": mapped, "missing_family_kmers": missing}
