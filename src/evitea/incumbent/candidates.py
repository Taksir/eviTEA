"""Build the candidate read pairs and loci from panel acquisition.

The candidate set is the working pairs plus the pairs added by the wide-anchor
step. A pair is localized by its host-mate record: the working-pair record
when one exists, otherwise the wide-anchor record, using the primary alignment
of *any* mapped status (unique, multiple or low-confidence). Its family set
comes from the resident-compatibility pair families, falling back to
wide-anchor then sparse-panel evidence. A pair with several families joins
every family bucket. Loci are single-link clusters of float host midpoints per
``(chrom, family)`` (gap 500 bp, >= 2 distinct pairs).

The ``source`` field is "G1F" for a working pair and "G1HR_RESCUE" for a
wide-anchor pair.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

CANDIDATE_CLUSTER_GAP = 500
CANDIDATE_MIN_PAIRS = 2

__all__ = [
    "CANDIDATE_CLUSTER_GAP",
    "CANDIDATE_MIN_PAIRS",
    "HostMateRecord",
    "PairPanelEvidence",
    "CandidatePair",
    "CandidateLocus",
    "split_families",
    "reconstruct_candidate_pairs",
    "cluster_candidate_loci",
]


def split_families(text: str | None) -> tuple[str, ...]:
    if not text:
        return ()
    values = []
    for token in str(text).replace(";", ",").split(","):
        token = token.strip()
        if token and token.lower() not in {"nan", "none"}:
            values.append(token)
    return tuple(sorted(set(values)))


@dataclass(frozen=True, slots=True)
class HostMateRecord:
    """One host-mate mapping row; empty ``chrom``/``midpoint`` means unmapped."""

    pair_id: int
    status: str
    chrom: str
    midpoint: str
    strand: str


@dataclass(frozen=True, slots=True)
class PairPanelEvidence:
    pair_id: int
    families: str
    mate1_hit: bool
    mate2_hit: bool


@dataclass(frozen=True, slots=True)
class CandidatePair:
    pair_id: int
    chrom: str
    midpoint: float
    host_status: str
    host_strand: str
    families: tuple[str, ...]
    mate1_hit: bool
    mate2_hit: bool
    source: str


@dataclass(frozen=True, slots=True)
class CandidateLocus:
    locus_id: str
    chrom: str
    family: str
    start: float
    end: float
    pair_ids: tuple[int, ...]

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


def reconstruct_candidate_pairs(
    *,
    final_pair_ids: Iterable[int],
    working_pair_ids: Iterable[int],
    primary_host: Mapping[int, HostMateRecord],
    rescue_host: Mapping[int, HostMateRecord],
    resident_families: Mapping[int, str],
    rescue_evidence: Mapping[int, PairPanelEvidence],
    panel_evidence: Mapping[int, PairPanelEvidence],
) -> tuple[tuple[CandidatePair, ...], dict[str, int]]:
    base = set(int(value) for value in working_pair_ids)
    final = sorted(set(int(value) for value in final_pair_ids))
    candidates: list[CandidatePair] = []
    missing_host = missing_family = missing_hitmask = 0
    for pair_id in final:
        host = primary_host.get(pair_id) or rescue_host.get(pair_id)
        if host is None or not host.chrom or not host.midpoint:
            missing_host += 1
            continue
        families = split_families(resident_families.get(pair_id, ""))
        if not families:
            evidence = rescue_evidence.get(pair_id) or panel_evidence.get(pair_id)
            families = split_families(evidence.families if evidence else "")
        if not families:
            missing_family += 1
            continue
        hit = rescue_evidence.get(pair_id) or panel_evidence.get(pair_id)
        if hit is None:
            missing_hitmask += 1
            mate1 = mate2 = False
        else:
            mate1, mate2 = hit.mate1_hit, hit.mate2_hit
        candidates.append(
            CandidatePair(
                pair_id=pair_id,
                chrom=host.chrom,
                midpoint=float(host.midpoint),
                host_status=host.status,
                host_strand=host.strand,
                families=families,
                mate1_hit=mate1,
                mate2_hit=mate2,
                source="G1F" if pair_id in base else "G1HR_RESCUE",
            )
        )
    meta = {
        "working_candidate_pairs": len(base),
        "final_candidate_pairs": len(final),
        "localized_candidate_pairs": len(candidates),
        "unlocalized_or_missing_host": missing_host,
        "missing_family": missing_family,
        "missing_hitmask": missing_hitmask,
    }
    return tuple(candidates), meta


def cluster_candidate_loci(
    candidates: Sequence[CandidatePair],
    gap: int = CANDIDATE_CLUSTER_GAP,
    min_pairs: int = CANDIDATE_MIN_PAIRS,
) -> tuple[CandidateLocus, ...]:
    buckets: dict[tuple[str, str], list[tuple[float, int]]] = defaultdict(list)
    for candidate in candidates:
        for family in candidate.families:
            buckets[(candidate.chrom, family)].append((candidate.midpoint, candidate.pair_id))
    loci: list[CandidateLocus] = []

    def emit(chrom: str, family: str, rows: list[tuple[float, int]]) -> None:
        pair_ids = tuple(sorted({pair_id for _, pair_id in rows}))
        if len(pair_ids) >= min_pairs:
            loci.append(
                CandidateLocus(f"L{len(loci) + 1:07d}", chrom, family, rows[0][0], rows[-1][0], pair_ids)
            )

    for (chrom, family), values in sorted(buckets.items()):
        values.sort()
        current: list[tuple[float, int]] = []
        last: float | None = None
        for value in values:
            if last is None or value[0] - last <= gap:
                current.append(value)
            else:
                emit(chrom, family, current)
                current = [value]
            last = value[0]
        if current:
            emit(chrom, family, current)
    return tuple(loci)
