"""Truth-blind construction and validation of the safe 21-mer language."""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType

from ..contracts import BaseSpan0, ProbeClass, SafeProbeRecordStatus, TEProbe
from ..errors import ArchitectureInvariantError, AssetCompatibilityError, InputContractError

PROBE_LENGTH = 21


@dataclass(frozen=True, slots=True)
class PreparedProbeLanguage:
    """Unique exact query units plus every source membership."""

    probes: tuple[TEProbe, ...]
    memberships: Mapping[str, tuple[TEProbe, ...]]

    def __post_init__(self) -> None:
        probe_ids = {probe.probe_id for probe in self.probes}
        if len(probe_ids) != len(self.probes) or probe_ids != set(self.memberships):
            raise InputContractError("prepared probe language identity mismatch")
        normalized = {
            probe_id: tuple(rows)
            for probe_id, rows in sorted(self.memberships.items())
        }
        if any(
            not rows or any(row.probe_id != probe_id for row in rows)
            for probe_id, rows in normalized.items()
        ):
            raise InputContractError("prepared probe memberships are malformed")
        object.__setattr__(self, "memberships", MappingProxyType(normalized))

    @property
    def codebook(self) -> Mapping[str, TEProbe]:
        return MappingProxyType({probe.probe_id: probe for probe in self.probes})


def reverse_complement(sequence: str) -> str:
    normalized = sequence.upper()
    if not normalized or any(base not in "ACGT" for base in normalized):
        raise InputContractError("probe sequence must contain only A/C/G/T bases")
    return normalized.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_probe_sequence(sequence: str) -> str:
    normalized = sequence.upper()
    return min(normalized, reverse_complement(normalized))


def load_probe_language_tsv(
    path: str | Path,
    probe_class: ProbeClass,
) -> PreparedProbeLanguage:
    """Load one prepared safe language and collapse shared exact query units."""

    if not isinstance(probe_class, ProbeClass):
        raise InputContractError("probe_class must be a ProbeClass")
    required = {
        "probe_class",
        "sequence",
        "query_sequences_both_strands",
        "family",
        "te_record_id",
        "source_start",
        "source_end",
        "source_strand",
        "independent_block_id",
        "compatible_families",
    }
    try:
        with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not required <= set(reader.fieldnames):
                raise AssetCompatibilityError("prepared codebook schema mismatch")
            rows = list(reader)
    except AssetCompatibilityError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise AssetCompatibilityError(f"cannot read prepared codebook: {exc}") from exc
    if not rows:
        raise AssetCompatibilityError("prepared codebook is empty")

    grouped: dict[str, list[dict[str, object]]] = {}
    blocks_by_record: dict[str, set[str]] = {}
    expected_class = probe_class.value
    for row in rows:
        try:
            sequence = canonical_probe_sequence(row["sequence"])
            if row["sequence"].upper() != sequence or len(sequence) != PROBE_LENGTH:
                raise ValueError("sequence is not a canonical 21-mer")
            if row["probe_class"] != expected_class:
                raise ValueError("probe class mismatch")
            oriented = tuple(
                value for value in row["query_sequences_both_strands"].split(";") if value
            )
            if set(oriented) != {sequence, reverse_complement(sequence)}:
                raise ValueError("oriented query set mismatch")
            source_start = int(row["source_start"])
            source_end = int(row["source_end"])
            if source_end - source_start != PROBE_LENGTH or source_start < 0:
                raise ValueError("source span mismatch")
            source_strand = row["source_strand"]
            if source_strand not in {"+", "-"}:
                raise ValueError("source strand mismatch")
            record_id = row["te_record_id"]
            family = row["family"]
            block_id = row["independent_block_id"]
            if not record_id or not family or not block_id:
                raise ValueError("empty codebook identity")
            compatible = json.loads(row.get("compatible_families", "[]"))
            if compatible != [family]:
                raise ValueError("codebook family compatibility mismatch")
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise AssetCompatibilityError("malformed prepared codebook row") from exc
        grouped.setdefault(sequence, []).append(
            {
                "family": family,
                "record_id": record_id,
                "span": BaseSpan0(source_start, source_end),
                "strand": source_strand,
                "block_id": block_id,
            }
        )
        blocks_by_record.setdefault(record_id, set()).add(block_id)

    probes: list[TEProbe] = []
    memberships: dict[str, tuple[TEProbe, ...]] = {}
    for sequence in sorted(grouped):
        source_rows = grouped[sequence]
        families = {str(row["family"]) for row in source_rows}
        if len(families) != 1:
            raise AssetCompatibilityError(
                "shared exact query unit maps to multiple TE families"
            )
        probe_id = (
            f"{probe_class.value.lower()}_unit_"
            f"{hashlib.sha256(sequence.encode('ascii')).hexdigest()[:20]}"
        )
        unique_memberships = {
            (
                str(row["record_id"]),
                row["span"],
                str(row["strand"]),
                str(row["block_id"]),
            )
            for row in source_rows
        }
        membership_rows: list[TEProbe] = []
        for record_id, span, strand, block_id in sorted(
            unique_memberships,
            key=lambda item: (item[0], item[1].start0, item[1].end0, item[2], item[3]),
        ):
            block_count = len(blocks_by_record[record_id])
            status = (
                SafeProbeRecordStatus.MULTI_INDEPENDENT_BLOCK
                if block_count >= 2
                else SafeProbeRecordStatus.SINGLE_INDEPENDENT_BLOCK
            )
            membership_rows.append(
                TEProbe(
                    probe_id=probe_id,
                    probe_class=probe_class,
                    sequence=sequence,
                    family=next(iter(families)),
                    te_record_id=record_id,
                    source_span0=span,
                    probe_source_strand=strand,
                    independent_block_id=block_id,
                    record_safe_status=status,
                )
            )
        memberships[probe_id] = tuple(membership_rows)
        probes.append(membership_rows[0])
    return PreparedProbeLanguage(tuple(probes), memberships)


def load_prepared_probe_language(
    prepared_te,
    language: str,
) -> PreparedProbeLanguage:
    """Load a hash-bound language from ``PreparedTEAssets``."""

    from ..assets.te import verify_te_assets

    verify_te_assets(prepared_te)
    try:
        entry = prepared_te.codebooks[language]
        relative_path = entry["path"]
    except (KeyError, TypeError) as exc:
        raise AssetCompatibilityError(
            f"prepared TE assets lack codebook language {language!r}"
        ) from exc
    if language == "terminal":
        probe_class = ProbeClass.TERMINAL_21MER
    elif language in {"internal_full", "internal_only"}:
        probe_class = ProbeClass.INTERNAL_21MER
    else:
        raise InputContractError("unknown prepared probe language")
    return load_probe_language_tsv(prepared_te.root / relative_path, probe_class)


@dataclass(frozen=True, slots=True)
class SafeProbeCandidate:
    """One exhaustively assessed TE 21-mer candidate.

    ``compatible_families`` is the complete family set in which the canonical
    sequence occurs. ``reference_occurrences`` is its exact host-reference
    count. Both are prepared without sample truth data.
    """

    candidate_id: str
    sequence: str
    family: str
    te_record_id: str
    source_span0: BaseSpan0
    probe_source_strand: str
    independent_block_id: str
    compatible_families: tuple[str, ...]
    reference_occurrences: int

    def __post_init__(self) -> None:
        for value, field in (
            (self.candidate_id, "candidate_id"),
            (self.family, "family"),
            (self.te_record_id, "te_record_id"),
            (self.independent_block_id, "independent_block_id"),
        ):
            if not isinstance(value, str) or not value:
                raise InputContractError(f"{field} must be non-empty text")
        if self.probe_source_strand not in {"+", "-"}:
            raise InputContractError("probe_source_strand must be '+' or '-'")
        if len(self.sequence) != PROBE_LENGTH or self.source_span0.length != PROBE_LENGTH:
            raise InputContractError("safe probe candidates must be exact 21-mers")
        if (
            isinstance(self.reference_occurrences, bool)
            or not isinstance(self.reference_occurrences, int)
            or self.reference_occurrences < 0
        ):
            raise InputContractError("reference_occurrences must be non-negative")


def classify_te_record_probe_status(
    probes: Iterable[TEProbe],
) -> SafeProbeRecordStatus:
    rows = tuple(probes)
    if not rows:
        return SafeProbeRecordStatus.UNSUPPORTED_TERMINAL_SAFE
    record_ids = {probe.te_record_id for probe in rows}
    if len(record_ids) != 1:
        raise InputContractError("record status classification mixes TE records")
    blocks = {probe.independent_block_id for probe in rows}
    if len(blocks) >= 2:
        return SafeProbeRecordStatus.MULTI_INDEPENDENT_BLOCK
    return SafeProbeRecordStatus.SINGLE_INDEPENDENT_BLOCK


def _overlaps(left: BaseSpan0, right: BaseSpan0) -> bool:
    return left.start0 < right.end0 and right.start0 < left.end0


def build_safe_probes(
    candidates: Iterable[SafeProbeCandidate],
    probe_class: ProbeClass,
) -> tuple[TEProbe, ...]:
    """Select deterministic family-unique, reference-absent independent probes."""

    if not isinstance(probe_class, ProbeClass):
        raise InputContractError("probe_class must be a ProbeClass")
    rows = tuple(candidates)
    declared_compatibility: dict[str, frozenset[str]] = {}
    for candidate in rows:
        canonical = canonical_probe_sequence(candidate.sequence)
        declared = frozenset(candidate.compatible_families)
        previous = declared_compatibility.setdefault(canonical, declared)
        if previous != declared:
            raise InputContractError(
                "canonical probe candidates disagree on exhaustive family compatibility"
            )
    safe = sorted(
        (
            candidate
            for candidate in rows
            if set(candidate.compatible_families) == {candidate.family}
            and candidate.reference_occurrences == 0
        ),
        key=lambda candidate: (
            candidate.te_record_id,
            candidate.source_span0.end0,
            candidate.source_span0.start0,
            candidate.probe_source_strand,
            candidate.candidate_id,
        ),
    )
    selected: list[TEProbe] = []
    canonical_seen: set[str] = set()
    spans_by_record: dict[str, list[BaseSpan0]] = {}
    for candidate in safe:
        canonical = canonical_probe_sequence(candidate.sequence)
        if canonical in canonical_seen:
            continue
        record_spans = spans_by_record.setdefault(candidate.te_record_id, [])
        if any(_overlaps(candidate.source_span0, span) for span in record_spans):
            continue
        canonical_seen.add(canonical)
        record_spans.append(candidate.source_span0)
        selected.append(
            TEProbe(
                probe_id=candidate.candidate_id,
                probe_class=probe_class,
                sequence=candidate.sequence.upper(),
                family=candidate.family,
                te_record_id=candidate.te_record_id,
                source_span0=candidate.source_span0,
                probe_source_strand=candidate.probe_source_strand,
                independent_block_id=candidate.independent_block_id,
                record_safe_status=SafeProbeRecordStatus.UNSUPPORTED_TERMINAL_SAFE,
            )
        )

    by_record: dict[str, list[TEProbe]] = {}
    for probe in selected:
        by_record.setdefault(probe.te_record_id, []).append(probe)
    finalized = tuple(
        replace(
            probe,
            record_safe_status=classify_te_record_probe_status(
                by_record[probe.te_record_id]
            ),
        )
        for probe in sorted(selected, key=lambda item: item.probe_id)
    )
    validate_probe_independence(finalized)
    return finalized


def validate_probe_independence(probes: Iterable[TEProbe]) -> None:
    """Raise on duplicate sequences, duplicate IDs, or overlapping source units."""

    rows = tuple(probes)
    probe_ids: set[str] = set()
    canonical_sequences: set[str] = set()
    spans_by_record: dict[str, list[BaseSpan0]] = {}
    for probe in rows:
        if probe.probe_id in probe_ids:
            raise ArchitectureInvariantError(f"duplicate probe_id: {probe.probe_id}")
        probe_ids.add(probe.probe_id)
        if len(probe.sequence) != PROBE_LENGTH or probe.source_span0.length != PROBE_LENGTH:
            raise ArchitectureInvariantError("safe codebook contains a non-21-mer")
        canonical = canonical_probe_sequence(probe.sequence)
        if canonical in canonical_sequences:
            raise ArchitectureInvariantError("safe codebook repeats a canonical probe")
        canonical_sequences.add(canonical)
        record_spans = spans_by_record.setdefault(probe.te_record_id, [])
        if any(_overlaps(probe.source_span0, span) for span in record_spans):
            raise ArchitectureInvariantError("safe codebook source intervals overlap")
        record_spans.append(probe.source_span0)
