"""Build truth-blind safe 21-mer codebooks from the hash-checked anchor catalogue archive."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..contracts import FastaRecord, ReferenceInput, RepeatRecord, RepeatContext, TELibraryInput
from ..errors import AssetCompatibilityError, InputContractError
from ..families import canonical_family

G1A_ARCHIVE_SHA256 = "f08f1c1936566b0ca469a1c583ac9a96164b2fbd94e6d82f4235b88c49c6b695"
FAMILY_CROSSWALK_SHA256 = "6616b6305e6eac2f97336620376859151083a4b9ed30ba3bd64d7eead7c70714"
CODEBOOK_SCHEMA = "evitea.frozen-g1a-safe-21mer.v1"
PROBE_LENGTH = 21

_CATALOGUE = "G1A_ANCHOR_CATALOGUE.tsv"
_PROVENANCE = "G1A_ANCHOR_PROVENANCE.tsv"
_RUN_MANIFEST = "G1A_RUN_MANIFEST.json"
_CATALOGUE_FIELDS = {
    "query_id", "k", "query_seq", "canonical_strand_key",
    "has_terminal_provenance", "has_internal_provenance", "low_complexity",
}
_PROVENANCE_FIELDS = {
    "query_id", "te_id", "family", "k", "start", "region",
    "source_orientation", "source_seq", "query_seq",
}
_CROSSWALK_FIELDS = {"FBte_ID", "FlyBase_Family"}
_ENCODING = [-1] * 256
for _base, _value in ((ord("A"), 0), (ord("C"), 1), (ord("G"), 2), (ord("T"), 3)):
    _ENCODING[_base] = _value
    _ENCODING[ord(chr(_base).lower())] = _value
_FBTE = re.compile(r"FBte\d+", re.IGNORECASE)
_FBTE_UNKNOWN = re.compile(r"FBte000XXXX", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class FrozenCodebookBuild:
    terminal_tsv: bytes
    internal_tsv: bytes
    internal_only_tsv: bytes
    terminal_rows: tuple[Mapping[str, str], ...]
    internal_rows: tuple[Mapping[str, str], ...]
    internal_only_rows: tuple[Mapping[str, str], ...]
    candidate_counts: Mapping[str, int]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise InputContractError(f"cannot read file {path}: {exc}") from exc
    return digest.hexdigest()


def _canonical(sequence: str) -> str:
    normalized = sequence.upper()
    if len(normalized) != PROBE_LENGTH or any(base not in "ACGT" for base in normalized):
        raise AssetCompatibilityError("anchor candidate is not a 21-mer of A/C/G/T bases only")
    reverse = normalized.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(normalized, reverse)


def _archive_rows(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if sha256_file(path) != G1A_ARCHIVE_SHA256:
        raise AssetCompatibilityError("anchor archive SHA-256 does not match the expected value")
    try:
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            if len(names) != len(set(names)):
                raise AssetCompatibilityError("anchor archive has duplicate members")
            if _RUN_MANIFEST not in names:
                raise AssetCompatibilityError("anchor archive is missing its run manifest")
            manifest = json.loads(archive.read(_RUN_MANIFEST).decode("utf-8-sig"))
            if not isinstance(manifest, dict) or set(manifest) != {
                "data_root", "files", "generated_at", "project_root", "python"
            }:
                raise AssetCompatibilityError("anchor archive run-manifest schema mismatch")
            entries = manifest["files"]
            if not isinstance(entries, list):
                raise AssetCompatibilityError("anchor archive member manifest is malformed")
            expected: dict[str, tuple[int, str]] = {}
            for entry in entries:
                if not isinstance(entry, dict) or set(entry) != {"bytes", "file", "sha256"}:
                    raise AssetCompatibilityError("anchor archive member entry schema mismatch")
                name = entry["file"]
                size = entry["bytes"]
                digest = entry["sha256"]
                if (
                    not isinstance(name, str)
                    or Path(name).name != name
                    or not isinstance(size, int)
                    or isinstance(size, bool)
                    or size < 0
                    or not isinstance(digest, str)
                    or not re.fullmatch(r"[0-9a-f]{64}", digest)
                    or name in expected
                ):
                    raise AssetCompatibilityError("anchor archive member manifest has invalid values")
                expected[name] = (size, digest)
            if not {_CATALOGUE, _PROVENANCE} <= set(expected):
                raise AssetCompatibilityError("anchor archive lacks catalogue/provenance")
            if set(names) != set(expected) | {_RUN_MANIFEST}:
                raise AssetCompatibilityError("anchor archive members do not match its manifest")
            for name, (size, digest) in expected.items():
                content = archive.read(name)
                if len(content) != size or hashlib.sha256(content).hexdigest() != digest:
                    raise AssetCompatibilityError(f"anchor archive member hash mismatch: {name}")

            def read_tsv(name: str, required: set[str]) -> list[dict[str, str]]:
                with archive.open(name, "r") as raw:
                    reader = csv.DictReader(
                        io.TextIOWrapper(raw, encoding="utf-8-sig", newline=""),
                        delimiter="\t",
                    )
                    if (
                        reader.fieldnames is None
                        or len(reader.fieldnames) != len(set(reader.fieldnames))
                        or not required <= set(reader.fieldnames)
                    ):
                        raise AssetCompatibilityError(f"anchor archive schema mismatch: {name}")
                    rows = list(reader)
                if not rows:
                    raise AssetCompatibilityError(f"anchor archive table is empty: {name}")
                return rows

            return (
                read_tsv(_CATALOGUE, _CATALOGUE_FIELDS),
                read_tsv(_PROVENANCE, _PROVENANCE_FIELDS),
            )
    except AssetCompatibilityError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError, zipfile.BadZipFile, KeyError) as exc:
        raise AssetCompatibilityError(f"cannot validate anchor archive: {exc}") from exc


def _crosswalk(path: Path) -> dict[str, str]:
    if sha256_file(path) != FAMILY_CROSSWALK_SHA256:
        raise AssetCompatibilityError("FBte family crosswalk SHA-256 does not match the expected value")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if (
                reader.fieldnames is None
                or len(reader.fieldnames) != len(set(reader.fieldnames))
                or not _CROSSWALK_FIELDS <= set(reader.fieldnames)
            ):
                raise AssetCompatibilityError("FBte family crosswalk schema mismatch")
            rows = list(reader)
    except AssetCompatibilityError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise AssetCompatibilityError(f"cannot read FBte family crosswalk: {exc}") from exc
    if not rows:
        raise AssetCompatibilityError("FBte family crosswalk is empty")
    mapping: dict[str, str] = {}
    for row in rows:
        te_id = row.get("FBte_ID", "").strip()
        family = row.get("FlyBase_Family", "").strip()
        if not te_id or not family:
            raise AssetCompatibilityError("FBte family crosswalk has an empty key/value")
        if te_id in mapping and mapping[te_id] != family:
            raise AssetCompatibilityError(f"conflicting crosswalk family for {te_id}")
        mapping[te_id] = family
    return mapping


def _parse_flag(row: Mapping[str, str], field: str) -> int:
    value = row.get(field, "")
    if value not in {"0", "1"}:
        raise AssetCompatibilityError(f"anchor archive field {field} is not a binary flag")
    return int(value)


def _parse_g1a_candidates(
    catalogue: Sequence[Mapping[str, str]],
    provenance: Sequence[Mapping[str, str]],
) -> tuple[dict[str, dict[str, str]], dict[str, str], dict[str, list[dict[str, Any]]], dict[str, str]]:
    by_query: dict[str, dict[str, str]] = {}
    for row in catalogue:
        query_id = row.get("query_id", "")
        if not query_id or query_id in by_query:
            raise AssetCompatibilityError("anchor catalogue has an empty/duplicate query_id")
        by_query[query_id] = dict(row)

    family_by_te: dict[str, str] = {}
    for row in provenance:
        te_id = row.get("te_id", "").strip()
        family = row.get("family", "").strip()
        if not te_id or not family:
            raise AssetCompatibilityError("anchor provenance has an empty TE/family label")
        previous = family_by_te.setdefault(te_id, family)
        if previous != family:
            raise AssetCompatibilityError(f"inconsistent anchor provenance family for {te_id}")
    if not family_by_te:
        raise AssetCompatibilityError("anchor provenance contains no TE records")

    candidate_sequences: dict[str, dict[str, str]] = {"terminal": {}, "internal": {}}
    for query_id, row in by_query.items():
        try:
            k = int(row["k"])
            low_complexity = _parse_flag(row, "low_complexity")
            terminal_flag = _parse_flag(row, "has_terminal_provenance")
            internal_flag = _parse_flag(row, "has_internal_provenance")
        except (KeyError, ValueError) as exc:
            raise AssetCompatibilityError("malformed anchor catalogue row") from exc
        if k != PROBE_LENGTH or low_complexity != 0:
            continue
        try:
            canonical = _canonical(row["query_seq"])
            declared_text = row["canonical_strand_key"].upper()
            declared = _canonical(declared_text)
        except KeyError as exc:
            raise AssetCompatibilityError("anchor catalogue candidate is incomplete") from exc
        if canonical != declared or declared_text != declared:
            raise AssetCompatibilityError("anchor catalogue strand-canonical key disagrees with query sequence")
        if terminal_flag:
            candidate_sequences["terminal"][query_id] = canonical
        if internal_flag:
            candidate_sequences["internal"][query_id] = canonical

    source_rows: dict[str, list[dict[str, Any]]] = {"terminal": [], "internal": []}
    allowed_regions = {
        "terminal": {"5p_terminal", "3p_terminal"},
        "internal": {"internal", "short_internal_fallback"},
    }
    for row in provenance:
        query_id = row.get("query_id", "")
        try:
            k = int(row["k"])
            start = int(row["start"])
        except (KeyError, ValueError) as exc:
            raise AssetCompatibilityError("malformed anchor provenance coordinates") from exc
        if k != PROBE_LENGTH or start < 0:
            continue
        for kind in ("terminal", "internal"):
            canonical = candidate_sequences[kind].get(query_id)
            if canonical is None or row.get("region") not in allowed_regions[kind]:
                continue
            source_seq = row.get("source_seq", "").upper()
            query_seq = row.get("query_seq", "").upper()
            if _canonical(source_seq) != canonical or _canonical(query_seq) != canonical:
                raise AssetCompatibilityError("anchor provenance sequence disagrees with candidate key")
            source_rows[kind].append(
                {
                    "canonical": canonical,
                    "te_id": row["te_id"].strip(),
                    "family": canonical_family(row["family"].strip()),
                    "start": start,
                    "end": start + PROBE_LENGTH,
                    "region": row["region"],
                    "source_orientation": row.get("source_orientation", ""),
                    "query_id": query_id,
                }
            )
    for kind in ("terminal", "internal"):
        if candidate_sequences[kind] and not source_rows[kind]:
            raise AssetCompatibilityError(f"anchor archive {kind} candidates lack matching provenance")
        qids_with_sources = {row["query_id"] for row in source_rows[kind]}
        if set(candidate_sequences[kind]) - qids_with_sources:
            raise AssetCompatibilityError(f"anchor archive {kind} catalogue/provenance mismatch")
        unique: dict[tuple[Any, ...], dict[str, Any]] = {}
        for row in source_rows[kind]:
            key = (
                row["canonical"], row["te_id"], row["family"], row["region"],
                row["start"], row["end"],
            )
            unique.setdefault(key, row)
        source_rows[kind] = sorted(
            unique.values(),
            key=lambda row: (
                row["te_id"], row["start"], row["region"], row["canonical"], row["query_id"]
            ),
        )
    return by_query, family_by_te, source_rows, {
        kind: len(set(candidate_sequences[kind].values())) for kind in candidate_sequences
    }


def _extract_te_id(record: FastaRecord) -> str:
    header = f"{record.record_id} {record.description}"
    if _FBTE_UNKNOWN.search(header):
        return "FBte000XXXX"
    match = _FBTE.search(header)
    if match:
        return "FBte" + re.sub(r"\D", "", match.group(0))
    return record.record_id.split("|", 1)[0]


def _encode(sequence: str) -> int:
    code = 0
    for base in sequence.encode("ascii"):
        value = _ENCODING[base] if base < 256 else -1
        if value < 0:
            raise ValueError(sequence)
        code = (code << 2) | value
    return code


def _query_map(candidates: Iterable[str]) -> dict[int, str]:
    mapping: dict[int, str] = {}
    for canonical in candidates:
        reverse = canonical.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        for oriented in {canonical, reverse}:
            code = _encode(oriented)
            previous = mapping.setdefault(code, canonical)
            if previous != canonical:
                raise AssetCompatibilityError("anchor query encoding collision")
    return mapping


def _scan(sequence: str, qmap: Mapping[int, str]) -> Iterable[tuple[int, int, str]]:
    mask = (1 << (2 * PROBE_LENGTH)) - 1
    code = 0
    valid = 0
    for end0, base in enumerate(sequence.encode("ascii"), start=1):
        value = _ENCODING[base] if base < 256 else -1
        if value < 0:
            code = 0
            valid = 0
            continue
        code = ((code << 2) | value) & mask
        valid += 1
        if valid >= PROBE_LENGTH:
            candidate = qmap.get(code)
            if candidate is not None:
                yield end0 - PROBE_LENGTH, end0, candidate


def _family_key(label: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", label.strip().lower())


def _annotation_te(name: str) -> str | None:
    if _FBTE_UNKNOWN.search(name or ""):
        return "FBte000XXXX"
    match = _FBTE.search(name or "")
    if not match:
        return None
    return "FBte" + re.sub(r"\D", "", match.group(0))


def _overlap_index(repeats: Sequence[RepeatRecord], crosswalk: Mapping[str, str]):
    by_contig: dict[str, list[tuple[int, int, str, str | None, str | None]]] = defaultdict(list)
    for repeat in repeats:
        te_id = _annotation_te(repeat.name)
        by_contig[repeat.chrom].append(
            (
                repeat.span0.start0, repeat.span0.end0, repeat.name, te_id,
                crosswalk.get(te_id) if te_id is not None else None,
            )
        )
    result = {}
    for contig, rows in by_contig.items():
        rows.sort(key=lambda row: (row[0], row[1], row[2]))
        maximum = -1
        prefix: list[int] = []
        for row in rows:
            maximum = max(maximum, row[1])
            prefix.append(maximum)
        result[contig] = (rows, [row[0] for row in rows], prefix)
    return result


def _overlaps(index, chrom: str, start: int, end: int):
    import bisect

    item = index.get(chrom)
    if item is None:
        return []
    rows, starts, prefix = item
    hi = bisect.bisect_left(starts, end)
    lo = hi
    while lo > 0 and prefix[lo - 1] > start:
        lo -= 1
    return [row for row in rows[lo:hi] if row[1] > start and row[0] < end]


def _reference_hit_class(
    target_family: str,
    target_te_id: str,
    start: int,
    end: int,
    overlaps: Sequence[tuple[int, int, str, str | None, str | None]],
) -> str:
    if not overlaps:
        return "NON_TE_OR_UNANNOTATED_HOST"
    full = [row for row in overlaps if row[0] <= start and row[1] >= end]
    partial = [row for row in overlaps if row not in full]
    if not full:
        return "TE_BOUNDARY_OR_PARTIAL"
    target_key = _family_key(target_family)
    target_is_record_label = target_family == target_te_id or target_family.startswith("FBte")
    safe_full = []
    unsafe_full = []
    for row in full:
        _, _, _, annotation_te, annotation_family = row
        if annotation_family and _family_key(annotation_family) == target_key:
            safe_full.append(row)
        elif target_is_record_label and annotation_te and annotation_te == target_te_id:
            safe_full.append(row)
        else:
            unsafe_full.append(row)
    unsafe_partial = []
    for row in partial:
        _, _, _, annotation_te, annotation_family = row
        if annotation_family and _family_key(annotation_family) == target_key:
            continue
        if target_is_record_label and annotation_te and annotation_te == target_te_id:
            continue
        unsafe_partial.append(row)
    if safe_full and not unsafe_full and not unsafe_partial:
        return "TARGET_CONTEXT_ANNOTATED_TE"
    if unsafe_full:
        if any(row[4] for row in unsafe_full):
            return "OTHER_FAMILY_ANNOTATED_TE"
        return "UNMAPPED_ANNOTATED_TE"
    return "MIXED_REFERENCE_CONTEXT"


def _overlap_span(start: int, end: int, chosen: Sequence[Mapping[str, Any]]) -> bool:
    return any(
        max(start, row["source_start"] if "source_start" in row else row["start"])
        < min(end, row["source_end"] if "source_end" in row else row["end"])
        for row in chosen
    )


def _safe_tsv(rows: Sequence[Mapping[str, str]]) -> bytes:
    fields = [
        "probe_id", "probe_class", "sequence", "query_sequences_both_strands", "family",
        "te_record_id", "source_start", "source_end", "source_strand", "source_region",
        "terminal_side", "independent_block_id", "compatible_families",
        "te_library_occurrence_count", "reference_occurrence_count", "reference_class",
        "source_te_distinct_occurrence_count", "source_exact_position_match",
        "positionally_interpretable", "terminal_distance_bp", "probe_rank_within_te_record",
        "library_role",
    ]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def build_frozen_candidate_codebooks(
    g1a_archive_path: str | Path,
    family_crosswalk_path: str | Path,
    te_library: TELibraryInput,
    family_assignments: Mapping[str, str],
    reference: ReferenceInput,
    repeat_context: RepeatContext,
) -> FrozenCodebookBuild:
    """Reassess only the anchor archive candidates against the prepared TE/reference assets."""

    archive_path = Path(g1a_archive_path)
    crosswalk_path = Path(family_crosswalk_path)
    catalogue, provenance = _archive_rows(archive_path)
    crosswalk = _crosswalk(crosswalk_path)
    _, family_by_te, sources, candidate_counts = _parse_g1a_candidates(catalogue, provenance)

    if reference.semantic_sha256 != repeat_context.reference_semantic_sha256:
        raise AssetCompatibilityError("RepeatMasker context is bound to another reference")
    if repeat_context is None:
        raise AssetCompatibilityError("prepared reference lacks RepeatMasker context")

    library_by_te: dict[str, FastaRecord] = {}
    library_family: dict[str, str] = {}
    for record in te_library.records:
        te_id = _extract_te_id(record)
        if te_id in library_by_te:
            raise AssetCompatibilityError(f"TE library has duplicate record identity: {te_id}")
        family = family_assignments.get(record.record_id, family_assignments.get(te_id))
        if family is None:
            raise AssetCompatibilityError(f"TE library family assignment is missing: {record.record_id}")
        library_by_te[te_id] = record
        library_family[te_id] = canonical_family(family)

    missing_g1a = sorted(set(family_by_te) - set(library_by_te))
    if missing_g1a:
        raise AssetCompatibilityError(f"TE library is missing anchor archive records: {missing_g1a[:20]}")
    for te_id, frozen_family in family_by_te.items():
        source_family = canonical_family(frozen_family)
        if _family_key(library_family[te_id]) != _family_key(source_family):
            raise AssetCompatibilityError(f"TE-family map disagrees with anchor provenance: {te_id}")
        mapped_family = crosswalk.get(te_id)
        if mapped_family is not None and _family_key(mapped_family) != _family_key(source_family):
            raise AssetCompatibilityError(f"FBte crosswalk disagrees with anchor provenance family: {te_id}")

    candidate_sequences = {
        kind: {row["canonical"] for row in source_rows}
        for kind, source_rows in sources.items()
    }
    all_queries = candidate_sequences["terminal"] | candidate_sequences["internal"]
    qmap = _query_map(all_queries)

    te_hits: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    positions: dict[str, dict[str, set[tuple[int, int]]]] = defaultdict(lambda: defaultdict(set))
    for te_id in sorted(library_by_te):
        record = library_by_te[te_id]
        family = library_family[te_id]
        for start, end, candidate in _scan(record.sequence, qmap):
            te_hits[candidate].append((te_id, family, start, end))
            positions[candidate][te_id].add((start, end))

    source_labels: dict[str, dict[str, set[tuple[str, str]]]] = {
        "terminal": defaultdict(set),
        "internal": defaultdict(set),
    }
    for kind, source_rows in sources.items():
        for row in source_rows:
            source_labels[kind][row["canonical"]].add((row["family"], row["te_id"]))
        for probe, labels in source_labels[kind].items():
            if any(
                _family_key(library_family.get(te_id, "")) != _family_key(family)
                for family, te_id in labels
            ):
                raise AssetCompatibilityError(f"candidate provenance family map mismatch: {probe}")

    hit_families = {
        probe: {family for _, family, _, _ in te_hits.get(probe, ())}
        for probe in all_queries
    }
    reference_index = _overlap_index(repeat_context.records, crosswalk)
    reference_counts: Counter[str] = Counter()
    reference_context_counts: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for record in reference.records:
        for start, end, candidate in _scan(record.sequence, qmap):
            reference_counts[candidate] += 1
            overlaps = _overlaps(reference_index, record.record_id, start, end)
            for kind in ("terminal", "internal"):
                labels = source_labels[kind].get(candidate, set())
                classes = [
                    _reference_hit_class(family, te_id, start, end, overlaps)
                    for family, te_id in sorted(labels)
                ]
                classification = (
                    "TARGET_CONTEXT_ANNOTATED_TE"
                    if classes and all(value == "TARGET_CONTEXT_ANNOTATED_TE" for value in classes)
                    else next(
                        (
                            value
                            for value in (
                                "OTHER_FAMILY_ANNOTATED_TE", "NON_TE_OR_UNANNOTATED_HOST",
                                "UNMAPPED_ANNOTATED_TE", "TE_BOUNDARY_OR_PARTIAL",
                                "MIXED_REFERENCE_CONTEXT",
                            )
                            if value in classes
                        ),
                        "NO_SOURCE_PROVENANCE",
                    )
                )
                reference_context_counts[(kind, candidate)][classification] += 1

    family_safe: dict[str, dict[str, bool]] = {"terminal": {}, "internal": {}}
    reference_safe: dict[str, dict[str, bool]] = {"terminal": {}, "internal": {}}
    acquisition_safe: dict[str, dict[str, bool]] = {"terminal": {}, "internal": {}}
    reference_class: dict[str, dict[str, str]] = {"terminal": {}, "internal": {}}
    for kind in ("terminal", "internal"):
        for probe in candidate_sequences[kind]:
            labels = source_labels[kind].get(probe, set())
            source_families = {family for family, _ in labels}
            compatible = (
                bool(labels)
                and len(source_families) == 1
                and len(hit_families[probe]) == 1
                and hit_families[probe] == source_families
            )
            context_safe = (
                reference_counts[probe] == 0
                or reference_context_counts[(kind, probe)]["TARGET_CONTEXT_ANNOTATED_TE"]
                == reference_counts[probe]
            )
            family_safe[kind][probe] = compatible
            reference_safe[kind][probe] = context_safe
            acquisition_safe[kind][probe] = compatible and context_safe
            reference_class[kind][probe] = (
                "REFERENCE_ABSENT"
                if reference_counts[probe] == 0
                else "TARGET_CONTEXT_ANNOTATED_TE_ONLY"
                if context_safe
                else "UNSAFE_OR_MIXED_REFERENCE_CONTEXT"
            )

    def assessed_sources(kind: str) -> list[dict[str, Any]]:
        assessed = []
        for source in sources[kind]:
            te_id = source["te_id"]
            record = library_by_te.get(te_id)
            if record is None or source["end"] > len(record.sequence):
                raise AssetCompatibilityError(f"anchor source interval lies outside TE record: {te_id}")
            distinct = sorted(positions[source["canonical"]].get(te_id, set()))
            source_match = (source["start"], source["end"]) in distinct
            if not source_match:
                continue
            row = dict(source)
            row.update(
                {
                    "reference_occurrence_count": reference_counts[source["canonical"]],
                    "reference_class": reference_class[kind][source["canonical"]],
                    "source_te_distinct_occurrence_count": len(distinct),
                    "source_exact_position_match": True,
                    "acquisition_safe": acquisition_safe[kind][source["canonical"]],
                    "te_record_length": len(record.sequence),
                    "te_library_occurrence_count": len(te_hits[source["canonical"]]),
                    "record_sequence": record.sequence,
                }
            )
            if kind == "terminal":
                side = "5p" if source["region"] == "5p_terminal" else "3p"
                distance = source["start"] if side == "5p" else len(record.sequence) - source["end"]
                opposite_recurrences = 0
                for start, end in distinct:
                    if (start, end) == (source["start"], source["end"]):
                        continue
                    if side == "5p" and len(record.sequence) - end <= distance:
                        opposite_recurrences += 1
                    elif side == "3p" and start <= distance:
                        opposite_recurrences += 1
                row.update(
                    {
                        "terminal_side": side,
                        "terminal_distance_bp": distance,
                        "positionally_interpretable": bool(
                            acquisition_safe[kind][source["canonical"]]
                            and len(distinct) == 1
                            and opposite_recurrences == 0
                        ),
                    }
                )
            else:
                row.update(
                    {
                        "terminal_side": "",
                        "terminal_distance_bp": "",
                        "positionally_interpretable": False,
                    }
                )
            assessed.append(row)
        return assessed

    terminal_sources = assessed_sources("terminal")
    eligible_terminal: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in terminal_sources:
        if row["acquisition_safe"]:
            eligible_terminal[row["te_id"]].append(row)

    selected_terminal: list[dict[str, str]] = []
    selected_terminal_sources: list[dict[str, Any]] = []
    for te_id in sorted(family_by_te):
        candidates = eligible_terminal.get(te_id, [])
        score = lambda row: (
            -int(row["positionally_interpretable"]),
            row["reference_occurrence_count"],
            row["terminal_distance_bp"],
            row["canonical"],
            row["start"],
        )
        by_side = {
            side: sorted((row for row in candidates if row["terminal_side"] == side), key=score)
            for side in ("5p", "3p")
        }
        chosen: list[dict[str, Any]] = []
        used_sequences: set[str] = set()
        for side in ("5p", "3p"):
            for row in by_side[side]:
                if sum(item["terminal_side"] == side for item in chosen) >= 2:
                    break
                if row["canonical"] in used_sequences or _overlap_span(row["start"], row["end"], chosen):
                    continue
                chosen.append(row)
                used_sequences.add(row["canonical"])
        for row in sorted(candidates, key=score):
            if len(chosen) >= 4:
                break
            if row["canonical"] in used_sequences or _overlap_span(row["start"], row["end"], chosen):
                continue
            chosen.append(row)
            used_sequences.add(row["canonical"])
        chosen.sort(key=lambda row: (row["terminal_side"], row["start"], 21, row["canonical"]))
        for rank, source in enumerate(chosen, start=1):
            selected_terminal_sources.append(source)
            selected_terminal.append(
                _codebook_row(source, "TERMINAL_21MER", rank, family_by_te[te_id])
            )

    terminal_sequences = {row["canonical"] for row in selected_terminal_sources}
    internal_sources = assessed_sources("internal")
    safe_internal_sources = sorted(
        (row for row in internal_sources if row["acquisition_safe"]),
        key=lambda row: (row["te_id"], row["start"], row["end"], row["canonical"], row["region"]),
    )
    internal_ranks: Counter[str] = Counter()
    selected_internal: list[dict[str, str]] = []
    selected_internal_only: list[dict[str, str]] = []
    for source in safe_internal_sources:
        internal_ranks[source["te_id"]] += 1
        output_row = _codebook_row(
            source,
            "INTERNAL_21MER",
            internal_ranks[source["te_id"]],
            family_by_te[source["te_id"]],
        )
        selected_internal.append(output_row)
        if source["canonical"] not in terminal_sequences:
            selected_internal_only.append(output_row)
    counts = {
        **{f"{kind}_candidate_sequences": count for kind, count in candidate_counts.items()},
        "terminal_selected": len(selected_terminal),
        "internal_full_sources": len(selected_internal),
        "internal_only_sources": len(selected_internal_only),
        "exhaustively_family_safe": sum(sum(values.values()) for values in family_safe.values()),
        "reference_context_safe": sum(sum(values.values()) for values in reference_safe.values()),
    }
    return FrozenCodebookBuild(
        terminal_tsv=_safe_tsv(selected_terminal),
        internal_tsv=_safe_tsv(selected_internal),
        internal_only_tsv=_safe_tsv(selected_internal_only),
        terminal_rows=tuple(selected_terminal),
        internal_rows=tuple(selected_internal),
        internal_only_rows=tuple(selected_internal_only),
        candidate_counts=counts,
    )


def _codebook_row(
    source: Mapping[str, Any],
    probe_class: str,
    rank: int,
    family: str,
) -> dict[str, str]:
    canonical = source["canonical"]
    reverse = canonical.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    record_sequence = source.get("record_sequence", "")
    source_sequence = record_sequence[source["start"] : source["end"]].upper()
    reverse_source = source_sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    strand = "+" if source_sequence == canonical else "-" if reverse_source == canonical else "+"
    identity = "|".join(
        (
            probe_class, source["te_id"], str(source["start"]), str(source["end"]),
            source["region"], canonical,
        )
    )
    probe_id = f"{probe_class.lower()}_{hashlib.sha256(identity.encode()).hexdigest()[:20]}"
    return {
        "probe_id": probe_id,
        "probe_class": probe_class,
        "sequence": canonical,
        "query_sequences_both_strands": canonical if canonical == reverse else f"{canonical};{reverse}",
        "family": family,
        "te_record_id": source["te_id"],
        "source_start": str(source["start"]),
        "source_end": str(source["end"]),
        "source_strand": strand,
        "source_region": source["region"],
        "terminal_side": str(source.get("terminal_side", "")),
        "independent_block_id": f"{source['te_id']}:{source['start'] // 125}",
        "compatible_families": json.dumps([family], separators=(",", ":")),
        "te_library_occurrence_count": str(source["te_library_occurrence_count"]),
        "reference_occurrence_count": str(source["reference_occurrence_count"]),
        "reference_class": str(source["reference_class"]),
        "source_te_distinct_occurrence_count": str(source["source_te_distinct_occurrence_count"]),
        "source_exact_position_match": "1",
        "positionally_interpretable": "1" if source["positionally_interpretable"] else "0",
        "terminal_distance_bp": str(source.get("terminal_distance_bp", "")),
        "probe_rank_within_te_record": str(rank),
        "library_role": "FROZEN_G1A_SAFE_ACQUISITION",
    }
