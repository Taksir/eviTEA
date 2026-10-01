"""Fail-closed raw FASTQ, FASTA, family-map, and repeat parsers."""

from __future__ import annotations

import gzip
import hashlib
import re
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

from .contracts import (
    BaseSpan0,
    FastaRecord,
    RawReadPair,
    ReadMate,
    ReadRecord,
    ReferenceInput,
    RepeatContext,
    RepeatRecord,
    TELibraryInput,
)
from .errors import ArchitectureInvariantError, InputContractError, PairSynchronizationError
from .families import FamilyMap, canonical_family

_FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_ASCII_HEADER_SPLIT = re.compile(r"[\t\f\v ]+")
_CASAVA_MATE = re.compile(r"^([12]):")


def _path(path: str | Path) -> Path:
    result = Path(path)
    if not result.is_file():
        raise InputContractError(f"input file does not exist: {result}")
    return result


def _has_suffix(path: Path, suffixes: tuple[str, ...]) -> bool:
    lower = path.name.casefold()
    return any(lower.endswith(suffix) for suffix in suffixes)


def _open_text(path: Path) -> TextIO:
    if path.name.casefold().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open("rt", encoding="ascii", newline="")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _strip_newline(line: str) -> str:
    return line[:-2] if line.endswith("\r\n") else line[:-1] if line.endswith("\n") else line


def _normalize_fastq_header(header: str, expected_mate: ReadMate) -> str:
    if not header.startswith("@"):
        raise InputContractError("FASTQ header must begin with '@'")
    content = header[1:]
    if not content:
        raise InputContractError("FASTQ header name is empty")
    fields = _ASCII_HEADER_SPLIT.split(content, maxsplit=1)
    token = fields[0]
    if not token:
        raise InputContractError("FASTQ header name is empty")

    terminal_mate: ReadMate | None = None
    if token.endswith("/1"):
        terminal_mate = ReadMate.READ1
        token = token[:-2]
    elif token.endswith("/2"):
        terminal_mate = ReadMate.READ2
        token = token[:-2]
    if terminal_mate is not None and terminal_mate is not expected_mate:
        raise PairSynchronizationError(
            f"FASTQ terminal mate marker does not match {expected_mate.name}"
        )
    if not token:
        raise InputContractError("normalized FASTQ pair name is empty")

    if len(fields) == 2:
        second_token = _ASCII_HEADER_SPLIT.split(fields[1], maxsplit=1)[0]
        match = _CASAVA_MATE.match(second_token)
        if match is not None and int(match.group(1)) != int(expected_mate):
            raise PairSynchronizationError(
                f"CASAVA mate marker does not match {expected_mate.name}"
            )
    return token


def _validate_sequence(sequence: str, context: str) -> str:
    if not sequence:
        raise InputContractError(f"{context} sequence is empty")
    upper = sequence.upper()
    invalid = sorted(set(upper) - _IUPAC_DNA)
    if invalid:
        raise InputContractError(f"{context} contains invalid IUPAC bases: {invalid}")
    return upper


def _iter_fastq_records(path: Path, mate: ReadMate) -> Iterator[ReadRecord]:
    try:
        with _open_text(path) as handle:
            record0 = 0
            while True:
                header_line = handle.readline()
                if header_line == "":
                    break
                sequence_line = handle.readline()
                plus_line = handle.readline()
                quality_line = handle.readline()
                if "" in (sequence_line, plus_line, quality_line):
                    raise InputContractError(
                        f"incomplete FASTQ record {record0} in {path}"
                    )
                header = _strip_newline(header_line)
                sequence = _validate_sequence(
                    _strip_newline(sequence_line), f"FASTQ record {record0}"
                )
                plus = _strip_newline(plus_line)
                quality = _strip_newline(quality_line)
                if not plus.startswith("+"):
                    raise InputContractError(
                        f"FASTQ record {record0} separator must begin with '+'"
                    )
                if len(sequence) != len(quality):
                    raise InputContractError(
                        f"FASTQ record {record0} sequence/quality lengths differ"
                    )
                if any(ord(character) < 33 or ord(character) > 126 for character in quality):
                    raise InputContractError(
                        f"FASTQ record {record0} quality is not printable ASCII"
                    )
                pair_name = _normalize_fastq_header(header, mate)
                yield ReadRecord(
                    header=header,
                    pair_name=pair_name,
                    sequence=sequence,
                    quality=quality,
                    mate=mate,
                )
                record0 += 1
    except (OSError, UnicodeError) as exc:
        raise InputContractError(f"cannot read FASTQ {path}: {exc}") from exc


def iter_paired_fastq(
    read1: str | Path,
    read2: str | Path,
) -> Iterator[RawReadPair]:
    """Yield synchronized read pairs without resynchronization or search-ahead."""

    read1_path = _path(read1)
    read2_path = _path(read2)
    if not _has_suffix(read1_path, _FASTQ_SUFFIXES):
        raise InputContractError(f"unsupported FASTQ suffix: {read1_path}")
    if not _has_suffix(read2_path, _FASTQ_SUFFIXES):
        raise InputContractError(f"unsupported FASTQ suffix: {read2_path}")

    records1 = _iter_fastq_records(read1_path, ReadMate.READ1)
    records2 = _iter_fastq_records(read2_path, ReadMate.READ2)
    ordinal0 = 0
    while True:
        try:
            record1 = next(records1)
            ended1 = False
        except StopIteration:
            record1 = None
            ended1 = True
        try:
            record2 = next(records2)
            ended2 = False
        except StopIteration:
            record2 = None
            ended2 = True

        if ended1 and ended2:
            return
        if ended1 != ended2:
            raise PairSynchronizationError(
                f"paired FASTQ streams have different lengths at ordinal {ordinal0}"
            )
        assert record1 is not None and record2 is not None
        if record1.pair_name != record2.pair_name:
            raise PairSynchronizationError(
                f"paired FASTQ names differ at ordinal {ordinal0}: "
                f"{record1.pair_name!r} != {record2.pair_name!r}"
            )
        yield RawReadPair(record1.pair_name, record1, record2)
        ordinal0 += 1


def _iter_fasta(path: Path) -> Iterator[FastaRecord]:
    try:
        with _open_text(path) as handle:
            record_id: str | None = None
            description = ""
            sequence_parts: list[str] = []
            rank = 0
            for line0, raw_line in enumerate(handle):
                line = _strip_newline(raw_line)
                if not line:
                    continue
                if line.startswith(">"):
                    if record_id is not None:
                        sequence = _validate_sequence(
                            "".join(sequence_parts), f"FASTA record {record_id!r}"
                        )
                        yield FastaRecord(record_id, description, sequence, rank)
                        rank += 1
                    description = line[1:]
                    if not description:
                        raise InputContractError(f"empty FASTA header at line {line0 + 1}")
                    record_id = description.split()[0]
                    sequence_parts = []
                    continue
                if record_id is None:
                    raise InputContractError(
                        f"FASTA sequence precedes first header at line {line0 + 1}"
                    )
                if any(character.isspace() for character in line):
                    raise InputContractError(
                        f"embedded whitespace in FASTA sequence at line {line0 + 1}"
                    )
                sequence_parts.append(line)
            if record_id is not None:
                sequence = _validate_sequence(
                    "".join(sequence_parts), f"FASTA record {record_id!r}"
                )
                yield FastaRecord(record_id, description, sequence, rank)
    except (OSError, UnicodeError) as exc:
        raise InputContractError(f"cannot read FASTA {path}: {exc}") from exc


def _validated_fasta(path: str | Path) -> tuple[Path, tuple[FastaRecord, ...], str, str]:
    fasta_path = _path(path)
    records = tuple(_iter_fasta(fasta_path))
    if not records:
        raise InputContractError(f"FASTA is empty: {fasta_path}")
    ids = [record.record_id for record in records]
    if len(ids) != len(set(ids)):
        duplicates = sorted({record_id for record_id in ids if ids.count(record_id) > 1})
        raise InputContractError(f"duplicate FASTA IDs: {duplicates}")
    semantic = hashlib.sha256()
    for record in records:
        semantic.update(f">{record.record_id}\n{record.sequence}\n".encode("ascii"))
    return fasta_path, records, _sha256_file(fasta_path), semantic.hexdigest()


def validate_reference_fasta(path: str | Path) -> ReferenceInput:
    """Validate a host FASTA and bind raw plus canonical semantic identities."""

    fasta_path, records, raw_sha256, semantic_sha256 = _validated_fasta(path)
    return ReferenceInput(
        path=fasta_path,
        records=records,
        raw_sha256=raw_sha256,
        semantic_sha256=semantic_sha256,
        contig_lengths={record.record_id: len(record.sequence) for record in records},
    )


def validate_te_fasta(path: str | Path) -> TELibraryInput:
    """Validate a TE FASTA and bind raw plus canonical semantic identities."""

    fasta_path, records, raw_sha256, semantic_sha256 = _validated_fasta(path)
    return TELibraryInput(fasta_path, records, raw_sha256, semantic_sha256)


FLYBASE_FAMILY_COLUMNS = ("FBte_ID", "FlyBase_Family")


def is_flybase_family_table(path: str | Path) -> bool:
    """True for a FlyBase-style table with ``FBte_ID`` and ``FlyBase_Family`` columns."""

    try:
        with _open_text(_path(path)) as handle:
            for raw_line in handle:
                line = _strip_newline(raw_line)
                if line and not line.startswith("#"):
                    fields = line.split("\t")
                    return all(name in fields for name in FLYBASE_FAMILY_COLUMNS)
    except (OSError, UnicodeError) as exc:
        raise InputContractError(f"cannot read family map {path}: {exc}") from exc
    return False


def _load_flybase_family_table(map_path: Path) -> FamilyMap:
    assignments: dict[str, str] = {}
    with _open_text(map_path) as handle:
        header: list[str] | None = None
        for line0, raw_line in enumerate(handle):
            line = _strip_newline(raw_line)
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                id_index = header.index(FLYBASE_FAMILY_COLUMNS[0])
                family_index = header.index(FLYBASE_FAMILY_COLUMNS[1])
                continue
            record_id = fields[id_index].strip() if id_index < len(fields) else ""
            family = fields[family_index].strip() if family_index < len(fields) else ""
            if not record_id or not family:
                continue
            canonical = canonical_family(family)
            if assignments.setdefault(record_id, canonical) != canonical:
                raise InputContractError(f"conflicting family assignments for {record_id!r}")
    if not assignments:
        raise InputContractError("family map is empty")
    return FamilyMap(assignments)


def load_family_map(path: str | Path) -> FamilyMap:
    """Load a TE-record/family relation with frozen alias handling.

    Accepts the two-column ``te_record_id<TAB>family`` form or a FlyBase-style
    table (``FBte_ID`` / ``FlyBase_Family`` columns, as the engines read it).
    """

    map_path = _path(path)
    if is_flybase_family_table(map_path):
        try:
            return _load_flybase_family_table(map_path)
        except (OSError, UnicodeError) as exc:
            raise InputContractError(f"cannot read family map {map_path}: {exc}") from exc
    assignments: dict[str, str] = {}
    saw_data_row = False
    try:
        with _open_text(map_path) as handle:
            for line0, raw_line in enumerate(handle):
                line = _strip_newline(raw_line)
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if not saw_data_row and fields == ["te_record_id", "family"]:
                    saw_data_row = True
                    continue
                saw_data_row = True
                if len(fields) != 2:
                    raise InputContractError(
                        f"family-map line {line0 + 1} must have exactly two tab fields"
                    )
                record_id, family = fields
                if not record_id:
                    raise InputContractError(f"empty TE ID at family-map line {line0 + 1}")
                canonical = canonical_family(family)
                previous = assignments.get(record_id)
                if previous is not None and previous != canonical:
                    raise InputContractError(
                        f"conflicting family assignments for {record_id!r}"
                    )
                assignments[record_id] = canonical
    except (OSError, UnicodeError) as exc:
        raise InputContractError(f"cannot read family map {map_path}: {exc}") from exc
    if not assignments:
        raise InputContractError("family map is empty")
    return FamilyMap(assignments)


def load_repeat_annotation(
    path: str | Path,
    reference: ReferenceInput,
) -> RepeatContext:
    """Load BED6-compatible repeat context against prepared reference bounds."""

    bed_path = _path(path)
    records: list[RepeatRecord] = []
    try:
        with _open_text(bed_path) as handle:
            for line0, raw_line in enumerate(handle):
                line = _strip_newline(raw_line)
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 6:
                    raise InputContractError(
                        f"repeat line {line0 + 1} must have at least six tab fields"
                    )
                chrom, start_text, end_text, name, score, strand = fields[:6]
                if chrom not in reference.contig_lengths:
                    raise InputContractError(
                        f"repeat line {line0 + 1} names unknown contig {chrom!r}"
                    )
                try:
                    start0 = int(start_text)
                    end0 = int(end_text)
                except ValueError as exc:
                    raise InputContractError(
                        f"repeat line {line0 + 1} has non-integer coordinates"
                    ) from exc
                try:
                    span0 = BaseSpan0(start0, end0)
                except ArchitectureInvariantError as exc:
                    raise InputContractError(
                        f"repeat line {line0 + 1} has invalid coordinates"
                    ) from exc
                if span0.length == 0 or span0.end0 > reference.contig_lengths[chrom]:
                    raise InputContractError(
                        f"repeat line {line0 + 1} lies outside reference bounds"
                    )
                if strand not in {"+", "-", "."}:
                    raise InputContractError(
                        f"repeat line {line0 + 1} has invalid strand {strand!r}"
                    )
                if not name:
                    raise InputContractError(f"repeat line {line0 + 1} has empty name")
                records.append(
                    RepeatRecord(chrom, span0, name, score, strand, len(records))
                )
    except (OSError, UnicodeError) as exc:
        raise InputContractError(f"cannot read repeat annotation {bed_path}: {exc}") from exc
    return RepeatContext(
        path=bed_path,
        records=tuple(records),
        raw_sha256=_sha256_file(bed_path),
        reference_semantic_sha256=reference.semantic_sha256,
    )
