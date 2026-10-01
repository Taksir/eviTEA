"""Shell-free Bowtie2 candidate wrapper and streaming SAM parser."""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..assets.reference import PreparedReference, verify_reference_asset
from ..contracts import ReadMate
from ..evidence.projection import SampleReadIndex
from ..errors import AssetCompatibilityError, BackendCompatibilityError, InputContractError
from .protocols import MapperRead, RawAlignment, ReadIdentity

CANDIDATE_STATUS = "CANDIDATE_UNAPPROVED"
FROZEN_REFERENCE_MAX_INSERT = 805
FROZEN_REFERENCE_MAX_REPORTED = 20


@dataclass(frozen=True, slots=True)
class PairedReferenceBatch:
    alignments: tuple[RawAlignment, ...]
    capped_pair_ids: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class Bowtie2CompatibilityProfile:
    """Explicit mapper identity; candidate status is not equivalence approval."""

    version: str
    options: tuple[str, ...] = ()
    max_reported: int = 2
    min_mapq: int = 20
    status: str = CANDIDATE_STATUS

    def __post_init__(self) -> None:
        if not self.version:
            raise InputContractError("Bowtie2 profile version must be non-empty")
        if self.max_reported != 2:
            raise InputContractError("Phase 3 Bowtie2 profile requires max_reported=2")
        if self.min_mapq != 20:
            raise InputContractError("Phase 3 Bowtie2 profile requires min_mapq=20")
        if self.status != CANDIDATE_STATUS:
            raise InputContractError("Bowtie2 profile cannot claim unvalidated approval")
        reserved = {"-a", "--all", "-k", "-x", "-U", "-S", "-f", "-q"}
        if any(option in reserved for option in self.options):
            raise InputContractError(
                "Bowtie2 profile options cannot override input/output or reporting cap"
            )

    def to_document(self) -> dict[str, object]:
        return {
            "backend": "bowtie2",
            "version": self.version,
            "options": list(self.options),
            "max_reported": self.max_reported,
            "min_mapq": self.min_mapq,
            "status": self.status,
        }


def candidate_bowtie2_2_5_5(*, options: Sequence[str] = ()) -> Bowtie2CompatibilityProfile:
    return Bowtie2CompatibilityProfile(version="2.5.5", options=tuple(options))


def _identity_for_sam(
    query_name: str,
    flag: int,
    identities: Mapping[str, ReadIdentity | Mapping[ReadMate, ReadIdentity]],
) -> ReadIdentity:
    try:
        identity_value = identities[query_name]
    except KeyError as exc:
        raise AssetCompatibilityError(
            f"SAM query name is absent from explicit identity map: {query_name}"
        ) from exc
    flag_mates: list[ReadMate] = []
    if flag & 0x40:
        flag_mates.append(ReadMate.READ1)
    if flag & 0x80:
        flag_mates.append(ReadMate.READ2)
    if len(flag_mates) > 1:
        raise AssetCompatibilityError(f"SAM row has conflicting mate flags: {query_name}")
    if isinstance(identity_value, Mapping):
        if len(flag_mates) != 1:
            raise AssetCompatibilityError(
                f"paired SAM row lacks an unambiguous mate flag: {query_name}"
            )
        try:
            identity = identity_value[flag_mates[0]]
        except KeyError as exc:
            raise AssetCompatibilityError(
                f"SAM mate is absent from explicit identity map: {query_name}"
            ) from exc
    else:
        identity = identity_value
    if flag_mates and flag_mates[0] is not identity.mate:
        raise AssetCompatibilityError(f"SAM mate flag conflicts with identity map: {query_name}")
    return identity


def parse_sam_line(
    line: str,
    identities: Mapping[str, ReadIdentity | Mapping[ReadMate, ReadIdentity]],
) -> RawAlignment | None:
    """Parse one SAM line; headers return ``None``."""

    if line.startswith("@"):
        return None
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11:
        raise BackendCompatibilityError("SAM row has fewer than 11 fields")
    query_name = fields[0]
    try:
        flag = int(fields[1])
        position1 = int(fields[3])
        mapq = int(fields[4])
        mate_position1 = int(fields[7])
        template_length = int(fields[8])
    except ValueError as exc:
        raise BackendCompatibilityError("SAM integer field is invalid") from exc
    if flag < 0 or mapq < 0:
        raise BackendCompatibilityError("SAM flag/MAPQ must be non-negative")
    identity = _identity_for_sam(query_name, flag, identities)
    unmapped = bool(flag & 0x4)
    if unmapped:
        if fields[2] != "*" or position1 != 0:
            raise BackendCompatibilityError("unmapped SAM row has a mapped coordinate")
        reference_name = None
        position0 = None
        cigar = None if fields[5] == "*" else fields[5]
    else:
        if fields[2] == "*" or position1 < 1 or fields[5] == "*":
            raise BackendCompatibilityError("mapped SAM row lacks a valid placement")
        reference_name = fields[2]
        position0 = position1 - 1
        cigar = fields[5]
    mate_reference_name = fields[6]
    if mate_reference_name == "=":
        mate_reference_name = reference_name
    elif mate_reference_name == "*":
        mate_reference_name = None
    mate_position0 = mate_position1 - 1 if mate_position1 > 0 else None
    return RawAlignment(
        query_name=query_name,
        pair_id=identity.pair_id,
        mate=identity.mate,
        flag=flag,
        reference_name=reference_name,
        position0=position0,
        mapq=mapq,
        cigar=cigar,
        sequence=fields[9],
        tags=tuple(fields[11:]),
        mate_reference_name=mate_reference_name,
        mate_position0=mate_position0,
        template_length=template_length,
    )


def iter_sam(
    lines: Iterable[str],
    identities: Mapping[str, ReadIdentity | Mapping[ReadMate, ReadIdentity]],
) -> Iterable[RawAlignment]:
    for line in lines:
        parsed = parse_sam_line(line, identities)
        if parsed is not None:
            yield parsed


def build_bowtie2_command(
    executable: str | Path,
    index_prefix: str | Path,
    profile: Bowtie2CompatibilityProfile,
    *,
    max_reported: int = 2,
) -> list[str]:
    if max_reported != profile.max_reported:
        raise BackendCompatibilityError("requested placement cap differs from profile")
    return [
        str(executable),
        "-x",
        str(index_prefix),
        "-f",
        "-U",
        "-",
        "-S",
        "-",
        "-k",
        str(max_reported),
        *profile.options,
    ]


def build_paired_reference_command(
    executable: str | Path,
    index_prefix: str | Path,
    read1_fasta: str | Path,
    read2_fasta: str | Path,
    output_sam: str | Path,
    *,
    threads: int = 1,
) -> list[str]:
    """Build the frozen evidence-complete paired-reference command."""

    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise InputContractError("threads must be a positive integer")
    return [
        str(executable),
        "-x",
        str(index_prefix),
        "-f",
        "-1",
        str(read1_fasta),
        "-2",
        str(read2_fasta),
        "--very-sensitive-local",
        "--no-mixed",
        "--no-discordant",
        "-I",
        "0",
        "-X",
        str(FROZEN_REFERENCE_MAX_INSERT),
        "-k",
        str(FROZEN_REFERENCE_MAX_REPORTED),
        "--no-unal",
        "--reorder",
        "-p",
        str(threads),
        "-S",
        str(output_sam),
    ]


def map_paired_reference_reads(
    executable: str | Path,
    sample_index: SampleReadIndex,
    pair_ids: Iterable[int],
    prepared_reference: PreparedReference,
    work_dir: str | Path,
    *,
    threads: int = 1,
) -> PairedReferenceBatch:
    """Run the frozen paired-reference producer over selected physical pairs."""

    verify_reference_asset(prepared_reference)
    index_value = prepared_reference.mapper_profile.get("index_prefix")
    if not isinstance(index_value, str) or not index_value:
        raise AssetCompatibilityError(
            "prepared reference mapper profile lacks index_prefix"
        )
    selected = tuple(sorted(set(pair_ids)))
    if any(isinstance(pair_id, bool) or not isinstance(pair_id, int) or pair_id < 0 for pair_id in selected):
        raise InputContractError("paired-reference pair IDs must be non-negative integers")
    root = Path(work_dir)
    root.mkdir(parents=True, exist_ok=True)
    read1 = root / "pairs_1.fa"
    read2 = root / "pairs_2.fa"
    sam = root / "local_reference.sam"
    if any(path.exists() for path in (read1, read2, sam)):
        raise InputContractError(
            "paired-reference work directory already contains stage outputs"
        )
    first_lines: list[str] = []
    second_lines: list[str] = []
    identities: dict[str, dict[ReadMate, ReadIdentity]] = {}
    backend_identity = 0
    for pair_id in selected:
        try:
            first = sample_index.sequences[(pair_id, ReadMate.READ1)]
            second = sample_index.sequences[(pair_id, ReadMate.READ2)]
        except KeyError as exc:
            raise InputContractError(
                f"paired-reference input lacks both mates for pair {pair_id}"
            ) from exc
        query_name = f"P{pair_id}"
        first_lines.extend((f">{query_name}", first))
        second_lines.extend((f">{query_name}", second))
        identities[query_name] = {
            ReadMate.READ1: ReadIdentity(
                backend_identity, pair_id, ReadMate.READ1
            ),
            ReadMate.READ2: ReadIdentity(
                backend_identity + 1, pair_id, ReadMate.READ2
            ),
        }
        backend_identity += 2
    read1.write_text(
        "\n".join(first_lines) + ("\n" if first_lines else ""), encoding="ascii"
    )
    read2.write_text(
        "\n".join(second_lines) + ("\n" if second_lines else ""), encoding="ascii"
    )
    command = build_paired_reference_command(
        executable,
        prepared_reference.root / index_value,
        read1,
        read2,
        sam,
        threads=threads,
    )
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            shell=False,
        )
    except OSError as exc:
        raise BackendCompatibilityError(
            "failed to start paired-reference Bowtie2"
        ) from exc
    if completed.returncode != 0:
        detail = completed.stderr[-2000:].strip()
        raise BackendCompatibilityError(
            f"paired-reference Bowtie2 exited with code {completed.returncode}: {detail}"
        )
    try:
        with sam.open("r", encoding="utf-8") as handle:
            alignments = tuple(iter_sam(handle, identities))
    except OSError as exc:
        raise BackendCompatibilityError(
            "paired-reference Bowtie2 did not produce readable SAM"
        ) from exc
    reported: dict[tuple[int, ReadMate], int] = {}
    for row in alignments:
        if row.is_mapped and not row.is_supplementary:
            key = (row.pair_id, row.mate)
            reported[key] = reported.get(key, 0) + 1
    capped = tuple(
        pair_id
        for pair_id in selected
        if any(
            reported.get((pair_id, mate), 0) >= FROZEN_REFERENCE_MAX_REPORTED
            for mate in (ReadMate.READ1, ReadMate.READ2)
        )
    )
    return PairedReferenceBatch(alignments, capped)


def _fasta_text(reads: Sequence[MapperRead]) -> tuple[str, dict[str, ReadIdentity]]:
    lines: list[str] = []
    identities: dict[str, ReadIdentity] = {}
    for ordinal0, read in enumerate(reads):
        if read.query_name in identities:
            raise InputContractError("mapper query names must be unique within a batch")
        if not read.sequence:
            raise InputContractError("mapper read sequence must be non-empty")
        identities[read.query_name] = ReadIdentity(ordinal0, read.pair_id, read.mate)
        lines.extend((f">{read.query_name}", read.sequence.upper()))
    return "\n".join(lines) + ("\n" if lines else ""), identities


class Bowtie2Backend:
    """Candidate external Bowtie2 runner with shell-free argument passing."""

    def __init__(
        self,
        executable: str | Path,
        profile: Bowtie2CompatibilityProfile,
    ) -> None:
        self.executable = Path(executable)
        self.profile = profile

    def map_reads(
        self,
        reads: Sequence[MapperRead],
        prepared_reference: PreparedReference,
        *,
        max_reported: int = 2,
    ) -> Iterable[RawAlignment]:
        verify_reference_asset(prepared_reference)
        index_value = prepared_reference.mapper_profile.get("index_prefix")
        if not isinstance(index_value, str) or not index_value:
            raise AssetCompatibilityError(
                "prepared reference mapper profile lacks index_prefix"
            )
        index_prefix = prepared_reference.root / index_value
        command = build_bowtie2_command(
            self.executable, index_prefix, self.profile, max_reported=max_reported
        )
        fasta, identities = _fasta_text(reads)
        with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stderr:
            try:
                process = subprocess.Popen(
                    command,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=stderr,
                    text=True,
                    shell=False,
                )
            except OSError as exc:
                raise BackendCompatibilityError("failed to start Bowtie2") from exc
            assert process.stdin is not None and process.stdout is not None
            try:
                process.stdin.write(fasta)
                process.stdin.close()
                yield from iter_sam(process.stdout, identities)
            except BaseException:
                process.kill()
                process.wait()
                raise
            finally:
                process.stdout.close()
            returncode = process.wait()
            if returncode != 0:
                stderr.seek(0)
                detail = stderr.read(2000).strip()
                raise BackendCompatibilityError(
                    f"Bowtie2 exited with code {returncode}: {detail}"
                )
