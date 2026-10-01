"""Structural checks for the synthetic backend observation contracts."""

from collections.abc import Iterable, Sequence
from pathlib import Path

from evitea.backends.bowtie2 import (
    build_paired_reference_command,
    parse_sam_line,
)
from evitea.backends.protocols import (
    ExactQueryBackend,
    HostMapperBackend,
    MapperRead,
    RawAlignment,
    RawReadHit,
    ReadIdentity,
)
from evitea.contracts import ReadMate


class _ExactBackend:
    def build(self, reads: Iterable[MapperRead], output_dir: str | Path) -> None:
        raise AssertionError("protocol checks must not build an index")

    def count_exact(self, sequence: str) -> int:
        return 0

    def locate_reads(self, sequence: str) -> Iterable[RawReadHit]:
        return ()


class _HostBackend:
    def map_reads(
        self,
        reads: Sequence[MapperRead],
        prepared_reference: object,
        *,
        max_reported: int = 2,
    ) -> Iterable[RawAlignment]:
        return ()


def test_backend_protocols_are_runtime_checkable_and_structural():
    assert isinstance(_ExactBackend(), ExactQueryBackend)
    assert isinstance(_HostBackend(), HostMapperBackend)


def test_raw_hit_carries_one_explicit_pair_mate_and_occurrence():
    hit = RawReadHit(
        sequence="ACGT",
        backend_row=23,
        backend_read_id=10,
        pair_id=0,
        mate=ReadMate.READ1,
        occurrence_ordinal0=2,
    )

    assert (hit.backend_row, hit.backend_read_id) == (23, 10)
    assert (hit.pair_id, hit.mate, hit.occurrence_ordinal0) == (0, ReadMate.READ1, 2)


def test_raw_alignment_exposes_sam_flags_without_collapsing_secondary_rows():
    primary = RawAlignment(
        query_name="toy-pair-0/1",
        pair_id=0,
        mate=ReadMate.READ1,
        flag=99,
        reference_name="chrToy",
        position0=0,
        mapq=42,
        cigar="8M",
        sequence="ACGTACGT",
    )
    secondary = RawAlignment(
        query_name="toy-pair-0/1",
        pair_id=0,
        mate=ReadMate.READ1,
        flag=355,
        reference_name="chrToy",
        position0=8,
        mapq=0,
        cigar="8M",
        sequence="ACGTACGT",
    )

    assert primary.is_mapped and not primary.is_secondary
    assert secondary.is_mapped and secondary.is_secondary
    assert primary.pair_id == secondary.pair_id == 0
    assert primary.mate is secondary.mate is ReadMate.READ1


def test_paired_sam_parser_resolves_shared_qname_and_mate_geometry():
    identities = {
        "P7": {
            ReadMate.READ1: ReadIdentity(10, 7, ReadMate.READ1),
            ReadMate.READ2: ReadIdentity(11, 7, ReadMate.READ2),
        }
    }
    row = parse_sam_line(
        "P7\t99\tchr2L\t101\t42\t10M\t=\t251\t160\tAAAAAAAAAA\t*\tNM:i:0\n",
        identities,
    )

    assert row is not None
    assert row.pair_id == 7
    assert row.mate is ReadMate.READ1
    assert row.mate_reference_name == "chr2L"
    assert row.mate_position0 == 250
    assert row.template_length == 160


def test_frozen_paired_reference_command_is_exact_and_ordered():
    command = build_paired_reference_command(
        "bowtie2", "ref/index", "r1.fa", "r2.fa", "out.sam", threads=3
    )

    assert command == [
        "bowtie2", "-x", "ref/index", "-f", "-1", "r1.fa", "-2", "r2.fa",
        "--very-sensitive-local", "--no-mixed", "--no-discordant",
        "-I", "0", "-X", "805", "-k", "20", "--no-unal", "--reorder",
        "-p", "3", "-S", "out.sam",
    ]

