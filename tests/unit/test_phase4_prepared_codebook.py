"""Runtime loading of prepared safe-language query units."""

import csv

import pytest

from evitea.contracts import (
    BaseSpan0,
    ProbeClass,
    ReadMate,
    TEEligibilityState,
    TEEventOrientation,
    TEProbeHit,
)
from evitea.errors import AssetCompatibilityError
from evitea.evidence.codebook import load_probe_language_tsv, reverse_complement
from evitea.evidence.fragments import decode_pair_family


FIELDS = (
    "probe_id",
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
)


def _write(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _row(sequence, record, start, block, family="FamilyA"):
    return {
        "probe_id": f"source-{record}-{start}",
        "probe_class": "TERMINAL_21MER",
        "sequence": sequence,
        "query_sequences_both_strands": f"{sequence};{reverse_complement(sequence)}",
        "family": family,
        "te_record_id": record,
        "source_start": str(start),
        "source_end": str(start + 21),
        "source_strand": "+",
        "independent_block_id": block,
        "compatible_families": f'["{family}"]',
    }


def test_loader_collapses_shared_query_but_preserves_all_record_memberships(tmp_path):
    first_raw = "AACCGGTTAACCGGTTAACCG"
    second_raw = "ACGTACGTACGTACGTACGTA"
    first = min(first_raw, reverse_complement(first_raw))
    second = min(second_raw, reverse_complement(second_raw))
    path = tmp_path / "terminal.tsv"
    _write(
        path,
        (
            _row(first, "TE_A", 0, "TE_A:0"),
            _row(first, "TE_B", 10, "TE_B:0"),
            _row(second, "TE_A", 150, "TE_A:1"),
        ),
    )

    language = load_probe_language_tsv(path, ProbeClass.TERMINAL_21MER)
    first_probe = next(probe for probe in language.probes if probe.sequence == first)
    second_probe = next(probe for probe in language.probes if probe.sequence == second)
    hits = (
        TEProbeHit(
            "sample", 3, ReadMate.READ1, first_probe.probe_id,
            BaseSpan0(0, 21), TEEventOrientation.FORWARD,
        ),
        TEProbeHit(
            "sample", 3, ReadMate.READ1, second_probe.probe_id,
            BaseSpan0(22, 43), TEEventOrientation.FORWARD,
        ),
    )

    evidence = decode_pair_family(
        "sample",
        3,
        hits,
        language.codebook,
        probe_memberships=language.memberships,
    )

    assert len(language.probes) == 2
    assert len(language.memberships[first_probe.probe_id]) == 2
    assert evidence.contributing_te_record_ids == ("TE_A", "TE_B")
    assert evidence.independent_block_ids == ("TE_A:0", "TE_A:1", "TE_B:0")
    assert evidence.eligibility_state is TEEligibilityState.ELIGIBLE


def test_loader_rejects_shared_query_with_cross_family_memberships(tmp_path):
    raw = "AACCGGTTAACCGGTTAACCG"
    sequence = min(raw, reverse_complement(raw))
    path = tmp_path / "conflict.tsv"
    _write(
        path,
        (
            _row(sequence, "TE_A", 0, "TE_A:0"),
            _row(sequence, "TE_B", 0, "TE_B:0", family="FamilyB"),
        ),
    )

    with pytest.raises(AssetCompatibilityError, match="multiple TE families"):
        load_probe_language_tsv(path, ProbeClass.TERMINAL_21MER)
