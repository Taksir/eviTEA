"""Synthetic CIGAR evidence-completeness and finite-alternative competition."""

import json
from pathlib import Path

from evitea.adjudication.reference import (
    ReferenceAlignmentObservation,
    aligned_query_spans,
    build_reference_alternatives,
    evaluate_reference_explanation,
    produce_reference_observation_state,
)
from evitea.backends.protocols import RawAlignment
from evitea.contracts import (
    BaseSpan0,
    CensorState,
    FamilyResolutionState,
    ReadMate,
    ReferenceExplanationCode,
    TEEligibilityState,
    TEFamilyEvidence,
)

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "phase4" / "F-REF-001"


def _evidence(pair_id):
    return TEFamilyEvidence(
        sample_id="synthetic-sample",
        pair_id=pair_id,
        family_state=FamilyResolutionState.RESOLVED_SINGLETON,
        eligibility_state=TEEligibilityState.ELIGIBLE,
        resolved_family="FamilyA",
        te_positive_mate=ReadMate.READ1,
        families_observed=("FamilyA",),
        probe_ids=("safe-a", "safe-b"),
        independent_block_ids=("block-a", "block-b"),
        contributing_te_record_ids=("TE_A",),
    )


def _observations(rows):
    return tuple(
        ReferenceAlignmentObservation(
            alternative_id=row["alternative_id"],
            placement_rank=row["placement_rank"],
            chrom=row["chrom"],
            te_alignment_span0=BaseSpan0(*row["te_alignment_span0"]),
            cigar=row["cigar"],
            is_local_proper=row["is_local_proper"],
            censor_state=CensorState[row.get("censor_state", "COMPLETE")],
        )
        for row in rows
    )


def test_cigar_only_counts_m_equals_x_query_bases_as_explained():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    expected = fixture["cigar_semantics"]["explained_query_spans0"]

    spans = aligned_query_spans(fixture["cigar_semantics"]["cigar"])

    assert [[span.start0, span.end0] for span in spans] == expected


def test_all_safe_occurrences_and_all_finite_alternatives_are_preserved():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    expected = json.loads((FIXTURE_DIR / "expected.json").read_text(encoding="utf-8"))
    required = tuple(BaseSpan0(*span) for span in fixture["required_query_spans0"])
    alternatives = build_reference_alternatives(
        _evidence(0), required, _observations(fixture["pair0_alternatives"])
    )
    state = evaluate_reference_explanation(_evidence(0), alternatives)

    target = expected["pair0"]
    assert state.code is ReferenceExplanationCode[target["code"]]
    assert state.finite_alternative_ids == tuple(target["finite_alternative_ids"])
    assert state.evidence_complete_alternative_ids == tuple(
        target["evidence_complete_alternative_ids"]
    )
    first = alternatives[0]
    assert first.chrom == "chrResidentToy"
    assert first.evidence_complete is False
    assert first.unexplained_query_spans0 == (BaseSpan0(7, 10),)


def test_one_missing_safe_occurrence_is_reference_unexplained():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    expected = json.loads((FIXTURE_DIR / "expected.json").read_text(encoding="utf-8"))
    required = tuple(BaseSpan0(*span) for span in fixture["required_query_spans0"])
    alternatives = build_reference_alternatives(
        _evidence(1), required, _observations(fixture["pair1_alternatives"])
    )

    state = evaluate_reference_explanation(_evidence(1), alternatives)

    assert state.code is ReferenceExplanationCode[expected["pair1"]["code"]]
    assert alternatives[0].unexplained_query_spans0 == (BaseSpan0(7, 10),)


def test_resource_censor_without_complete_alternative_stays_unresolved():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    expected = json.loads((FIXTURE_DIR / "expected.json").read_text(encoding="utf-8"))
    required = tuple(BaseSpan0(*span) for span in fixture["required_query_spans0"])
    alternatives = build_reference_alternatives(
        _evidence(2), required, _observations(fixture["pair2_alternatives"])
    )

    state = evaluate_reference_explanation(
        _evidence(2), alternatives, search_censor_state=CensorState.REFERENCE_RESOURCE_CENSORED
    )

    assert state.code is ReferenceExplanationCode[expected["pair2"]["code"]]
    assert state.censor_state is CensorState[expected["pair2"]["censor_state"]]
    assert state.evidence_complete_alternative_ids == ()


def test_resource_censor_overrides_an_observed_complete_alternative():
    fixture = json.loads((FIXTURE_DIR / "input.json").read_text(encoding="utf-8"))
    required = tuple(BaseSpan0(*span) for span in fixture["required_query_spans0"])
    alternatives = build_reference_alternatives(
        _evidence(0), required, _observations(fixture["pair0_alternatives"])
    )
    complete_ids = tuple(
        alternative.alternative_id
        for alternative in alternatives
        if alternative.evidence_complete
    )

    state = evaluate_reference_explanation(
        _evidence(0),
        alternatives,
        search_censor_state=CensorState.REFERENCE_RESOURCE_CENSORED,
    )

    assert state.code is ReferenceExplanationCode.INDETERMINATE_CENSORED
    assert state.evidence_complete_alternative_ids == complete_ids


def _paired_row(*, mate, flag, position0, mate_position0, cigar="10M"):
    return RawAlignment(
        query_name="P9",
        pair_id=9,
        mate=mate,
        flag=flag,
        reference_name="chr2L",
        position0=position0,
        mapq=30,
        cigar=cigar,
        sequence="A" * 10,
        mate_reference_name="chr2L",
        mate_position0=mate_position0,
        template_length=150,
    )


def test_s30_producer_retains_secondary_te_rows_in_complete_pair_signatures():
    evidence = _evidence(9)
    rows = (
        _paired_row(
            mate=ReadMate.READ1,
            flag=0x1 | 0x2 | 0x40,
            position0=99,
            mate_position0=199,
        ),
        _paired_row(
            mate=ReadMate.READ1,
            flag=0x1 | 0x2 | 0x40 | 0x100,
            position0=99,
            mate_position0=199,
            cigar="2S8M",
        ),
        _paired_row(
            mate=ReadMate.READ2,
            flag=0x1 | 0x2 | 0x80,
            position0=199,
            mate_position0=99,
        ),
        _paired_row(
            mate=ReadMate.READ1,
            flag=0x1 | 0x2 | 0x40,
            position0=399,
            mate_position0=499,
        ),
    )

    result = produce_reference_observation_state(
        evidence,
        (BaseSpan0(0, 2),),
        rows,
    )

    assert len(result.observations) == 2
    assert [row.placement_rank for row in result.observations] == [1, 2]
    assert result.alternatives[0].evidence_complete is True
    assert result.alternatives[1].unexplained_query_spans0 == (BaseSpan0(0, 2),)
    assert result.state.code is ReferenceExplanationCode.EVIDENCE_COMPLETE


def test_s30_producer_requires_both_mates_and_propagates_search_censor():
    evidence = _evidence(9)
    singleton = (
        _paired_row(
            mate=ReadMate.READ1,
            flag=0x1 | 0x2 | 0x40,
            position0=99,
            mate_position0=199,
        ),
    )

    complete = produce_reference_observation_state(
        evidence,
        (BaseSpan0(0, 2),),
        singleton,
    )
    censored = produce_reference_observation_state(
        evidence,
        (BaseSpan0(0, 2),),
        (),
        search_censor_state=CensorState.REFERENCE_RESOURCE_CENSORED,
    )

    assert complete.observations == ()
    assert complete.state.code is ReferenceExplanationCode.NO_LOCAL_PROPER
    assert censored.state.code is ReferenceExplanationCode.INDETERMINATE_CENSORED

