from __future__ import annotations

import hashlib
import json
from pathlib import Path

from evitea.contracts import BreakpointSide, ClosedBreakpointInterval0
from evitea.incumbent.repeat_aware import (
    RepeatAwareGroupObservation,
    RepeatAwareLocusObservation,
    RepeatAwareObservations,
    RepeatAwarePreparedAssets,
    RepeatAwareSampleContext,
    RepeatAwareVariant,
    alignment_variant_ops,
    run_repeat_aware_incumbent,
)


FIXTURE_ROOT = (
    Path(__file__).parents[1]
    / "fixtures"
    / "phase5"
    / "F-REPEAT-AWARE-001"
)


def _fixture() -> tuple[dict, dict, dict]:
    manifest = json.loads((FIXTURE_ROOT / "fixture.json").read_text(encoding="utf-8"))
    input_doc = json.loads((FIXTURE_ROOT / manifest["input"]).read_text(encoding="utf-8"))
    expected = json.loads((FIXTURE_ROOT / manifest["expected"]).read_text(encoding="utf-8"))
    return manifest, input_doc, expected


def _op(row: list) -> RepeatAwareVariant:
    return RepeatAwareVariant(
        kind=row[0], ref_pos0=row[1], ref=row[2], alt=row[3]
    )


def _run_case(input_doc: dict, case_name: str):
    case = input_doc["cases"][case_name]
    locus_row = input_doc["locus"]
    locus = RepeatAwareLocusObservation(
        sample_id=locus_row["sample_id"],
        locus_id=locus_row["locus_id"],
        legacy_event_id=locus_row["legacy_event_id"],
        family=locus_row["family"],
        chrom=locus_row["chrom"],
        interval0=ClosedBreakpointInterval0(*locus_row["interval0"]),
        left_coord0=locus_row["left_coord0"],
        right_coord0=locus_row["right_coord0"],
    )
    groups = tuple(
        RepeatAwareGroupObservation(
            group_id=row["id"],
            locus_id=locus.locus_id,
            side=BreakpointSide(row["side"]),
            candidate_family="FamilyA",
            representative_pair_id=row["pair"],
            te_facing_length=row["length"],
            family_edit_costs=row["costs"],
            te_record_id=row.get("record"),
            alignment_nm=row.get("nm"),
            operations=tuple(_op(op) for op in row.get("ops", [])),
        )
        for row in case["groups"]
    )
    return run_repeat_aware_incumbent(
        RepeatAwareSampleContext(locus.sample_id, (locus,)),
        RepeatAwarePreparedAssets({"te1": case["te_sequence"]}),
        RepeatAwareObservations(groups),
    )


def test_fixture_manifest_binds_only_synthetic_case_bytes():
    manifest, input_doc, _expected = _fixture()
    assert manifest["fixture_id"] == "F-REPEAT-AWARE-001"
    assert manifest["synthetic"] is True
    assert "no truth" in manifest["provenance"].lower()
    assert input_doc["synthetic"] is True
    assert manifest["case_count"] == len(input_doc["cases"])
    assert set(manifest["sha256"]) == {manifest["input"], manifest["expected"]}
    for relative_path, expected_hash in manifest["sha256"].items():
        actual_hash = hashlib.sha256((FIXTURE_ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash


def test_alignment_cigar_regression_reconstructs_sub_ins_del_operations():
    _manifest, input_doc, expected = _fixture()
    case = input_doc["cigar_regression"]
    operations = alignment_variant_ops(
        case["query"],
        case["reference"],
        case["reference_start0"],
        case["cigar"],
    )
    assert [list(operation.exact_key) for operation in operations] == expected[
        "cigar_operations"
    ]


def test_exact_three_over_thirty_five_boundary_and_tied_best_family():
    _manifest, input_doc, expected = _fixture()
    boundary = expected["rate_boundary"]
    passing = RepeatAwareGroupObservation(
        "pass", "loc", BreakpointSide.LEFT, "FamilyA", 1, boundary["length"],
        {"FamilyA": boundary["pass_nm"], "FamilyB": boundary["pass_nm"] + 1},
    )
    failing = RepeatAwareGroupObservation(
        "fail", "loc", BreakpointSide.LEFT, "FamilyA", 2, boundary["length"],
        {"FamilyA": boundary["fail_nm"], "FamilyB": boundary["fail_nm"] + 1},
    )
    assert passing.baseline_pass is True
    assert failing.baseline_pass is False

    result = _run_case(input_doc, "best_family_tie")
    target = next(row for row in result.group_decisions if row.group_id == "target")
    expected_target = expected["decisions"]["best_family_tie"]
    assert target.status == expected_target["status"]
    assert target.eligible_rate_only_failure is expected_target[
        "eligible_rate_only_failure"
    ]
    assert target.crossvalidated_pass is expected_target["crossvalidated_pass"]
    assert target.final_pass is False


def test_leave_one_pair_out_requires_two_distinct_other_pair_ids():
    _manifest, input_doc, expected = _fixture()
    duplicate_result = _run_case(input_doc, "duplicate_pair_support")
    duplicate = next(
        row for row in duplicate_result.group_decisions if row.group_id == "target"
    )
    duplicate_expected = expected["decisions"]["duplicate_pair_support"]
    assert duplicate.status == duplicate_expected["status"]
    assert duplicate.adjusted_nm == duplicate_expected["adjusted_nm"]
    assert duplicate.crossvalidated_pass is duplicate_expected["crossvalidated_pass"]
    assert list(duplicate.training_pair_ids) == duplicate_expected["training_pair_ids"]

    same_pair_doc = json.loads(json.dumps(input_doc))
    same_pair_groups = same_pair_doc["cases"]["duplicate_pair_support"]["groups"]
    same_pair_groups[0]["pair"] = 1
    same_pair_groups[1]["pair"] = 1
    same_pair_groups[2]["pair"] = 20
    same_pair_result = _run_case(same_pair_doc, "duplicate_pair_support")
    same_pair_target = next(
        row for row in same_pair_result.group_decisions if row.group_id == "target"
    )
    assert same_pair_target.status == "NO_RECURRENT_VARIANT_IN_HELDOUT"
    assert same_pair_target.training_pair_ids == ()

    distinct_result = _run_case(input_doc, "distinct_pair_support")
    distinct = next(
        row for row in distinct_result.group_decisions if row.group_id == "target"
    )
    distinct_expected = expected["decisions"]["distinct_pair_support"]
    assert distinct.status == distinct_expected["status"]
    assert distinct.adjusted_nm == distinct_expected["adjusted_nm"]
    assert distinct.crossvalidated_pass is distinct_expected["crossvalidated_pass"]
    assert list(distinct.training_pair_ids) == distinct_expected["training_pair_ids"]


def test_exact_allele_support_tie_fails_closed():
    _manifest, input_doc, expected = _fixture()
    result = _run_case(input_doc, "exact_allele_support_tie")
    target = next(row for row in result.group_decisions if row.group_id == "target")
    expected_target = expected["decisions"]["exact_allele_support_tie"]
    assert target.status == expected_target["status"]
    assert target.adjusted_nm == expected_target["adjusted_nm"]
    assert target.crossvalidated_pass is expected_target["crossvalidated_pass"]
    assert target.credited_variant_labels == ()


def test_repeat_equivalent_indel_requires_linked_recurrent_substitution():
    _manifest, input_doc, expected = _fixture()
    unanchored_result = _run_case(input_doc, "repeat_unanchored_indel")
    unanchored = next(
        row for row in unanchored_result.group_decisions if row.group_id == "target"
    )
    unanchored_expected = expected["decisions"]["repeat_unanchored_indel"]
    assert unanchored.status == unanchored_expected["status"]
    assert unanchored.adjusted_nm == unanchored_expected["adjusted_nm"]
    assert list(unanchored.credited_variant_labels) == unanchored_expected[
        "credited_variant_labels"
    ]
    assert unanchored.crossvalidated_pass is unanchored_expected[
        "crossvalidated_pass"
    ]

    linked_result = _run_case(input_doc, "repeat_linked_substitution")
    linked = next(row for row in linked_result.group_decisions if row.group_id == "target")
    linked_expected = expected["decisions"]["repeat_linked_substitution"]
    assert linked.status == linked_expected["status"]
    assert linked.adjusted_nm == linked_expected["adjusted_nm"]
    assert set(linked.credited_variant_labels) == {
        "INS@3:>A",
        "SUB@1:A>C",
    }
    assert linked.crossvalidated_pass is linked_expected["crossvalidated_pass"]


def test_unique_indel_placement_is_credited():
    _manifest, input_doc, expected = _fixture()
    result = _run_case(input_doc, "unique_indel_placement")
    target = next(row for row in result.group_decisions if row.group_id == "target")
    expected_target = expected["decisions"]["unique_indel_placement"]
    assert target.status == expected_target["status"]
    assert target.adjusted_nm == expected_target["adjusted_nm"]
    assert target.credited_variant_labels == ("INS@3:>A",)
    assert target.crossvalidated_pass is expected_target["crossvalidated_pass"]


def test_event_requires_two_distinct_pairs_on_both_sides_and_preserves_closed_interval():
    _manifest, input_doc, expected = _fixture()
    case_expected = expected["events"]["bilateral_event"]
    result = _run_case(input_doc, "bilateral_event")
    assert result.locus_decisions[0].passes_bilateral_support is True
    assert result.locus_decisions[0].left_pair_ids == (1, 2)
    assert result.locus_decisions[0].right_pair_ids == (3, 4)
    assert len(result.events) == 1
    event = result.events[0]
    assert event.interval0 == ClosedBreakpointInterval0(*case_expected["interval0"])
    assert event.left_support_pairs == case_expected["left_support_pairs"]
    assert event.right_support_pairs == case_expected["right_support_pairs"]
    assert event.total_support_pairs == 4
    assert event.fragment_pair_ids == (1, 2, 3, 4)

    insufficient = _run_case(input_doc, "one_right_event")
    insufficient_expected = expected["events"]["one_right_event"]
    assert insufficient.locus_decisions[0].passes_bilateral_support is False
    assert len(insufficient.events) == 0
    assert insufficient_expected["emitted"] is False
