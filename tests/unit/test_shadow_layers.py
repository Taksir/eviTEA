"""Policies for appending the extension layers to the library-layer call set."""

from __future__ import annotations

import pytest

from evitea.errors import InputContractError
from evitea.shadow_layers import (
    AnchorPECandidate,
    ExtensionLayerInputs,
    FullMemoryEvent,
    Internal21Instance,
    Internal21Locus,
    ShadowCallset,
    ShadowLabels,
    append_asymmetric_layer,
    append_bilateral_layer,
    append_full_memory_layer,
    append_internal21_layer,
    compose_frozen_shadow,
    paired_end_state,
)

FIELDS = (
    "event_id\tchr\tstart\tend\tfamily\tfamilies\tleft_coord\tright_coord\tdelta\t"
    "left_support\tright_support\ttotal_support\tevidence_type\tscore\tfragment_ids"
)
LABELS = ShadowLabels(bilateral_prefix="", asymmetric="S", full_memory="S", internal21="S")


def base(*rows: str) -> ShadowCallset:
    return ShadowCallset.from_bytes(("\n".join((FIELDS, *rows)) + "\n").encode())


def r9e(event_id: str, chrom: str, start: int, end: int, family: str) -> str:
    return (
        f"{event_id}\t{chrom}\t{start}\t{end}\t{family}\t['{family}']\t-1\t{end}\t0\t0\t3\t3\t"
        f"one_sided_right\t3.0\t1,2,3"
    )


def cand(cluster: int, family: str, bp: int, left=(), right=(), side="LEFT", chrom="c1"):
    return AnchorPECandidate(cluster, family, chrom, bp, side, (10, 11), tuple(left), tuple(right))


BOTH = dict(left=(0, 1), right=(2, 3))


def test_paired_end_state_needs_two_pairs_across_two_folds():
    assert paired_end_state((0, 5), (1, 2)) == "RIGHT_ONLY"  # 0 and 5 share fold 0
    assert paired_end_state((0, 1), (2, 3)) == "BILATERAL"
    assert paired_end_state((), (1,)) == "NEITHER"


def test_channel_a_replication_is_required():
    with pytest.raises(InputContractError):
        AnchorPECandidate(1, "F", "c1", 100, "LEFT", (0, 5), (), ())


def test_bilateral_appends_after_exact_base_prefix_in_native_schema():
    start = base(r9e("E1", "c1", 5000, 5000, "F"))
    result = append_bilateral_layer(start, [cand(7, "F", 100, **BOTH)], LABELS)
    assert result.callset.data.startswith(start.data)
    new = result.callset.data[len(start.data):].decode()
    assert new == (
        "V3ANCHORPE_0007_001\tc1\t100\t100\tF\t['F']\t100\t-1\t0\t2\t2\t4\t"
        "v3_anchor_pe_bilateral\t4\t0,1,2,3,10,11\n"
    )


def test_bilateral_absorbs_same_family_within_500_and_ignores_other_families():
    start = base(r9e("E1", "c1", 600, 700, "F"), r9e("E2", "c1", 100, 100, "G"))
    result = append_bilateral_layer(
        start, [cand(1, "F", 100, **BOTH), cand(2, "F", 1201, **BOTH)], LABELS
    )
    statuses = {row.event_id: row.status for row in result.ledger}
    assert statuses["V3ANCHORPE_0001_001"] == "ABSORBED_EXISTING_SAME_FAMILY"  # 500 bp away
    assert statuses["V3ANCHORPE_0002_002"] == "ADDED"  # 501 bp away
    assert result.unresolved == ()


def test_asymmetric_admits_only_one_sided_candidates():
    start = base(r9e("E1", "c1", 9000, 9000, "F"))
    cands = [
        cand(1, "F", 100, **BOTH),
        cand(2, "F", 2000, left=(0, 1)),
        cand(3, "F", 4000),
    ]
    result = append_asymmetric_layer(start, cands, LABELS)
    assert result.appended_event_ids == ("V3ANCHORPE_ASYM_S_0002_001",)


def test_full_memory_uses_aliases_and_skips_neither():
    start = base(r9e("E1", "c1", 100, 100, "roo"))
    events = [
        FullMemoryEvent("S_FULL_000001", "rooA", "c1", 300, "RIGHT", (1, 2), (3, 4), (), "LEFT_ONLY"),
        FullMemoryEvent("S_FULL_000002", "F", "c1", 5000, "RIGHT", (1, 2), (), (), "NEITHER"),
    ]
    result = append_full_memory_layer(start, events, LABELS, sample_id="S")
    assert [r.status for r in result.ledger] == ["ABSORBED_EXISTING_SAME_FAMILY"]
    assert result.appended_event_ids == ()


def test_full_memory_new_same_family_multiplicity_fails_closed():
    start = base(r9e("E1", "c9", 100, 100, "roo"))
    events = [
        FullMemoryEvent("S_FULL_000001", "F", "c1", 1000, "LEFT", (1, 2), (3, 4), (), "LEFT_ONLY"),
        FullMemoryEvent("S_FULL_000002", "F", "c1", 1400, "LEFT", (5, 6), (7, 8), (), "LEFT_ONLY"),
        FullMemoryEvent("S_FULL_000003", "G", "c1", 1200, "LEFT", (5, 6), (7, 8), (), "LEFT_ONLY"),
    ]
    result = append_full_memory_layer(start, events, LABELS, sample_id="S")
    assert result.appended_event_ids == ("V3FULLMEM_AP_S_000003",)
    assert {u.reason for u in result.unresolved} == {"UNRESOLVED_NEW_SAME_FAMILY_MULTIPLICITY"}


def inst(rotation, bp, side="LEFT"):
    return Internal21Instance(rotation, bp, side, (1, 2), (3, 4), ())


def test_internal21_representative_rule():
    locus = Internal21Locus("F", "c1", 100, 120, (inst("H12", 110, "RIGHT"), inst("H01", 105), inst("H02", 105)))
    chosen, recurrence = locus.representative()
    assert (chosen.breakpoint0, chosen.rotation, recurrence) == (105, "H01", 2)


def test_internal21_different_family_nearby_is_unresolved_and_new_ids_are_sequential():
    start = base(r9e("E1", "c1", 1000, 1000, "G"))
    loci = [
        Internal21Locus("F", "c1", 1400, 1450, (inst("H01", 1400),)),
        Internal21Locus("F", "c2", 10, 20, (inst("H01", 15),)),
        Internal21Locus("F", "c3", 10, 20, (inst("H01", 15),)),
    ]
    result = append_internal21_layer(start, loci, LABELS)
    assert [u.reason for u in result.unresolved] == ["UNRESOLVED_NEARBY_DIFFERENT_FAMILY_CALL"]
    assert result.appended_event_ids == ("V3INTA_AP_S_000001", "V3INTA_AP_S_000002")


def test_internal21_new_cross_family_competition_fails_closed():
    start = base(r9e("E1", "c9", 1000, 1000, "G"))
    loci = [
        Internal21Locus("F", "c1", 100, 100, (inst("H01", 100),)),
        Internal21Locus("H", "c1", 550, 560, (inst("H01", 555),)),
    ]
    result = append_internal21_layer(start, loci, LABELS)
    assert result.appended_event_ids == ()
    assert {u.reason for u in result.unresolved} == {"UNRESOLVED_NEW_CROSS_FAMILY_COMPETITION"}


def test_compose_runs_layers_in_frozen_order():
    start = base(r9e("E1", "c9", 1000, 1000, "G"))
    inputs = ExtensionLayerInputs(
        bilateral_candidates=(cand(1, "F", 100, **BOTH), cand(2, "F", 100, left=(0, 1), chrom="c2")),
        full_memory_events=(),
        internal21_loci=(),
    )
    layers = compose_frozen_shadow(start.data, inputs, LABELS, sample_id="S")
    assert [layer.layer for layer in layers] == ["BILATERAL", "ASYMMETRIC", "FULL_MEMORY", "INTERNAL21"]
    assert layers[-1].callset.row_count == 3
