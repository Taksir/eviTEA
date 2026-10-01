"""Phase 1 checks for frozen scientific constants and family aliases."""

from dataclasses import FrozenInstanceError

import pytest

from evitea.constants import (
    ARCHITECTURE_ID,
    EVENT_EQUIV_DISTANCE,
    MIN_FOLDS,
    MIN_PAIRS,
    N_FOLDS,
    Q99,
    SCIENCE_PROFILE,
    repeat_aware_edit_rate_passes,
)
from evitea.errors import InputContractError
from evitea.families import canonical_family, validate_family_map


def test_frozen_science_constants_and_profile_are_present():
    assert ARCHITECTURE_ID
    assert (N_FOLDS, Q99, MIN_PAIRS, MIN_FOLDS, EVENT_EQUIV_DISTANCE) == (
        5,
        805,
        2,
        2,
        500,
    )
    assert SCIENCE_PROFILE is not None
    assert SCIENCE_PROFILE.n_folds == N_FOLDS
    assert SCIENCE_PROFILE.q99 == Q99
    assert SCIENCE_PROFILE.min_pairs == MIN_PAIRS
    assert SCIENCE_PROFILE.min_folds == MIN_FOLDS
    assert SCIENCE_PROFILE.event_equiv_distance == EVENT_EQUIV_DISTANCE
    assert SCIENCE_PROFILE.family_aliases == (("rooA", "roo"), ("Rt1a", "Rt1b"))
    assert SCIENCE_PROFILE.certifying_paired_end_states == (
        "BILATERAL",
        "LEFT_ONLY",
        "RIGHT_ONLY",
    )
    with pytest.raises(FrozenInstanceError):
        SCIENCE_PROFILE.q99 = 806


def test_repeat_aware_edit_rate_uses_exact_integer_boundary():
    assert repeat_aware_edit_rate_passes(3, 35)
    assert repeat_aware_edit_rate_passes(6, 70)
    assert not repeat_aware_edit_rate_passes(4, 35)


def test_only_the_two_frozen_family_aliases_are_canonicalized():
    assert canonical_family("rooA") == "roo"
    assert canonical_family("Rt1a") == "Rt1b"
    assert canonical_family("roo") == "roo"
    assert canonical_family("Rt1b") == "Rt1b"
    assert canonical_family("P7_candidate") == "P7_candidate"


def test_family_map_canonicalizes_aliases_and_requires_one_known_row_per_record():
    resolved = validate_family_map(
        [("te_roo", "rooA"), ("te_rt1", "Rt1a")],
        {"te_roo", "te_rt1"},
    )

    assert dict(resolved) == {"te_roo": "roo", "te_rt1": "Rt1b"}


@pytest.mark.parametrize(
    "records, te_ids",
    [
        ([("te1", "roo"), ("te1", "Rt1b")], {"te1"}),
        ([("te1", "roo")], {"te1", "te2"}),
        ([("extra", "roo")], {"te1"}),
    ],
)
def test_family_map_rejects_conflicts_missing_or_unknown_te_ids(records, te_ids):
    with pytest.raises(InputContractError):
        validate_family_map(records, te_ids)
