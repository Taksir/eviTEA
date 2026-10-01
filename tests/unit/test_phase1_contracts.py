"""Phase 1 checks for coordinate and stable-serialization contracts."""

import json

import pytest

from evitea.contracts import (
    BaseSpan0,
    PairedEndSupportState,
    ClosedBreakpointInterval0,
    InterbasePoint0,
    PairIdentityMode,
    ReadMate,
    stable_json_dumps,
    stable_to_primitive,
)
from evitea.errors import ArchitectureInvariantError


@pytest.mark.parametrize(
    "factory, args",
    [
        (BaseSpan0, (-1, 0)),
        (BaseSpan0, (3, 2)),
        (InterbasePoint0, (-1,)),
        (ClosedBreakpointInterval0, (-1, 0)),
        (ClosedBreakpointInterval0, (5, 4)),
    ],
)
def test_coordinate_types_reject_negative_or_inverted_values(factory, args):
    with pytest.raises(ArchitectureInvariantError):
        factory(*args)


def test_interbase_endpoint_is_valid_only_within_contig_bounds():
    endpoint = InterbasePoint0(10)
    assert endpoint.require_within(10) is endpoint
    with pytest.raises(ArchitectureInvariantError):
        InterbasePoint0(11).require_within(10)


def test_closed_interval_distance_uses_inclusive_endpoints():
    left = ClosedBreakpointInterval0(2, 5)
    right = ClosedBreakpointInterval0(9, 12)
    touching = ClosedBreakpointInterval0(5, 8)

    assert left.distance_to(right) == 4
    assert left.distance_to(touching) == 0
    assert ClosedBreakpointInterval0(7, 7).distance_to(InterbasePoint0(10)) == 3


def test_frozen_core_state_enums_have_only_documented_values():
    assert {state.value for state in PairedEndSupportState} == {
        "BILATERAL",
        "LEFT_ONLY",
        "RIGHT_ONLY",
        "NEITHER",
    }
    assert {mode.value for mode in PairIdentityMode} == {
        "EVITEA_ORDINAL0",
        "LEGACY_FROZEN",
    }
    assert {mate.name: mate.value for mate in ReadMate} == {"READ1": 1, "READ2": 2}


def test_stable_serialization_is_independent_of_mapping_insertion_order():
    first = {
        "z_state": PairedEndSupportState.NEITHER,
        "a_span": BaseSpan0(1, 4),
        "folds": {4, 1, 3},
    }
    second = {
        "folds": {3, 4, 1},
        "a_span": BaseSpan0(1, 4),
        "z_state": PairedEndSupportState.NEITHER,
    }

    primitive = stable_to_primitive(first)
    encoded = stable_json_dumps(first)

    assert primitive == stable_to_primitive(second)
    assert encoded == stable_json_dumps(second)
    assert json.loads(encoded) == primitive
    assert encoded.index('"a_span"') < encoded.index('"folds"') < encoded.index('"z_state"')
