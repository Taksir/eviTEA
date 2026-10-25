"""Architecture identifier and scientific constants."""

from __future__ import annotations

from dataclasses import dataclass

from .errors import ArchitectureInvariantError

ARCHITECTURE_ID = "EVITEA_V1_REPEAT_AWARE_EVIDENCE_ARCHITECTURE_20260928_REVISION_1"
N_FOLDS = 5
Q99 = 805  # default fragment-size 99th percentile (bp). A run uses the sample's own value (``sample_geometry``)
MIN_PAIRS = 2
MIN_FOLDS = 2
EVENT_EQUIV_DISTANCE = 500
FROZEN_FAMILY_ALIASES = (("rooA", "roo"), ("Rt1a", "Rt1b"))
CERTIFYING_PAIRED_END_STATES = ("BILATERAL", "LEFT_ONLY", "RIGHT_ONLY")


def _require_nonnegative_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ArchitectureInvariantError(f"{name} must be a non-negative integer")


def repeat_aware_edit_rate_passes(adjusted_nm: int, te_facing_length: int) -> bool:
    """Return True when the edit rate is at most ``3/35``, using integer arithmetic."""

    _require_nonnegative_integer(adjusted_nm, "adjusted_nm")
    if (
        isinstance(te_facing_length, bool)
        or not isinstance(te_facing_length, int)
        or te_facing_length <= 0
    ):
        raise ArchitectureInvariantError("te_facing_length must be a positive integer")
    return 35 * adjusted_nm <= 3 * te_facing_length


@dataclass(frozen=True, slots=True)
class ScienceProfile:
    """Immutable scientific parameters tied to an architecture identifier."""

    architecture_id: str = ARCHITECTURE_ID
    n_folds: int = N_FOLDS
    q99: int = Q99
    min_pairs: int = MIN_PAIRS
    min_folds: int = MIN_FOLDS
    event_equiv_distance: int = EVENT_EQUIV_DISTANCE
    family_aliases: tuple[tuple[str, str], ...] = FROZEN_FAMILY_ALIASES
    certifying_paired_end_states: tuple[str, ...] = CERTIFYING_PAIRED_END_STATES

    def __post_init__(self) -> None:
        # q99 (the fragment-size 99th percentile) is not checked: each run estimates it from the sample (default 805)
        # and records the value used in the run manifest. See ``sample_geometry``.
        frozen_values = (
            self.n_folds,
            self.min_pairs,
            self.min_folds,
            self.event_equiv_distance,
            self.family_aliases,
            self.certifying_paired_end_states,
        )
        expected = (
            N_FOLDS,
            MIN_PAIRS,
            MIN_FOLDS,
            EVENT_EQUIV_DISTANCE,
            FROZEN_FAMILY_ALIASES,
            CERTIFYING_PAIRED_END_STATES,
        )
        if self.architecture_id == ARCHITECTURE_ID and frozen_values != expected:
            raise ArchitectureInvariantError(
                "scientific values other than the defaults require a different architecture identity"
            )


SCIENCE_PROFILE = ScienceProfile()
