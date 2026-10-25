"""TE-family labels and family-map validation."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from types import MappingProxyType

from .constants import FROZEN_FAMILY_ALIASES
from .errors import InputContractError

_ALIASES = MappingProxyType(dict(FROZEN_FAMILY_ALIASES))


def canonical_family(label: str) -> str:
    """Return the family label after applying the two built-in aliases only."""

    if not isinstance(label, str):
        raise InputContractError("family label must be text")
    normalized = label.strip()
    if not normalized:
        raise InputContractError("family label must be non-empty")
    return _ALIASES.get(normalized, normalized)


@dataclass(frozen=True, slots=True)
class FamilyMap(Mapping[str, str]):
    """Immutable mapping from TE record IDs to families (aliases applied)."""

    assignments: Mapping[str, str]

    def __post_init__(self) -> None:
        canonical: dict[str, str] = {}
        for record_id, family in self.assignments.items():
            if not isinstance(record_id, str) or not record_id.strip():
                raise InputContractError("TE record IDs must be non-empty text")
            if record_id != record_id.strip():
                raise InputContractError("TE record IDs may not have surrounding whitespace")
            canonical[record_id] = canonical_family(family)
        object.__setattr__(
            self, "assignments", MappingProxyType(dict(sorted(canonical.items())))
        )

    def __getitem__(self, key: str) -> str:
        return self.assignments[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.assignments)

    def __len__(self) -> int:
        return len(self.assignments)


def validate_family_map(
    records: Mapping[str, str] | Iterable[tuple[str, str]],
    te_record_ids: Iterable[str],
    *,
    allow_unassigned: bool = False,
) -> FamilyMap:
    """Validate that every TE record ID has exactly one family (aliases applied).

    ``allow_unassigned`` accepts a species-wide map (for example the FlyBase
    ``FBte_ID -> FlyBase_Family`` table).  IDs absent from the TE library are
    dropped, and library records without a family stay unassigned.
    """

    rows = records.items() if isinstance(records, Mapping) else records
    assignments: dict[str, str] = {}
    for row in rows:
        try:
            record_id, family = row
        except (TypeError, ValueError) as exc:
            raise InputContractError("family-map rows must contain exactly two fields") from exc
        canonical = canonical_family(family)
        previous = assignments.get(record_id)
        if previous is not None and previous != canonical:
            raise InputContractError(
                f"conflicting family assignments for TE record {record_id!r}"
            )
        assignments[record_id] = canonical

    expected = tuple(te_record_ids)
    if len(expected) != len(set(expected)):
        raise InputContractError("TE record IDs must be unique")
    expected_set = set(expected)
    observed_set = set(assignments)
    unknown = sorted(observed_set - expected_set)
    missing = sorted(expected_set - observed_set)
    if allow_unassigned:
        return FamilyMap({key: value for key, value in assignments.items() if key in expected_set})
    if unknown:
        raise InputContractError(f"family map contains unknown TE IDs: {unknown}")
    if missing:
        raise InputContractError(f"family map is missing TE IDs: {missing}")
    return FamilyMap(assignments)
