"""Immutable physical-pair identity assignment."""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from .constants import N_FOLDS
from .contracts import PairIdentityMode, PhysicalPair, RawReadPair
from .errors import InputContractError


def _validate_pair_id(pair_id: int) -> None:
    if isinstance(pair_id, bool) or not isinstance(pair_id, int) or pair_id < 0:
        raise InputContractError("pair_id must be a non-negative integer")


def fold_id(pair_id: int) -> int:
    """Return the five-fold label for a physical-pair ID."""

    _validate_pair_id(pair_id)
    return pair_id % N_FOLDS


def assign_pair_identity(
    raw_pair: RawReadPair,
    ordinal0: int,
    sample_id: str = "sample",
) -> PhysicalPair:
    """Assign the ``pair_id = synchronized ordinal0`` identity."""

    _validate_pair_id(ordinal0)
    if not isinstance(sample_id, str) or not sample_id:
        raise InputContractError("sample_id must be non-empty text")
    return PhysicalPair(
        sample_id=sample_id,
        raw_pair_ordinal0=ordinal0,
        pair_id=ordinal0,
        fold_id=fold_id(ordinal0),
        pair_name=raw_pair.pair_name,
        read1=raw_pair.read1,
        read2=raw_pair.read2,
        identity_mode=PairIdentityMode.EVITEA_ORDINAL0,
    )


def iter_physical_pairs(
    raw_pairs: Iterable[RawReadPair],
    *,
    sample_id: str = "sample",
) -> Iterator[PhysicalPair]:
    """Assign monotonically increasing ordinal-zero identities to a stream."""

    for ordinal0, raw_pair in enumerate(raw_pairs):
        yield assign_pair_identity(raw_pair, ordinal0, sample_id)
