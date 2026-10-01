"""Shared helpers for the full-memory branch."""

from __future__ import annotations

from ..constants import FROZEN_FAMILY_ALIASES, N_FOLDS

PRIMARY_CHROMS = ("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
ROTATIONS = tuple((a, b) for a in range(N_FOLDS) for b in range(a + 1, N_FOLDS))
_ALIASES = dict(FROZEN_FAMILY_ALIASES)


def normfam(name: object) -> str:
    text = str(name or "").strip()
    return _ALIASES.get(text, text)


def as_int(value: object, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def rotation_label(heldout: tuple[int, int]) -> str:
    return f"H{heldout[0]}{heldout[1]}"
