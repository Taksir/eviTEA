"""Read length of the sample, resolved once per run.

The library layer measures TE coverage in read-length windows.  The length comes, in this order, from

1. ``--read-length`` (``runtime.read_length``), when given;
2. the prepared sample (``backend_profile.uniform_length``, recorded by ``prepare-sample``);
3. the msBWT index itself: a uniform-length index of ``n`` reads holds ``n`` end markers and ``n * (L + 1)``
   symbols (``totalCounts.npy``), so ``L = total / n - 1``.

``DEFAULT_READ_LENGTH`` (125 bp) is the length of the development samples; it is used only when none of the above
is available (for example by direct callers of the library functions).

"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import InputContractError

__all__ = [
    "DEFAULT_READ_LENGTH",
    "MIN_READ_LENGTH",
    "ReadLength",
    "read_length_from_index",
    "resolve_read_length",
    "validate_read_length",
]

DEFAULT_READ_LENGTH = 125
# Shortest read the window geometry accepts: the library layer's longest exact word is 21 bp.
MIN_READ_LENGTH = 21

@dataclass(frozen=True, slots=True)
class ReadLength:
    length: int
    source: str  # "option", "prepared-sample", "msbwt-index" or "default"


def validate_read_length(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < MIN_READ_LENGTH:
        raise InputContractError(f"read length must be an integer of at least {MIN_READ_LENGTH} bp, got {value!r}")
    return value


def read_length_from_index(index_dir: str | Path) -> int:
    """Read length of a uniform-length msBWT index, from its symbol counts."""

    import numpy as np

    path = Path(index_dir) / "totalCounts.npy"
    try:
        counts = np.load(path, allow_pickle=False)
    except (OSError, ValueError) as exc:
        raise InputContractError(f"cannot read the msBWT symbol counts: {path}") from exc
    reads, total = int(counts[0]), int(counts.sum())
    if reads <= 0 or total % reads:
        raise InputContractError(
            f"the msBWT index {index_dir} does not hold reads of one length; pass --read-length"
        )
    return validate_read_length(total // reads - 1)


def resolve_read_length(
    option: int | None,
    backend_profile: Mapping[str, Any] | None,
    msbwt_index: str | Path | None,
) -> ReadLength:
    if option is not None:
        return ReadLength(validate_read_length(option), "option")
    if backend_profile is not None and backend_profile.get("uniform_length") is not None:
        return ReadLength(validate_read_length(int(backend_profile["uniform_length"])), "prepared-sample")
    if msbwt_index is not None and (Path(msbwt_index) / "totalCounts.npy").is_file():
        return ReadLength(read_length_from_index(msbwt_index), "msbwt-index")
    return ReadLength(DEFAULT_READ_LENGTH, "default")

