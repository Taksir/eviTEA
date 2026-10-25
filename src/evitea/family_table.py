"""The TE family map: one reader for the two accepted formats.

1. The FlyBase table: a header row with ``FBte_ID`` and ``FlyBase_Family`` columns (other columns are ignored).  It is read
   as before, so FlyBase ids (``FBte`` + digits) keep their exact parsing.
2. A plain two-column table, ``record<TAB>family``, with any record names (``B1_Mm``, ``L1MdA_I``, ``Ty1``, ...).  An optional
   first row naming the columns (``te_record_id<TAB>family``, ``record<TAB>family``, ``id<TAB>family``, ...) is skipped.

Blank lines and lines starting with ``#`` are ignored in both.  Record names in the plain table are taken as written; they must
match the first word of the TE library's FASTA headers.  Labels are returned as written (no alias or case folding); each
reader applies its own normalisation.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

from .errors import InputContractError

__all__ = [
    "FLYBASE_COLUMNS",
    "TWO_COLUMN",
    "FLYBASE",
    "OTHER",
    "family_pairs",
    "family_table_kind",
    "is_two_column_table",
    "read_rows",
]

FLYBASE_COLUMNS = ("FBte_ID", "FlyBase_Family")
FLYBASE, TWO_COLUMN, OTHER = "flybase", "two-column", "other"
_HEADER_IDS = frozenset({"te_record_id", "record", "record_id", "te_id", "id", "te", "fbte_id", "name"})
_HEADER_FAMILIES = frozenset({"family", "te_family", "flybase_family"})


def read_rows(path: str | Path) -> list[list[str]]:
    """Tab-separated rows, without blank lines and ``#`` comments."""

    path = Path(path)
    try:
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise InputContractError(f"cannot read family map {path}: {exc}") from exc
    return [row for row in rows if any(field.strip() for field in row) and not row[0].startswith("#")]


def _kind(rows: list[list[str]]) -> str:
    if not rows:
        return OTHER
    if all(name in rows[0] for name in FLYBASE_COLUMNS):
        return FLYBASE
    return TWO_COLUMN if len(rows[0]) == 2 else OTHER


def family_table_kind(path: str | Path) -> str:
    """``flybase``, ``two-column`` (``record<TAB>family``) or ``other`` (any other layout, read by each reader as before)."""

    return _kind(read_rows(path))


def is_two_column_table(path: str | Path) -> bool:
    return family_table_kind(path) == TWO_COLUMN


def _is_header(row: list[str]) -> bool:
    return row[0].strip().lower() in _HEADER_IDS and row[1].strip().lower() in _HEADER_FAMILIES


def family_pairs(path: str | Path) -> list[tuple[str, str]]:
    """``(record, family)`` pairs in file order, both non-empty and stripped, from either format.

    A table in neither format is refused, as is a table that yields no pair.
    """

    rows = read_rows(path)
    kind = _kind(rows)
    if kind == FLYBASE:
        i, j = rows[0].index(FLYBASE_COLUMNS[0]), rows[0].index(FLYBASE_COLUMNS[1])
        pairs = [(r[i].strip(), r[j].strip()) for r in rows[1:] if len(r) > max(i, j)]
    elif kind == TWO_COLUMN:
        body = rows[1:] if _is_header(rows[0]) else rows
        for number, row in enumerate(body, 1):
            if len(row) != 2:
                raise InputContractError(f"family map {path}: row {number} must have two tab-separated fields: {row!r}")
        pairs = [(r[0].strip(), r[1].strip()) for r in body]
    else:
        raise InputContractError(
            f"family map {path} is neither a FlyBase table ({FLYBASE_COLUMNS[0]}, {FLYBASE_COLUMNS[1]} columns) nor a "
            "two-column record<TAB>family table"
        )
    pairs = [(record, family) for record, family in pairs if record and family]
    if not pairs:
        raise InputContractError(f"family map {path} is empty")
    return pairs
