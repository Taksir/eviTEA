"""Own-family resident-copy filter (``evitea filter-resident-copies``): a post-filter on a finished call table.

A core call is dropped when a reference-genome copy of the call's own TE family touches or overlaps it.  It is a
separate subcommand over finished tables: the ``call`` command and
the library-free layer never run it, so a run's bytes cannot change because this module exists.

Rule:

* the copy is a RepeatMasker ``.out`` line on the same chromosome (chromosome names normalized as ``norm_chrom``:
  a bare name such as ``2L`` or ``1`` gets ``chr`` when the reference has that ``chr`` contig), copy begin taken 0-based (``int(f[5]) - 1``), copy end
  ``int(f[6])``, copy id the record id (FBte id) in column 10, divergence column 2; lines with fewer than 15 fields, or whose first
  field is not a number (the header), are skipped;
* gap = max(0, copy_begin - call_end, call_start - copy_end); the call is dropped iff a copy of its own family has gap 0
  (touching counts; no strand, similarity or distance condition);
* the call's own family is its ``family`` column.  A family that is itself a record id (an FBte id, or a record of the map)
  stands for that id; a family name stands for every record id the family map lists under it (HeT-A and a few others have several); a family with neither is
  never dropped and is counted as ``unmapped_family``;
* only core calls are judged: an ``event_id`` starting ``LF_`` (library-free layer) passes through unchanged;
* kept rows are written byte for byte, in input order; dropped calls go to a side table with the copy that decided it
  (the smallest ``(begin, end, id, divergence)`` among the matching copies).
"""

from __future__ import annotations

import bisect
import csv
import hashlib
import io
from collections import defaultdict
from collections.abc import Collection, Iterable, Mapping
from dataclasses import astuple, dataclass, fields
from pathlib import Path

from .contigs import match_contig_name
from .errors import InputContractError
from .family_table import family_pairs

__all__ = [
    "CopyIndex",
    "DROPPED_COLUMNS",
    "DroppedRow",
    "FamilyIds",
    "FilterResult",
    "LISTED_COPY_TABLES",
    "dropped_tsv_bytes",
    "filter_calls",
    "load_copy_table",
    "load_family_ids",
    "parse_copy_table",
    "verify_copy_table",
]

LIBRARY_FREE_PREFIX = "LF_"
REQUIRED_COLUMNS = ("event_id", "chr", "start", "end", "family")

# sha256 of the reference copy tables this filter may read without ``--allow-unlisted-copies``
LISTED_COPY_TABLES = {
    "b141b76497d5446e30e3d9a67da68db08983c9ec6e2aedddf0c15cff9008516d":
        "dm6_TEMP2_compatible.fa.out (RepeatMasker, 43,971 lines; the table of the wrong-call anatomy)",
}


def norm_chrom(value: str, contigs: Collection[str]) -> str:
    """Return ``chr`` + name for a bare name when that contig exists in the reference (``contigs``); see ``evitea.contigs``."""

    return match_contig_name(str(value).strip(), contigs)


class CopyIndex:
    """Reference copies per chromosome, sorted by ``(begin, end, id, divergence)``.

    ``contigs``: the reference's contig names, against which bare copy-table names are matched (``norm_chrom``).
    """

    def __init__(self, records: Iterable[tuple[str, int, int, str, str]], contigs: Collection[str]):
        self.contigs = frozenset(contigs)
        by_chrom: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
        for chrom, begin, end, fbte, divergence in records:
            by_chrom[norm_chrom(chrom, self.contigs)].append((int(begin), int(end), fbte, divergence))
        self.by_chrom = {chrom: sorted(v) for chrom, v in by_chrom.items()}
        self.starts = {chrom: [c[0] for c in v] for chrom, v in self.by_chrom.items()}
        self.max_len = {chrom: max(c[1] - c[0] for c in v) for chrom, v in self.by_chrom.items()}

    def __len__(self) -> int:
        return sum(len(v) for v in self.by_chrom.values())

    def match(self, chrom: str, start: int, end: int, fbte_ids: frozenset[str]) -> tuple[str, int, int, str, str] | None:
        """Return the first copy (sorted order) of any id in ``fbte_ids`` at gap 0 from ``[start, end]``, or None."""

        chrom = norm_chrom(chrom, self.contigs)
        copies = self.by_chrom.get(chrom)
        if not copies:
            return None
        starts = self.starts[chrom]
        lo = bisect.bisect_left(starts, start - self.max_len[chrom])
        hi = bisect.bisect_right(starts, end)  # begin <= call end
        for begin, c_end, fbte, divergence in copies[lo:hi]:
            if c_end >= start and fbte in fbte_ids:  # gap 0: begin <= end and call_start <= copy_end
                return (chrom, begin, c_end, fbte, divergence)
        return None


def parse_copy_table(lines: Iterable[str]) -> list[tuple[str, int, int, str, str]]:
    """Return (chrom, 0-based begin, end, FBte id, divergence) for each RepeatMasker ``.out`` line."""

    records = []
    for line in lines:
        f = line.split()
        if len(f) < 15 or not f[0].isdigit():
            continue
        records.append((f[4], int(f[5]) - 1, int(f[6]), f[9], f[1]))
    return records


def load_copy_table(path: str | Path) -> list[tuple[str, int, int, str, str]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return parse_copy_table(handle)


def verify_copy_table(path: str | Path, *, allow_unlisted: bool) -> tuple[str, bool]:
    """Return (sha256, listed).  A table whose hash is not in ``LISTED_COPY_TABLES`` is refused unless allowed."""

    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    listed = digest in LISTED_COPY_TABLES
    if not listed and not allow_unlisted:
        raise InputContractError(
            f"reference TE copy table {path} (sha256 {digest}) is not a listed table; "
            "pass --allow-unlisted-copies to use it"
        )
    return digest, listed


class FamilyIds:
    """Map a family name to the record ids (FBte ids, or the record names of a two-column map) that stand for it."""

    def __init__(self, ids_by_name: Mapping[str, Iterable[str]]):
        self._ids = {name: frozenset(ids) for name, ids in ids_by_name.items()}
        self._records = frozenset(record for ids in self._ids.values() for record in ids)

    def of(self, family: str) -> frozenset[str] | None:
        if family.startswith("FBte") or family in self._records:
            return frozenset({family})
        return self._ids.get(family) or None


def load_family_ids(path: str | Path) -> FamilyIds:
    """Read the family map: the FlyBase table (``FBte_ID``, ``FlyBase_Family``) or two columns ``record<TAB>family``."""

    ids: dict[str, set[str]] = defaultdict(set)
    for record, name in family_pairs(path):
        ids[name].add(record)
    return FamilyIds(ids)


@dataclass(frozen=True, slots=True)
class DroppedRow:
    event_id: str
    chrom: str
    start: str
    end: str
    family: str
    family_fbte: str
    copy_chr: str
    copy_begin: int
    copy_end: int
    copy_divergence: str


_COLUMN_NAMES = {"chrom": "chr"}
DROPPED_COLUMNS = tuple(_COLUMN_NAMES.get(f.name, f.name) for f in fields(DroppedRow))


@dataclass(frozen=True, slots=True)
class FilterResult:
    filtered: bytes
    dropped: tuple[DroppedRow, ...]
    kept: int
    unmapped_family: int
    library_free_passed: int


def _integer(value: str, name: str, event_id: str) -> int:
    try:
        return int(value.strip())
    except ValueError:
        raise InputContractError(f"call {event_id}: {name} {value!r} is not an integer") from None


def call_contigs(calls: bytes) -> frozenset[str]:
    """Return the contig names of a call table's ``chr`` column (reference names), used when the reference is not given."""

    lines = calls.decode("utf-8-sig").splitlines()
    if not lines:
        return frozenset()
    header = lines[0].split("\t")
    if "chr" not in header:
        return frozenset()
    column = header.index("chr")
    return frozenset(f[column] for f in (line.split("\t") for line in lines[1:]) if len(f) > column and f[column])


def filter_calls(calls: bytes, index: CopyIndex, family_ids: FamilyIds) -> FilterResult:
    """Drop the core calls inside a copy of their own family; every other row is kept byte for byte."""

    if not calls:
        raise InputContractError("call table is empty")
    lines = calls.decode("utf-8-sig").splitlines(keepends=True)
    if calls.startswith(b"\xef\xbb\xbf"):
        raise InputContractError("call table starts with a byte-order mark")
    header = lines[0].rstrip("\r\n").split("\t")
    missing = [name for name in REQUIRED_COLUMNS if name not in header]
    if missing:
        raise InputContractError(f"call table lacks {missing}")
    col = {name: header.index(name) for name in REQUIRED_COLUMNS}
    kept_lines = [lines[0]]
    dropped: list[DroppedRow] = []
    unmapped = passed = kept = 0
    for line in lines[1:]:
        text = line.rstrip("\r\n")
        if not text.strip():
            kept_lines.append(line)
            continue
        f = text.split("\t")
        if len(f) <= max(col.values()):
            raise InputContractError(f"call table row has {len(f)} fields, fewer than its header: {text[:80]!r}")
        event_id = f[col["event_id"]]
        if event_id.startswith(LIBRARY_FREE_PREFIX):
            passed += 1
            kept += 1
            kept_lines.append(line)
            continue
        start = _integer(f[col["start"]], "start", event_id)
        end = _integer(f[col["end"]], "end", event_id)
        ids = family_ids.of(f[col["family"]])
        if ids is None:
            unmapped += 1
            kept += 1
            kept_lines.append(line)
            continue
        hit = index.match(f[col["chr"]], start, end, ids)
        if hit is None:
            kept += 1
            kept_lines.append(line)
            continue
        chrom, begin, c_end, fbte, divergence = hit
        dropped.append(DroppedRow(event_id, f[col["chr"]], f[col["start"]], f[col["end"]], f[col["family"]], fbte, chrom,
                                  begin, c_end, divergence))
    return FilterResult("".join(kept_lines).encode("utf-8"), tuple(dropped), kept, unmapped, passed)


def dropped_tsv_bytes(rows: Iterable[DroppedRow]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(DROPPED_COLUMNS)
    for row in rows:
        writer.writerow(astuple(row))
    return handle.getvalue().encode("utf-8")
