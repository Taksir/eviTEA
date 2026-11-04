"""The floor-1 add-on of the library-free step (opt-in, ``--floor1-addon``; never run by default).

When the weak-site pass's allele depth is at or under its cut, a second library-free step runs with the extension floor at 1
(``<output>/floor1_addon/``), and the LIBRARY_FREE rows of its ``calls.combined.tsv`` that lie more than 500 bp from every row of the
main ``calls.combined.tsv`` (any source, same chromosome) are added to the main output.  The rule is the desk rule of
``scripts/union_desk.py`` (setting U1): the add-on rows are compared with the main rows only, not with each other.

Added rows are appended at the end of the main ``calls.combined.tsv`` and ``calls.library_free.tsv`` in the add-on file's order, with
``event_id`` = ``F1_`` + the add-on's id.  The records of ``calls.combined.vcf`` are merged in VCF order (a VCF is position-sorted;
the header lists the families of all records).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence

__all__ = ["ID_PREFIX", "NEAR_BP", "SOURCE", "RULE", "added_rows", "read_rows", "append_tsv_rows", "merge_vcf"]

ID_PREFIX = "F1_"
NEAR_BP = 500
SOURCE = "LIBRARY_FREE"
RULE = "add-on LIBRARY_FREE rows of floor 1 not within 500 bp of any main row"


def read_rows(data: bytes) -> list[dict[str, str]]:
    """Rows of a tab-separated table with a header line."""

    lines = data.decode("utf-8").splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:] if line.strip()]


def added_rows(main_rows: Iterable[Mapping[str, str]], addon_rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    """The add-on rows to add: ``source`` LIBRARY_FREE, and no main row of any source on the same ``chr`` with
    ``main.start - 500 <= addon.end and addon.start <= main.end + 500`` (a gap of exactly 500 is near, 501 is not).  The rows come
    back in the add-on's order with ``event_id`` prefixed ``F1_``."""

    have: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in main_rows:
        have[row["chr"]].append((int(row["start"]), int(row["end"])))
    out = []
    for row in addon_rows:
        if row.get("source") != SOURCE:
            continue
        start, end = int(row["start"]), int(row["end"])
        if any(a - NEAR_BP <= end and start <= b + NEAR_BP for a, b in have[row["chr"]]):
            continue
        out.append({**row, "event_id": ID_PREFIX + row["event_id"]})
    return out


def combined_lines(rows: Sequence[Mapping[str, str]]) -> list[str]:
    return [f"{r['event_id']}\t{r['chr']}\t{r['start']}\t{r['end']}\t{r['family']}\t{SOURCE}" for r in rows]


def append_tsv_rows(data: bytes, lines: Sequence[str]) -> bytes:
    """``data`` (text ending in a newline) with ``lines`` added at the end."""

    text = data.decode("utf-8")
    if text and not text.endswith("\n"):
        text += "\n"
    return (text + "".join(line + "\n" for line in lines)).encode("utf-8")


def library_free_lines(addon_lf_tsv: bytes, ids: Sequence[str]) -> tuple[list[str], str]:
    """(lines of the add-on's ``calls.library_free.tsv`` for the add-on event ids ``ids``, with the ``F1_`` prefix, in the
    add-on file's order; the header line)."""

    lines = addon_lf_tsv.decode("utf-8").splitlines()
    wanted = set(ids)
    out = []
    for line in lines[1:]:
        first = line.split("\t", 1)[0]
        if first in wanted:
            out.append(ID_PREFIX + line)
    return out, (lines[0] if lines else "")


def merge_vcf(existing: bytes, new: bytes) -> bytes:
    """One VCF from two: the header of ``existing`` with the ALT lines of the families of both, the records of both in the order
    ``vcf.calls_vcf_bytes`` writes them (contig order, position, id).  ``new`` was written with the same contigs."""

    def parts(data: bytes):
        header, records = [], []
        for line in data.decode("utf-8").splitlines():
            (header if line.startswith("#") else records).append(line)
        return header, records

    head_a, rec_a = parts(existing)
    head_b, rec_b = parts(new)
    contig_rank = {line.split("ID=", 1)[1].split(",", 1)[0]: n for n, line in enumerate(h for h in head_a if h.startswith("##contig="))}
    family_prefix = "##ALT=<ID=INS:ME:"
    families = {h for h in (*head_a, *head_b) if h.startswith(family_prefix)}
    out_header = []
    placed = False
    for line in head_a:
        if line.startswith(family_prefix):
            if not placed:
                out_header.extend(sorted(families, key=lambda h: h[len(family_prefix):].split(",", 1)[0]))
                placed = True
            continue
        out_header.append(line)
    if not placed:
        raise ValueError("the VCF has no ALT line for a family")

    def key(record: str):
        fields = record.split("\t")
        return contig_rank[fields[0]], int(fields[1]), fields[2]

    records = sorted([*rec_a, *rec_b], key=key)
    return ("\n".join([*out_header, *records]) + "\n").encode("utf-8")
