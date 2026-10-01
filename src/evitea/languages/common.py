"""Shared primitives of the v3 language builders (ported, semantics frozen).

Ported from the frozen Route-B v3 research builders
(``route_b_v3_rebuild_safe_terminal21_codebook.py`` and
``route_b_v3_fill_underpowered_terminal_codebook.py``).  Where the two research
scripts carried slightly different helpers, both variants are kept and named
after the stage that used them; outputs must stay byte-identical to the
frozen tables.

Exact k-mer scanning is vectorized with NumPy but emits exactly the research
scanner's hits, in the same order: every window of ``k`` A/C/G/T bases whose
2-bit code is a query, by increasing start (per ``k`` in query-map order).
"""

from __future__ import annotations

import bisect
import csv
import io
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "build_interval_index",
    "canonical_family",
    "canonical_seq",
    "encode",
    "extract_fbte_rebuild",
    "extract_fbte_gapfill",
    "extract_te_id_rebuild",
    "extract_te_id_gapfill",
    "load_repeatmasker_rebuild",
    "load_repeatmasker_gapfill",
    "oriented_code_map",
    "overlap_query",
    "parse_family_map_rebuild",
    "parse_family_map_gapfill",
    "parse_fasta",
    "read_tsv_bytes",
    "reference_hit_class",
    "revcomp",
    "scan_exact",
    "tsv_bytes",
]

_TRANS = str.maketrans("ACGT", "TGCA")
_LUT = np.full(256, 255, dtype=np.uint8)
for _c, _v in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("a", 0), ("c", 1), ("g", 2), ("t", 3)):
    _LUT[ord(_c)] = _v
_BLOCK = 1 << 22


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """``(header, uppercase sequence)``; header is the line after '>' stripped."""

    name = None
    chunks: list[str] = []
    with Path(path).open("rt", encoding="ascii", errors="strict") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def revcomp(seq: str) -> str:
    return seq.translate(_TRANS)[::-1]


def canonical_seq(seq: str) -> str:
    return min(seq, revcomp(seq))


def encode(seq: str) -> int:
    value = 0
    for character in seq.encode("ascii"):
        code = int(_LUT[character])
        if code == 255:
            raise ValueError(seq)
        value = (value << 2) | code
    return value


def oriented_code_map(canonical_probes: Iterable[str]) -> dict[int, str]:
    """Both orientations of every probe -> the canonical probe (collision-checked)."""

    qmap: dict[int, str] = {}
    for probe in canonical_probes:
        for query in (probe, revcomp(probe)):
            code = encode(query)
            old = qmap.get(code)
            if old is not None and old != probe:
                raise RuntimeError("2-bit query collision")
            qmap[code] = probe
    return qmap


def _window_codes(values: np.ndarray, k: int, start: int, stop: int) -> np.ndarray:
    """2-bit codes of the last ``min(k, 32)`` bases of windows starting in [start, stop)."""

    width = min(k, 32)
    offset = k - width
    codes = np.zeros(stop - start, dtype=np.uint64)
    for j in range(width):
        codes <<= np.uint64(2)
        codes |= values[start + offset + j: stop + offset + j].astype(np.uint64)
    return codes


def scan_exact(seq: str, qmap: Mapping[int, str], k: int) -> Iterator[tuple[int, int, str]]:
    """Exact hits ``(start0, end0, canonical probe)`` of the research scanner, in order."""

    n = len(seq)
    if not qmap or n < k:
        return
    raw = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    values = _LUT[raw]
    invalid = values == 255
    values = np.where(invalid, 0, values).astype(np.uint8)
    bad = np.concatenate(([0], np.cumsum(invalid, dtype=np.int64)))
    low_mask = (1 << 64) - 1
    low_keys = np.array(sorted({code & low_mask for code in qmap}), dtype=np.uint64)
    windows = n - k + 1
    for block in range(0, windows, _BLOCK):
        stop = min(windows, block + _BLOCK)
        codes = _window_codes(values, k, block, stop)
        ok = (bad[block + k: stop + k] - bad[block: stop]) == 0
        hits = np.flatnonzero(ok & np.isin(codes, low_keys))
        for index in hits.tolist():
            start = block + index
            code = encode(seq[start: start + k]) if k > 32 else int(codes[index])
            probe = qmap.get(code)
            if probe is not None:
                yield start, start + k, probe


# ---------------------------------------------------------------------------
# Identifiers and families (per-stage research variants)
# ---------------------------------------------------------------------------


def extract_fbte_rebuild(text: str | None) -> str | None:
    if re.search(r"\bFBte000XXXX\b", text or "", re.I):
        return "FBte000XXXX"
    match = re.search(r"FBte\d+", text or "", re.I)
    if not match:
        return None
    return "FBte" + re.sub(r"\D", "", match.group(0))


def extract_fbte_gapfill(text: str | None) -> str | None:
    if re.search(r"\bFBte000XXXX\b", text or "", re.I):
        return "FBte000XXXX"
    match = re.search(r"\bFBte\d+\b", text or "", re.I)
    if not match:
        return None
    return "FBte" + re.sub(r"\D", "", match.group(0))


def extract_te_id_rebuild(header: str) -> str:
    return extract_fbte_rebuild(header) or header.split()[0].split("|")[0]


def extract_te_id_gapfill(header: str) -> str:
    return extract_fbte_gapfill(header) or header.split()[0].split("|")[0]


def canonical_family(text: str | None) -> str:
    return re.sub(r"[^a-z0-9]+", "", (text or "").strip().lower())


def _family_map_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _parse_family_map(path: str | Path, extract) -> dict[str, str]:
    rows = _family_map_rows(path)
    if not rows:
        raise RuntimeError("Empty family map.")
    headers = list(rows[0])
    id_column = next((h for h in headers if "fbte" in h.lower()), None)
    family_column = next((h for h in headers if "family" in h.lower()), None)
    if not id_column or not family_column:
        raise RuntimeError(f"Could not identify family-map columns: {headers}")
    out: dict[str, str] = {}
    for row in rows:
        te = extract(row.get(id_column, ""))
        family = (row.get(family_column, "") or "").strip()
        if te and family:
            out[te] = family
    if not out:
        raise RuntimeError("No family-map rows parsed.")
    return out


def parse_family_map_rebuild(path: str | Path) -> dict[str, str]:
    return _parse_family_map(path, extract_fbte_rebuild)


def parse_family_map_gapfill(path: str | Path) -> dict[str, str]:
    return _parse_family_map(path, extract_fbte_gapfill)


# ---------------------------------------------------------------------------
# RepeatMasker context
# ---------------------------------------------------------------------------

Annotation = tuple[int, int, str, str | None, str | None]


def _bed_rows(path: str | Path) -> Iterator[tuple[str, int, int, str]]:
    with Path(path).open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) < 3:
                continue
            yield columns[0], int(columns[1]), int(columns[2]), columns[3] if len(columns) >= 4 else ""


def _sorted(intervals: dict[str, list[Annotation]]) -> dict[str, list[Annotation]]:
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: (x[0], x[1], x[2]))
    return intervals


def load_repeatmasker_rebuild(path: str | Path, fbte_to_family: Mapping[str, str]) -> dict[str, list[Annotation]]:
    intervals: dict[str, list[Annotation]] = defaultdict(list)
    for chrom, start, end, name in _bed_rows(path):
        te = extract_fbte_rebuild(name)
        family = fbte_to_family.get(te) if te else None
        intervals[chrom].append((start, end, name, te, family))
    return _sorted(intervals)


def load_repeatmasker_gapfill(path: str | Path, fbte_to_family: Mapping[str, str]) -> dict[str, list[Annotation]]:
    by_norm: dict[str, str] = {}
    for family in fbte_to_family.values():
        by_norm.setdefault(canonical_family(family), family)
    intervals: dict[str, list[Annotation]] = defaultdict(list)
    for chrom, start, end, name in _bed_rows(path):
        te = extract_fbte_gapfill(name)
        if te:
            family = fbte_to_family.get(te)
        else:
            te, family = None, None
            raw = (name or "").strip()
            for token in (raw, raw.split("#")[0], raw.split("|")[0], raw.split(":")[0]):
                found = by_norm.get(canonical_family(token))
                if found:
                    family = found
                    break
        intervals[chrom].append((start, end, name, te, family))
    return _sorted(intervals)


def build_interval_index(intervals: Mapping[str, list[Annotation]]):
    out = {}
    for chrom, values in intervals.items():
        starts = [x[0] for x in values]
        prefix: list[int] = []
        running = -1
        for x in values:
            running = max(running, x[1])
            prefix.append(running)
        out[chrom] = (values, starts, prefix)
    return out


def overlap_query(index, chrom: str, start: int, end: int) -> list[Annotation]:
    item = index.get(chrom)
    if item is None:
        return []
    values, starts, prefix = item
    hi = bisect.bisect_left(starts, end)
    lo = hi
    while lo > 0 and prefix[lo - 1] > start:
        lo -= 1
    return [x for x in values[lo:hi] if x[1] > start and x[0] < end]


def reference_hit_class(target_family: str, target_te_id: str, start: int, end: int,
                        overlaps: Sequence[Annotation]) -> str:
    """Context class of one exact reference hit for one (family, record) target."""

    if not overlaps:
        return "NON_TE_OR_UNANNOTATED_HOST"
    full = [x for x in overlaps if x[0] <= start and x[1] >= end]
    partial = [x for x in overlaps if x not in full]
    if not full:
        return "TE_BOUNDARY_OR_PARTIAL"
    target_norm = canonical_family(target_family)
    record_label = target_family == target_te_id or target_family.startswith("FBte")
    safe_full, unsafe_full = [], []
    for x in full:
        _, _, _, annotation_te, annotation_family = x
        if annotation_family and canonical_family(annotation_family) == target_norm:
            safe_full.append(x)
        elif record_label and annotation_te and annotation_te == target_te_id:
            safe_full.append(x)
        else:
            unsafe_full.append(x)
    unsafe_partial = []
    for x in partial:
        _, _, _, annotation_te, annotation_family = x
        if annotation_family and canonical_family(annotation_family) == target_norm:
            continue
        if record_label and annotation_te and annotation_te == target_te_id:
            continue
        unsafe_partial.append(x)
    if safe_full and not unsafe_full and not unsafe_partial:
        return "TARGET_CONTEXT_ANNOTATED_TE"
    if unsafe_full:
        if any(x[4] for x in unsafe_full):
            return "OTHER_FAMILY_ANNOTATED_TE"
        return "UNMAPPED_ANNOTATED_TE"
    return "MIXED_REFERENCE_CONTEXT"


# ---------------------------------------------------------------------------
# TSV (the research writers' exact text)
# ---------------------------------------------------------------------------


def _norm_rebuild(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (list, tuple, set)):
        return ";".join(map(str, value))
    return str(value)


def _norm_internal(value: Any) -> Any:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (list, tuple, set)):
        return ";".join(map(str, value))
    return value


def tsv_bytes(rows: Sequence[Mapping[str, Any]], fields: Sequence[str] | None = None, *, style: str = "rebuild") -> bytes:
    """Serialize like the research ``write_tsv`` of the given stage.

    ``rebuild``/``gapfill``: values normalized (bool -> 1/0, lists ';'-joined);
    ``internal``: bool -> int, lists ';'-joined; ``raw``: values as given.
    """

    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    handle = io.StringIO(newline="")
    if not fields:
        return b""
    writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        if style in ("rebuild", "gapfill"):
            writer.writerow({k: _norm_rebuild(row.get(k, "")) for k in fields})
        elif style == "internal":
            writer.writerow({k: _norm_internal(row.get(k, "")) for k in fields})
        elif style == "raw":
            writer.writerow({k: row.get(k, "") for k in fields})
        else:
            raise ValueError(style)
    return handle.getvalue().encode("utf-8")


def read_tsv_bytes(data: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""), delimiter="\t"))
