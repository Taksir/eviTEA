"""VCF view of a call set (sites only, VCF 4.2).

Each call becomes one record with a symbolic insertion allele ``<INS:ME:family>``:

* ``POS`` is the base before the call: a point call at interbase position ``p`` (0-based) gives ``POS = p`` (the 1-based
  number of base ``p - 1``); a range call ``[start, end)`` (0-based, half-open, as in ``calls.tsv``) gives ``POS = start``
  and ``END = end``, with ``CIPOS=0,<end - start>`` and the ``IMPRECISE`` flag: the insertion lies somewhere in the range.
* ``REF`` is the reference base at ``POS`` (``N`` when the reference is not given).
* ``INFO`` carries ``SVTYPE=INS``, the family (``FAMILY``, and ``FAMILIES`` when several were listed), the call source
  (``SOURCE=library`` or ``library-free``), the stage that made the call (``ORIGIN``), its support counts and, for library-free
  calls, the target-site duplication length (``TSD``).
* ``FILTER`` is ``PASS``: the call set holds only reported calls.  No genotype is given (sites-only VCF).

The file is deterministic: no date and no absolute paths are written; records follow the reference's contig order, then
``POS``, then ``ID``.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ._version import __version__
from .errors import InputContractError

__all__ = [
    "COMBINED_VCF_OUTPUT",
    "VCF_OUTPUT",
    "VcfCall",
    "calls_vcf_bytes",
    "fasta_contigs",
    "library_free_vcf_calls",
    "native_row_vcf_calls",
    "native_vcf_calls",
    "reference_bases",
    "vcf_position",
]

VCF_OUTPUT = "calls.vcf"
COMBINED_VCF_OUTPUT = "calls.combined.vcf"
_UNSAFE = re.compile(r"[^A-Za-z0-9_.\-+]")


@dataclass(frozen=True, slots=True)
class VcfCall:
    call_id: str
    chrom: str
    start: int  # 0-based
    end: int  # 0-based, exclusive; equal to ``start`` for a point call
    family: str
    source: str  # "library" or "library-free"
    families: tuple[str, ...] = ()
    point: bool = False
    info: Mapping[str, object] = field(default_factory=dict)


def _allele_name(family: str) -> str:
    return _UNSAFE.sub("_", family) or "unknown"


def _info_value(value: object) -> str:
    text = str(value)
    return text.replace(" ", "_").replace(";", ",").replace("=", ":")


def native_vcf_calls(views: Iterable[Mapping[str, object]]) -> list[VcfCall]:
    """Return the calls of the library layer, from the ``calls.tsv`` views (``callset._sorted_views``)."""

    out = []
    for view in views:
        families = view.get("families") or ()
        if isinstance(families, str):
            families = tuple(f for f in re.split(r"[\[\]\"', ]+", families) if f)
        point = view.get("coordinate_kind") == "POINT"
        start = int(view["start"])
        out.append(VcfCall(
            call_id=str(view["call_id"]), chrom=str(view["chrom"]), start=start,
            end=start if point else int(view["end"]), family=str(view["family"]), source="library",
            families=tuple(str(f) for f in families), point=point,
            info={"ORIGIN": view.get("call_origin", ""), "EVIDENCE": view.get("evidence_type", ""),
                  "LEFT_SUPPORT": view.get("left_support", ""), "RIGHT_SUPPORT": view.get("right_support", ""),
                  "SUPPORT": view.get("total_support", "")},
        ))
    return out


def native_row_vcf_calls(rows: Iterable[Mapping[str, str]]) -> list[VcfCall]:
    """Return the calls of a native call table (``calls.native.tsv`` rows); a call with ``start == end`` is a point.  Only
    ``event_id``, ``chr``, ``start``, ``end`` and ``family`` are required; the other columns are used when present."""

    from .callset import CALL_ORIGINS

    views = []
    for row in rows:
        start, end = int(row["start"]), int(row["end"])
        views.append({
            "call_id": row["event_id"], "chrom": row["chr"], "start": start, "end": end, "family": row["family"],
            "families": row.get("families", ""), "coordinate_kind": "POINT" if start == end else "RANGE",
            "call_origin": CALL_ORIGINS.get(row.get("evidence_type", ""), ""), "evidence_type": row.get("evidence_type", ""),
            "left_support": row.get("left_support", ""), "right_support": row.get("right_support", ""),
            "total_support": row.get("total_support", ""),
        })
    return native_vcf_calls(views)


def fasta_contigs(fasta: str | Path) -> list[tuple[str, int]]:
    """Return (name, length) of every FASTA record, in file order."""

    contigs: list[tuple[str, int]] = []
    with Path(fasta).open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                contigs.append((line[1:].split()[0], 0))
            elif contigs:
                contigs[-1] = (contigs[-1][0], contigs[-1][1] + len(line.strip()))
    return contigs


def library_free_vcf_calls(added: bytes) -> list[VcfCall]:
    """Return the calls of ``calls.library_free.tsv`` (columns event_id, chr, start, end, family, ..., tsd)."""

    lines = added.decode("utf-8").splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    out = []
    for line in lines[1:]:
        if not line.strip():
            continue
        row = dict(zip(head, line.split("\t")))
        info: dict[str, object] = {"ORIGIN": "LIBRARY_FREE"}
        if row.get("tsd", ""):
            info["TSD"] = row["tsd"]
        out.append(VcfCall(row["event_id"], row["chr"], int(row["start"]), int(row["end"]), row["family"],
                           "library-free", info=info))
    return out


_HEADER_INFO = (
    ("SVTYPE", "1", "String", "Type of structural variant"),
    ("END", "1", "Integer", "End of the call range (1-based, inclusive)"),
    ("CIPOS", "2", "Integer", "Confidence interval of the insertion point: it lies within POS..POS+CIPOS[2]"),
    ("IMPRECISE", "0", "Flag", "The insertion point lies somewhere in POS..END"),
    ("FAMILY", "1", "String", "TE family of the inserted element"),
    ("FAMILIES", ".", "String", "All TE families listed for the call, when more than one"),
    ("SOURCE", "1", "String", "Call source: library (TE library layer) or library-free"),
    ("ORIGIN", "1", "String", "Stage that made the call (call_origin of calls.tsv, or LIBRARY_FREE)"),
    ("EVIDENCE", "1", "String", "Native evidence type of a library call"),
    ("LEFT_SUPPORT", "1", "Integer", "Read pairs supporting the left junction"),
    ("RIGHT_SUPPORT", "1", "Integer", "Read pairs supporting the right junction"),
    ("SUPPORT", "1", "Integer", "Read pairs supporting the call"),
    ("TSD", "1", "Integer", "Target-site duplication length of a library-free call (bp)"),
)


def vcf_position(call: VcfCall) -> int:
    """Return ``POS`` of a call (1-based): the base before the point or range."""

    return max(1, call.start)


def reference_bases(fasta: str | Path, calls: Iterable[VcfCall]) -> dict[tuple[str, int], str]:
    """Return the REF base of every call, read from the FASTA in one pass without holding whole contigs."""

    wanted: dict[str, set[int]] = {}
    for call in calls:
        wanted.setdefault(call.chrom, set()).add(vcf_position(call))
    bases: dict[tuple[str, int], str] = {}
    contig, offset, targets = None, 0, ()
    with Path(fasta).open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                contig, offset = line[1:].split()[0], 0
                targets = sorted(wanted.get(contig, ()))
                continue
            text = line.strip()
            if targets:
                for pos in targets:
                    if offset < pos <= offset + len(text):
                        bases[(contig, pos)] = text[pos - offset - 1].upper()
            offset += len(text)
    return bases


def calls_vcf_bytes(
    calls: Sequence[VcfCall],
    contigs: Sequence[tuple[str, int]],
    *,
    sample_id: str,
    bases: Mapping[tuple[str, int], str] | None = None,
    reference_name: str | None = None,
) -> bytes:
    """Return the VCF 4.2 text of ``calls``.  ``contigs`` is (name, length) in reference order.  ``bases`` gives the REF
    bases (``reference_bases``; ``N`` without it)."""

    order = {name: rank for rank, (name, _) in enumerate(contigs)}
    lengths = dict(contigs)
    records = []
    for call in calls:
        if call.chrom not in order:
            raise InputContractError(f"call {call.call_id} is on {call.chrom!r}, which is not a reference contig")
        if not 0 <= call.start <= call.end <= lengths[call.chrom]:
            raise InputContractError(f"call {call.call_id} lies outside {call.chrom}")
        pos = vcf_position(call)
        ref = (bases or {}).get((call.chrom, pos), "N")
        if ref not in ("A", "C", "G", "T", "N"):
            ref = "N"
        info = ["SVTYPE=INS"]
        if not call.point:
            info += [f"END={max(pos, call.end)}", f"CIPOS=0,{max(0, call.end - pos)}", "IMPRECISE"]
        info.append(f"FAMILY={_info_value(call.family)}")
        if len(call.families) > 1:
            info.append("FAMILIES=" + ",".join(_info_value(f) for f in call.families))
        info.append(f"SOURCE={call.source}")
        for key, value in call.info.items():
            if value not in ("", None):
                info.append(f"{key}={_info_value(value)}")
        records.append(((order[call.chrom], pos, call.call_id), "\t".join((
            call.chrom, str(pos), call.call_id, ref, f"<INS:ME:{_allele_name(call.family)}>", ".", "PASS", ";".join(info),
        ))))
    records.sort(key=lambda item: item[0])
    families = sorted({_allele_name(call.family) for call in calls})
    header = ["##fileformat=VCFv4.2", f"##source=EviTEA {__version__}"]
    if reference_name:
        header.append(f"##reference={reference_name}")
    header += [f"##contig=<ID={name},length={length}>" for name, length in contigs]
    header.append('##ALT=<ID=INS:ME,Description="Insertion of a mobile element">')
    header += [f'##ALT=<ID=INS:ME:{name},Description="Insertion of a {name} element">' for name in families]
    header += [f'##INFO=<ID={key},Number={number},Type={kind},Description="{text}">' for key, number, kind, text in _HEADER_INFO]
    header.append('##FILTER=<ID=PASS,Description="All filters passed">')
    header.append(f"##SAMPLE=<ID={sample_id}>")
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return ("\n".join(header + [line for _, line in records]) + "\n").encode("utf-8")
