"""Strict whole-read TE-record fragment decoder with masked single-probe
source competition (the relaxed-provisional novelty mask input).

For every candidate TE record of the TE-positive read, the decoder aligns the
read and its reverse complement by infix (``HW``) edit distance.  A pair is
*strictly fragment-decodable* when one record is uniquely best, its
orientation is resolved, every assigned probe occurs in the read orientation
implied by the alignment, one best location spans every positionally
interpretable probe, and the alignment's terminal side agrees with the
codebook.  Strict single-probe pairs outside the primary set additionally get
an exact masked-source competition: the decoder record must be among the
records with the most exact bases adjacent to the (masked) probe.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path

from ..errors import BackendCompatibilityError, InputContractError
from .codebook import CodeUnit, SafeCodebook, revcomp

__all__ = ["DECODER_FIELDS", "decode_fragments", "decode_fragments_external", "edlib_available", "load_te_library", "load_record_families"]

DECODER_FIELDS = (
    "pair_id", "old_primary", "te_mate", "host_mate", "all_te_query_unit_count", "unit_ids",
    "candidate_records", "decoder_class", "strict_fragment_decodable", "family",
    "qualifying_te_records", "decoder_best_edit", "decoder_orientation", "decoder_terminal_side",
    "full_read_et3", "masked_source_class", "masked_family_confirmed", "masked_record_confirmed",
    "masked_best_outside_exact_bases",
)


def edlib_available() -> bool:
    try:
        import edlib  # type: ignore  # noqa: F401
    except ImportError:
        return False
    return True


def _edlib():
    try:
        import edlib  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise BackendCompatibilityError("edlib is required for the paired-end anchor fragment decoder") from exc
    return edlib


def load_te_library(path: str | Path) -> dict[str, str]:
    """``first header token -> upper-case sequence``."""

    out: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with Path(path).open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out[name.split()[0]] = "".join(chunks).upper()
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
        if name is not None:
            out[name.split()[0]] = "".join(chunks).upper()
    if not out:
        raise InputContractError("TE library is empty")
    return out


def load_record_families(codebook: SafeCodebook, family_map_rows: Iterable[Mapping[str, str]]) -> dict[str, str]:
    """Codebook record families, filled for absent records from the family map."""

    families = codebook.record_families()
    rows = list(family_map_rows)
    if rows:
        header = list(rows[0])
        fb = next((h for h in header if "fbte" in h.lower() or h == "te_record_id"), None)
        fam = next((h for h in header if "family" in h.lower()), None)
        if fb and fam:
            for row in rows:
                te_id, family = str(row.get(fb, "")).strip(), str(row.get(fam, "")).strip()
                if te_id and family:
                    families.setdefault(te_id, family)
    return families


def _alignment_summary(edlib, read: str, te_seq: str) -> dict:
    fwd = edlib.align(read, te_seq, mode="HW", task="locations")
    rev = edlib.align(revcomp(read).upper(), te_seq, mode="HW", task="locations")
    ef, er = int(fwd.get("editDistance", -1)), int(rev.get("editDistance", -1))
    if ef < 0 or er < 0:
        best = min(x for x in (ef, er) if x >= 0) if max(ef, er) >= 0 else -1
        return {"orientation": "UNRESOLVED", "best_edit": best, "locations": []}
    if ef < er:
        return {"orientation": "+", "best_edit": ef, "locations": [tuple(map(int, x)) for x in fwd.get("locations") or []]}
    if er < ef:
        return {"orientation": "-", "best_edit": er, "locations": [tuple(map(int, x)) for x in rev.get("locations") or []]}
    locs = list({tuple(x) for x in (fwd.get("locations") or []) + (rev.get("locations") or [])})
    return {"orientation": "TIE", "best_edit": ef, "locations": [tuple(map(int, x)) for x in locs]}


def _terminal_side(locs: list[tuple[int, int]], te_len: int) -> str:
    if not locs:
        return "UNRESOLVED"
    sides = set()
    for start, end in locs:
        d5, d3 = max(0, start), max(0, (te_len - 1) - end)
        sides.add("5p" if d5 < d3 else "3p" if d3 < d5 else "MIDPOINT_TIE")
    return next(iter(sides)) if len(sides) == 1 else "MULTI_LOCATION_SIDE_CONFLICT"


def _probe_consistency(read: str, queries: list[str], te_id: str, orientation: str,
                       locs: list[tuple[int, int]], index) -> tuple[bool, bool, str]:
    """(orientation consistent, source span supported, codebook terminal side)."""

    if orientation not in {"+", "-"}:
        return False, False, "UNRESOLVED"
    conflicts = unresolved = 0
    metas: list[CodeUnit] = []
    for q0 in queries:
        q = q0.upper()
        hits = list(index.get((te_id, q), [])) or list(index.get((te_id, revcomp(q).upper()), []))
        if not hits:
            unresolved += 1
            continue
        canonical = hits[0].canonical
        expected = canonical if orientation == "+" else revcomp(canonical).upper()
        opposite = revcomp(canonical).upper() if orientation == "+" else canonical
        expected_hit, opposite_hit = expected in read, opposite in read
        conflicts += int(not expected_hit)
        unresolved += int(expected_hit and opposite_hit)
        metas.extend(hits)
    interp = [m for m in metas if m.positionally_interpretable and m.source_start >= 0 and m.source_end >= 0]
    sides = {m.terminal_side for m in interp if m.terminal_side in {"5p", "3p"}}
    side = next(iter(sides)) if len(sides) == 1 else ("MIXED" if len(sides) > 1 else "UNRESOLVED")
    span = bool(interp and locs) and any(
        all(st <= m.source_start and en >= (m.source_end - 1) for m in interp) for st, en in locs
    )
    return conflicts == 0 and unresolved == 0 and bool(queries), span, side


def _find_all(text: str, pattern: str):
    start = 0
    while True:
        j = text.find(pattern, start)
        if j < 0:
            return
        yield j
        start = j + 1


def _outside_score(read: str, rp: int, k: int, source: str, sp: int) -> tuple[int, int]:
    left = 0
    while rp - 1 - left >= 0 and sp - 1 - left >= 0 and read[rp - 1 - left] == source[sp - 1 - left]:
        left += 1
    right = 0
    while rp + k + right < len(read) and sp + k + right < len(source) and read[rp + k + right] == source[sp + k + right]:
        right += 1
    return left, right


def decode_fragments(
    inputs: Iterable[Mapping[str, object]],
    codebook: SafeCodebook,
    te_records: Mapping[str, str],
    record_families: Mapping[str, str],
) -> list[dict[str, object]]:
    """Decode ``{pair_id, te_mate, host_mate, unit_ids, old_primary, te_read}`` rows."""

    edlib = _edlib()
    units = {u.unit_id: u for u in codebook.units}
    index = codebook.probe_index()
    seed_cache: dict[str, list[tuple[str, str, str, int]]] = {}

    def seed_occurrences(query: str):
        hit = seed_cache.get(query)
        if hit is None:
            hit = []
            for tid, seq in te_records.items():
                for orient, text in (("+", seq), ("-", revcomp(seq))):
                    for sp in _find_all(text, query):
                        hit.append((tid, orient, text, sp))
            seed_cache[query] = hit
        return hit

    out: list[dict[str, object]] = []
    for row in inputs:
        pid = int(row["pair_id"])
        read = str(row["te_read"]).upper()
        uids = [int(x) for x in row["unit_ids"]]  # type: ignore[union-attr]
        if not uids:
            raise InputContractError(f"pair {pid}: no codebook units")
        old_primary = bool(row["old_primary"])
        records = sorted({units[u].te_id for u in uids if units[u].te_id})
        candidates = []
        for tid in records:
            if tid not in te_records:
                raise InputContractError(f"TE record {tid} missing from the TE library")
            queries = [units[u].canonical for u in uids if units[u].te_id == tid]
            summary = _alignment_summary(edlib, read, te_records[tid])
            side = _terminal_side(summary["locations"], len(te_records[tid]))
            consistent, span, cb_side = _probe_consistency(read, queries, tid, summary["orientation"], summary["locations"], index)
            candidates.append({"te_id": tid, "family": record_families.get(tid, ""), **summary,
                               "alignment_terminal_side": side, "orientation_consistent": consistent,
                               "span_supported": span, "codebook_terminal_side": cb_side})
        valid = [c for c in candidates if int(c["best_edit"]) >= 0]
        best_edit = min((int(c["best_edit"]) for c in valid), default=-1)
        bests = [c for c in valid if int(c["best_edit"]) == best_edit]
        state = "UNIQUE_BEST_RECORD" if len(bests) == 1 else "BEST_RECORD_TIE" if len(bests) > 1 else "NO_RECORD_ALIGNMENT"
        b = bests[0] if len(bests) == 1 else (
            min(candidates, key=lambda c: (999999 if int(c["best_edit"]) < 0 else int(c["best_edit"]), c["te_id"]))
            if candidates else None
        )
        side_agree = bool(b and b["codebook_terminal_side"] in {"5p", "3p"} and b["alignment_terminal_side"] == b["codebook_terminal_side"])
        strict = bool(b and state == "UNIQUE_BEST_RECORD" and b["orientation"] in {"+", "-"}
                      and b["orientation_consistent"] and b["span_supported"] and side_agree)
        if strict:
            cls = "FRAGMENT_DECODABLE_TERMINAL_RECORD"
        elif state != "UNIQUE_BEST_RECORD":
            cls = state
        elif b and b["orientation"] not in {"+", "-"}:
            cls = "WHOLE_READ_ORIENTATION_UNRESOLVED"
        elif b and not b["orientation_consistent"]:
            cls = "PROBE_ORIENTATION_CONFLICT_WITH_RECORD_ALIGNMENT"
        elif b and not b["span_supported"]:
            cls = "PROBE_SOURCE_SPAN_INCOMPATIBLE_WITH_RECORD_ALIGNMENT"
        elif b and not side_agree:
            cls = "CODEBOOK_TERMINAL_SIDE_DISAGREES_WITH_RECORD_ALIGNMENT"
        else:
            cls = "UNRESOLVED_FRAGMENT_DECODER_STATE"
        result: dict[str, object] = {
            "pair_id": pid, "old_primary": old_primary, "te_mate": row["te_mate"], "host_mate": row["host_mate"],
            "all_te_query_unit_count": len(uids), "unit_ids": ";".join(map(str, uids)),
            "candidate_records": ";".join(records), "decoder_class": cls, "strict_fragment_decodable": int(strict),
            "family": "" if b is None else b["family"], "qualifying_te_records": b["te_id"] if strict else "",
            "decoder_best_edit": "" if b is None else b["best_edit"],
            "decoder_orientation": "" if b is None else b["orientation"],
            "decoder_terminal_side": "" if b is None else b["alignment_terminal_side"],
            "full_read_et3": int(bool(b and 0 <= int(b["best_edit"]) <= 3)),
            "masked_source_class": "NOT_SINGLE_PROBE_NEW_STRICT", "masked_family_confirmed": 0,
            "masked_record_confirmed": 0, "masked_best_outside_exact_bases": "",
        }
        if strict and not old_primary and len(uids) == 1:
            unit = units[uids[0]]
            k = unit.k or len(unit.canonical)
            strands = list(dict.fromkeys(unit.sequences)) or [unit.canonical, revcomp(unit.canonical)]
            seeds = [(q, rp) for q in strands for rp in _find_all(read, q)]
            if not seeds:
                raise InputContractError(f"pair {pid}: exact seed absent from the TE read")
            per_record: dict[str, tuple] = {}
            for q, rp in seeds:
                for tid, orient, source, sp in seed_occurrences(q):
                    left, right = _outside_score(read, rp, k, source, sp)
                    cand = (left + right, left, right, orient, rp, sp, q)
                    prev = per_record.get(tid)
                    if prev is None or cand[:3] > prev[:3]:
                        per_record[tid] = cand
            decoded, dfam = b["te_id"], b["family"]  # type: ignore[index]
            scores = {tid: v[0] for tid, v in per_record.items()}
            best = max(scores.values()) if scores else 0
            best_records = {tid for tid, s in scores.items() if s == best}
            best_fams = {record_families.get(tid) or f"UNMAPPED::{tid}" for tid in best_records}
            if decoded not in scores:
                mcls = "DECODER_RECORD_NOT_IN_GLOBAL_SEED_SOURCE_SET"
            elif best <= 0:
                mcls = "MASKED_NO_OUTSIDE_EXACT_SUPPORT"
            elif decoded not in best_records:
                mcls = "MASKED_COMPETITION_FAVORS_OTHER_SOURCE"
            elif len(best_fams) != 1 or dfam not in best_fams:
                mcls = "MASKED_BEST_FAMILY_AMBIGUOUS_OR_DISAGREES"
            elif len(best_records) == 1:
                mcls = "MASKED_UNIQUE_RECORD_CONFIRMS_DECODER"
            else:
                mcls = "MASKED_UNIQUE_FAMILY_RECORD_AMBIGUOUS"
            result.update({
                "masked_source_class": mcls,
                "masked_family_confirmed": int(mcls in {"MASKED_UNIQUE_RECORD_CONFIRMS_DECODER", "MASKED_UNIQUE_FAMILY_RECORD_AMBIGUOUS"}),
                "masked_record_confirmed": int(mcls == "MASKED_UNIQUE_RECORD_CONFIRMS_DECODER"),
                "masked_best_outside_exact_bases": best,
            })
        out.append(result)
    return out


def decode_fragments_external(
    python: Path,
    rows: list[Mapping[str, object]],
    *,
    codebook_path: Path,
    te_fasta: Path,
    family_map: Path,
    work: Path,
) -> list[dict[str, object]]:
    """Run the same decoder under another interpreter that has ``edlib``.

    For an environment without edlib (edlib publishes no
    CPython 3.14 Windows wheel).  The code and results are identical.
    """

    work.mkdir(parents=True, exist_ok=True)
    inp, out = work / "decoder_input.jsonl", work / "decoder_output.jsonl"
    with inp.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps({**row, "unit_ids": list(row["unit_ids"])}) + "\n")  # type: ignore[arg-type]
    src = Path(__file__).resolve().parents[2]
    command = [str(python), "-B", "-c",
               f"import sys; sys.path.insert(0, {str(src)!r}); from evitea.anchor_pe.decoder import _worker; _worker(sys.argv[1:])",
               str(codebook_path), str(te_fasta), str(family_map), str(inp), str(out)]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise BackendCompatibilityError(f"external decoder failed: {completed.stderr[-3000:]}")
    with out.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _worker(argv: list[str]) -> None:
    import csv

    from .codebook import load_safe_codebook

    codebook_path, te_fasta, family_map, inp, out = map(Path, argv)
    codebook = load_safe_codebook(codebook_path)
    with family_map.open(encoding="utf-8-sig", newline="") as handle:
        families = load_record_families(codebook, list(csv.DictReader(handle, delimiter="\t")))
    with inp.open(encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    decoded = decode_fragments(rows, codebook, load_te_library(te_fasta), families)
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        for row in decoded:
            handle.write(json.dumps(row) + "\n")
    sys.stdout.flush()
