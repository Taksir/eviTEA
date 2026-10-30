"""Genotype output: ``genotypes.tsv``, one row per call of ``calls.combined.tsv`` of a library-free output folder.

The step only adds a file. It reads ``calls.combined.tsv``, ``calls.library_free.tsv`` and ``library_free_family_calls.tsv`` and writes
``genotypes.tsv`` and ``genotype_summary.json`` into an output folder; no call file is opened for writing.

``ref_share`` is the share of the read evidence at the break that reads straight through the reference (the empty allele): 0 = every read
carries the insertion. ``genotype`` is hom at ``ref_share <= HOM_MAX``, het up to ``HET_MAX``, low above it (``_p4_lib.genotype``), and
``.`` when no count could be made.

Step A (``method`` lf_walk), library-free calls: the reference-path fraction the detection recorded for the call's place,
``calling.e_value(fracL, fracR, rstat)`` (the minimum of the two sides when the right walk broke, else the left side). It is the number the
places table's ``geno`` column is made from, read here from the 3-decimal ``fracL`` / ``fracR`` of the family-call table. A call whose
share falls within 0.0005 of a threshold can therefore fall on the other side of it than the unrounded ``geno`` column. ``n_ref`` and
``n_total`` are not kept by the detection: ``.``.

Step B (``method`` ref_walk), library calls: no empty-allele count exists, so the reference-path walk of the library-free layer
(``_p4_lib.walk``: ``frac = nref / T``) is run on the msBWT at the call's break. The left walk starts at the first sampled word before the
call's ``start`` and crosses the break to the right; the mirror walk starts after ``end``. ``ref_share`` is the smaller of the two fractions
when both walks reach a break, else the one that does; ``n_ref`` / ``n_total`` are the counts of that walk. The depth ``D`` is
``_p4_lib.local_D`` at the middle of the call (the median of the sampled word counts within 2 kb). A call whose anchor words are
outside the contig or hold a letter other than A, C, G, T, or are repeated more than twice the local depth, gets ``.``; so does a call
where neither walk finds a divergence (the insertion allele is below the walk's branch support). No repeat map is read.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import time
from collections import Counter
from pathlib import Path
from typing import Callable

from .errors import InputContractError

__all__ = ["COLUMNS", "InputFiles", "run_genotype", "genotype_tsv_bytes", "step_a_rows", "walk_at_break"]

COLUMNS = ("event_id", "chr", "start", "end", "family", "source", "ref_share", "n_ref", "n_total", "genotype", "method")
COMBINED, ADDED, FAMILY_CALLS = "calls.combined.tsv", "calls.library_free.tsv", "library_free_family_calls.tsv"
OUT_TSV, OUT_SUMMARY = "genotypes.tsv", "genotype_summary.json"
InputFiles = (COMBINED, ADDED, FAMILY_CALLS)
MISSING = "."


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    text = Path(path).read_bytes().decode("utf-8-sig")
    return [dict(row) for row in csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")]


def _share_text(value: float) -> str:
    return MISSING if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.3f}"


def _label(share: float | None) -> str:
    from .library_free import _p4_lib as P

    return MISSING if share is None or math.isnan(share) else P.genotype(share)


def step_a_rows(added: list[dict[str, str]], family_calls: list[dict[str, str]]) -> dict[str, tuple[str, str, str, str]]:
    """event_id -> (ref_share, n_ref, n_total, genotype) of every added library-free call.

    ``calls.library_free.tsv`` lists the ADDED rows of ``library_free_family_calls.tsv`` in order (``assemble.added_calls_tsv_bytes``), so
    the two are read side by side and each pair must agree on the contig and on ``bL``."""
    from .library_free import calling as CL

    rows = [r for r in family_calls if r.get("status") == "ADDED"]
    if len(rows) != len(added):
        raise InputContractError(
            f"{ADDED} has {len(added)} calls but {FAMILY_CALLS} has {len(rows)} ADDED rows: the two files are not from one run")
    out: dict[str, tuple[str, str, str, str]] = {}
    for call, row in zip(added, rows):
        if call["chr"] != row["chrom"] or str(call["bL"]) != str(row["bL"]):
            raise InputContractError(
                f"{ADDED} call {call['event_id']} ({call['chr']}:{call['bL']}) does not match the ADDED row of {FAMILY_CALLS} "
                f"({row['chrom']}:{row['bL']})")
        share = CL.e_value(row["fracL"], row["fracR"], row["rstat"])
        out[call["event_id"]] = (_share_text(share), MISSING, MISSING, _label(share))
    return out


def walk_at_break(C, ref: str, start: int, end: int):
    """Step B for one library call: ``(ref_share, n_ref, n_total, why)``; ``ref_share`` None when no count could be made.

    ``C`` counts both strands (``_p4_lib.Counter`` or ``fastcount.FastCounter``); ``ref`` is the upper-case contig."""
    from .library_free import _p4_lib as P

    D = P.local_D(C, ref, (start + end) // 2)
    if not D == D or D <= 0:
        return None, None, None, "no_depth"
    first = start - P.K + 1
    cand = dict(first=first, last=max(first, end - 1), D=D)
    lw, a = P.left_walk(C, ref, cand)
    rw, b, _, _ = P.right_walk(C, ref, cand)
    for walk_result, anchor in ((lw, a), (rw, b)):
        if walk_result["status"] == "edge":
            continue
        word = ref[anchor:anchor + P.K]
        if len(word) < P.K or not P.is_acgt(word) or C.cnt(word) > P.ANCHOR_RC * D:
            walk_result["status"] = "edge"
    sides = [w for w in (lw, rw) if w["status"] == "break"]
    if not sides:
        reasons = sorted({w["status"] for w in (lw, rw)})
        return None, None, None, "+".join(reasons)
    best = min(sides, key=lambda w: w["frac"])
    return float(best["frac"]), int(best["nref"]), int(best["T"]), "break"


def genotype_tsv_bytes(rows: list[dict[str, str]]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(COLUMNS)
    for row in rows:
        writer.writerow([row[c] for c in COLUMNS])
    return handle.getvalue().encode("utf-8")


def run_genotype(*, library_free_dir: Path, output: Path, reference_fasta: Path | None, msbwt_index: Path | None,
                 progress: Callable[[str], None] | None = None, command_line=(), source_identity: dict | None = None,
                 counter: str = "fast") -> dict:
    """Write ``genotypes.tsv`` and ``genotype_summary.json`` into ``output`` (made if absent; may be ``library_free_dir``).

    Without ``msbwt_index`` (or ``reference_fasta``) only Step A is made and the library calls are left out of the table."""

    say = progress or (lambda message: None)
    started = time.time()
    lf = Path(library_free_dir)
    paths = {name: lf / name for name in InputFiles}
    for name, path in paths.items():
        if not path.is_file():
            raise InputContractError(f"library-free folder has no {name}: {lf}")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for name in (OUT_TSV, OUT_SUMMARY):
        if (output / name).exists():
            raise InputContractError(f"{output / name} exists; the genotype step does not overwrite")
    before = {name: _sha256(path) for name, path in paths.items()}
    combined = _read_tsv(paths[COMBINED])
    a_rows = step_a_rows(_read_tsv(paths[ADDED]), _read_tsv(paths[FAMILY_CALLS]))
    walk_possible = msbwt_index is not None and reference_fasta is not None
    lib_calls = [r for r in combined if r["source"] == "EVITEA"]
    counts: Counter = Counter()
    walk_seconds = 0.0
    walk_reasons: Counter = Counter()
    b_rows: dict[str, tuple[str, str, str, str]] = {}
    if walk_possible and lib_calls:
        from .backends.msbwt import load_msbwt
        from .breakpoints.runtime import read_fasta
        from .library_free import _p4_lib as P
        from .library_free import fastcount as FC

        if counter not in ("plain", "fast"):
            raise InputContractError(f"--counter must be plain or fast, not {counter}")
        wanted = {r["chr"] for r in lib_calls}
        genome = read_fasta(Path(reference_fasta))
        genome = {name: seq for name, seq in genome.items() if name in wanted}
        bwt = load_msbwt(msbwt_index)
        C = FC.FastCounter(bwt, use=True) if counter == "fast" else P.Counter(bwt, cache=True)
        upper: dict[str, str] = {}
        walk_started = time.time()
        for n, call in enumerate(lib_calls, 1):
            chrom = call["chr"]
            if chrom not in genome:
                b_rows[call["event_id"]] = (MISSING, MISSING, MISSING, MISSING)
                walk_reasons["contig_not_in_reference"] += 1
                continue
            ref = upper.setdefault(chrom, genome[chrom].upper())
            share, n_ref, n_total, why = walk_at_break(C, ref, int(call["start"]), int(call["end"]))
            walk_reasons[why] += 1
            b_rows[call["event_id"]] = (_share_text(share), MISSING if n_ref is None else str(n_ref),
                                        MISSING if n_total is None else str(n_total), _label(share))
            if len(getattr(C, "_counts", getattr(C, "cache", ()))) > 4_000_000:
                C.clear()
            if n % 200 == 0:
                say(f"genotype: {n} of {len(lib_calls)} library calls walked")
        walk_seconds = time.time() - walk_started
    rows = []
    for call in combined:
        if call["source"] == "LIBRARY_FREE":
            if call["event_id"] not in a_rows:
                raise InputContractError(f"library-free call {call['event_id']} of {COMBINED} is not an ADDED call of {ADDED}")
            share, n_ref, n_total, geno = a_rows[call["event_id"]]
            method = "lf_walk"
        elif call["source"] == "EVITEA":
            if not walk_possible:
                counts["library_calls_left_out"] += 1
                continue
            share, n_ref, n_total, geno = b_rows[call["event_id"]]
            method = "ref_walk"
        else:
            raise InputContractError(f"{COMBINED}: unknown source {call['source']!r}")
        rows.append({**{k: call[k] for k in ("event_id", "chr", "start", "end", "family", "source")},
                     "ref_share": share, "n_ref": n_ref, "n_total": n_total, "genotype": geno, "method": method})
        counts[f"{call['source']}:{geno}"] += 1
    (output / OUT_TSV).write_bytes(genotype_tsv_bytes(rows))
    after = {name: _sha256(path) for name, path in paths.items()}
    if after != before:
        raise InputContractError("an input call file changed while the genotype step ran")
    summary = {
        "library_free_dir": str(lf),
        "rows": len(rows),
        "calls_combined_rows": len(combined),
        "genotype_counts": dict(sorted(counts.items())),
        "step_b_run": walk_possible,
        "step_b_outcomes": dict(sorted(walk_reasons.items())),
        "counter": counter,
        "seconds": {"step_b_walks": round(walk_seconds, 2), "total": round(time.time() - started, 2)},
        "input_sha256": before,
        "genotypes_sha256": hashlib.sha256((output / OUT_TSV).read_bytes()).hexdigest(),
        "thresholds": _thresholds(),
        "msbwt_index": None if msbwt_index is None else str(msbwt_index),
        "command_line": list(command_line),
        "source_identity": source_identity or {},
    }
    (output / OUT_SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    say(f"genotype: {len(rows)} rows written to {output / OUT_TSV}")
    return summary


def _thresholds() -> dict:
    from .library_free import _p4_lib as P

    return {"hom_max": P.HOM_MAX, "het_max": P.HET_MAX, "ref_ok": P.REF_OK, "anchor_rc": P.ANCHOR_RC}
