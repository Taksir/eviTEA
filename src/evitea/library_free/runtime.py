"""``evitea library-free``: detection over a sample msBWT, family calls, and the merge with a finished run.

Detection runs one contig at a time (resumable: each contig's places are kept
in the work directory with the identity they were computed from), then one
grouping over all scanned contigs in the given order. The family-call layer
reads the written places and segments tables, so a detection can be reused
with ``--detection-dir``. The source run directory is only read.
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import sys
import time
from collections import Counter as Tally
from collections.abc import Callable, Sequence
from pathlib import Path

from ..contigs import MITO_CONTIG_NAMES, find_mito_contig
from .. import vcf as V
from ..errors import InputContractError, OperationalResourceError
from ..outputs import write_bytes, write_run_summary
from ..shadow_layers import ShadowCallset
from . import NEAR_EVITEA_BP, RULE_C_MAX, T_E
from .options import Options

__all__ = ["MIN_DEFAULT_CONTIG", "combined_calls_tsv_bytes", "run_library_free"]

MIN_DEFAULT_CONTIG = 1_000_000
PLACES = "library_free_places.tsv"
SEGMENTS = "library_free_segments.fa"
SECOND_PLACES = "library_free_second_places.tsv"
SECOND_SEGMENTS = "library_free_second_segments.fa"
FAMILY_CALLS = "library_free_family_calls.tsv"
ADDED = "calls.library_free.tsv"
COMBINED = "calls.combined.tsv"
SUMMARY = "library_free_summary.json"
PIECES_READONLY_ENV = "EVITEA_PIECES_READONLY"
PIECE_SCHEMA = "evitea.library_free.piece.v2"
CODE_FILES = ("_p4_lib.py", "detector.py", "place.py", "scan.py", "walks.py", "extension.py", "priority.py", "options.py")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def _depth_cap_report():
    """The absolute upper bound of the local-depth window (None = no upper bound)."""

    from .. import depth_window

    value = depth_window.upper_bound()
    return None if value == float("inf") else value


def _peak_rss_gb() -> float | None:
    try:
        import resource
    except ImportError:
        return None
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024, 3)


def combined_calls_tsv_bytes(native: bytes, added: bytes) -> bytes:
    """EviTEA's native calls, then the added library-free calls: event_id, chr, start, end, family, source."""

    lines = ["event_id\tchr\tstart\tend\tfamily\tsource"]
    for row in ShadowCallset.from_bytes(native).rows:
        lines.append(f"{row['event_id']}\t{row['chr']}\t{row['start']}\t{row['end']}\t{row['family']}\tEVITEA")
    text = added.decode("utf-8").splitlines()
    if text:
        head = text[0].split("\t")
        for line in text[1:]:
            r = dict(zip(head, line.split("\t")))
            lines.append(f"{r['event_id']}\t{r['chr']}\t{r['start']}\t{r['end']}\t{r['family']}\tLIBRARY_FREE")
    return ("\n".join(lines) + "\n").encode("utf-8")


def combined_calls_vcf_bytes(native: bytes, added: bytes, reference, sample_id: str) -> bytes:
    """``calls.combined.vcf``: the library calls (as in ``calls.tsv``) and the added library-free calls, as VCF records."""

    calls = [*V.native_row_vcf_calls(ShadowCallset.from_bytes(native).rows), *V.library_free_vcf_calls(added)]
    fasta = Path(reference.reference_fasta_path)
    contigs = getattr(reference, "contigs", None) or V.fasta_contigs(fasta)
    return V.calls_vcf_bytes(calls, contigs, sample_id=sample_id, bases=V.reference_bases(fasta, calls),
                             reference_name=fasta.name)


def _code_sha256() -> str:
    """One hash over the detection code: it is part of the identity of a saved scan piece."""

    here = Path(__file__).parent
    digest = hashlib.sha256()
    for name in CODE_FILES:
        digest.update(name.encode() + b"\0" + _sha256(here / name).encode() + b"\0")
    return digest.hexdigest()


def _detect(*, reference_fasta: Path, genome: dict, msbwt_index: Path, contigs: Sequence[str] | None, work: Path,
            progress: Callable[[str], None], opts: Options, mito_contig: str | None = None):
    """Scan the contigs (resumable), then group each pass over all scanned contigs. Return the tables of the first pass
    (places bytes, segments bytes), of the second scan (or None) and an info dict."""

    from ..backends.msbwt import load_msbwt
    from . import _p4_lib as P
    from . import detector as D
    from . import priority as PR
    from . import tables as TB
    from .scan import scan_contig_passes
    from .walks import RefWords

    try:
        import mappy  # noqa: F401
    except ImportError as exc:
        raise OperationalResourceError("library-free detection needs mappy (minimap2); run on Linux or WSL") from exc
    if contigs is None:
        contigs = tuple(name for name, seq in genome.items() if len(seq) >= MIN_DEFAULT_CONTIG)
    missing = [c for c in contigs if c not in genome]
    if missing or len(set(contigs)) != len(contigs):
        raise InputContractError(f"contigs must be distinct reference contigs; not found: {missing}")
    mito = find_mito_contig(genome, mito_contig)
    progress(f"library-free: mitochondrial contig {mito}" if mito else
             "library-free: no mitochondrial contig found (" + ", ".join(MITO_CONTIG_NAMES)
             + "); the mitochondrial segment test is off (--mito-contig names one)")
    about = Path(msbwt_index) / "about.npy"
    identity_base = {
        "schema": PIECE_SCHEMA,
        "code_sha256": _code_sha256(),
        "reference_sha256": _sha256(reference_fasta),
        "msbwt_about_sha256": _sha256(about) if about.is_file() else None,
        "msbwt_index": str(msbwt_index),
        "depth_cap": _depth_cap_report(),
        "detection": opts.detection_identity,
        "mito_contig": mito,
    }
    pieces_dir = work / "pieces"
    pieces_dir.mkdir(parents=True, exist_ok=True)
    rm = words = None
    timings: dict[str, float] = {}
    first_all: list[dict] = []
    second_all: list[dict] = []
    stats: dict[str, dict] = {}
    for chrom in contigs:
        identity = {**identity_base, "contig": chrom}
        data, marker = pieces_dir / f"{chrom}.pkl", pieces_dir / f"{chrom}.complete.json"
        piece = None
        if data.is_file() and marker.is_file():
            recorded = json.loads(marker.read_text(encoding="utf-8"))
            if recorded.get("identity") == identity and _sha256(data) == recorded.get("sha256"):
                with data.open("rb") as handle:
                    piece = pickle.load(handle)
                progress(f"library-free: {chrom} reused from the work directory ({len(piece['first'])} places)")
        if piece is None and os.environ.get(PIECES_READONLY_ENV):
            raise InputContractError(
                f"{PIECES_READONLY_ENV} is set and the saved scan piece of {chrom} in {pieces_dir} does not match this run's "
                "detection switches, code or inputs; this run would scan again and replace it. Use a work directory of its own.")
        if piece is None:
            if rm is None:
                progress("library-free: reference multi-copy words")
                rm = D.reference_multi(genome)
            if opts.walk_out and words is None:
                progress("library-free: reference words for the walk-out")
                words = RefWords(genome)
            started = time.time()
            progress(f"library-free: scanning {chrom}")
            counter = P.Counter(load_msbwt(msbwt_index), cache=False)
            firewall = P.MappyFirewall(str(reference_fasta), genome[mito] if mito else None)
            if opts.is_previous:
                first, second, info = D.scan_contig(counter, genome[chrom], chrom, rm, firewall), None, {}
            else:
                first, second, info = scan_contig_passes(counter, genome[chrom], chrom, rm, firewall, opts, words)
            info = {k: v for k, v in info.items() if k not in ("walk_log", "cross_log")}
            piece = dict(first=first, second=second, info=info)
            payload = pickle.dumps(piece, protocol=pickle.HIGHEST_PROTOCOL)
            tmp = data.with_suffix(".tmp")
            tmp.write_bytes(payload)
            tmp.replace(data)
            marker.write_text(json.dumps({"identity": identity, "sha256": hashlib.sha256(payload).hexdigest()}, indent=2),
                              encoding="utf-8")
            timings[chrom] = round(time.time() - started, 1)
            progress(f"library-free: {chrom} done, {len(first)} places, {timings[chrom]} s, peak RSS {_peak_rss_gb()} GB")
        stats[chrom] = piece.get("info", {})
        first_all.extend(piece["first"])
        if piece.get("second") is not None:
            second_all.extend(piece["second"])
    progress("library-free: grouping")
    grouped = D.group_and_call(first_all) if opts.is_previous else PR.group_and_call(first_all, opts)
    tables = {"first": (TB.places_tsv_bytes(first_all, opts), TB.segments_fa_bytes(grouped))}
    if opts.het_scan != "off":
        grouped2 = PR.group_and_call(second_all, opts)
        tables["second"] = (TB.places_tsv_bytes(second_all, opts, second=True), TB.segments_fa_bytes(grouped2))
    stages = Tally(p["stage"] for p in first_all)
    info = {"contigs": list(contigs), "scan_seconds": timings, "stages": dict(sorted(stages.items())),
            "grouped_places": len(grouped)}
    if "second" in tables:
        info["second_scan"] = {"stages": dict(sorted(Tally(p["stage"] for p in second_all).items())),
                               "origins": dict(sorted(Tally(p.get("origin", "") for p in second_all if p["stage"] == "S4_pass").items())),
                               "grouped_places": len(grouped2)}
    if any(stats.values()):
        info["scan_counts"] = stats
    if opts.relative_depth_floor:
        info["depth_floors"] = {chrom: stats[chrom].get("depth_floor") for chrom in contigs}
    info["mito_contig"] = mito
    return tables, info


def _load_copies(opts: Options, reference, repeat_annotation, genome, progress, repeat_map=None, library_calls=None):
    """(CopyIndex or None, bases of the reference) for the resident-copy test. The repeat map stands in for the annotation of a run
    in which the library layer made no call (``library_calls`` == 0), unless ``repeat_annotation`` is given."""

    if opts.resident_guard == "off":
        return None, 0
    from .resident import CopyIndex

    path = repeat_annotation
    if path is None and repeat_map is not None and library_calls == 0:
        progress("library-free: no library call; the resident-copy test reads the repeat map")
        path = repeat_map
    if path is None:
        try:
            path = reference.repeat_annotation_path
        except Exception:  # the prepared reference binds no annotation
            path = None
    if path is None:
        progress("library-free: no repeat annotation; the resident-copy test is off")
        return None, 0
    copies = CopyIndex.from_file(path, contigs=set(genome))
    return copies, sum(len(v) for v in genome.values())


def run_library_free(*, run_dir: Path, reference, msbwt_index: Path | None, output: Path, work_dir: Path | None,
                     contigs: Sequence[str] | None, rule_c: bool, detection_dir: Path | None,
                     command_line: Sequence[str] = (), source_identity: dict | None = None,
                     progress: Callable[[str], None] | None = None, options: Options | None = None,
                     repeat_annotation: Path | None = None, repeat_map: Path | None = None,
                     mito_contig: str | None = None) -> dict:
    from ..breakpoints.runtime import read_fasta
    from ..pipeline import prepare_output_directory
    from . import assemble as AS
    from . import calling as CL
    from . import detector as D

    opts = options if options is not None else Options.previous()
    say = progress or (lambda message: print(message, file=sys.stderr, flush=True))
    started = time.time()
    run_dir = Path(run_dir)
    native = (run_dir / "calls.native.tsv").read_bytes()
    try:
        sample_id = json.loads((run_dir / "run_summary.json").read_text(encoding="utf-8"))["sample_id"]
    except (OSError, json.JSONDecodeError, KeyError) as exc:
        raise InputContractError(f"run directory has no readable run_summary.json: {run_dir}") from exc
    output = prepare_output_directory(output, resume=False)
    work = Path(work_dir) if work_dir is not None else output / ".evitea-work"
    reference_fasta = Path(reference.reference_fasta_path)
    genome = None
    if detection_dir is not None:
        detection_dir = Path(detection_dir)
        tables = {"first": ((detection_dir / PLACES).read_bytes(), (detection_dir / SEGMENTS).read_bytes())}
        if (detection_dir / SECOND_PLACES).is_file():
            tables["second"] = ((detection_dir / SECOND_PLACES).read_bytes(), (detection_dir / SECOND_SEGMENTS).read_bytes())
        info = {"reused_detection": str(detection_dir)}
    else:
        if msbwt_index is None:
            raise InputContractError("library-free detection needs an msBWT index")
        genome = read_fasta(reference_fasta)
        tables, info = _detect(reference_fasta=reference_fasta, genome=genome, msbwt_index=Path(msbwt_index), contigs=contigs,
                               work=work, progress=say, opts=opts, mito_contig=mito_contig)
    if ("second" in tables) != (opts.het_scan != "off") and detection_dir is not None:
        raise InputContractError(
            "the detection directory " + ("holds" if "second" in tables else "has no") + " second-scan tables but --het-scan is "
            + opts.het_scan)
    write_bytes(output / PLACES, tables["first"][0])
    write_bytes(output / SEGMENTS, tables["first"][1])
    if "second" in tables:
        write_bytes(output / SECOND_PLACES, tables["second"][0])
        write_bytes(output / SECOND_SEGMENTS, tables["second"][1])
    say("library-free: family calls (v2 filter" + (", rule c" if rule_c else "") + ")")
    first = (D.parse_places(tables["first"][0]), D.parse_segments(tables["first"][1]))
    second = (D.parse_places(tables["second"][0]), D.parse_segments(tables["second"][1])) if "second" in tables else None
    if genome is None and (opts.check_scope != "none" or opts.resident_guard != "off"):
        genome = read_fasta(reference_fasta)
    native_rows = ShadowCallset.from_bytes(native).rows
    copies, genome_bp = (_load_copies(opts, reference, repeat_annotation, genome, say, repeat_map, len(native_rows))
                         if genome is not None else (None, 0))
    if repeat_map is not None and (repeat_annotation is not None or native_rows):
        say("library-free: --repeat-map not used: " + ("--repeat-annotation is given" if repeat_annotation is not None
                                                      else f"the library layer made {len(native_rows)} calls"))
    c_started = time.time()
    index = CL.CallIndex.from_native(native_rows)
    rows, assembly = AS.assemble(first, second, index, opts, identity_filter=rule_c, reference=genome, copies=copies, genome_bp=genome_bp)
    c_seconds = round(time.time() - c_started, 1)
    extended = not opts.is_previous
    added = AS.added_calls_tsv_bytes(rows)
    write_bytes(output / FAMILY_CALLS, AS.family_calls_tsv_bytes(rows, extended=extended))
    write_bytes(output / ADDED, added)
    write_bytes(output / COMBINED, combined_calls_tsv_bytes(native, added))
    write_bytes(output / V.COMBINED_VCF_OUTPUT, combined_calls_vcf_bytes(native, added, reference, sample_id))
    statuses = Tally(r.status for r in rows)
    names = [PLACES, SEGMENTS, FAMILY_CALLS, ADDED, COMBINED, V.COMBINED_VCF_OUTPUT] + ([SECOND_PLACES, SECOND_SEGMENTS] if second is not None else [])
    summary = {
        "sample_id": sample_id,
        "source_run": str(run_dir),
        "source_native_callset_sha256": hashlib.sha256(native).hexdigest(),
        "command_line": list(command_line),
        "source_identity": source_identity or {},
        "parameters": {"v2_reference_path_fraction_max": T_E, "rule_c": bool(rule_c), "rule_c_similarity_max": RULE_C_MAX,
                       "near_evitea_bp": NEAR_EVITEA_BP},
        "options": opts.as_dict(),
        "detection": info,
        "assembly": assembly,
        "depth_cap": _depth_cap_report(),
        "family_calls": len(rows),
        "family_call_statuses": dict(sorted(statuses.items())),
        "families_with_similarity": assembly["families_with_similarity"],
        "rule_c_seconds": c_seconds,
        "seconds": round(time.time() - started, 1),
        "peak_rss_gb": _peak_rss_gb(),
        "output_hashes": {name: _sha256(output / name) for name in names},
    }
    if repeat_map is not None:
        summary["repeat_map"] = {"path": str(repeat_map), "used": repeat_annotation is None and not native_rows and copies is not None,
                                 "library_calls": len(native_rows)}
    write_run_summary(output / SUMMARY, summary)
    say(f"library-free: {len(rows)} family calls, {statuses.get('ADDED', 0)} added")
    return summary
