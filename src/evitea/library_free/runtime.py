"""``evitea library-free``: detection over a sample msBWT, family calls, and the merge with a finished run.

Detection runs one contig at a time (resumable: each contig's places are kept
in the work directory with the identity they were computed from), then one
grouping over all scanned contigs in the given order. Before the grouping the
adapter test (``adapter.py``) sets aside places whose segments hold a
sequencing-adapter core; it reads the places only, so saved scan pieces are
reused with it on or off. The family-call layer reads the written places and
segments tables, so a detection can be reused with ``--detection-dir`` (its
grouping, and so its adapter test, stands as it was). The single-copy rule
(``single_copy.py``) acts on the family calls; its contig depth profile is kept
in the work directory with its identity. The source run directory is only read.
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import sys
import time
from collections import Counter as Tally
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

from ..contigs import MITO_CONTIG_NAMES, default_primary_contigs, find_mito_contig
from .. import vcf as V
from ..errors import InputContractError, OperationalResourceError
from ..outputs import write_bytes, write_run_summary
from ..resident_filter import RESIDENT_DROPPED_OUTPUT, RESIDENT_SUMMARY_KEY, dropped_tsv_bytes
from ..shadow_layers import ShadowCallset
from . import NEAR_EVITEA_BP, RULE_C_MAX, T_E
from .options import Options

__all__ = ["MIN_DEFAULT_CONTIG", "combined_calls_tsv_bytes", "default_contigs", "run_library_free"]

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
CODE_FILES = ("_p4_lib.py", "detector.py", "place.py", "scan.py", "walks.py", "extension.py", "priority.py", "options.py", "cross.py")
MERGE_VIEW = "merge_view"
REPEAT_PROFILE = "repeat_profile.json"
COPY_PROFILE = "copy_profile.json"
COPY_SCHEMA = "evitea.library_free.copy_profile.v1"


def default_contigs(lengths: Mapping[str, int]) -> tuple[str, ...]:
    """The contigs scanned without ``--contigs`` (``{name: length}`` in reference order), in reference order: every primary
    chromosome (``contigs.default_primary_contigs``: at least a tenth of the longest contig, not an organelle, not an unplaced
    sequence) and every contig of at least ``MIN_DEFAULT_CONTIG``. The list differs from the contigs of at least 1 Mb only when
    the longest contig is under 10 Mb, so a genome of short chromosomes (yeast) is scanned whole and the lists of dm6 and mm10
    are those of the contigs of at least 1 Mb."""

    primary = set(default_primary_contigs(lengths))
    return tuple(name for name, n in lengths.items() if n >= MIN_DEFAULT_CONTIG or name in primary)


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


def _repeat_cutoff(opts: Options, genome: dict, rm, work: Path, identity: dict, progress):
    """(cutoff, summary) of ``--repeat-cutoff`` (``cross.py``); the profile is kept in the work directory with its identity."""
    from . import _p4_lib as P

    if opts.repeat_cutoff == "fixed":
        return P.REP_MAX, {"mode": "fixed", "cutoff": P.REP_MAX}
    from . import cross as X

    path = work / REPEAT_PROFILE
    if path.is_file():
        saved = json.loads(path.read_text(encoding="utf-8"))
        if saved.get("identity") == identity:
            return saved["summary"]["cutoff"], saved["summary"]
    if rm is None:
        return None, None
    progress("library-free: repeat content of the reference windows (repeat cutoff)")
    cut, summary = X.repeat_cutoff(X.repeat_profile(genome, rm))
    summary = {"mode": "genome-share", **summary}
    path.write_text(json.dumps({"identity": identity, "summary": summary}, indent=2), encoding="utf-8")
    return cut, summary


def _detect(*, reference_fasta: Path, genome: dict, msbwt_index: Path, contigs: Sequence[str] | None, work: Path,
            progress: Callable[[str], None], opts: Options, mito_contig: str | None = None, adapter_test: bool = True):
    """Scan the contigs (resumable), then group each pass over all scanned contigs (after the adapter test, unless
    ``adapter_test`` is False or every rule switch is off). Return the tables of the first pass
    (places bytes, segments bytes), of the second scan (or None), with a cross-species rule on also of the merge's first pass
    (``merge_first``), and an info dict."""

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
        contigs = default_contigs({name: len(seq) for name, seq in genome.items()})
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
    cutoff = cutoff_info = None
    if opts.any_cross:
        from . import cross as X

        cut_identity = {"code_sha256": identity_base["code_sha256"], "reference_sha256": identity_base["reference_sha256"],
                        "anchor": X.REPEAT_ANCHOR, "min_contig": X.MIN_PROFILE_CONTIG}
        identity_base["repeat_cutoff"] = {"mode": opts.repeat_cutoff, **cut_identity}
        cutoff, cutoff_info = _repeat_cutoff(opts, genome, None, work, cut_identity, progress)
    timings: dict[str, float] = {}
    first_all: list[dict] = []
    first_merge_all: list[dict] = []
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
            if opts.any_cross and cutoff is None:
                cutoff, cutoff_info = _repeat_cutoff(opts, genome, rm, work, cut_identity, progress)
                progress(f"library-free: repeat cutoff {cutoff}")
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
                first, second, info = scan_contig_passes(counter, genome[chrom], chrom, rm, firewall, opts, words, cutoff)
            first_merge = info.pop("first_merge", None)
            info = {k: v for k, v in info.items() if k not in ("walk_log", "cross_log")}
            piece = dict(first=first, second=second, info=info)
            if first_merge is not None:
                piece["first_merge"] = first_merge
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
        if piece.get("first_merge") is not None:
            first_merge_all.extend(piece["first_merge"])
        if piece.get("second") is not None:
            second_all.extend(piece["second"])
    adapter = None
    if adapter_test and not opts.is_previous:
        from .adapter import mark_adapter_places

        adapter = {"first": mark_adapter_places(first_all)}
        if first_merge_all:
            adapter["merge_first"] = mark_adapter_places(first_merge_all)
        if second_all:
            adapter["second"] = mark_adapter_places(second_all)
        progress("library-free: adapter test: " + ", ".join(f"{k} {n}" for k, n in adapter.items()) + " places set aside")
    progress("library-free: grouping")
    if opts.is_previous:
        grouped = D.group_and_call(first_all)
    else:
        from .same_break import set_aside_same_break

        same = {"first": set_aside_same_break(first_all, opts)}
        grouped = PR.group_and_call(first_all, opts)
    tables = {"first": (TB.places_tsv_bytes(first_all, opts), TB.segments_fa_bytes(grouped))}
    if opts.any_cross and opts.merge_view:
        om = opts.without_cross()
        same["merge_first"] = set_aside_same_break(first_merge_all, om)
        grouped_m = PR.group_and_call(first_merge_all, om)
        tables["merge_first"] = (TB.places_tsv_bytes(first_merge_all, om), TB.segments_fa_bytes(grouped_m))
    if opts.het_scan != "off":
        same["second"] = set_aside_same_break(second_all, opts)
        grouped2 = PR.group_and_call(second_all, opts, first=False)
        tables["second"] = (TB.places_tsv_bytes(second_all, opts, second=True), TB.segments_fa_bytes(grouped2))
    if not opts.is_previous:
        progress("library-free: one place per break: " + ", ".join(f"{k} {n}" for k, n in same.items()) + " places set aside")
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
    if "merge_first" in tables:
        info["merge_first"] = {"stages": dict(sorted(Tally(p["stage"] for p in first_merge_all).items())),
                               "grouped_places": len(grouped_m)}
    if opts.any_cross:
        info["repeat_cutoff"] = cutoff_info
        info["cross_rules"] = {
            "places_by_rule": {f"{r}:{st}": n for (r, st), n in sorted(Tally((p.get("cross_rule", ""), p["stage"])
                                                                               for p in first_all if p.get("cross_rule")).items())},
            "per_contig": {chrom: stats[chrom].get("cross_rules") for chrom in contigs}}
    info["mito_contig"] = mito
    if not opts.is_previous:
        info["adapter_test"] = {"places": adapter} if adapter is not None else "off"
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


def _single_copy(choice, genome: dict, reference_fasta: Path, msbwt_index: Path | None, contigs: Sequence[str] | None, work: Path,
                 progress: Callable[[str], None]):
    """``single_copy.Decision`` of ``--single-copy-contigs``. For auto, the depth profile of the default contigs and the run's
    contigs, kept in the work directory (``copy_profile.json``) with the identity of the reference and the msBWT."""
    from . import single_copy as SC

    if choice.mode == "list":
        missing = [c for c in choice.names if c not in genome]
        if missing:
            raise InputContractError(f"--single-copy-contigs names contigs that are not in the reference: {missing}")
    if choice.mode != "auto":
        return SC.Decision(choice)
    if msbwt_index is None:
        return SC.Decision(choice, reason="no msBWT index: the contig depth is not measured")
    wanted = set(default_contigs({name: len(seq) for name, seq in genome.items()})) | set(contigs or ())
    measured = [c for c in genome if c in wanted]
    about = Path(msbwt_index) / "about.npy"
    identity = {"schema": COPY_SCHEMA, "reference_sha256": _sha256(reference_fasta),
                "msbwt_about_sha256": _sha256(about) if about.is_file() else None, "msbwt_index": str(msbwt_index),
                "words": SC.N_WORDS, "word_length": SC.K, "contigs": measured}
    path = work / COPY_PROFILE
    if path.is_file():
        saved = json.loads(path.read_text(encoding="utf-8"))
        if saved.get("identity") == identity:
            progress("library-free: single-copy contigs: depth profile reused from the work directory")
            return SC.Decision(choice, profile=saved["profile"])
    from ..backends.msbwt import load_msbwt
    from . import _p4_lib as P

    started = time.time()
    progress(f"library-free: single-copy contigs: depth of {len(measured)} contigs")
    counter = P.Counter(load_msbwt(msbwt_index), cache=False)
    profile = SC.measure(genome, counter.cnt, measured, progress=progress)
    profile["seconds"] = round(time.time() - started, 1)
    work.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps({"identity": identity, "profile": profile}, indent=2), encoding="utf-8")
    tmp.replace(path)
    low = [c for c, m in profile["contigs"].items() if m["ratio"] is not None and m["words"] >= SC.MIN_WORDS
           and m["ratio"] <= SC.DEPTH_RATIO_MAX]
    progress(f"library-free: single-copy contigs: genome depth {profile['genome_depth']}, low depth: {', '.join(low) or 'none'}"
             f" ({profile['seconds']} s)")
    return SC.Decision(choice, profile=profile)


def run_library_free(*, run_dir: Path, reference, msbwt_index: Path | None, output: Path, work_dir: Path | None,
                     contigs: Sequence[str] | None, rule_c: bool, detection_dir: Path | None,
                     command_line: Sequence[str] = (), source_identity: dict | None = None,
                     progress: Callable[[str], None] | None = None, options: Options | None = None,
                     repeat_annotation: Path | None = None, repeat_map: Path | None = None,
                     mito_contig: str | None = None, resident_filter=None, adapter_test: bool = True,
                     single_copy="auto") -> dict:
    """``resident_filter`` (``resident_filter.OwnFamilyCopyFilter``): drop the library calls inside a reference copy of their own
    family from ``calls.combined.tsv`` and ``calls.combined.vcf`` only.  Detection and assembly read the unfiltered
    ``calls.native.tsv`` of ``run_dir``, so the family calls and the added calls do not depend on it; ``None``: no filter.
    ``adapter_test``: set aside the places whose segments hold a sequencing-adapter core before the grouping (``adapter.py``);
    not run with ``detection_dir``, whose grouping stands. ``single_copy``: ``--single-copy-contigs`` (text or
    ``single_copy.Choice``); the merge view takes the main view's single-copy contigs."""

    from ..breakpoints.runtime import read_fasta
    from ..pipeline import prepare_output_directory

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
        if opts.any_cross and opts.merge_view and (detection_dir / MERGE_VIEW / PLACES).is_file():
            tables["merge_first"] = ((detection_dir / MERGE_VIEW / PLACES).read_bytes(), (detection_dir / MERGE_VIEW / SEGMENTS).read_bytes())
        info = {"reused_detection": str(detection_dir)}
        if not opts.is_previous:
            info["adapter_test"] = "as in the reused detection"
    else:
        if msbwt_index is None:
            raise InputContractError("library-free detection needs an msBWT index")
        genome = read_fasta(reference_fasta)
        tables, info = _detect(reference_fasta=reference_fasta, genome=genome, msbwt_index=Path(msbwt_index), contigs=contigs,
                               work=work, progress=say, opts=opts, mito_contig=mito_contig, adapter_test=adapter_test)
    if ("second" in tables) != (opts.het_scan != "off") and detection_dir is not None:
        raise InputContractError(
            "the detection directory " + ("holds" if "second" in tables else "has no") + " second-scan tables but --het-scan is "
            + opts.het_scan)
    from . import single_copy as SC

    choice = single_copy if isinstance(single_copy, SC.Choice) else SC.Choice.parse(single_copy)
    rule_on = not opts.is_previous and choice.mode != "none"
    if genome is None and (opts.check_scope != "none" or opts.resident_guard != "off" or rule_on):
        genome = read_fasta(reference_fasta)
    native_rows = ShadowCallset.from_bytes(native).rows
    copies, genome_bp = (_load_copies(opts, reference, repeat_annotation, genome, say, repeat_map, len(native_rows))
                         if genome is not None else (None, 0))
    if repeat_map is not None and (repeat_annotation is not None or native_rows):
        say("library-free: --repeat-map not used: " + ("--repeat-annotation is given" if repeat_annotation is not None
                                                      else f"the library layer made {len(native_rows)} calls"))
    decision = None
    if not opts.is_previous:
        decision = _single_copy(choice, genome, reference_fasta, msbwt_index, contigs, work, say) if rule_on else SC.Decision(choice)
    common = dict(native=native, native_rows=native_rows, sample_id=sample_id, reference=reference, genome=genome, copies=copies,
                  genome_bp=genome_bp, rule_c=rule_c, run_dir=run_dir, command_line=command_line, source_identity=source_identity,
                  repeat_map=repeat_map, repeat_annotation=repeat_annotation, say=say, started=started,
                  resident_filter=resident_filter)
    views = [(output, tables["first"], opts, info)]
    if "merge_first" in tables:
        merge_dir = output / MERGE_VIEW
        merge_dir.mkdir(parents=True, exist_ok=True)
        views.append((merge_dir, tables["merge_first"], opts.without_cross(), _merge_view_info(info)))
        info = {k: v for k, v in info.items() if k != "merge_first"}
        views[0] = (output, tables["first"], opts, info)
    summary = None
    for out_dir, first_tables, view_opts, view_info in views:
        result = _assemble_and_write(out_dir, first_tables, tables.get("second"), view_opts, view_info, single_copy=decision,
                                     **common)
        if summary is None:
            summary = result
            if decision is not None:
                decision = SC.Decision.from_view(choice, result["assembly"]["single_copy"].get("flagged", ()), "main view")
        else:
            summary["merge_view"] = {"output": str(out_dir), "family_calls": result["family_calls"],
                                     "family_call_statuses": result["family_call_statuses"]}
            write_run_summary(output / SUMMARY, summary)
    return summary


def _merge_view_info(info: dict) -> dict:
    """The detection summary of the merge view: the merge's first-pass stages and grouped places, no cross-species entries."""
    out = {k: v for k, v in info.items() if k not in ("repeat_cutoff", "cross_rules", "merge_first")}
    out.update(info.get("merge_first", {}))
    if "scan_counts" in out:
        out["scan_counts"] = {c: {k: v for k, v in st.items() if k != "cross_rules"} for c, st in out["scan_counts"].items()}
    return out


def _assemble_and_write(output: Path, first_tables, second_tables, opts: Options, info: dict, *, native, native_rows, sample_id,
                        reference, genome, copies, genome_bp, rule_c, run_dir, command_line, source_identity, repeat_map,
                        repeat_annotation, say, started, resident_filter=None, single_copy=None) -> dict:
    """Family calls of one set of rules from the detection tables, and the files of ``library-free`` in ``output``."""
    from . import assemble as AS
    from . import calling as CL
    from . import detector as D

    write_bytes(output / PLACES, first_tables[0])
    write_bytes(output / SEGMENTS, first_tables[1])
    if second_tables is not None:
        write_bytes(output / SECOND_PLACES, second_tables[0])
        write_bytes(output / SECOND_SEGMENTS, second_tables[1])
    say("library-free: family calls (reference-path filter" + (", family-identity filter" if rule_c else "") + ")"
        + ("" if output.name != MERGE_VIEW else ", merge rules"))
    first = (D.parse_places(first_tables[0]), D.parse_segments(first_tables[1]))
    second = (D.parse_places(second_tables[0]), D.parse_segments(second_tables[1])) if second_tables is not None else None
    c_started = time.time()
    index = CL.CallIndex.from_native(native_rows)
    rows, assembly = AS.assemble(first, second, index, opts, identity_filter=rule_c, reference=genome, copies=copies, genome_bp=genome_bp,
                                 single_copy=single_copy)
    c_seconds = round(time.time() - c_started, 1)
    extended = not opts.is_previous
    added = AS.added_calls_tsv_bytes(rows)
    write_bytes(output / FAMILY_CALLS, AS.family_calls_tsv_bytes(rows, extended=extended))
    write_bytes(output / ADDED, added)
    shown, resident = native, None
    if resident_filter is not None:  # final output only: `native` (the full table) is what the assembly above read
        resident = resident_filter(native)
        shown = resident.filtered
        write_bytes(output / RESIDENT_DROPPED_OUTPUT, dropped_tsv_bytes(resident.dropped))
    write_bytes(output / COMBINED, combined_calls_tsv_bytes(shown, added))
    write_bytes(output / V.COMBINED_VCF_OUTPUT, combined_calls_vcf_bytes(shown, added, reference, sample_id))
    statuses = Tally(r.status for r in rows)
    names = [PLACES, SEGMENTS, FAMILY_CALLS, ADDED, COMBINED, V.COMBINED_VCF_OUTPUT] + ([SECOND_PLACES, SECOND_SEGMENTS] if second is not None else [])
    if resident is not None:
        names.append(RESIDENT_DROPPED_OUTPUT)
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
    if resident is not None:
        summary[RESIDENT_SUMMARY_KEY] = resident_filter.summary(resident)
    if repeat_map is not None:
        summary["repeat_map"] = {"path": str(repeat_map), "used": repeat_annotation is None and not native_rows and copies is not None,
                                 "library_calls": len(native_rows)}
    write_run_summary(output / SUMMARY, summary)
    say(f"library-free: {len(rows)} family calls, {statuses.get('ADDED', 0)} added")
    return summary
