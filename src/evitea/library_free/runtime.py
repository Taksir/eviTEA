"""``evitea library-free``: detection over a sample msBWT, family calls, and the merge with a finished run.

Detection runs one contig at a time (resumable: each contig's places are kept
in the work directory with the identity they were computed from), then one
grouping over all scanned contigs in the given order.  The family-call layer
reads the written places and segments tables, so a detection can be reused
with ``--detection-dir``.  The source run directory is only read.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import sys
import time
from collections import Counter as Tally
from collections.abc import Callable, Sequence
from pathlib import Path

from ..errors import InputContractError, OperationalResourceError
from ..outputs import write_bytes, write_run_summary
from ..shadow_layers import ShadowCallset
from . import NEAR_EVITEA_BP, RULE_C_MAX, T_E

__all__ = ["MIN_DEFAULT_CONTIG", "combined_calls_tsv_bytes", "run_library_free"]

MIN_DEFAULT_CONTIG = 1_000_000
PLACES = "library_free_places.tsv"
SEGMENTS = "library_free_segments.fa"
FAMILY_CALLS = "library_free_family_calls.tsv"
ADDED = "calls.library_free.tsv"
COMBINED = "calls.combined.tsv"
SUMMARY = "library_free_summary.json"
PIECE_SCHEMA = "evitea.library_free.piece.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


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


def _detect(*, reference_fasta: Path, msbwt_index: Path, contigs: Sequence[str] | None, work: Path,
            progress: Callable[[str], None]) -> tuple[bytes, bytes, dict]:
    from ..backends.msbwt import load_msbwt
    from ..breakpoints.runtime import read_fasta
    from . import _p4_lib as P
    from . import detector as D

    try:
        import mappy  # noqa: F401
    except ImportError as exc:
        raise OperationalResourceError("library-free detection needs mappy (minimap2); run on Linux or WSL") from exc
    genome = read_fasta(reference_fasta)
    if contigs is None:
        contigs = tuple(name for name, seq in genome.items() if len(seq) >= MIN_DEFAULT_CONTIG)
    missing = [c for c in contigs if c not in genome]
    if missing or len(set(contigs)) != len(contigs):
        raise InputContractError(f"contigs must be distinct reference contigs; not found: {missing}")
    about = Path(msbwt_index) / "about.npy"
    identity_base = {
        "schema": PIECE_SCHEMA,
        "detector_sha256": _sha256(Path(P.__file__)),
        "reference_sha256": _sha256(reference_fasta),
        "msbwt_about_sha256": _sha256(about) if about.is_file() else None,
        "msbwt_index": str(msbwt_index),
    }
    pieces_dir = work / "pieces"
    pieces_dir.mkdir(parents=True, exist_ok=True)
    rm = None
    timings: dict[str, float] = {}
    all_places: list[dict] = []
    for chrom in contigs:
        identity = {**identity_base, "contig": chrom}
        data, marker = pieces_dir / f"{chrom}.pkl", pieces_dir / f"{chrom}.complete.json"
        places = None
        if data.is_file() and marker.is_file():
            recorded = json.loads(marker.read_text(encoding="utf-8"))
            if recorded.get("identity") == identity and _sha256(data) == recorded.get("sha256"):
                with data.open("rb") as handle:
                    places = pickle.load(handle)
                progress(f"library-free: {chrom} reused from the work directory ({len(places)} places)")
        if places is None:
            if rm is None:
                progress("library-free: reference multi-copy words")
                rm = D.reference_multi(genome)
            started = time.time()
            progress(f"library-free: scanning {chrom}")
            counter = P.Counter(load_msbwt(msbwt_index), cache=False)
            firewall = P.MappyFirewall(str(reference_fasta), genome.get("chrM"))
            places = D.scan_contig(counter, genome[chrom], chrom, rm, firewall)
            payload = pickle.dumps(places, protocol=pickle.HIGHEST_PROTOCOL)
            tmp = data.with_suffix(".tmp")
            tmp.write_bytes(payload)
            tmp.replace(data)
            marker.write_text(json.dumps({"identity": identity, "sha256": hashlib.sha256(payload).hexdigest()},
                                         indent=2), encoding="utf-8")
            timings[chrom] = round(time.time() - started, 1)
            progress(f"library-free: {chrom} done, {len(places)} places, {timings[chrom]} s, peak RSS {_peak_rss_gb()} GB")
        all_places.extend(places)
    progress("library-free: grouping")
    grouped = D.group_and_call(all_places)
    stages = Tally(p["stage"] for p in all_places)
    info = {"contigs": list(contigs), "scan_seconds": timings, "stages": dict(sorted(stages.items())),
            "grouped_places": len(grouped)}
    return D.places_tsv_bytes(all_places), D.segments_fa_bytes(grouped), info


def run_library_free(*, run_dir: Path, reference, msbwt_index: Path | None, output: Path, work_dir: Path | None,
                     contigs: Sequence[str] | None, rule_c: bool, detection_dir: Path | None,
                     command_line: Sequence[str] = (), source_identity: dict | None = None,
                     progress: Callable[[str], None] | None = None) -> dict:
    from ..pipeline import prepare_output_directory
    from . import calling as CL
    from . import detector as D

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
    if detection_dir is not None:
        detection_dir = Path(detection_dir)
        places_bytes = (detection_dir / PLACES).read_bytes()
        segments_bytes = (detection_dir / SEGMENTS).read_bytes()
        info = {"reused_detection": str(detection_dir)}
    else:
        if msbwt_index is None:
            raise InputContractError("library-free detection needs an msBWT index")
        places_bytes, segments_bytes, info = _detect(reference_fasta=Path(reference.reference_fasta_path),
                                                     msbwt_index=Path(msbwt_index), contigs=contigs, work=work,
                                                     progress=say)
    write_bytes(output / PLACES, places_bytes)
    write_bytes(output / SEGMENTS, segments_bytes)
    say("library-free: family calls (v2 filter" + (", rule c" if rule_c else "") + ")")
    places = D.parse_places(places_bytes)
    segments = D.parse_segments(segments_bytes)
    c_started = time.time()
    c_by_cluster = CL.family_similarity(places, segments)
    c_seconds = round(time.time() - c_started, 1)
    index = CL.CallIndex.from_native(ShadowCallset.from_bytes(native).rows)
    rows = CL.family_calls(places, index, c_by_cluster, rule_c=rule_c)
    added = CL.added_calls_tsv_bytes(rows, sample_id=sample_id)
    write_bytes(output / FAMILY_CALLS, CL.family_calls_tsv_bytes(rows))
    write_bytes(output / ADDED, added)
    write_bytes(output / COMBINED, combined_calls_tsv_bytes(native, added))
    statuses = Tally(r.status for r in rows)
    summary = {
        "sample_id": sample_id,
        "source_run": str(run_dir),
        "source_native_callset_sha256": hashlib.sha256(native).hexdigest(),
        "command_line": list(command_line),
        "source_identity": source_identity or {},
        "parameters": {"v2_reference_path_fraction_max": T_E, "rule_c": bool(rule_c), "rule_c_similarity_max": RULE_C_MAX,
                       "near_evitea_bp": NEAR_EVITEA_BP},
        "detection": info,
        "family_calls": len(rows),
        "family_call_statuses": dict(sorted(statuses.items())),
        "families_with_similarity": len(c_by_cluster),
        "rule_c_seconds": c_seconds,
        "seconds": round(time.time() - started, 1),
        "peak_rss_gb": _peak_rss_gb(),
        "output_hashes": {name: _sha256(output / name) for name in (PLACES, SEGMENTS, FAMILY_CALLS, ADDED, COMBINED)},
    }
    write_run_summary(output / SUMMARY, summary)
    say(f"library-free: {len(rows)} family calls, {statuses.get('ADDED', 0)} added")
    return summary
