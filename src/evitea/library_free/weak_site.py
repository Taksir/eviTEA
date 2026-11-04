"""Weak-site pass of the library-free step (opt-in: ``--weak-site-pass`` with ``--weak-site-allele-depth-cut DA``).

The pass of ``incumbent/weak_sites.py`` (pooled reads on both sides, or two split reads on one side, at a site no call covers) is
switched here by the depth of the *inserted allele*, not by the depth of the sample: ``Da = D x a`` with ``D`` the single-copy depth
(``single_copy.measure``, ``genome_depth``) and ``a = 1 - median(e)`` over the confirmed places of the main view (the places of the
local heterozygous guard, ``het_guard.confirmed``: first pass, status ADDED or NEAR_EVITEA, both junctions placed, TSD 2 to 20, a
measured e), with ``a = 1`` when there are fewer than ``het_guard.NMIN`` (10) of them. The pass runs only when ``Da`` is at or under the
cut. A pooled sample or a heterozygous one carries the inserted allele in a share of the reads, which ``D`` alone does not see.

The pass reads the library layer's saved stage tables (``split_junction``, ``mate_rescue``, ``junction_groups``) from the work
directory of the ``call`` run that made the run folder (named in its ``run_manifest.json``, else ``<run>/.evitea-work``), and the
results of the scan only. It is not part of the identity of a saved scan piece. It runs after the library-free calls and statuses are
final and before ``calls.combined.tsv`` and ``calls.combined.vcf`` are written. The 500 bp cover is every row the combined table
will hold (the library rows after the own-family copy filter, and the added library-free calls). The new rows go through the same
own-family copy filter and are kept as library rows (source EVITEA, ``evidence_type`` ``weak_site_pooled_floor``, tier 4). Breakpoint
refinement is not run on them: they keep the span of their evidence. Main view only; the merge view is unchanged.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import statistics
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..errors import InputContractError

__all__ = ["STAGES", "WeakSite", "allele_depth", "apply", "load_stage_tables", "stage_directory"]

STAGES = ("split_junction", "mate_rescue", "junction_groups")
ENGINE_CHECKPOINTS = Path("engine") / "work" / "incumbent" / "checkpoints"


def stage_directory(run_dir: Path) -> Path:
    """The checkpoint folder of the library layer of the ``call`` run in ``run_dir`` (from its run manifest)."""

    run_dir = Path(run_dir)
    manifest = run_dir / "run_manifest.json"
    try:
        values = json.loads(manifest.read_text(encoding="utf-8"))["values"]
    except (OSError, ValueError, KeyError) as exc:
        raise InputContractError(f"--weak-site-pass needs the run manifest of the library run: {manifest} is missing or malformed") from exc
    work = ((values.get("normalized_config") or {}).get("runtime") or {}).get("work_dir")
    return (Path(work) if work else run_dir / ".evitea-work") / ENGINE_CHECKPOINTS


def load_stage_tables(run_dir: Path):
    """(split_junction, mate_rescue, junction_groups) results saved by the library layer; an error naming what is missing."""

    folder = stage_directory(run_dir)
    out = []
    for stage in STAGES:
        data, marker = folder / f"{stage}.pkl", folder / f"{stage}.complete.json"
        if not data.is_file() or not marker.is_file():
            raise InputContractError(
                f"--weak-site-pass needs the library layer's saved stage table {stage}: {data} (and its .complete.json) is missing; "
                "the call run's work directory must be kept (not --no-library, not a deleted work directory)")
        payload = data.read_bytes()
        if hashlib.sha256(payload).hexdigest() != json.loads(marker.read_text(encoding="utf-8")).get("sha256"):
            raise InputContractError(f"--weak-site-pass: the saved stage table {data} does not match its completion marker")
        out.append(pickle.loads(payload))
    return tuple(out)


def allele_depth(depth: float, confirmed_e: Sequence[float], minimum: int) -> dict:
    """D, the number of confirmed places, their median e, a = 1 - median e (1 under ``minimum`` places) and Da = D x a."""

    n = len(confirmed_e)
    median = float(statistics.median(confirmed_e)) if n else None
    a = 1.0 if n < minimum else 1.0 - median
    return {"single_copy_depth": float(depth), "confirmed_places": n, "median_e": median, "a": a,
            "allele_depth": float(depth) * a, "a_rule": f"1 - median e; 1 under {minimum} confirmed places"}


@dataclass
class WeakSite:
    """What the main view needs: the cut, the single-copy depth, the saved stage tables and the sample id."""

    cut: float
    depth: float
    tables: tuple
    confirmed_e: list = field(default_factory=list)      # filled by ``assemble`` (the guard's confirmed places of the view)


def apply(weak: WeakSite, shown: bytes, added: bytes, resident_filter, sample_id: str):
    """(table the combined files are made from, dropped rows of the filter, summary entry).

    ``shown``: the library rows of the combined table (after the own-family copy filter); ``added``: ``calls.library_free.tsv``."""

    from ..incumbent.events import legacy_tsv_bytes
    from ..incumbent.weak_sites import weak_site_calls
    from ..shadow_layers import ShadowCallset
    from . import het_guard as HG

    info = allele_depth(weak.depth, weak.confirmed_e, HG.NMIN)
    info["cut"] = weak.cut
    ran = info["allele_depth"] <= weak.cut
    info["ran"] = ran
    if not ran:
        return shown, (), info
    base = ShadowCallset.from_bytes(shown)
    cover = [(row["chr"], int(row["start"]), int(row["end"])) for row in base.rows]
    lines = added.decode("utf-8").splitlines()
    if lines:
        head = lines[0].split("\t")
        for line in lines[1:]:
            if line.strip():
                r = dict(zip(head, line.split("\t")))
                cover.append((r["chr"], int(r["start"]), int(r["end"])))
    split, mates, groups = weak.tables
    found = weak_site_calls(sample_id, split, mates, groups, cover)
    rows = ShadowCallset.from_bytes(legacy_tsv_bytes(found.calls)).rows if found.calls else ()
    dropped = ()
    if found.calls and resident_filter is not None:
        result = resident_filter(legacy_tsv_bytes(found.calls))
        dropped = result.dropped
        gone = {d.event_id for d in dropped}
        rows = tuple(r for r in rows if r["event_id"] not in gone)
    info.update(sites=found.sites, uncovered_sites=found.uncovered, passing_rule=len(found.calls), dropped_by_copy_filter=len(dropped),
                rows_added=len(rows), cover_rows=len(cover))
    return (base.append(rows).data if rows else shown), dropped, info
