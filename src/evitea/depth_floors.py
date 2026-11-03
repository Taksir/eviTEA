"""Opt-in depth-aware floors of the library layer (``--weak-site-pass``, ``--host-cluster-support``, ``--candidate-min-pairs``).

All off by default; with none given nothing here runs and every output is the freeze's.  The depth D is the sample's single-copy
depth of the library-free layer (``library_free.single_copy.measure``, ``genome_depth``: the median sample count of words that
occur once in the reference, over the default contigs), measured the same way, so one number decides the library-free and the
library-layer floors.  It is computed only when a depth cut is given, and kept in the work directory.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .config import DepthFloorConfig
from .errors import InputContractError

__all__ = ["DEFAULT_SUPPORT", "ResolvedFloors", "resolve_floors", "sample_depth"]

DEFAULT_SUPPORT = 2     # host_mate_gate.HOST_CLUSTER_SUPPORT and candidates.CANDIDATE_MIN_PAIRS


@dataclass(frozen=True, slots=True)
class ResolvedFloors:
    """The floors this run uses, after the depth cuts."""

    depth: float | None
    host_cluster_support: int
    candidate_min_pairs: int
    weak_site_pass: bool        # the pass is on and D is at or under its cut

    def record(self, config: DepthFloorConfig) -> dict[str, object]:
        out: dict[str, object] = {"single_copy_depth": self.depth}
        if config.weak_site_pass:
            out["weak_site_pass"] = {"depth_cut": config.weak_site_depth_cut, "ran": self.weak_site_pass}
        if (config.host_cluster_support, config.candidate_min_pairs) != (DEFAULT_SUPPORT, DEFAULT_SUPPORT) \
                or config.candidate_floor_depth_cut is not None:
            out["candidate_floor"] = {
                "host_cluster_support_requested": config.host_cluster_support,
                "candidate_min_pairs_requested": config.candidate_min_pairs,
                "depth_cut": config.candidate_floor_depth_cut,
                "host_cluster_support_used": self.host_cluster_support,
                "candidate_min_pairs_used": self.candidate_min_pairs,
            }
        return out


def sample_depth(reference_fasta: Path, msbwt_index: Path, work: Path, progress: Callable[[str], None]) -> float:
    """``genome_depth`` of the sample, as the library-free layer measures it (default contigs; profile kept in ``work``)."""

    from .breakpoints.runtime import read_fasta
    from .library_free import single_copy as SC
    from .library_free.runtime import _single_copy

    decision = _single_copy(SC.Choice("auto"), read_fasta(Path(reference_fasta)), Path(reference_fasta), Path(msbwt_index),
                            None, Path(work), progress)
    depth = None if decision.profile is None else decision.profile.get("genome_depth")
    if depth is None:
        raise InputContractError("a depth cut needs the sample's single-copy depth, which could not be measured ("
                                 + (decision.reason or "no single-copy words counted") + ")")
    return float(depth)


def resolve_floors(config: DepthFloorConfig, reference_fasta: Path, msbwt_index: Path, work: Path,
                   progress: Callable[[str], None]) -> ResolvedFloors:
    depth = sample_depth(reference_fasta, msbwt_index, work / "depth_floors", progress) if config.needs_depth else None
    support, pairs = config.host_cluster_support, config.candidate_min_pairs
    if config.candidate_floor_depth_cut is not None and depth > config.candidate_floor_depth_cut:
        support = pairs = DEFAULT_SUPPORT
    weak = config.weak_site_pass and depth <= config.weak_site_depth_cut
    progress(f"depth floors: single-copy depth {depth}; host-cluster support {support}, minimum pairs per place {pairs}"
             + (f"; weak-site pass {'runs' if weak else 'skipped (depth above the cut)'}" if config.weak_site_pass else ""))
    return ResolvedFloors(depth, support, pairs, weak)
