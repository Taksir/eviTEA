"""Increment in observable exact junction evidence from internal 21-mers, under fixed proposals."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from ..constants import EVENT_EQUIV_DISTANCE, N_FOLDS
from ..errors import InputContractError
from ..full_memory.common import ROTATIONS, as_int, normfam, rotation_label
from ..full_memory.strict_a import BIndex, a_groups, exact_anchors, score_b
from ..shadow_layers import Internal21Instance, Internal21Locus

__all__ = ["IncrementalLocus", "observe_increment", "novel_unconfounded_loci", "INSTANCE_FIELDS", "LOCUS_FIELDS"]

INSTANCE_FIELDS = ("sample", "analysis_arm", "rotation", "family", "chrom", "breakpoint0", "anchor_host_side",
                   "proposal_id", "B_state", "A_pair_ids", "internal_A_pair_ids", "left_B_pair_ids", "right_B_pair_ids")
LOCUS_FIELDS = ("sample", "analysis_arm", "family", "chrom", "locus_lo", "locus_hi", "exact_key_count", "rotation_count",
                "rotations", "B_states", "A_pair_count", "internal_A_pair_count", "B_pair_count", "classification")


@dataclass
class IncrementalLocus:
    sample: str
    family: str
    chrom: str
    lo: int
    hi: int
    exact_keys: set = field(default_factory=set)
    rotations: set = field(default_factory=set)
    states: set = field(default_factory=set)
    a_pair_ids: set = field(default_factory=set)
    b_pair_ids: set = field(default_factory=set)
    internal_a_pair_ids: set = field(default_factory=set)
    classification: str = ""
    instances: list = field(default_factory=list)

    def ledger_row(self, arm: str) -> dict:
        return {"sample": self.sample, "analysis_arm": arm, "family": self.family, "chrom": self.chrom,
                "locus_lo": self.lo, "locus_hi": self.hi, "exact_key_count": len(self.exact_keys),
                "rotation_count": len(self.rotations), "rotations": ";".join(sorted(self.rotations)),
                "B_states": ";".join(sorted(self.states)), "A_pair_count": len(self.a_pair_ids),
                "internal_A_pair_count": len(self.internal_a_pair_ids), "B_pair_count": len(self.b_pair_ids),
                "classification": self.classification}


def _key(e: Mapping) -> tuple[str, str, int, str]:
    return (normfam(e["family"]), str(e["chrom"]), int(e["breakpoint0"]), str(e["anchor_host_side"]))


def _cover(proposals: list[Mapping], event: Mapping) -> Mapping | None:
    family, chrom, bp = normfam(event["family"]), event["chrom"], int(event["breakpoint0"])
    hits = [p for p in proposals if normfam(p["family"]) == family and p["chrom"] == chrom
            and as_int(p["candidate_lo"]) <= bp <= as_int(p["candidate_hi"])]
    if not hits:
        return None
    return min(hits, key=lambda p: (as_int(p["candidate_hi"]) - as_int(p["candidate_lo"]), as_int(p["candidate_lo"]), p["candidate_id"]))


def _certify(witnesses: list[Mapping], index: BIndex, heldout: set[int], proposals: list[Mapping],
             internal_ids: set[int]) -> list[dict]:
    out = []
    for anchor in exact_anchors(a_groups(witnesses), heldout):
        b = score_b(anchor, index, heldout, exclude_a=True)
        if b["B_state"] == "NEITHER":
            continue
        event = {**anchor, **b}
        proposal = _cover(proposals, event)
        if proposal is None:
            continue
        a_ids = set(anchor["A_pair_ids"])
        b_ids = set(b["left_B_pair_ids"]) | set(b["right_B_pair_ids"])
        discovery = {as_int(x) for x in str(proposal["discovery_pair_ids"]).split(";") if x}
        if a_ids & b_ids:
            raise InputContractError("exact junction evidence and paired-end anchor support share a physical pair")
        if discovery & (a_ids | b_ids):
            raise InputContractError("pairs that discovered the proposal also appear in its confirming evidence")
        out.append({**event, "proposal_id": proposal["candidate_id"], "internal_A_pair_ids": sorted(a_ids & internal_ids)})
    return out


def _components(sample: str, events: list[dict]) -> list[IncrementalLocus]:
    by: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for e in events:
        by[(normfam(e["family"]), e["chrom"])].append(e)
    loci: list[IncrementalLocus] = []

    def make(family: str, chrom: str, xs: list[dict]) -> IncrementalLocus:
        bps = [int(e["breakpoint0"]) for e in xs]
        return IncrementalLocus(
            sample, family, chrom, min(bps), max(bps), {_key(e) for e in xs}, {e["rotation"] for e in xs},
            {e["B_state"] for e in xs}, set().union(*(set(e["A_pair_ids"]) for e in xs)),
            set().union(*((set(e["left_B_pair_ids"]) | set(e["right_B_pair_ids"])) for e in xs)),
            set().union(*(set(e["internal_A_pair_ids"]) for e in xs)), "", list(xs))

    for (family, chrom), xs in by.items():
        xs.sort(key=lambda e: int(e["breakpoint0"]))
        current: list[dict] = []
        for e in xs:
            if not current or int(e["breakpoint0"]) - max(int(z["breakpoint0"]) for z in current) <= EVENT_EQUIV_DISTANCE:
                current.append(e)
            else:
                loci.append(make(family, chrom, current))
                current = [e]
        if current:
            loci.append(make(family, chrom, current))
    return loci


def _classify(loci: list[IncrementalLocus], prior: list[tuple[str, str, int]]) -> None:
    for x in loci:
        represented = any(f == x.family and c == x.chrom and x.lo - EVENT_EQUIV_DISTANCE <= bp <= x.hi + EVENT_EQUIV_DISTANCE
                          for f, c, bp in prior)
        x.classification = "PRIOR_SAME_FAMILY_REPRESENTATION" if represented else "NOVEL_UNCONFOUNDED"
    novel = [x for x in loci if x.classification == "NOVEL_UNCONFOUNDED"]
    for i, a in enumerate(novel):
        for b in novel[i + 1:]:
            if a.chrom != b.chrom or a.family == b.family:
                continue
            if a.lo <= b.hi + EVENT_EQUIV_DISTANCE and b.lo <= a.hi + EVENT_EQUIV_DISTANCE:
                a.classification = b.classification = "UNRESOLVED_CROSS_FAMILY_LOCAL_COMPETITION"


def observe_increment(
    sample: str,
    rotation_proposals: Iterable[Mapping],
    old_witnesses: list[Mapping],
    internal_witnesses: list[Mapping],
    index: BIndex,
    prior_events: Iterable[Mapping],
) -> tuple[list[dict], list[IncrementalLocus]]:
    """Incremental exact events per rotation and their classified loci.

    ``prior_events`` is the full-memory event table; its validated events
    (``validated_anchorpe == 1``) are the prior same-family representation.
    """

    by_rotation: dict[str, list[Mapping]] = defaultdict(list)
    for p in rotation_proposals:
        if p["sample"] == sample:
            by_rotation[p["rotation"]].append(p)
    internal_ids = {as_int(w["pair_id"]) for w in internal_witnesses}
    union = list(old_witnesses) + list(internal_witnesses)
    events: list[dict] = []
    for held in ROTATIONS:
        label, folds = rotation_label(held), set(held)
        props = by_rotation[label]
        old = {_key(e) for e in _certify(list(old_witnesses), index, folds, props, set())}
        for e in _certify(union, index, folds, props, internal_ids):
            if _key(e) not in old:
                events.append({**e, "sample": sample, "rotation": label})
    prior = [(normfam(r["family"]), r["chrom"], as_int(r["breakpoint0"])) for r in prior_events
             if r.get("sample") == sample and as_int(r.get("validated_anchorpe")) == 1]
    loci = _components(sample, events)
    _classify(loci, prior)
    return events, loci


def novel_unconfounded_loci(loci: Iterable[IncrementalLocus]) -> tuple[Internal21Locus, ...]:
    """Typed integration input: every NOVEL_UNCONFOUNDED locus with the
    incremental instances of its family and chromosome."""

    loci = list(loci)
    out = []
    for x in loci:
        if x.classification != "NOVEL_UNCONFOUNDED":
            continue
        members = [e for y in loci if y.family == x.family and y.chrom == x.chrom for e in y.instances]
        out.append(Internal21Locus(x.family, x.chrom, x.lo, x.hi, tuple(
            Internal21Instance(e["rotation"], int(e["breakpoint0"]), e["anchor_host_side"], tuple(e["A_pair_ids"]),
                               tuple(e["left_B_pair_ids"]), tuple(e["right_B_pair_ids"])) for e in members)))
    return tuple(out)
