"""Exact soft-clip anchors ("A") and q99 paired-end geometry support ("B")."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping

from ..constants import MIN_FOLDS, MIN_PAIRS, Q99
from ..shadow_layers import AnchorPECandidate
from .alignments import family_occurrences, terminal_softclips
from .codebook import CodeUnit, fold_of

__all__ = ["channel_a", "channel_b_universe", "score_candidates"]


def channel_a(
    targets: list[dict],
    pair_to_cluster: Mapping[int, int],
    ledger: Mapping[int, dict],
    te_records: Mapping[int, list[dict]],
    by_family: Mapping[str, list[CodeUnit]],
    q99: int = Q99,
) -> tuple[list[dict], list[dict]]:
    """Witnesses and one replicated anchor per target cluster.

    Only ``LOCAL_PROPER_TE_UNEXPLAINED`` pairs of a target cluster contribute.
    A terminal soft clip of a proper TE-mate record, on the cluster's
    chromosome, within ``[start - q99, end + q99]`` (q99: the fragment-size
    99th percentile, 805 bp by default), that overlaps a
    candidate-family safe-probe occurrence by >= 1 query base, witnesses the
    interbase breakpoint at its clip edge (leading clip gives RIGHT at record
    start; trailing clip gives LEFT at record end).  Anchors need >= 2 pairs
    over >= 2 folds.  Per cluster the anchor with most pairs, then folds, then
    lower breakpoint, then LEFT is retained.
    """

    target_by = {int(c["cluster_id"]): c for c in targets}
    witnesses: list[dict] = []
    for pid, meta in ledger.items():
        cid = pair_to_cluster.get(pid)
        if cid not in target_by or meta["reference_state"] != "LOCAL_PROPER_TE_UNEXPLAINED":
            continue
        cluster = target_by[cid]
        family = meta["family"]
        lo, hi = max(0, int(cluster["start"]) - q99), int(cluster["end"]) + q99
        seen: set[tuple] = set()
        for rec in te_records.get(pid, []):
            if rec["chrom"] != cluster["chrom"]:
                continue
            seq = rec["seq"]
            occ = family_occurrences(seq, by_family[family])
            left, right = terminal_softclips(rec["cigar"], len(seq))
            cases = []
            if left is not None:
                cases.append((left, "LEADING", "RIGHT", int(rec["start"])))
            if right is not None:
                cases.append((right, "TRAILING", "LEFT", int(rec["end"])))
            for (a, b), clip_side, host_side, bp in cases:
                if not lo <= bp <= hi:
                    continue
                units = sorted({uid for s, e, uid in occ if max(s, a) < min(e, b)})
                if not units:
                    continue
                key = (cid, rec["chrom"], bp, family, host_side, pid)
                if key in seen:
                    continue
                seen.add(key)
                witnesses.append({
                    "cluster_id": cid, "pair_id": pid, "fold": fold_of(pid), "family": family,
                    "chrom": rec["chrom"], "breakpoint0": bp, "anchor_host_side": host_side,
                    "cigar": rec["cigar"], "alignment_strand": rec["strand"], "te_softclip_side": clip_side,
                    "te_units_in_softclip": ";".join(map(str, units)),
                })
    groups: dict[tuple, set[int]] = defaultdict(set)
    for w in witnesses:
        groups[(w["cluster_id"], w["chrom"], w["breakpoint0"], w["family"], w["anchor_host_side"])].add(int(w["pair_id"]))
    by_cluster: dict[int, list[dict]] = defaultdict(list)
    for (cid, chrom, bp, family, side), pids in groups.items():
        folds = sorted({fold_of(p) for p in pids})
        if len(pids) >= MIN_PAIRS and len(folds) >= MIN_FOLDS:
            by_cluster[cid].append({
                "cluster_id": cid, "chrom": chrom, "breakpoint0": bp, "family": family,
                "anchor_host_side": side, "channel_A_pairs": len(pids), "channel_A_folds": len(folds),
                "channel_A_pair_ids": ";".join(map(str, sorted(pids))),
                "channel_A_fold_ids": ";".join(map(str, folds)),
            })
    anchors = []
    for cid in sorted(by_cluster):
        rows = sorted(by_cluster[cid], key=lambda r: (-int(r["channel_A_pairs"]), -int(r["channel_A_folds"]),
                                                      int(r["breakpoint0"]), 0 if r["anchor_host_side"] == "LEFT" else 1))
        anchors.append(rows[0])
    return witnesses, anchors


def channel_b_universe(ledger: Mapping[int, dict], host: Mapping[int, dict]) -> list[dict]:
    """``NO_LOCAL_PROPER`` fragments with a ``UNIQUE_HIGH`` opposite host mate."""

    rows = []
    for pid, meta in ledger.items():
        h = host.get(pid)
        if meta["reference_state"] == "NO_LOCAL_PROPER" and h and h["host_class"] == "UNIQUE_HIGH":
            rows.append({"pair_id": pid, "fold": fold_of(pid), "family": meta["family"], "chrom": h["chrom"],
                         "start": h["start"], "end": h["end"], "strand": h["strand"], "mapq": h["mapq"]})
    return rows


def score_candidates(
    anchors: list[dict], b_rows: list[dict], q99: int = Q99
) -> tuple[list[dict], tuple[AnchorPECandidate, ...]]:
    """Directional q99 completion at each exact anchor (A).

    ``+`` host mate: LEFT side when ``end <= bp <= end + q99``; ``-``: RIGHT
    side when ``start - q99 <= bp <= start`` (q99 805 bp by default).  The
    union of all A pairs is excluded from B.
    """

    a_union = {int(x) for a in anchors for x in str(a["channel_A_pair_ids"]).split(";") if x}
    index: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in b_rows:
        if int(r["pair_id"]) not in a_union:
            index[(r["family"], r["chrom"])].append(r)
    rows: list[dict] = []
    typed: list[AnchorPECandidate] = []
    for a in anchors:
        bp = int(a["breakpoint0"])
        left, right = [], []
        for r in index[(a["family"], a["chrom"])]:
            s, e = int(r["start"]), int(r["end"])
            if r["strand"] == "+" and e <= bp <= e + q99:
                left.append(r)
            elif r["strand"] == "-" and s - q99 <= bp <= s:
                right.append(r)
        lp = sorted({int(r["pair_id"]) for r in left})
        rp = sorted({int(r["pair_id"]) for r in right})
        lf, rf = sorted({fold_of(p) for p in lp}), sorted({fold_of(p) for p in rp})
        bilateral = len(lp) >= MIN_PAIRS and len(lf) >= MIN_FOLDS and len(rp) >= MIN_PAIRS and len(rf) >= MIN_FOLDS
        rows.append({**a, "left_B_pairs": len(lp), "left_B_folds": len(lf), "right_B_pairs": len(rp),
                     "right_B_folds": len(rf), "distinct_B_pairs": len(set(lp) | set(rp)),
                     "left_B_pair_ids": ";".join(map(str, lp)), "right_B_pair_ids": ";".join(map(str, rp)),
                     "bilateral_B_2p2fold": int(bilateral)})
        typed.append(AnchorPECandidate(
            cluster_id=int(a["cluster_id"]), family=a["family"], chrom=a["chrom"], breakpoint0=bp,
            anchor_host_side=a["anchor_host_side"],
            a_pair_ids=tuple(int(x) for x in str(a["channel_A_pair_ids"]).split(";") if x),
            left_b_pair_ids=tuple(lp), right_b_pair_ids=tuple(rp),
        ))
    return rows, tuple(typed)
