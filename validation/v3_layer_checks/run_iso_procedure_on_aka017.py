"""Validation-only probe: does the frozen ISO-1 AnchorPE transfer procedure,
run unchanged on AKA-017, reproduce AKA-017's frozen AnchorPE candidates?

Why: AKA-017's 91-row ``ANCHOR_PE_CANDIDATES.tsv`` (the direct input of the
bilateral and asymmetric layers) came from a long historical chain, whereas
ISO-1's came from one self-contained raw-to-candidates script
(``route_b_v3_anchor_pe_iso1_transfer_freeze.py``, also reused unchanged for the
GIM-024 blind run).  If the self-contained procedure reproduces the AKA table,
it is the implementable S20-S70 specification for those layers.

The frozen ISO helper is imported by SHA-256 (as the GIM-024 script does) and
only its sample bindings (pair count, read length, paths) change.  Every
scientific function and constant is used as frozen.  Outputs go to ``--work``;
the research tree is only read.  No truth/comparator/GIM file is opened.

Heavy: codebook projection on the 31.45 M-pair msBWT, one FASTQ pass per mate,
edlib whole-read decoding, three Bowtie2 runs.  Resumable: rerun the same
command; finished steps are reused from ``--work``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import pickle
import sys
import time
from collections import defaultdict
import os
from pathlib import Path


def _env_path(name):
    """Path from environment variable ``name``; None when it is not set."""
    value = os.environ.get(name)
    return Path(value) if value else None

import numpy as np

ISO_HELPER_SHA = "8c29f990c6f007a2a905d76b7a72bd0d8607f9cd0149415538020492c9f6d6b9"
PAIR_COUNT = 31_451_315
READ_LENGTH = 125
MSBWT_CORE = {
    "about.npy": "2d754c34e1b025a2baf22b9ecf7fd47fcee85d806630a788d6cad773b1276af0",
    "comp_msbwt.npy": "488071ca7df84aca41271ac74a6f78e0d9d227bac039bf8e9aea7cfa87100bc4",
    "totalCounts.npy": "ef001cf6ac38fbedf2a74eb4e8886461b8e92b8e2feabaaaea6f2b8bbdeb3f2b",
}
SHARED = {
    "reference": "eeb896ecb227ec5b2c779faa46f064cf781adddaae509c7f50b9d6503c194397",
    "te_library": "8dff7f6c25b357d64184f232eff7e533058492ca7b118a06dac7d54bce13a922",
    "repeatmasker": "29ddba016d1afa6fa9b7c0fd45882ac5a3931408779580a946053f9bbab11e04",
    "family_map": "6616b6305e6eac2f97336620376859151083a4b9ed30ba3bd64d7eead7c70714",
    "codebook": "d631f16a0252e4820fdfbe4f2822165595ebf726018618468602fea5dbfb3e00",
}
FROZEN_AKA = {
    "candidates": "reports/V3_OBSERVED_ANCHOR_BILATERAL_PE_COMPLETION_FALSIFICATION_20260920/ANCHOR_PE_CANDIDATES.tsv",
    "clusters": "reports/V3_EVIDENCE_COMPLETE_WHOLE_ACQUISITION_RERUN_20260919/EVIDENCE_COMPLETE_UNEXPLAINED_CLUSTERS.tsv",
    "witnesses": "reports/V3_CORRECTED_CROSS_CHANNEL_BREAKPOINT_FALSIFICATION_20260920/CROSS_CHANNEL_BREAKPOINT_WITNESSES.tsv",
    "ledger": "reports/V3_EVIDENCE_COMPLETE_WHOLE_ACQUISITION_RERUN_20260919/EVIDENCE_COMPLETE_WHOLE_PAIR_LEDGER.tsv",
    "b_universe": "reports/V3_OBSERVED_ANCHOR_BILATERAL_PE_COMPLETION_FALSIFICATION_20260920/ANCHOR_PE_B_UNIVERSE.tsv",
}
FORBIDDEN = ("truth", "teforest", "temp2_calls", "gim024", "gim-024", "quarantined", "srr9951106")


def say(msg: str) -> None:
    print(f"[ISO->AKA {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def guard(path: Path) -> Path:
    if any(t in str(path).lower() for t in FORBIDDEN):
        raise RuntimeError(f"truth/comparator firewall blocked {path}")
    return path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with guard(path).open("rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def require(path: Path, expected: str, label: str) -> None:
    got = sha256(path)
    if got != expected:
        raise RuntimeError(f"{label} SHA mismatch: {path}\n got={got}\n expected={expected}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with guard(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def import_iso(path: Path):
    require(path, ISO_HELPER_SHA, "frozen ISO transfer helper")
    spec = importlib.util.spec_from_file_location(f"iso_helper_{time.time_ns()}", str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def done(path: Path) -> bool:
    return path.is_file() and Path(str(path) + ".done").is_file()


def mark(path: Path) -> None:
    Path(str(path) + ".done").write_text("complete\n", encoding="ascii")


def ids(raw: str) -> tuple[int, ...]:
    return tuple(sorted(int(x) for x in str(raw or "").replace(",", ";").split(";") if x.strip()))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--route-b-v3", type=Path, required=True)
    ap.add_argument("--data", type=Path, default=_env_path("EVITEA_DATA"), required=_env_path("EVITEA_DATA") is None)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--edlib-python", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--host-read-names", choices=("iso", "aka-historical"), default="iso",
                    help="Bowtie2 query names for single-mate host mapping: ISO/GIM transfer "
                         "'P{pid}' or the AKA historical chain's 'p{pid}|h{host_mate}'. Bowtie2 seeds "
                         "its per-read pseudo-random choices from the read name, so this changes "
                         "multi-mapping outcomes.")
    args = ap.parse_args(argv)

    v3, data, out = args.route_b_v3.resolve(), args.data.resolve(), args.work.resolve()
    out.mkdir(parents=True, exist_ok=True)
    iso_script = v3 / "route_b_v3_anchor_pe_iso1_transfer_freeze.py"
    iso = import_iso(iso_script)
    iso.PAIR_COUNT, iso.READ_LENGTH = PAIR_COUNT, READ_LENGTH
    assert (iso.Q99, iso.MIN_PAIRS, iso.MIN_FOLDS, iso.CLUSTER_GAP) == (805, 2, 2, 500)

    shared = data / "TEMP2_ISO1"
    msbwt = data / "AKA017" / "msbwt_SRR9951105"
    reference = shared / "Reference" / "dm6_TEMP2_compatible.fa"
    te_library = shared / "TE consensus library" / "transposon_sequence_set.fa"
    repeatmasker = shared / "dm6_repeatmasker_TE_TEMP2.bed6"
    codebook = v3 / "reports" / "V3_FINAL_SAFE_TERMINAL_CODEBOOK_20260918" / "FINAL_SAFE_TERMINAL_CODEBOOK.tsv"
    fammap = v3 / "fbte_to_flybase_family_FROZEN.tsv"
    fq1 = data / "AKA017" / "fastq" / "SRR9951105_1.fastq"
    fq2 = data / "AKA017" / "fastq" / "SRR9951105_2.fastq"

    binding = out / "INPUT_BINDING.json"
    if not binding.is_file():
        for p, h, label in [(reference, SHARED["reference"], "reference"), (te_library, SHARED["te_library"], "te_library"),
                            (repeatmasker, SHARED["repeatmasker"], "repeatmasker"), (fammap, SHARED["family_map"], "family_map"),
                            (codebook, SHARED["codebook"], "codebook")]:
            require(p, h, label)
        for name, h in MSBWT_CORE.items():
            require(msbwt / name, h, f"AKA msBWT {name}")
        say("hashing AKA FASTQs (recorded, not historically frozen)")
        binding.write_text(json.dumps({
            "fastq_r1_sha256": sha256(fq1), "fastq_r2_sha256": sha256(fq2),
            "msbwt_core": MSBWT_CORE, "shared": SHARED, "iso_helper_sha256": ISO_HELPER_SHA,
            "pair_count": PAIR_COUNT, "read_length": READ_LENGTH,
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    about = np.load(msbwt / "about.npy", mmap_mode="r", allow_pickle=False)
    units, by_family, _ = iso.load_codebook(codebook)

    hits_pkl = out / "pair_hits.pkl"
    if done(hits_pkl):
        pair_hits = pickle.loads(hits_pkl.read_bytes())
        say(f"reusing codebook projection: {len(pair_hits):,} pairs")
    else:
        bwt = iso.load_msbwt(msbwt)
        pair_hits = iso.project_codebook(bwt, about, units, out / "AKA017_SAFE_CODEBOOK_OBSERVATION.tsv")
        hits_pkl.write_bytes(pickle.dumps({k: {1: set(v[1]), 2: set(v[2])} for k, v in pair_hits.items()}))
        mark(hits_pkl)
    pair_rows, eligible = iso.classify_pairs(pair_hits, units)
    one_mate_meta = {}
    for pid, hits in pair_hits.items():
        nonempty = [m for m in (1, 2) if hits[m]]
        if len(nonempty) != 1:
            continue
        tm = nonempty[0]
        one_mate_meta[pid] = {"pair_id": pid, "te_mate": tm, "host_mate": 2 if tm == 1 else 1,
                              "unit_ids": sorted(hits[tm]), "old_primary": bool(pid in eligible and eligible[pid]["primary"])}
    say(f"hit pairs={len(pair_rows):,}; family-resolved one-mate={len(eligible):,}; decoder universe={len(one_mate_meta):,}")

    seqs_pkl = out / "selected_reads.pkl"
    if done(seqs_pkl):
        seqs = pickle.loads(seqs_pkl.read_bytes())
    else:
        say("streaming AKA FASTQs once per mate")
        seqs = iso.read_selected_pairs(fq1, fq2, set(one_mate_meta))
        seqs_pkl.write_bytes(pickle.dumps(seqs))
        mark(seqs_pkl)
    p1, p2 = out / "AKA017_SELECTED_R1.fa", out / "AKA017_SELECTED_R2.fa"
    iso.write_pair_fastas(p1, p2, eligible, seqs)

    dec_in, dec_out = out / "AKA017_RELAXED_FRAGMENT_DECODER_INPUT.tsv", out / "AKA017_RELAXED_FRAGMENT_DECODER.tsv"
    if not done(dec_out):
        iso.write_tsv(dec_in, [{"pair_id": pid, **{k: one_mate_meta[pid][k] for k in ("te_mate", "host_mate", "unit_ids", "old_primary")},
                                "te_read": seqs[pid][one_mate_meta[pid]["te_mate"] - 1]} for pid in sorted(one_mate_meta)],
                      ["pair_id", "te_mate", "host_mate", "unit_ids", "old_primary", "te_read"])
        iso.run_edlib_subprocess(args.edlib_python.resolve(), iso_script, v3, codebook, te_library, fammap, dec_in, dec_out)
        mark(dec_out)
    decoder_by = {int(r["pair_id"]): r for r in read_tsv(dec_out)}
    strict_new = {pid: r for pid, r in decoder_by.items()
                  if str(r.get("strict_fragment_decodable", "0")) == "1" and str(r.get("old_primary", "0")) != "1"}
    say(f"strict decoded outside old primary: {len(strict_new):,}")

    bt, bb = iso.find_bowtie2(shared)
    prefix = iso.ensure_bowtie_index(shared, reference, out, bb)

    def bowtie(sam: Path, bt_args: list[str]) -> None:
        if not done(sam):
            iso.run_checked(bt, bt_args + ["-S", str(sam), "-p", str(args.threads)])
            mark(sam)

    suffix = "" if args.host_read_names == "iso" else "_AKANAMES"

    def host_map(stem: str, meta: dict) -> dict:
        """Single-mate host mapping; qnames rewritten to P{pid} before the frozen parser."""
        sam = out / f"{stem}{suffix}.sam"
        if args.host_read_names == "iso":
            fa = out / f"{stem}.fa"
            iso.write_host_fasta(fa, meta, seqs)
            bowtie(sam, ["-x", str(prefix), "-f", "-U", str(fa), "-k", "2", "--very-sensitive", "--no-unal", "--no-hd", "--no-sq"])
            return iso.parse_host_sam(sam)
        fa = out / f"{stem}{suffix}.fa"
        with fa.open("w", encoding="ascii", newline="\n") as f:
            for pid in sorted(meta):
                hm = meta[pid]["host_mate"]
                f.write(f">p{pid}|h{hm}\n{seqs[pid][hm - 1]}\n")
        raw = out / f"{stem}{suffix}.raw.sam"
        bowtie(raw, ["-x", str(prefix), "-f", "-U", str(fa), "-k", "2", "--very-sensitive", "--no-unal", "--no-hd", "--no-sq"])
        if not done(sam):
            with raw.open(encoding="ascii") as src, sam.open("w", encoding="ascii", newline="\n") as dst:
                for line in src:
                    name, rest = line.split("\t", 1)
                    dst.write("P" + name[1:].split("|", 1)[0] + "\t" + rest)
            mark(sam)
        return iso.parse_host_sam(sam)

    host = host_map("AKA017_HOST_MATES", eligible)
    ref_sam = out / "AKA017_EVIDENCE_COMPLETE_REFERENCE.sam"
    bowtie(ref_sam, ["-x", str(prefix), "-f", "-1", str(p1), "-2", str(p2), "--very-sensitive-local", "--no-mixed", "--no-discordant",
                     "-I", "0", "-X", str(iso.Q99), "-k", "20", "--no-unal", "--reorder", "--no-hd", "--no-sq"])
    ledger, te_records = iso.parse_reference_sam(ref_sam, eligible, seqs, by_family)

    unique = [{"pair_id": pid, "family": m["family"], **{k: host[pid][k] for k in ("chrom", "start", "end", "strand")},
               "primary": m["primary"], "reference_state": m["reference_state"]}
              for pid, m in ledger.items() if host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"]
    decoded_host = {}
    if strict_new:
        decoded_host = host_map("AKA017_RELAXED_DECODED_HOST_MATES", {pid: {"host_mate": int(r["host_mate"])} for pid, r in strict_new.items()})

    # ---- the ISO main() clustering/anchoring block, verbatim in logic ----
    old_primary_unique = [r for r in unique if r["primary"]]
    primary_clusters = iso.cluster_mappings(old_primary_unique, require_folds=False, prefix="PRIMARY")
    decoded_unique = [{"pair_id": pid, "family": r.get("family", ""), **{k: decoded_host[pid][k] for k in ("chrom", "start", "end", "strand")},
                       "primary": False, "reference_state": "FRAGMENT_DECODED_EXPANSION"}
                      for pid, r in strict_new.items() if decoded_host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"]
    combined = iso.cluster_mappings(old_primary_unique + decoded_unique, require_folds=False, prefix="FRAGMENT_DECODED")
    pidx = defaultdict(list)
    for c in primary_clusters:
        pidx[(c["family"], c["chrom"])].append(c)
    relaxed_clusters = []
    for c in combined:
        members = sorted({int(x) for x in str(c["pair_ids"]).split(";") if x} & set(strict_new))
        if not members:
            continue
        d = min((iso.interval_gap(int(c["start"]), int(c["end"]), int(b["start"]), int(b["end"])) for b in pidx[(c["family"], c["chrom"])]), default=10**12)
        if d <= iso.CLUSTER_GAP:
            continue
        confirmed = [pid for pid in members if str(strict_new[pid].get("masked_family_confirmed", "0")) == "1"
                     and decoded_host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"]
        if len(confirmed) < iso.MIN_PAIRS:
            continue
        folds = sorted({iso.fold_of(p) for p in confirmed})
        relaxed_clusters.append({"cluster_id": c["cluster_id"], "cluster_label": f"RELAXED_{int(c['cluster_id']):05d}", "family": c["family"],
                                 "chrom": c["chrom"], "start": c["start"], "end": c["end"], "pair_count": len(confirmed),
                                 "fold_count": len(folds), "pair_ids": ";".join(map(str, confirmed)), "folds": ";".join(map(str, folds))})
    corrected = iso.cluster_mappings([r for r in unique if r["reference_state"] != "EVIDENCE_COMPLETE_CLEAN_REFERENCE"],
                                     require_folds=True, prefix="CORRECTED", keep_all=True)
    fbmap, known = iso.load_family_map(fammap)
    corrected = iso.annotate_corrected_clusters(corrected, primary_clusters + relaxed_clusters, iso.parse_repeatmasker(repeatmasker, fbmap, known))
    iso.write_tsv(out / f"AKA017_EVIDENCE_COMPLETE_UNEXPLAINED_CLUSTERS{suffix}.tsv", corrected,
                  ["cluster_id", "cluster_label", "family", "chrom", "start", "end", "pair_count", "fold_count", "event_2p2fold", "pair_ids", "folds",
                   "nearest_baseline_same_family_bp", "nearest_repeatmasker_same_family_bp", "novel", "resident_compatible", "novel_nonresident"])
    targets = [c for c in corrected if c["novel_nonresident"] and int(c.get("event_2p2fold", 0)) == 1]
    witnesses, anchors = iso.build_a_witnesses(targets, iso.assign_pairs_to_clusters(corrected, {}), ledger, te_records, by_family)
    b_rows = [{"pair_id": pid, "fold": iso.fold_of(pid), "family": m["family"], **{k: host[pid][k] for k in ("chrom", "start", "end", "strand", "mapq")}}
              for pid, m in ledger.items() if m["reference_state"] == "NO_LOCAL_PROPER" and host.get(pid, {}).get("host_class") == "UNIQUE_HIGH"]
    a_union = {int(x) for a in anchors for x in str(a["channel_A_pair_ids"]).split(";") if x}
    candidates = iso.score_b_completion(anchors, b_rows, a_union)
    cand_fields = ["cluster_id", "chrom", "breakpoint0", "family", "anchor_host_side", "channel_A_pairs", "channel_A_folds", "channel_A_pair_ids",
                   "channel_A_fold_ids", "left_B_pairs", "left_B_folds", "right_B_pairs", "right_B_folds", "distinct_B_pairs",
                   "left_B_pair_ids", "right_B_pair_ids", "bilateral_B_2p2fold"]
    iso.write_tsv(out / f"AKA017_ANCHOR_PE_CANDIDATES{suffix}.tsv", candidates, cand_fields)

    # ---- comparison with the frozen AKA chain ----
    def cand_key(r, with_cluster: bool):
        base = (r["family"], r["chrom"], int(r["breakpoint0"]), r["anchor_host_side"], ids(r["channel_A_pair_ids"]),
                ids(r["left_B_pair_ids"]), ids(r["right_B_pair_ids"]))
        return ((int(r["cluster_id"]),) + base) if with_cluster else base

    frozen_c = read_tsv(v3 / FROZEN_AKA["candidates"])
    frozen_clusters = read_tsv(v3 / FROZEN_AKA["clusters"])
    fz_targets = [r for r in frozen_clusters if str(r.get("novel_nonresident_event")).strip() in {"1", "True", "true"}]
    comparison = {
        "counts": {
            "corrected_clusters": [len(corrected), len(frozen_clusters)],
            "novel_nonresident_targets": [len(targets), len(fz_targets)],
            "replicated_A_anchors": [len(anchors), len(frozen_c)],
            "B_universe_unique_high": [len(b_rows), 8783],
            "bilateral": [sum(int(r["bilateral_B_2p2fold"]) for r in candidates), sum(int(r["bilateral_B_2p2fold"]) for r in frozen_c)],
        },
        "target_pair_sets_equal": sorted(ids(c["pair_ids"]) for c in targets) == sorted(ids(c["pair_ids"]) for c in fz_targets),
        "candidates_equal_ignoring_cluster_id": sorted(cand_key(r, False) for r in candidates) == sorted(cand_key(r, False) for r in frozen_c),
        "candidates_equal_including_cluster_id": sorted(cand_key(r, True) for r in candidates) == sorted(cand_key(r, True) for r in frozen_c),
    }
    mine = {cand_key(r, False) for r in candidates}
    theirs = {cand_key(r, False) for r in frozen_c}
    comparison["only_iso_procedure"] = [list(map(str, k[:4])) for k in sorted(mine - theirs)[:20]]
    comparison["only_frozen_chain"] = [list(map(str, k[:4])) for k in sorted(theirs - mine)[:20]]
    comparison["verdict"] = (
        "ISO_PROCEDURE_REPRODUCES_AKA017_CANDIDATES_EXACTLY" if comparison["candidates_equal_including_cluster_id"]
        else "ISO_PROCEDURE_REPRODUCES_AKA017_CANDIDATES_UP_TO_CLUSTER_IDS" if comparison["candidates_equal_ignoring_cluster_id"]
        else "ISO_PROCEDURE_DIVERGES_FROM_AKA017_FROZEN_CHAIN"
    )
    comparison["truth_opened"] = False
    comparison["host_read_names"] = args.host_read_names
    comparison["target_cluster_ids_equal"] = sorted((int(c["cluster_id"]), ids(c["pair_ids"])) for c in targets) == sorted((int(c["cluster_id"]), ids(c["pair_ids"])) for c in fz_targets)
    (out / f"ISO_PROCEDURE_AKA017_COMPARISON{suffix}.json").write_text(json.dumps(comparison, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    say(json.dumps(comparison["counts"]))
    say(comparison["verdict"])
    return 0 if comparison["verdict"].startswith("ISO_PROCEDURE_REPRODUCES") else 1


if __name__ == "__main__":
    raise SystemExit(main())
