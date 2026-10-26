"""VM-only: ISO1 acquisition of candidate loci, stage by stage.

Runs the ``evitea.incumbent`` stage functions on the ISO1 msBWT in short
resumable steps (each step caches its result under ``--cache``) and compares
every stage with the recorded ISO1 stage tables
(``TE_locator_route_b/reports/G1*`` and ``TE_locator_route_b_v2/reports``).
Bowtie2 host mappings are split into chunks (per-read results do not depend on
batch composition; parsers are order-independent).  Archive members whose names
contain ``TRUTH`` are never read.

    python iso1_acquisition_vm.py g1a|g1b|g1c|g1d|g1e|g1f-fasta|g1f-map N|g1f|
                                  g1hr-plan|g1hr-map N|g1hr|loci|all-status
"""

from __future__ import annotations

import csv
import io
import json
import os
import pickle
import sys
import time
import zipfile
from pathlib import Path

HOME = Path(os.path.expanduser("~"))
sys.path.insert(0, str(Path(os.environ.get("EVITEA_REPO", ".")) / "src"))
sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

from evitea.incumbent import MsbwtSample  # noqa: E402
from evitea.incumbent.acquisition import _host_record  # noqa: E402
from evitea.incumbent.anchor_census import (  # noqa: E402
    build_anchor_catalogue,
    census_anchor_intervals,
    load_te_family_map,
    parse_te_consensus,
)
from evitea.incumbent.candidates import PairPanelEvidence, cluster_candidate_loci, reconstruct_candidate_pairs  # noqa: E402
from evitea.incumbent.coverage_panel import (  # noqa: E402
    baseline_anchor_loci,
    hard_gap_h1_loci,
    optimize_rescue_panel,
    short_exact_census,
)
from evitea.incumbent.external import run_bowtie2, write_fasta  # noqa: E402
from evitea.incumbent.host_mate_gate import (  # noqa: E402
    PRIMARY_HOST_OPTIONS,
    one_mate_host_requests,
    parse_primary_host_sam,
    resident_pair_families,
    supported_host_cluster_pairs,
    working_candidate_pairs,
)
from evitea.incumbent.pair_aware_panel import select_pair_aware_panel  # noqa: E402
from evitea.incumbent.recall_repair import (  # noqa: E402
    RESCUE_HOST_OPTIONS,
    enumerate_rescue_pairs,
    parse_rescue_host_sam,
    plan_rescue_anchors,
    rescue_host_requests,
    rescue_query_name,
    supported_rescue_pairs,
)
from evitea.incumbent.sparse_panel import route_sparse_panel, select_sparse_panel  # noqa: E402

ASSET = Path(os.environ["EVITEA_DATA"]) / "TEMP2_ISO1"
TE_FASTA = ASSET / "TE consensus library" / "transposon_sequence_set.fa"
FAMMAP = ASSET / "fbte_to_flybase_family.tsv"
G1 = Path(os.environ["EVITEA_ROUTE_B"]) / "reports"
V2 = Path(os.environ["EVITEA_ROUTE_B_V2"]) / "reports"
MSBWT = Path(os.environ.get("ISO1_MSBWT", HOME / "msbwt_iso"))
CACHE = Path(os.environ.get("ISO1_CACHE", HOME / "scratch/iso1"))
BOWTIE2 = Path(os.environ.get("BOWTIE2", HOME / "bowtie2-2.5.5-linux-x86_64/bowtie2"))
HOST_INDEX = Path(os.environ.get("HOST_INDEX", ASSET / "Reference/bowtie2_index/dm6"))
CHUNK = 60_000
ZIPS = {s: f"{d}/TE_LOCATOR_ROUTE_B_{s.upper()}_RESULTS_FOR_CHATGPT.zip" for s, d in {
    "g1a": "G1A_native_results", "g1b": "G1B_sparse_panel_routing", "g1c": "G1C_coverage_rescue",
    "g1d": "G1D_optimized_rescue_routing", "g1e": "G1E_pair_aware_rescue", "g1f": "G1F_host_mate_locus_gate",
    "g1g": "G1G_reference_resident_compatibility", "g1hr": "G1HR_candidate_recall_repair"}.items()}


def _member(stage: str, member: str):
    if "TRUTH" in member.upper():
        raise SystemExit("truth-bearing member refused")
    return zipfile.ZipFile(G1 / ZIPS[stage]).open(member)


def ztsv(stage, member):
    with _member(stage, member) as raw:
        return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8", newline=""), delimiter="\t"))


def znpy(stage, member):
    with _member(stage, member) as raw:
        return np.load(io.BytesIO(raw.read()), allow_pickle=False)


def load(name):
    return pickle.loads((CACHE / f"{name}.pkl").read_bytes())


def save(name, value):
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = CACHE / f"{name}.pkl.tmp"
    tmp.write_bytes(pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL))
    os.replace(tmp, CACHE / f"{name}.pkl")


def record(stage, checks):
    path = CACHE / "ISO1_ACQUISITION_CHECKS.json"
    doc = json.loads(path.read_text()) if path.exists() else {}
    doc[stage] = checks
    path.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    print(stage, json.dumps(checks), flush=True)


def index():
    return MsbwtSample(MSBWT)


def s(value):
    return "" if value is None else str(value)


def g(text):
    return text[:-2] if text.endswith(".0") else text


def step_g1a():
    cat = build_anchor_catalogue(parse_te_consensus(TE_FASTA), load_te_family_map(FAMMAP))
    census_anchor_intervals(cat, index())
    hist = ztsv("g1a", "G1A_ANCHOR_INTERVALS.tsv")
    cat_bad = sum(
        (q.query_id, str(q.k), q.query_seq, ",".join(q.source_te_ids), ",".join(q.source_families))
        != (h["query_id"], h["k"], h["query_seq"], h["source_te_ids"], h["source_families"])
        for q, h in zip(cat.queries, hist))
    int_bad = sum((q.interval_lo, q.interval_hi) != (int(h["interval_lo"]), int(h["interval_hi"]))
                  for q, h in zip(cat.queries, hist))
    prov = ztsv("g1a", "G1A_ANCHOR_PROVENANCE.tsv")
    prov_bad = sum((x.query_id, x.te_id, x.family, str(x.k), str(x.start), x.region, x.orientation)
                   != (h["query_id"], h["te_id"], h["family"], h["k"], h["start"], h["region"], h["source_orientation"])
                   for x, h in zip(cat.sources, prov))
    save("g1a", cat)
    record("g1a", {"queries": len(cat.queries), "catalogue_equal": cat_bad == 0 and len(cat.queries) == len(hist),
                   "intervals_equal": int_bad == 0, "provenance_equal": prov_bad == 0 and len(prov) == len(cat.sources)})


def step_g1b():
    cat = load("g1a")
    panel = select_sparse_panel(cat)
    hist = ztsv("g1b", "G1B_SELECTED_PANEL.tsv")
    panel_eq = [(p.query_id, p.query_seq, str(p.width), str(p.channel_mask), ",".join(sorted(p.source_families))) for p in panel] == \
        [(h["query_id"], h["query_seq"], h["interval_width"], h["channel_mask"], h["source_families"]) for h in hist]
    pairs = route_sparse_panel(panel, index())
    hp = ztsv("g1b", "G1B_CANDIDATE_PAIRS.tsv")
    pairs_eq = [(str(p.pair_id), str(int(p.mate1_hit)), str(int(p.mate2_hit)), str(p.channel_mask), ",".join(p.families)) for p in pairs] == \
        [(h["pair_id"], h["mate1_hit"], h["mate2_hit"], h["channel_mask"], h["families"]) for h in hp]
    save("g1b", (panel, pairs))
    record("g1b", {"panel": len(panel), "panel_equal": panel_eq, "pairs": len(pairs), "pairs_equal": pairs_eq})


def step_g1c():
    cat = load("g1a")
    panel, _ = load("g1b")
    fam = cat.family_by_te()
    short = short_exact_census(cat.records, fam, index())
    hist = ztsv("g1c", "G1C_SHORT_EXACT_LOCUS_COUNTS.tsv")
    short_eq = [(x.te_id, x.family, x.orientation, str(x.start), str(x.k), x.query_seq, str(x.width)) for x in short] == \
        [(h["te_id"], h["family"], h["orientation"], h["start"], h["k"], h["query_seq"], h["interval_width"]) for h in hist]
    h1 = hard_gap_h1_loci(cat.records, fam, baseline_anchor_loci(cat, panel), short)
    hh = ztsv("g1c", "G1C_H1_LOCUS_COUNTS.tsv")
    h1_eq = [(x.te_id, x.family, x.orientation, str(x.start), x.source_21mer) for x in h1] == \
        [(h["te_id"], h["family"], h["orientation"], h["start"], h["source_21mer"]) for h in hh]
    save("g1c", (short, h1))
    record("g1c", {"short": len(short), "short_equal": short_eq, "h1": len(h1), "h1_equal": h1_eq})


def step_g1d():
    cat = load("g1a")
    panel, _ = load("g1b")
    short, h1 = load("g1c")
    rp = optimize_rescue_panel(cat, panel, short, h1, index())
    hist = [r for r in ztsv("g1d", "G1D_SELECTED_ACTIONS.tsv") if r["B"] == "1024"]
    actions_eq = [a.action_id for a in rp.selected_actions] == [h["action_id"] for h in hist]
    gc = [r for r in ztsv("g1d", "G1D_GROUP_COVERAGE.tsv") if r["B"] == "1024"]
    mine = {k: (repr(v.selected_coverage), str(int(v.selected_terminal5)), str(int(v.selected_terminal3)), str(v.window_count))
            for k, v in rp.group_targets.items()}
    his = {(r["te_id"], r["orientation"]): (r["selected_coverage"], r["selected_terminal5"], r["selected_terminal3"], r["window_count"]) for r in gc}
    cs = ztsv("g1d", "G1D_CHOSEN_PANEL_SEQUENCES.tsv")
    seqs_eq = [(q, str(w), ",".join(sorted(rp.selected_sources[q]))) for q, w in sorted(rp.selected_sequences.items())] == \
        [(r["sequence"], r["interval_width"], r["sources"]) for r in cs]
    save("g1d", rp)
    record("g1d", {"actions": len(rp.selected_actions), "actions_equal": actions_eq, "groups_equal": mine == his,
                   "sequences_equal": seqs_eq})


def step_g1e():
    cat = load("g1a")
    panel, _ = load("g1b")
    short, h1 = load("g1c")
    rp = load("g1d")
    pa = select_pair_aware_panel(cat, panel, short, rp, index())
    hp, hc = znpy("g1e", "G1E_CANDIDATE_PAIR_IDS.npy"), znpy("g1e", "G1E_PAIR_READ_COUNTS.npy")
    hs = ztsv("g1e", "G1E_SELECTED_PANEL_SEQUENCES.tsv")
    hg = ztsv("g1e", "G1E_GROUP_COVERAGE.tsv")
    groups_eq = [(x.te_id, x.family, x.orientation, str(x.window_count), repr(x.selected_coverage)) for x in pa.groups] == \
        [(r["te_id"], r["family"], r["orientation"], r["window_count"], r["g1e_selected_coverage"]) for r in hg]
    save("g1e", pa)
    record("g1e", {"pairs": int(len(pa.pair_ids)), "pairs_equal": bool(np.array_equal(hp, pa.pair_ids)),
                   "read_counts_equal": bool(np.array_equal(hc, pa.pair_read_counts)),
                   "sequences_equal": sorted(pa.selected_sequences) == [r["sequence"] for r in hs], "groups_equal": groups_eq})


def _write_chunks(stem, records):
    for old in CACHE.glob(f"{stem}_*.fa"):
        old.unlink()
    n = 0
    for n, start in enumerate(range(0, len(records), CHUNK)):
        write_fasta(CACHE / f"{stem}_{n:03d}.fa", records[start:start + CHUNK])
    return n + 1 if records else 0


def _map_chunk(stem, options, n, threads=2):
    fa, sam = CACHE / f"{stem}_{n:03d}.fa", CACHE / f"{stem}_{n:03d}.sam"
    tmp = sam.with_suffix(".sam.tmp")
    run_bowtie2(BOWTIE2, ["-x", str(HOST_INDEX), "-f", "-U", str(fa), "-S", str(tmp), *options, "-p", str(threads)])
    os.replace(tmp, sam)
    print("mapped", sam.name, flush=True)


def _concat(stem, chunks):
    out = CACHE / f"{stem}.sam"
    with out.open("w", encoding="utf-8") as handle:
        for n in range(chunks):
            handle.write((CACHE / f"{stem}_{n:03d}.sam").read_text(encoding="utf-8"))
    return out


def step_g1f_fasta():
    pa = load("g1e")
    idx = index()
    requests = one_mate_host_requests(pa, idx)
    dollars = idx.mate_dollars(requests)
    records = [(f"p{pid}|h{requests[pid] + 1}", idx.recover(dollars[pid])) for pid in sorted(requests)]
    chunks = _write_chunks("g1f_host", records)
    save("g1f_requests", (requests, chunks))
    print("requests", len(requests), "chunks", chunks)


def step_g1f():
    requests, chunks = load("g1f_requests")
    primary = parse_primary_host_sam(_concat("g1f_host", chunks), requests)
    hist = ztsv("g1f", "G1F_HOST_MATE_MAPPING.tsv")
    mine = [(str(r.pair_id), r.status, str(r.alignment_count), r.chrom, s(r.pos), s(r.midpoint), s(r.mapq), r.strand, r.cigar) for r in primary.values()]
    his = [(h["pair_id"], h["status"], h["alignment_count"], h["chrom"], h["pos"], h["midpoint"], h["mapq"], h["strand"], h["cigar"]) for h in hist]
    _, panel_pairs = load("g1b")
    host = supported_host_cluster_pairs(primary)
    working = working_candidate_pairs(host, panel_pairs)
    families = resident_pair_families(host, load("g1e"))
    hf = ztsv("g1g", "G1G_HOST_CLUSTER_PAIR_FAMILIES.tsv")
    save("g1f", (primary, host, working, families))
    diff = [(a, b) for a, b in zip(mine, his) if a != b]
    record("g1f_g1g", {"mappings": len(mine), "mapping_equal": mine == his, "mapping_diffs": len(diff),
                       "first_diffs": [list(map(list, d)) for d in diff[:3]],
                       "working": int(len(working)),
                       "working_equal": bool(np.array_equal(working, znpy("g1f", "G1F_WORKING_CANDIDATE_PAIR_IDS.npy"))),
                       "families_equal": [(str(k), v) for k, v in families.items()] == [(h["pair_id"], h["families"]) for h in hf]})


def step_g1hr_plan():
    cat = load("g1a")
    pa = load("g1e")
    anchors = plan_rescue_anchors(cat, pa.groups)
    h = ztsv("g1hr", "G1HR_RESCUE_ANCHORS.tsv")
    anchors_eq = [(a.query_id, a.query_seq, str(a.k), str(a.width), a.family) for a in anchors] == \
        [(r["query_id"], r["query_seq"], r["k"], r["interval_width"], r["family"]) for r in h]
    idx = index()
    evidence = enumerate_rescue_pairs(anchors, pa.pair_ids, idx)
    he = ztsv("g1hr", "G1HR_NEW_PAIR_EVIDENCE.tsv")
    evidence_eq = [(str(e.pair_id), str(int(e.mate1_hit)), str(int(e.mate2_hit)), ",".join(e.families)) for e in evidence.values()] == \
        [(r["pair_id"], r["mate1_hit"], r["mate2_hit"], r["families"]) for r in he]
    requests = rescue_host_requests(evidence)
    dollars = idx.mate_dollars(requests)
    # ISO1 profile: the recorded run names these queries {pid}.
    chunks = _write_chunks("g1hr_host", [(rescue_query_name(p, "numeric"), idx.recover(dollars[p])) for p in sorted(requests)])
    save("g1hr_plan", (anchors, evidence, requests, chunks))
    record("g1hr_plan", {"anchors": len(anchors), "anchors_equal": anchors_eq, "evidence": len(evidence),
                         "evidence_equal": evidence_eq, "requests": len(requests), "chunks": chunks})


def step_g1hr():
    anchors, evidence, requests, chunks = load("g1hr_plan")
    mapping = parse_rescue_host_sam(_concat("g1hr_host", chunks), requests)
    h = ztsv("g1hr", "G1HR_NEW_HOST_MAPPING.tsv")
    mine = [(str(r.pair_id), r.status, str(r.alignment_count), r.chrom, s(r.pos), s(r.midpoint), s(r.mapq), r.strand, r.cigar) for r in mapping.values()]
    his = [(r["pair_id"], r["status"], r["alignment_count"], r["chrom"], g(r["pos"]), g(r["midpoint"]), g(r["mapq"]), r["strand"], r["cigar"]) for r in h]
    primary, host, working, families = load("g1f")
    repaired = np.union1d(working, np.asarray(sorted(supported_rescue_pairs(mapping)), dtype=np.uint64))
    save("g1hr", (mapping, repaired))
    diff = [(a, b) for a, b in zip(mine, his) if a != b]
    record("g1hr", {"mappings": len(mine), "mapping_equal": mine == his, "mapping_diffs": len(diff),
                    "first_diffs": [list(map(list, d)) for d in diff[:3]], "repaired": int(len(repaired)),
                    "repaired_equal": bool(np.array_equal(repaired, znpy("g1hr", "G1HR_REPAIRED_CANDIDATE_PAIR_IDS.npy")))})


def step_loci():
    primary, host, working, families = load("g1f")
    anchors, evidence, requests, chunks = load("g1hr_plan")
    mapping, repaired = load("g1hr")
    _, panel_pairs = load("g1b")
    candidates, meta = reconstruct_candidate_pairs(
        final_pair_ids=repaired.tolist(), working_pair_ids=working.tolist(),
        primary_host={p: _host_record(r) for p, r in primary.items()},
        rescue_host={p: _host_record(r) for p, r in mapping.items()}, resident_families=families,
        rescue_evidence={p: PairPanelEvidence(p, ",".join(r.families), r.mate1_hit, r.mate2_hit) for p, r in evidence.items()},
        panel_evidence={r.pair_id: PairPanelEvidence(r.pair_id, ",".join(r.families), r.mate1_hit, r.mate2_hit) for r in panel_pairs})
    loci = cluster_candidate_loci(candidates)
    with (V2 / "V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1" / "V2FCRC1_CANDIDATE_LOCI.tsv").open(encoding="utf-8", newline="") as handle:
        hist = list(csv.DictReader(handle, delimiter="\t"))
    mine = [{"locus_id": L.locus_id, "chr": L.chrom, "candidate_family": L.family, "candidate_start_midpoint": f"{L.start:.1f}",
             "candidate_end_midpoint": f"{L.end:.1f}", "candidate_midpoint": f"{L.midpoint:.1f}",
             "candidate_pair_count": str(len(L.pair_ids))} for L in loci]
    equal = mine == [{k: h[k] for k in mine[0]} for h in hist] if mine else not hist
    save("loci", loci)
    record("candidate_loci", {"loci": len(loci), "frozen": len(hist), "equal": equal, "meta": {k: int(v) for k, v in meta.items()}})


def main(argv):
    started = time.time()
    step = argv[0]
    table = {"g1a": step_g1a, "g1b": step_g1b, "g1c": step_g1c, "g1d": step_g1d, "g1e": step_g1e,
             "g1f-fasta": step_g1f_fasta, "g1f": step_g1f, "g1hr-plan": step_g1hr_plan, "g1hr": step_g1hr,
             "loci": step_loci}
    if step == "g1f-map":
        _map_chunk("g1f_host", PRIMARY_HOST_OPTIONS, int(argv[1]))
    elif step == "g1hr-map":
        _map_chunk("g1hr_host", RESCUE_HOST_OPTIONS, int(argv[1]))
    elif step == "all-status":
        print((CACHE / "ISO1_ACQUISITION_CHECKS.json").read_text())
    else:
        table[step]()
    print(f"[{step}] {time.time() - started:.1f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
