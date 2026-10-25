"""VM-only: the ISO1 stages after candidate acquisition, in resumable chunks.

The stages run from split-junction events to repeat-aware calls.  The script
follows ``evitea.incumbent.engine.run_incumbent_caller`` step for step, using the
same functions, on the ISO1 candidate loci produced (and checked) by
``iso1_acquisition_vm.py``.  Per-locus layers run in locus chunks through
``map_locus_batches`` (per-locus results are independent; rows are
concatenated in locus order).  Each stage is compared with the recorded ISO1
stage tables; the final callset with the recorded ISO1 callset bytes
(192 rows, b3e48681...).

    python iso1_downstream_vm.py prep | sj N | sj-final | mr N | mr-final |
                                 route | strict | groups | final | status
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

HOME = Path(os.path.expanduser("~"))
sys.path.insert(0, str(Path(os.environ.get("EVITEA_REPO", ".")) / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.dont_write_bytecode = True

import iso1_acquisition_vm as acq  # noqa: E402

from evitea.incumbent.engine import FAMILY_KMER, FAMILY_KMER_MIN_ENTROPY, _te_index_prefix  # noqa: E402
from evitea.incumbent.events import legacy_tsv_bytes  # noqa: E402
from evitea.incumbent.junction_groups import run_junction_groups  # noqa: E402
from evitea.incumbent.mate_rescue import (  # noqa: E402
    certify_strict_pairs,
    finalize_mate_rescue,
    guard_candidate_loci,
    one_sided_targets,
    route_one_sided_targets,
)
from evitea.incumbent.parallel import map_locus_batches  # noqa: E402
from evitea.incumbent.repeat_rescue import crossvalidate, reaggregate, rescue_calls  # noqa: E402
from evitea.incumbent.sequence import fasta_records_full_header, load_reference_contigs  # noqa: E402
from evitea.incumbent.split_junction import finalize_split_junctions  # noqa: E402
from evitea.incumbent.te_certification import align_segments, certify_groups, full_te_segment, used_groups  # noqa: E402
from evitea.incumbent.te_library import (  # noqa: E402
    LocalTEIndex,
    load_family_kmer_panels,
    load_heuristic_family_map,
    load_record_family_map,
    load_te_records,
    normalize_family,
)

SAMPLE = "ISO1"
# Loci per resumable chunk (per layer; mate-anchored search has heavier loci).
CHUNKS = {"split_junction": 400, "mate_rescue": 100}
V2 = acq.V2
FROZEN_R9E = "b3e48681568d616142fe534995d3eeb9e5163907992e9d0893919aea42e6ebc5"
load, save, record = acq.load, acq.save, acq.record


def state():
    return load("ds_state")


def step_prep():
    loci = load("loci")
    families = {locus.family for locus in loci}
    family_kmers, _ = load_family_kmer_panels(acq.TE_FASTA, acq.FAMMAP, families, FAMILY_KMER, FAMILY_KMER_MIN_ENTROPY)
    te_index = LocalTEIndex(load_te_records(acq.TE_FASTA, load_record_family_map(acq.FAMMAP)))
    reference = load_reference_contigs(acq.ASSET / "Reference" / "dm6_TEMP2_compatible.fa", {locus.chrom for locus in loci})
    reads = acq.index().pair_reads(pid for locus in loci for pid in locus.pair_ids)
    save("ds_state", {"reads": reads, "reference": reference, "family_kmers": family_kmers, "te_index": te_index})
    print("loci", len(loci), "pairs", len(reads))


def _chunk(layer, loci, n):
    CHUNK = CHUNKS[layer]
    part = loci[n * CHUNK:(n + 1) * CHUNK]
    rows, counts = map_locus_batches(layer, part, state(), workers=int(os.environ.get("ISO1_WORKERS", "2")))
    save(f"{layer}_{n:03d}", (rows, counts))
    print(layer, n, len(part), "loci", len(rows), "rows")


def _auto(layer, loci, started, budget=float(os.environ.get("ISO1_START_BUDGET", "50"))):
    """Run missing chunks while less than ``budget`` seconds have elapsed."""

    CHUNK = CHUNKS[layer]
    total = (len(loci) + CHUNK - 1) // CHUNK
    for n in range(total):
        if (acq.CACHE / f"{layer}_{n:03d}.pkl").exists():
            continue
        if time.time() - started > budget:
            break
        _chunk(layer, loci, n)
    done = sum((acq.CACHE / f"{layer}_{n:03d}.pkl").exists() for n in range(total))
    print(f"{layer}: {done}/{total} chunks")


def _gather(layer, loci):
    rows, counts = [], {}
    CHUNK = CHUNKS[layer]
    for n in range((len(loci) + CHUNK - 1) // CHUNK):
        part_rows, part_counts = load(f"{layer}_{n:03d}")
        rows.extend(part_rows)
        for key, value in part_counts.items():
            counts[key] = counts.get(key, 0) + value
    return rows, counts


def step_sj_final():
    rows, counts = _gather("split_junction", load("loci"))
    split = finalize_split_junctions(SAMPLE, rows, counts)
    fc = V2 / "V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"
    save("split", split)
    record("fcrc1", {"events": len(split.all_events), "calls": len(split.calls),
                     "all_events_bytes_equal": legacy_tsv_bytes(split.all_events) == (fc / "V2FCRC1_ALL_EVENTS.tsv").read_bytes(),
                     "calls_bytes_equal": legacy_tsv_bytes(split.calls) == (fc / "V2FCRC1_CALLS_SCORE4.tsv").read_bytes()})


def eligible():
    split = load("split")
    return guard_candidate_loci(load("loci"), split.calls)


def step_mr_final():
    elig, guarded = eligible()
    rows, counts = _gather("mate_rescue", elig)
    mates = finalize_mate_rescue(SAMPLE, elig, guarded, rows, counts)
    split = load("split")
    r4 = V2 / "V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4" / "V2R4_COMBINED_CALLS.tsv"
    save("mates", mates)
    record("r4", {"calls": len(mates.calls), "ledger": len(mates.ledger),
                  "combined_bytes_equal": legacy_tsv_bytes((*split.calls, *mates.calls)) == r4.read_bytes()})


def step_route():
    split, mates = load("split"), load("mates")
    targets = one_sided_targets(mates.ledger)
    reserved = {pid for call in (*split.calls, *mates.calls) for pid in call.fragment_pair_ids}
    routing = route_one_sided_targets(targets=targets, reserved_pair_ids=reserved,
                                      reference=state()["reference"], locator=acq.index())
    routed_reads = acq.index().pair_reads(routing.pair_to_loci)
    save("routing", (targets, reserved, routing, routed_reads))
    record("r6_route", {"targets": len(targets), "routed_pairs": len(routing.pair_to_loci)})


def step_strict():
    targets, reserved, routing, routed_reads = load("routing")
    st = state()
    strict = certify_strict_pairs(routing=routing, reads=routed_reads, reference=st["reference"],
                                  family_kmers=st["family_kmers"], te_index=st["te_index"], bowtie2=acq.BOWTIE2,
                                  host_index_prefix=acq.HOST_INDEX, work_directory=acq.CACHE / "strict_pairs")
    save("strict", strict)
    record("r6_strict", {"winners": len(strict.winners)})


def step_groups():
    targets, reserved, routing, routed_reads = load("routing")
    strict = load("strict")
    st = state()
    idx = acq.index()
    groups = run_junction_groups(targets=targets, routing=routing, strict_winners=strict.winners,
                                 reserved_call_pair_ids=reserved, reads=routed_reads, reference=st["reference"],
                                 family_kmers=st["family_kmers"], locator=idx, bowtie2=acq.BOWTIE2,
                                 host_index_prefix=acq.HOST_INDEX, work_directory=acq.CACHE / "junction_groups",
                                 pair_reader=idx.pair_reads)
    save("groups", groups)
    record("r7", {"winners": len(groups.winners), "ledger": len(groups.ledger)})


def step_final():
    split, mates, groups = load("split"), load("mates"), load("groups")
    targets, reserved, routing, routed_reads = load("routing")
    used = used_groups(groups.ledger, groups.winners, groups.accepted_proposals)
    idx = acq.index()
    source_reads = {pid: routed_reads[pid] for pid in {x.source_pair_id for x in used} if pid in routed_reads}
    missing = {x.source_pair_id for x in used} - set(source_reads)
    if missing:
        source_reads.update(idx.pair_reads(missing))
    segments = {x.group_hash: full_te_segment(source_reads[x.source_pair_id][x.source_mate - 1], x.read_view,
                                              x.junction_semantic, x.split) for x in used}

    class _Inputs:
        te_index_prefix = None
        bowtie2_build = acq.BOWTIE2.with_name("bowtie2-build")
        te_fasta = acq.TE_FASTA

    alignments = align_segments(segments, bowtie2=acq.BOWTIE2, te_index_prefix=_te_index_prefix(_Inputs, acq.CACHE),
                                work_directory=acq.CACHE / "te_segments")
    record_families = {k: normalize_family(v) for k, v in load_heuristic_family_map(acq.FAMMAP).items()}
    certifications = certify_groups(used, segments, alignments, record_families)
    te_sequences = {h.split()[0]: s for h, s in fasta_records_full_header(acq.TE_FASTA)}
    decisions = crossvalidate(certifications, te_sequences, record_families)
    group_pass = {row.group.group_hash: row.passes or decisions[row.group.group_hash].crossvalidated_pass for row in certifications}
    rescued = rescue_calls(SAMPLE, groups.ledger, reaggregate(groups.ledger, group_pass))
    calls = (*split.calls, *mates.calls, *rescued)
    data = legacy_tsv_bytes(calls)
    (acq.CACHE / "ISO1_INCUMBENT_COMBINED_CALLS.tsv").write_bytes(data)
    frozen = V2 / "V2_R9E_REPEAT_AWARE_DIVERGENCE" / "R9E_COMBINED_CALLS.tsv"
    sha = hashlib.sha256(data).hexdigest()
    record("r9e_final", {"calls": len(calls), "used_groups": len(used), "rescued": len(rescued), "sha256": sha,
                         "bytes_equal_frozen_file": data == frozen.read_bytes(), "sha256_equal_frozen": sha == FROZEN_R9E})


def main(argv):
    started = time.time()
    step = argv[0]
    if step == "sj":
        _chunk("split_junction", load("loci"), int(argv[1]))
    elif step == "mr":
        _chunk("mate_rescue", eligible()[0], int(argv[1]))
    elif step == "sj-auto":
        _auto("split_junction", load("loci"), started)
    elif step == "mr-auto":
        _auto("mate_rescue", eligible()[0], started)
    elif step == "status":
        print(json.dumps(json.loads((acq.CACHE / "ISO1_ACQUISITION_CHECKS.json").read_text()), indent=1))
    else:
        {"prep": step_prep, "sj-final": step_sj_final, "mr-final": step_mr_final, "route": step_route,
         "strict": step_strict, "groups": step_groups, "final": step_final}[step]()
    print(f"[{step}] {time.time() - started:.1f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
