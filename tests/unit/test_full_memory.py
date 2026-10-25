"""Full-memory layer behavior on small synthetic inputs."""

from __future__ import annotations

from evitea.full_memory.finite_state import finite_host_alternatives, pair_accounting
from evitea.full_memory.proposals import ProposalRecord, discover_regions, full_data_proposals
from evitea.full_memory.selection import full_memory_events, select_anchors
from evitea.full_memory.strict_a import a_groups, allowed_pairs, build_b_index, derive_witnesses, exact_anchors, full_data_events


def test_host_cap_censoring_and_classes(tmp_path):
    sam = tmp_path / "h.sam"
    rows = [("P1", 0, "c", 100, 42), ("P2", 0, "c", 100, 42), ("P2", 256, "c", 900, 1), ("P3", 0, "c", 5, 3)]
    sam.write_text("".join(f"{q}\t{f}\t{c}\t{p}\t{m}\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n" for q, f, c, p, m in rows))
    host, summary = finite_host_alternatives(sam, {1, 2, 3, 4})
    assert summary[1] == {"host_alt_count": 1, "host_reported_alignment_count": 1, "host_cap_censored": 0, "host_class": "UNIQUE_HIGH"}
    assert summary[2]["host_cap_censored"] == 1 and summary[2]["host_class"] == "MULTI_OR_LOW"
    assert summary[4]["host_class"] == "UNMAPPED"
    acc = pair_accounting(
        [{"pair_id": p, "family": "F", "families": "F", "te_mate": 1, "host_mate": 2, "eligibility_state": "NO_MULTI_RECORD_GE2_BLOCKS"} for p in (1, 3, 4)]
        + [{"pair_id": 9, "family": "", "families": "F;G", "te_mate": 0, "host_mate": 0, "eligibility_state": "FAMILY_CONFLICT"}],
        summary, {p: {"proper_signature_count": 0, "complete_signature_count": 0, "reference_state": "REFERENCE_UNEXPLAINED",
                      "reference_cap_censored": 0} for p in (1, 3, 4)})
    states = {r["pair_id"]: (r["experiment_state"], r["exclusion_reason"]) for r in acc}
    assert states == {1: ("ELIGIBLE", ""), 3: ("HOST_UNRESOLVED", "LOW_CONFIDENCE_SINGLE_HOST_PLACEMENT"),
                      4: ("UNMAPPED", "HOST_UNMAPPED"), 9: ("TE_SOURCE_UNRESOLVED", "FAMILY_CONFLICT")}
    assert {h["pair_id"] for h in host} == {1, 2, 3}


def test_proposals_primary_chromosomes_and_two_fold_segments():
    recs = [ProposalRecord(0, 0, "F", "chr2L", 100, 905), ProposalRecord(6, 1, "F", "chr2L", 500, 1305),
            ProposalRecord(5, 0, "F", "chr2L", 600, 1405),  # fold 0 again
            ProposalRecord(1, 1, "F", "chr4", 100, 905), ProposalRecord(2, 2, "F", "chr4", 100, 905)]
    (region,) = discover_regions(recs, {0, 1, 2, 3, 4})
    assert (region["lo"], region["hi"], sorted(region["pair_ids"])) == (500, 1305, [0, 5, 6])
    assert full_data_proposals("S", recs)[0]["candidate_id"] == "S_FULL_000001"


def test_strict_a_selection_and_full_memory_event():
    accounting = [{"pair_id": p, "fold": p % 5, "family": "F", "experiment_state": "ELIGIBLE",
                   "reference_state": "REFERENCE_UNEXPLAINED", "reference_cap_censored": 0} for p in (0, 1, 7)]
    refs = [{"pair_id": p, "family": "F", "reference_cap_censored": 0, "evidence_complete": 0, "chrom": "chr2L",
             "te_cigar": "40S85M", "te_start": 1000, "te_end": 1085, "unexplained_query_intervals": "30-51"} for p in (0, 1, 7)]
    witnesses = derive_witnesses(allowed_pairs(accounting), refs)
    assert {(w["breakpoint0"], w["anchor_host_side"]) for w in witnesses} == {(1000, "RIGHT")}
    (anchor,) = exact_anchors(a_groups(witnesses), {0, 1, 2, 3, 4})
    assert anchor["A_pair_ids"] == [0, 1, 7]
    b_rows = [{"pair_id": p, "fold": p % 5, "family": "F", "chrom": "chr2L", "start": s, "end": s + 100, "strand": st}
              for p, s, st in ((10, 500, "+"), (21, 600, "+"), (12, 1200, "-"), (33, 1300, "-"))]
    index, _ = build_b_index(b_rows)
    props = [{"candidate_id": "S_FULL_000001", "family": "F", "chrom": "chr2L", "candidate_lo": 900, "candidate_hi": 1100}]
    ledger = full_data_events("S", props, a_groups(witnesses), index, old_keys=set())
    assert ledger[0]["proposal_only_A"] == 1 and ledger[0]["B_state"] == "BILATERAL"
    (event,) = full_memory_events(select_anchors(ledger, "S"), index)
    assert event.candidate_id == "S_FULL_000001" and event.b_state == "BILATERAL"
    assert full_data_events("S", props, a_groups(witnesses), index, old_keys={("F", "chr2L", 1000, "RIGHT")})[0]["proposal_only_A"] == 0
