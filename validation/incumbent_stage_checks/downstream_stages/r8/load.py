import sys, os, csv, json
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from collections import defaultdict
from evitea.incumbent import mate_rescue as mr, junction_groups as jg, te_certification as tc
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017")); R=A/"route_b_r9e_validation_full55x"/"shadow_v2"/"reports"
def rd(p): return list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"))
def ids(s): return tuple(int(x) for x in s.split(",") if x)
def load():
    hl=rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_LOCUS_LEDGER.tsv")
    ledger=tuple(jg.JunctionLocus(r["locus_id"],r["chr"],r["family"],float(r["candidate_start"]),float(r["candidate_end"]),ids(r["existing_r4_left_ids"]),ids(r["existing_r4_right_ids"]),ids(r["strict_r6_new_left_ids"]),ids(r["strict_r6_new_right_ids"]),tuple(x for x in r["r7_new_junction_left_group_hashes"].split(",") if x),tuple(x for x in r["r7_new_junction_right_group_hashes"].split(",") if x),ids(r["r7_new_junction_left_representative_ids"]),ids(r["r7_new_junction_right_representative_ids"])) for r in hl)
    winners={}
    for r in rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_JUNCTION_FRAGMENT_COMPETITION.tsv"):
        if r["kept"]=="True":
            winners[r["group_hash"]]=jg.GroupHypothesis(r["group_hash"],int(r["representative_pair_id"]),ids(r["group_pair_ids"]),r["locus_id"],r["side"],r["family"],r["canonical_signature"],int(r["source_pair_id"]),float(r["proposal_score"]),int(r["te_seed_hits"]),float(r["te_seed_fraction"]),int(r["raw_msbwt_rows"]),float(r["host_candidate_distance"]))
    acc=[]
    for r in rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_DISCOVERY_JUNCTIONS.tsv"):
        if r["prequery_status"]!="R7_INSERTION_SPECIFIC_JUNCTION_ACCEPTED_PREQUERY": continue
        acc.append(jg.JunctionProposal(r["qid"],r["locus_id"],r["chrom"],r["family"],r["missing_side"],float(r["candidate_start"]),float(r["candidate_end"]),float(r["candidate_midpoint"]),int(r["source_pair_id"]),int(r["source_mate"]),r["read_view"],r["junction_semantic"],int(r["split"]),int(r["te_bp"]),int(r["te_seed_hits"]),float(r["te_seed_fraction"]),int(r["host_seed_matches_in_scope"]),int(r["host_seed_first_pos"]),r["signature"],r["canonical_signature"],r["host_half"],float(r["proposal_score"])))
    return ledger,winners,acc
if __name__=="__main__":
    ledger,winners,acc=load()
    groups=tc.used_groups(ledger,winners,acc)
    r8c={r["group_hash"]:r for r in rd(R/"V2_R8C_FULL_TE_SEGMENT_FAMILY_CERTIFICATION"/"V2R8C_GROUP_FAMILY_CERTIFICATION.tsv")}
    ok=[g.group_hash for g in groups]==list(r8c) and all((str(g.source_pair_id),str(g.source_mate),g.read_view,g.junction_semantic,str(g.split),str(g.te_bp),g.candidate_family)==(r8c[g.group_hash]["source_pair_id"],r8c[g.group_hash]["source_mate"],r8c[g.group_hash]["read_view"],r8c[g.group_hash]["junction_semantic"],r8c[g.group_hash]["split"],r8c[g.group_hash]["recorded_te_bp"],r8c[g.group_hash]["candidate_family"]) for g in groups)
    print("used groups equal",ok,len(groups))
    json.dump(sorted({g.source_pair_id for g in groups}),open(os.path.expanduser("~/scratch/r8/pairs.json"),"w"))
