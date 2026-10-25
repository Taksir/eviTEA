import sys, os, csv, json
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from collections import defaultdict
from evitea.incumbent import mate_rescue as mr, junction_groups as jg
from validation.legacy_adapters.fastq_random_access import pair_reads
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017")); R=A/"route_b_r9e_validation_full55x"/"shadow_v2"/"reports"
def rd(p): return list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"))
def ids(s): return tuple(int(x) for x in s.split(",") if x)
ledger=[mr.MateRescueLocus(r["locus_id"],r["chr"],r["family"],float(r["candidate_start"]),float(r["candidate_end"]),float(r["candidate_midpoint"]),ids(r["left_fragment_ids"]),ids(r["right_fragment_ids"]),r["pass_two_sided_rescue"]=="True") for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_LOCUS_LEDGER.tsv")]
targets=mr.one_sided_targets(ledger)
strict=[mr.StrictPairWinner(r["locus_id"],r["side"],int(r["pair_id"])) for r in rd(R/"V2_REFERENCE_DISCRIMINATIVE_RESCUE_R6"/"V2R6_FRAGMENT_COMPETITION.tsv") if r["kept"]=="True"]
hyp=defaultdict(list)
H=rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_JUNCTION_GROUP_HYPOTHESES.tsv")
for r in H:
    hyp[r["group_hash"]].append(jg.GroupHypothesis(r["group_hash"],int(r["representative_pair_id"]),ids(r["group_pair_ids"]),r["locus_id"],r["side"],r["family"],r["canonical_signature"],int(r["source_pair_id"]),float(r["proposal_score"]),int(r["te_seed_hits"]),float(r["te_seed_fraction"]),int(r["raw_msbwt_rows"]),float(r["host_candidate_distance"])))
counts={}
winners,led=jg.compete_and_link(hyp,targets,strict,counts)
kept={r["group_hash"]:(r["locus_id"],r["side"],r["source_pair_id"],r["canonical_signature"]) for r in rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_JUNCTION_FRAGMENT_COMPETITION.tsv") if r["kept"]=="True"}
print("winners equal",{g:(w.locus_id,w.side,str(w.source_pair_id),w.canonical_signature) for g,w in winners.items()}==kept,len(kept))
hl=rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_LOCUS_LEDGER.tsv")
j=lambda t:",".join(map(str,t))
mine=[(x.locus_id,x.chrom,x.family,j(x.existing_left_ids),j(x.existing_right_ids),j(x.strict_left_ids),j(x.strict_right_ids),",".join(x.left_group_hashes),",".join(x.right_group_hashes),j(x.left_representatives),j(x.right_representatives),str(x.candidate_start),str(x.candidate_end)) for x in led]
hist=[(r["locus_id"],r["chr"],r["family"],r["existing_r4_left_ids"],r["existing_r4_right_ids"],r["strict_r6_new_left_ids"],r["strict_r6_new_right_ids"],r["r7_new_junction_left_group_hashes"],r["r7_new_junction_right_group_hashes"],r["r7_new_junction_left_representative_ids"],r["r7_new_junction_right_representative_ids"],r["candidate_start"],r["candidate_end"]) for r in hl]
print("ledger equal",mine==hist,len(mine))
if mine!=hist:
    for a,b in zip(mine,hist):
        if a!=b: print(a);print(b);break
# check the group hash on a sample of pairs
sample=[(r["group_hash"],ids(r["group_pair_ids"])[0]) for r in H[::2000]]
rs=pair_reads(A/"fastq"/"SRR9951105_1.fastq",A/"fastq"/"SRR9951105_2.fastq",[p for _,p in sample])
print("group hash", all(jg.pair_group_hash(*rs[p])==g for g,p in sample), len(sample))
