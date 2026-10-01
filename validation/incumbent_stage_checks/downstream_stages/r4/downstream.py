import sys, os, json, csv, hashlib, math
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from collections import defaultdict
from evitea.incumbent import mate_rescue as mr
from evitea.incumbent.events import IncumbentCall, IncumbentCallClass, legacy_tsv_bytes
from validation.legacy_adapters.route_b_artifacts import load_candidate_universe_inputs
from evitea.incumbent.candidates import reconstruct_candidate_pairs, cluster_candidate_loci
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017/route_b_r9e_validation_full55x"))
R=A/"shadow_v2"/"reports"
def rd(p): return list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"))
fc=rd(R/"V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"/"V2FCRC1_CALLS_SCORE4.tsv")
fcalls=[IncumbentCall("AKA017",r["event_id"],IncumbentCallClass.SPLIT_JUNCTION,r["chr"],int(r["start"]),int(r["end"]),r["family"],int(r["left_coord"]),int(r["right_coord"]),int(r["delta"]),int(r["left_support"]),int(r["right_support"]),r["evidence_type"],float(r["score"]),tuple(int(x) for x in r["fragment_ids"].split(","))) for r in fc]
print("fcrc1 roundtrip", legacy_tsv_bytes(fcalls)==(R/"V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"/"V2FCRC1_CALLS_SCORE4.tsv").read_bytes())
c,_=reconstruct_candidate_pairs(**load_candidate_universe_inputs(A/"shadow_v1"/"reports")); loci=cluster_candidate_loci(c)
guard=[L.locus_id for L in loci if mr.locus_redundant_with_calls(L,fcalls,50)]
hist_guard=[r["locus_id"] for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_CANDIDATE_LOCI.tsv") if r["guarded_by_fcrc1"]=="True"]
print("guard equal",guard==hist_guard,len(guard))
ev=defaultdict(list)
for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_CERTIFIED_FRAGMENT_EVIDENCE.tsv"):
    ev[int(r["pair_id"])].append(mr.MateEvidence(r["locus_id"],r["chrom"],r["candidate_family"],int(r["pair_id"]),r["side"],int(r["te_mate"]),int(r["host_mate"]),r["host_view"],int(r["host_start"]),int(r["host_end"]),float(r["host_midpoint"]),int(r["host_edits"]),float(r["host_candidate_distance"]),int(r["te_seed_hits"]),float(r["te_seed_fraction"]),int(r["te_edits"])))
winners=defaultdict(list)
for p in sorted(ev):
    w,_=mr.resolve_pair_evidence(ev[p])
    if w: winners[w.locus_id].append(w)
# reuse run_mate_rescue tail by constructing ledger same way
elig=[L for L in loci if L.locus_id not in set(guard)]
ledger=[];calls=[]
for L in sorted(elig,key=lambda x:x.locus_id):
    rows=winners.get(L.locus_id,[]); left=tuple(sorted({x.pair_id for x in rows if x.side=="LEFT"})); right=tuple(sorted({x.pair_id for x in rows if x.side=="RIGHT"}))
    ok=len(left)>=2 and len(right)>=2
    la=sorted(x.host_midpoint for x in rows if x.side=="LEFT"); ra=sorted(x.host_midpoint for x in rows if x.side=="RIGHT")
    ledger.append({"locus_id":L.locus_id,"chr":L.chrom,"family":L.family,"candidate_start":f"{L.start:.1f}","candidate_end":f"{L.end:.1f}","candidate_midpoint":f"{L.midpoint:.1f}","candidate_pair_count":str(len(L.pair_ids)),"left_fragments":str(len(left)),"right_fragments":str(len(right)),"total_distinct_fragments":str(len(set(left)|set(right))),"min_left_required":"2","min_right_required":"2","pass_two_sided_rescue":str(ok),"left_anchor_min":"" if not la else f"{min(la):.3f}","left_anchor_max":"" if not la else f"{max(la):.3f}","right_anchor_min":"" if not ra else f"{min(ra):.3f}","right_anchor_max":"" if not ra else f"{max(ra):.3f}","left_fragment_ids":",".join(map(str,left)),"right_fragment_ids":",".join(map(str,right))})
    if ok:
        u=tuple(sorted(set(left)|set(right)))
        calls.append(IncumbentCall("AKA017",f"V2R4_MATE_{len(calls)+1:06d}",IncumbentCallClass.MATE_RESCUE,L.chrom,int(math.floor(min(L.start,L.end))),int(math.ceil(max(L.start,L.end)))+1,L.family,None,None,None,len(left),len(right),"two_sided_mate_cluster_rescue",float(len(u)),u))
print("ledger equal", ledger==rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_LOCUS_LEDGER.tsv"))
print("combined equal", legacy_tsv_bytes(fcalls+calls)==(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_COMBINED_CALLS.tsv").read_bytes(), len(calls))
