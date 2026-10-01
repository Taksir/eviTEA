import sys, os, csv, hashlib
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from collections import defaultdict
from evitea.incumbent import split_junction as sj
from evitea.incumbent.te_library import QualifiedObservation, cluster_observations
from evitea.incumbent.events import legacy_tsv_bytes
R=Path(os.path.expandvars("$EVITEA_DATA/AKA017/route_b_r9e_validation_full55x/shadow_v2/reports/V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"))
loci={r["locus_id"]:r for r in csv.DictReader(open(R/"V2FCRC1_CANDIDATE_LOCI.tsv"),delimiter="\t")}
ctx=defaultdict(list)
for r in csv.DictReader(open(os.path.expanduser("~/scratch/fcrc1/pass_rows.tsv")),delimiter="\t"):
    ref=r["ref_full_edits"]; cont=r["continuation_edits"]
    cj=sj.CertifiedJunction(r["locus_id"],loci[r["locus_id"]]["chr"],int(r["pair_id"]),int(r["mate"]),r["view"],"L" if r["side_kind"]=="HOST_LEFT_TE_RIGHT" else "R",int(r["breakpoint"]),r["candidate_family"],r["verified_family"],"?","?",int(r["te_bp"]),int(r["host_bp"]),int(r["host_edits"]),int(r["te_edits"]),int(r["split_cost"]),None if ref=="" else int(ref),int(r["ref_advantage"]),None if cont=="" else int(cont),float(r["candidate_distance"]),0.0)
    ctx[cj.context_key].append(cj)
winners=defaultdict(list)
for rows in ctx.values():
    w,_=sj.resolve_context(rows)
    if w: winners[w.pair_id].append(w)
final=[]
for p in sorted(winners):
    k,_=sj.resolve_fragment(winners[p]); final.extend(k)
obs=[QualifiedObservation(str(x.pair_id),x.chrom,x.event_side,x.breakpoint,x.te_bp,x.verified_family,x.fbte_id,x.te_edits,x.te_orientation) for x in final]
events=sj.pair_side_clusters(cluster_observations(obs),"AKA017")
calls=[e for e in events if e.score>=4.0]
print(len(ctx),len(final),len(events),len(calls))
for name,rows in (("V2FCRC1_ALL_EVENTS.tsv",events),("V2FCRC1_CALLS_SCORE4.tsv",calls)):
    mine=legacy_tsv_bytes(rows); hist=(R/name).read_bytes()
    print(name, hashlib.sha256(mine).hexdigest()==hashlib.sha256(hist).hexdigest(), len(mine), len(hist))
