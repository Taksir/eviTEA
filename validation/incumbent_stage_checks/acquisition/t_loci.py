from stages import *
from evitea.incumbent.acquisition import _host_record
from evitea.incumbent.candidates import *
m,host,work,fams,ppairs=pickle.load(open(WORK/"g1fg.pkl","rb")); ev=pickle.load(open(WORK/"g1hr_ev.pkl","rb")); rm,rep=pickle.load(open(WORK/"g1hr_map.pkl","rb"))
c,meta=reconstruct_candidate_pairs(final_pair_ids=rep.tolist(),working_pair_ids=work.tolist(),primary_host={p:_host_record(r) for p,r in m.items()},rescue_host={p:_host_record(r) for p,r in rm.items()},resident_families=fams,rescue_evidence={p:PairPanelEvidence(p,",".join(r.families),r.mate1_hit,r.mate2_hit) for p,r in ev.items()},panel_evidence={r.pair_id:PairPanelEvidence(r.pair_id,",".join(r.families),r.mate1_hit,r.mate2_hit) for r in ppairs})
loci=cluster_candidate_loci(c); print(meta, len(loci))
V2=Path(os.path.expandvars("$EVITEA_DATA/AKA017/route_b_r9e_validation_full55x/shadow_v2/reports"))
hist=list(csv.DictReader(open(V2/"V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"/"V2FCRC1_CANDIDATE_LOCI.tsv"),delimiter="\t"))
mine=[{"locus_id":L.locus_id,"chr":L.chrom,"candidate_family":L.family,"candidate_start_midpoint":f"{L.start:.1f}","candidate_end_midpoint":f"{L.end:.1f}","candidate_midpoint":f"{L.midpoint:.1f}","candidate_pair_count":str(len(L.pair_ids))} for L in loci]
print("loci equal", mine==[{k:h[k] for k in mine[0]} for h in hist], len(hist))
pickle.dump(loci,open(WORK/"loci.pkl","wb"))
