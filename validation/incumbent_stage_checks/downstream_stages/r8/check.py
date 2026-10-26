import sys, os, csv, json, hashlib, time
sys.path.insert(0,"src"); sys.path.insert(0,"."); sys.path.insert(0,os.path.expanduser("~/scratch/r8"))
from pathlib import Path
from load import load, rd, R
from evitea.incumbent import te_certification as tc, repeat_rescue as rr
from evitea.incumbent.events import IncumbentCall, IncumbentCallClass, legacy_tsv_bytes
from evitea.incumbent.te_library import load_heuristic_family_map, normalize_family
from evitea.incumbent.sequence import fasta_records_full_header
D=Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1"))
ledger,winners,acc=load(); groups=tc.used_groups(ledger,winners,acc)
reads={int(k):tuple(v) for k,v in json.load(open(os.path.expanduser("~/scratch/r8/reads.json"))).items()}
segments={g.group_hash:tc.full_te_segment(reads[g.source_pair_id][g.source_mate-1],g.read_view,g.junction_semantic,g.split) for g in groups}
print("segment lengths ok",all(len(segments[g.group_hash])==g.te_bp for g in groups))
t=time.time()
al=tc.align_segments(segments,bowtie2=os.path.expanduser("~/bowtie2-2.5.5-linux-x86_64/bowtie2"),te_index_prefix=os.path.expanduser("~/scratch/teidx/flybase_te"),work_directory=Path(os.path.expanduser("~/scratch/r8/work")))
print("align secs",round(time.time()-t,1))
# compare with recorded T0 rows
hist=set()
for line in open(R/"V2_R8D_BATCHED_TE_CERTIFICATION_REFINEMENT"/"V2R8D_TRIM_VARIANTS.sam",encoding="utf-8"):
    if line.startswith("@") or not line.startswith("T0__"): continue
    p=line.rstrip("\r\n").split("\t"); fl=int(p[1])
    if fl&4 or fl&0x800: continue
    nm=[x for x in p[11:] if x.startswith("NM:i:")]
    if not nm: continue
    hist.add((p[0][4:],p[2],bool(fl&16),int(p[3]),p[5],p[9].upper(),int(nm[0][5:])))
mine={(g,a.record_id,a.reverse,a.pos1,a.cigar,a.sequence,a.nm) for g,rows in al.items() for a in rows}
print("T0 rows equal",mine==hist,len(mine),len(hist))
fmap={k:normalize_family(v) for k,v in load_heuristic_family_map(D/"fbte_to_flybase_family.tsv").items()}
certs=tc.certify_groups(groups,segments,al,fmap)
h8={r["group_hash"]:r for r in rd(R/"V2_R8E_NORMALIZED_FAMILY_CERTIFICATION"/"V2R8E_GROUP_NORMALIZED_CERTIFICATION.tsv")}
print("R8E equal",all((str(c.passes),c.status,"" if c.best_nm is None else str(c.best_nm))==(h8[c.group.group_hash]["normalized_untrimmed_pass"],h8[c.group.group_hash]["normalized_untrimmed_status"],h8[c.group.group_hash]["normalized_untrimmed_best_nm"]) for c in certs))
te={h.split()[0]:s for h,s in fasta_records_full_header(D/"TE consensus library"/"transposon_sequence_set.fa")}
dec=rr.crossvalidate(certs,te,fmap)
cv={r["group_hash"]:r for r in rd(R/"V2_R9E_REPEAT_AWARE_DIVERGENCE"/"R9E_CROSSVALIDATION_GROUP_LEDGER.tsv")}
print("CV equal",all((str(d.crossvalidated_pass),d.status if d.status!="NOT_ELIGIBLE" else "NOT_ELIGIBLE")==(cv[g]["crossvalidated_divergence_pass"],cv[g]["cv_status"]) for g,d in dec.items()), sum(d.crossvalidated_pass for d in dec.values()))
gp={c.group.group_hash:(c.passes or dec[c.group.group_hash].crossvalidated_pass) for c in certs}
agg=rr.reaggregate(ledger,gp); resc=rr.rescue_calls("AKA017",ledger,agg)
core=[r for r in rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_COMBINED_CALLS.tsv") if not r["event_id"].startswith("V2R7_MIXED_")]
def ids(s): return tuple(int(x) for x in s.split(",") if x)
corec=[]
for r in core:
    sj=r["event_id"].startswith("V2FCRC1")
    corec.append(IncumbentCall("AKA017",r["event_id"],IncumbentCallClass.SPLIT_JUNCTION if sj else IncumbentCallClass.MATE_RESCUE,r["chr"],int(r["start"]),int(r["end"]),r["family"],int(r["left_coord"]) if sj else None,int(r["right_coord"]) if sj else None,int(r["delta"]) if sj else None,int(r["left_support"]),int(r["right_support"]),r["evidence_type"],float(r["score"]),ids(r["fragment_ids"])))
b=legacy_tsv_bytes(corec+list(resc))
print("R9E rescues",len(resc),"final sha",hashlib.sha256(b).hexdigest())
