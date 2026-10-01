import sys, os, json, time, csv
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from evitea.incumbent.candidates import CandidateLocus
from evitea.incumbent.te_library import load_te_records, load_record_family_map, LocalTEIndex, load_family_kmer_panels
from evitea.incumbent.sequence import load_reference_contigs
from evitea.incumbent import mate_rescue as mr
D=Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1")); out=Path(os.path.expanduser("~/scratch/r4"))
R=Path(os.path.expandvars("$EVITEA_DATA/AKA017/route_b_r9e_validation_full55x/shadow_v2/reports/V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"))
loci=[CandidateLocus(l[0],l[1],l[2],l[3],l[4],tuple(l[5])) for l in json.loads((out/"loci.json").read_text())]
reads={int(k):tuple(v) for k,v in json.loads((out/"reads.json").read_text()).items()}
te=D/"TE consensus library"/"transposon_sequence_set.fa"; fmap=D/"fbte_to_flybase_family.tsv"
panels,_=load_family_kmer_panels(te,fmap,{L.family for L in loci},21,1.10)
index=LocalTEIndex(load_te_records(te,load_record_family_map(fmap)))
ref=load_reference_contigs(D/"Reference"/"dm6_TEMP2_compatible.fa",{L.chrom for L in loci})
got=[]
t=time.time()
for L in loci:
    contig=ref[L.chrom]; c=int(round(L.midpoint)); lo=max(0,c-800); scope=contig[lo:min(len(contig),c+801)]
    k=panels.get(L.family,set())
    for pid in L.pair_ids:
        m1,m2=reads[pid]
        for te_mate,te_read,host_mate,host_read in ((1,m1,2,m2),(2,m2,1,m1)):
            h,_=mr.best_local_host_anchor(host_read,scope,lo,L.start,L.end)
            if h is None: continue
            side=mr.classify_side(h.view,h.midpoint,L.midpoint)
            if side is None: continue
            w,_=mr.best_te_mate_window(te_read,L.family,k,index,scope)
            if w is None: continue
            got.append((L.locus_id,str(pid),side,str(te_mate),str(host_mate),h.view,str(h.start),str(h.end),str(h.midpoint),str(h.edits),str(h.candidate_distance),str(w.start),str(w.end),str(w.seed_hits),str(w.seed_total),str(w.seed_fraction),w.family,w.fbte_id,str(w.te_edits),w.te_orientation))
print("secs",round(time.time()-t,1))
want=set(out.joinpath("loci.txt").read_text().split())
hist=[(r["locus_id"],r["pair_id"],r["side"],r["te_mate"],r["host_mate"],r["host_view"],r["host_start"],r["host_end"],r["host_midpoint"],r["host_edits"],r["host_candidate_distance"],r["te_window_start"],r["te_window_end"],r["te_seed_hits"],r["te_seed_total"],r["te_seed_fraction"],r["verified_family"],r["fbte_id"],r["te_edits"],r["te_orientation"]) for r in csv.DictReader(open(R/"V2R4_CERTIFIED_FRAGMENT_EVIDENCE.tsv"),delimiter="\t") if r["locus_id"] in want]
print("equal",got==hist,len(got),len(hist))
if got!=hist:
    print(set(got)-set(hist)); print(set(hist)-set(got))
