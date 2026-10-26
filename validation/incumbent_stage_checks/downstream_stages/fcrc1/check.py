import sys, os, json, time, csv
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from dataclasses import asdict
from evitea.incumbent.candidates import CandidateLocus
from evitea.incumbent.te_library import load_te_records, load_record_family_map, LocalTEIndex, load_family_kmer_panels
from evitea.incumbent.sequence import load_reference_contigs
from evitea.incumbent import split_junction as sj
D=Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1"))
out=Path(os.path.expanduser("~/scratch/fcrc1"))
loci=[CandidateLocus(l[0],l[1],l[2],l[3],l[4],tuple(l[5])) for l in json.loads((out/"loci.json").read_text())]
reads={int(k):tuple(v) for k,v in json.loads((out/"reads.json").read_text()).items()}
te=D/"TE consensus library"/"transposon_sequence_set.fa"; fmap=D/"fbte_to_flybase_family.tsv"
t=time.time()
panels,_=load_family_kmer_panels(te,fmap,{L.family for L in loci},21,1.10)
index=LocalTEIndex(load_te_records(te,load_record_family_map(fmap)))
ref=load_reference_contigs(D/"Reference"/"dm6_TEMP2_compatible.fa",{L.chrom for L in loci})
print("setup",round(time.time()-t,1))
rows=[]
def sink(n,reason,cj,diag):
    r={**{k:v for k,v in asdict(n).items()},"status":"PASS" if cj else "REJECT","reason":reason,"breakpoint":"","verified_family":"","host_edits":"","te_edits":"","split_cost":"","ref_full_edits":"","ref_advantage":"","continuation_edits":"","candidate_distance":""}
    if cj is not None:
        r.update({"breakpoint":cj.breakpoint,"verified_family":cj.verified_family,"host_edits":cj.host_edits,"te_edits":cj.te_edits,"split_cost":cj.split_cost,"ref_full_edits":"" if cj.ref_full_edits is None else cj.ref_full_edits,"ref_advantage":cj.ref_advantage,"continuation_edits":"" if cj.continuation_edits is None else cj.continuation_edits,"candidate_distance":f"{cj.candidate_distance:.3f}"})
    else: r.update(diag)
    rows.append(r)
t=time.time()
res=sj.run_split_junction_layer(sample_id="AKA017",loci=loci,reads=reads,reference=ref,family_kmers=panels,te_index=index,hypothesis_sink=sink)
print("run",round(time.time()-t,1),"hyp",len(rows))
fields=open(out/"ledger_subset.tsv").readline().rstrip("\n").split("\t")
hist=[r for r in csv.DictReader(open(out/"ledger_subset.tsv"),delimiter="\t")]
def norm(r): return tuple("" if r.get(f) is None else str(r.get(f)) for f in fields if f!="chr")
# the order of the loci in the recorded table is not checked
A=[norm(r) for r in rows]; B=[norm(r) for r in hist]
print("equal",A==B,len(A),len(B))
if A!=B:
    import itertools
    for i,(x,y) in enumerate(itertools.zip_longest(A,B)):
        if x!=y: print(i);print(x);print(y);break
