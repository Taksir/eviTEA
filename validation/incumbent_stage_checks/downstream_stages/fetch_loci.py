import sys, os, json, time
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from validation.legacy_adapters.route_b_artifacts import load_candidate_universe_inputs
from validation.legacy_adapters.fastq_random_access import SraFastqIndex
from evitea.incumbent.candidates import reconstruct_candidate_pairs, cluster_candidate_loci
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017/route_b_r9e_validation_full55x"))
out=Path(sys.argv[1])
wanted=set(out.joinpath("loci.txt").read_text().split())
if not out.joinpath("loci.json").exists():
    c,_=reconstruct_candidate_pairs(**load_candidate_universe_inputs(A/"shadow_v1"/"reports"))
    loci=[L for L in cluster_candidate_loci(c) if L.locus_id in wanted]
    out.joinpath("loci.json").write_text(json.dumps([[L.locus_id,L.chrom,L.family,L.start,L.end,list(L.pair_ids)] for L in loci]))
loci=json.loads(out.joinpath("loci.json").read_text())
cache=out/"reads.json"; reads=json.loads(cache.read_text()) if cache.exists() else {}
F=Path(os.path.expandvars("$EVITEA_DATA/AKA017/fastq"))
a,b=SraFastqIndex(F/"SRR9951105_1.fastq"),SraFastqIndex(F/"SRR9951105_2.fastq")
t=time.time()
for L in loci:
    for p in L[5]:
        if str(p) in reads: continue
        reads[str(p)]=[a.sequence(p),b.sequence(p)]
        if time.time()-t>150: break
    if time.time()-t>150: break
cache.write_text(json.dumps(reads)); need=sum(len(L[5]) for L in loci); print(len(reads),"/",need)
