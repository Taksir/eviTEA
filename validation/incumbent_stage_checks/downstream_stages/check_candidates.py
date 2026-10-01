import sys, os, time
sys.path.insert(0, os.path.expanduser("~/work/EviTEA/src")); sys.path.insert(0, os.path.expanduser("~/work/EviTEA"))
from pathlib import Path
from validation.legacy_adapters.route_b_artifacts import load_candidate_universe_inputs, read_tsv
from evitea.incumbent.candidates import reconstruct_candidate_pairs, cluster_candidate_loci
root = Path(sys.argv[1]); v2 = Path(sys.argv[2])
t=time.time(); inputs = load_candidate_universe_inputs(root)
cands, meta = reconstruct_candidate_pairs(**inputs); loci = cluster_candidate_loci(cands)
print(meta, len(loci), round(time.time()-t,1))
hist = read_tsv(v2/"V2_FIXED_CANDIDATE_REFERENCE_COMPETITION_R1"/"V2FCRC1_CANDIDATE_LOCI.tsv")
mine = [{"locus_id":L.locus_id,"chr":L.chrom,"candidate_family":L.family,"candidate_start_midpoint":f"{L.start:.1f}","candidate_end_midpoint":f"{L.end:.1f}","candidate_midpoint":f"{L.midpoint:.1f}","candidate_pair_count":str(len(L.pair_ids))} for L in loci]
print("loci equal:", mine == hist, len(hist))
