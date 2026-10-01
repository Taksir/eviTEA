from stages import *
from evitea.incumbent.host_mate_gate import *
import numpy as np
pa=pickle.load(open(WORK/"g1e.pkl","rb")); idx=CachedIndex()
req=one_mate_host_requests(pa, idx)
m=parse_primary_host_sam(WORK/"g1f_host.sam", req)
hist=ztsv("g1f","G1F_HOST_MATE_MAPPING.tsv")
def s(x): return "" if x is None else str(x)
mine=[(str(r.pair_id),r.status,str(r.alignment_count),r.chrom,s(r.pos),s(r.midpoint),s(r.mapq),r.strand,r.cigar) for r in m.values()]
his=[(h["pair_id"],h["status"],h["alignment_count"],h["chrom"],h["pos"],h["midpoint"],h["mapq"],h["strand"],h["cigar"]) for h in hist]
diff=[(a,b) for a,b in zip(mine,his) if a!=b]
print("mapping equal", mine==his, len(mine), len(his), len(diff), diff[:3])
cat=catalogue(); panel=select_sparse_panel(cat); ppairs=route_sparse_panel(panel, CachedIndex(with_dollars=True))
host=supported_host_cluster_pairs(m); work=working_candidate_pairs(host, ppairs)
print("working equal", np.array_equal(work, znpy("g1f","G1F_WORKING_CANDIDATE_PAIR_IDS.npy")), len(work))
fams=resident_pair_families(host, pa)
hf=ztsv("g1g","G1G_HOST_CLUSTER_PAIR_FAMILIES.tsv")
print("families equal", [(str(k),v) for k,v in fams.items()]==[(h["pair_id"],h["families"]) for h in hf], len(fams), len(hf))
pickle.dump((m,host,work,fams,ppairs),open(WORK/"g1fg.pkl","wb"))
