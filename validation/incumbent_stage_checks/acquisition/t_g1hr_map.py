from stages import *
from evitea.incumbent.recall_repair import *
import numpy as np
ev=pickle.load(open(WORK/"g1hr_ev.pkl","rb")); req=rescue_host_requests(ev)
m=parse_rescue_host_sam(WORK/"g1hr_host.sam", req)
h=ztsv("g1hr","G1HR_NEW_HOST_MAPPING.tsv")
def f(x): return "" if x is None else str(x)
def g(x): return x[:-2] if x.endswith(".0") else x
mine=[(str(r.pair_id),r.status,str(r.alignment_count),r.chrom,f(r.pos),f(r.midpoint),f(r.mapq),r.strand,r.cigar) for r in m.values()]
his=[(r["pair_id"],r["status"],r["alignment_count"],r["chrom"],g(r["pos"]),g(r["midpoint"]),g(r["mapq"]),r["strand"],r["cigar"]) for r in h]
d=[(a,b) for a,b in zip(mine,his) if a!=b]; print("mapping equal", mine==his, len(mine), len(his), d[:3])
sup=supported_rescue_pairs(m); m1,host,work,fams,ppairs=pickle.load(open(WORK/"g1fg.pkl","rb"))
rep=np.union1d(work, np.asarray(sorted(sup),dtype=np.uint64))
print("repaired equal", np.array_equal(rep, znpy("g1hr","G1HR_REPAIRED_CANDIDATE_PAIR_IDS.npy")), len(rep))
pickle.dump((m,rep),open(WORK/"g1hr_map.pkl","wb"))
