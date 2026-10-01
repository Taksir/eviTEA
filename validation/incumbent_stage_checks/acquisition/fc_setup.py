from stages import *
loci=pickle.load(open(WORK/"loci.pkl","rb"))
idx=CachedIndex()
t=time.time(); reads=idx.pair_reads(p for L in loci for p in L.pair_ids); print(len(reads), time.time()-t)
pickle.dump(reads,open(WORK/"loci_reads.pkl","wb"))
