from stages import *
from evitea.incumbent.recall_repair import *
from evitea.incumbent.external import write_fasta
import numpy as np
anchors=pickle.load(open(WORK/"g1hr_anchors.pkl","rb")); pa=pickle.load(open(WORK/"g1e.pkl","rb"))
idx=CachedIndex(with_dollars=True)
t=time.time(); ev=enumerate_rescue_pairs(anchors, pa.pair_ids, idx); print("enum", len(ev), time.time()-t)
h=ztsv("g1hr","G1HR_NEW_PAIR_EVIDENCE.tsv")
print("evidence equal", [(str(e.pair_id),str(int(e.mate1_hit)),str(int(e.mate2_hit)),",".join(e.families)) for e in ev.values()]==[(r["pair_id"],r["mate1_hit"],r["mate2_hit"],r["families"]) for r in h], len(h))
req=rescue_host_requests(ev); dol=idx.mate_dollars(req)
write_fasta(WORK/"g1hr_host.fa", [(f"pair_{p}", idx.recover(dol[p])) for p in sorted(req)])
pickle.dump(ev,open(WORK/"g1hr_ev.pkl","wb")); print("done", time.time()-t)
