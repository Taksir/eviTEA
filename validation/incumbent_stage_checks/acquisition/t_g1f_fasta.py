from stages import *
from evitea.incumbent.host_mate_gate import *
from evitea.incumbent.external import write_fasta
pa=pickle.load(open(WORK/"g1e.pkl","rb")); idx=CachedIndex()
t=time.time(); req=one_mate_host_requests(pa, idx); print("req", len(req), time.time()-t)
dol=idx.mate_dollars(req); print("dollars", time.time()-t)
recs=[]
for pid in sorted(req):
    recs.append((f"p{pid}|h{req[pid]+1}", idx.recover(dol[pid])))
    if time.time()-t>160: print("TIMEOUT"); break
print("recovered", len(recs), time.time()-t)
write_fasta(WORK/"g1f_host.fa", recs)
