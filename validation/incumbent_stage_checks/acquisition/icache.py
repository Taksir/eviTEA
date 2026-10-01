"""Persistent interval / dollar caches for chunked VM runs."""
from common import *
import numpy as np
from MUSCython import MultiStringBWTCython as M
BWT_DIR=os.path.expanduser("~/msbwt_aka")
_bwt=None
def bwt():
    global _bwt
    if _bwt is None: _bwt=M.loadBWT(BWT_DIR, useMemmap=True)
    return _bwt
IPATH=WORK/"intervals.pkl"; DPATH=WORK/"dollars"
def load_intervals():
    return pickle.loads(IPATH.read_bytes()) if IPATH.exists() else {}
def save_intervals(d):
    tmp=IPATH.with_suffix(".tmp"); tmp.write_bytes(pickle.dumps(d)); tmp.replace(IPATH)
def fill_intervals(seqs, budget=150):
    d=load_intervals(); b=bwt(); t=time.time(); n=0
    todo=[s for s in seqs if s not in d]
    for s in todo:
        lo,hi=b.findIndicesOfStr(s); d[s]=(int(lo),int(hi)); n+=1
        if time.time()-t>budget: break
    save_intervals(d)
    left=sum(1 for s in seqs if s not in d)
    print(f"filled {n}, remaining {left}")
    return left
def load_dollars(): 
    DPATH.mkdir(exist_ok=True); out={}
    for p in sorted(DPATH.glob("*.pkl")): out.update(pickle.loads(p.read_bytes()))
    return out
def fill_dollars(seqs, budget=150):
    have=load_dollars(); d=load_intervals(); b=bwt(); t=time.time(); new={}
    for s in seqs:
        if s in have or s in new: continue
        lo,hi=d[s]
        new[s]=np.fromiter((int(b.getSequenceDollarID(r)) for r in range(lo,hi)),dtype=np.uint64,count=hi-lo)
        if time.time()-t>budget: break
    if new: (DPATH/f"{int(time.time()*1000)}.pkl").write_bytes(pickle.dumps(new))
    left=sum(1 for s in seqs if s not in have and s not in new)
    print(f"dollar filled {len(new)}, remaining {left}")
    return left
class CachedIndex:
    """MsbwtSample with prefilled caches (misses go to the real BWT)."""
    def __new__(cls, with_dollars=False):
        from evitea.incumbent.msbwt_queries import MsbwtSample
        return MsbwtSample(BWT_DIR, bwt=bwt(), interval_cache=load_intervals(), dollar_cache=load_dollars() if with_dollars else None)
