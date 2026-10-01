from stages import *
from evitea.incumbent.pair_aware_panel import *
import numpy as np
cat=catalogue(); panel=select_sparse_panel(cat); idx=CachedIndex(with_dollars=True)
short,h1=pickle.load(open(WORK/"g1c.pkl","rb")); rp=pickle.load(open(WORK/"g1d.pkl","rb"))
t=time.time(); pa=select_pair_aware_panel(cat, panel, short, rp, idx); print("g1e", time.time()-t)
hp=znpy("g1e","G1E_CANDIDATE_PAIR_IDS.npy"); hc=znpy("g1e","G1E_PAIR_READ_COUNTS.npy")
print("pairs equal", np.array_equal(hp, pa.pair_ids), np.array_equal(hc, pa.pair_read_counts), len(hp))
hs=ztsv("g1e","G1E_SELECTED_PANEL_SEQUENCES.tsv"); print("seqs equal", sorted(pa.selected_sequences)==[r["sequence"] for r in hs])
hg=ztsv("g1e","G1E_GROUP_COVERAGE.tsv")
print("groups equal", [(g.te_id,g.family,g.orientation,str(g.window_count),repr(g.selected_coverage)) for g in pa.groups]==[(r["te_id"],r["family"],r["orientation"],r["window_count"],r["g1e_selected_coverage"]) for r in hg])
pickle.dump(pa,open(WORK/"g1e.pkl","wb"))
