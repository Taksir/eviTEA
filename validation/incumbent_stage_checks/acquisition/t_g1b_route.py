from stages import *
cat=catalogue(); panel=select_sparse_panel(cat)
idx=CachedIndex(with_dollars=True)
t=time.time(); pairs=route_sparse_panel(panel, idx); print("route", time.time()-t, len(pairs))
hist=ztsv("g1b","G1B_CANDIDATE_PAIRS.tsv")
mine=[(str(p.pair_id),str(int(p.mate1_hit)),str(int(p.mate2_hit)),str(p.channel_mask),",".join(p.families)) for p in pairs]
his=[(h["pair_id"],h["mate1_hit"],h["mate2_hit"],h["channel_mask"],h["families"]) for h in hist]
print("pairs equal", mine==his, len(his))
