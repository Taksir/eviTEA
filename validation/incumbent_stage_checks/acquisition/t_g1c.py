from stages import *
cat=catalogue(); panel=select_sparse_panel(cat); idx=CachedIndex()
fam=cat.family_by_te()
t=time.time(); short=short_exact_census(cat.records, fam, idx); print("short", len(short), time.time()-t)
hist=ztsv("g1c","G1C_SHORT_EXACT_LOCUS_COUNTS.tsv")
mine=[(l.te_id,l.family,l.orientation,str(l.start),str(l.k),l.query_seq,str(l.width)) for l in short]
his=[(h["te_id"],h["family"],h["orientation"],h["start"],h["k"],h["query_seq"],h["interval_width"]) for h in hist]
print("short equal", mine==his, len(his))
base=baseline_anchor_loci(cat, panel)
h1=hard_gap_h1_loci(cat.records, fam, base, short)
hh=ztsv("g1c","G1C_H1_LOCUS_COUNTS.tsv")
print("h1 equal", [(l.te_id,l.family,l.orientation,str(l.start),l.source_21mer) for l in h1]==[(h["te_id"],h["family"],h["orientation"],h["start"],h["source_21mer"]) for h in hh], len(h1), len(hh))
pickle.dump((short,h1),open(WORK/"g1c.pkl","wb"))
mem=set(m for l in h1 for m in _h1_family(l.source_21mer)); print("h1 members", len(mem))
