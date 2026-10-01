from stages import *
cat=catalogue(); panel=select_sparse_panel(cat); idx=CachedIndex()
short,h1=pickle.load(open(WORK/"g1c.pkl","rb"))
t=time.time(); rp=optimize_rescue_panel(cat, panel, short, h1, idx); print("opt", time.time()-t, rp.cap)
hist=[r for r in ztsv("g1d","G1D_SELECTED_ACTIONS.tsv") if r["B"]=="1024"]
mine=[a.action_id for a in rp.selected_actions]
print("actions equal", mine==[h["action_id"] for h in hist], len(mine), len(hist))
gc=[r for r in ztsv("g1d","G1D_GROUP_COVERAGE.tsv") if r["B"]=="1024"]
m2={k:(repr(v.selected_coverage),str(int(v.selected_terminal5)),str(int(v.selected_terminal3)),str(v.window_count)) for k,v in rp.group_targets.items()}
h2={(r["te_id"],r["orientation"]):(r["selected_coverage"],r["selected_terminal5"],r["selected_terminal3"],r["window_count"]) for r in gc}
print("groups equal", m2==h2, [ (k,m2[k],h2.get(k)) for k in m2 if m2[k]!=h2.get(k)][:3])
cs=ztsv("g1d","G1D_CHOSEN_PANEL_SEQUENCES.tsv")
print("seqs equal", [(s,str(w),",".join(sorted(rp.selected_sources[s]))) for s,w in sorted(rp.selected_sequences.items())]==[(r["sequence"],r["interval_width"],r["sources"]) for r in cs], len(cs))
pickle.dump(rp,open(WORK/"g1d.pkl","wb"))
