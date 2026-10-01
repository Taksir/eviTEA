from stages import *
cat=catalogue(); panel=select_sparse_panel(cat)
hist=ztsv("g1b","G1B_SELECTED_PANEL.tsv")
mine=[(p.query_id,p.query_seq,str(p.width),str(p.channel_mask),",".join(sorted(p.source_families))) for p in panel]
his=[(h["query_id"],h["query_seq"],h["interval_width"],h["channel_mask"],h["source_families"]) for h in hist]
print("panel equal", mine==his, len(mine), len(his), sum(p.width for p in panel))
import sys
if len(sys.argv)>1: fill_dollars([p.query_seq for p in panel], budget=150)
