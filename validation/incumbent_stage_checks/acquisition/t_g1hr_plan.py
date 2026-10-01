from stages import *
from evitea.incumbent.recall_repair import *
cat=catalogue(); pa=pickle.load(open(WORK/"g1e.pkl","rb"))
anchors=plan_rescue_anchors(cat, pa.groups)
h=ztsv("g1hr","G1HR_RESCUE_ANCHORS.tsv")
print("anchors equal", [(a.query_id,a.query_seq,str(a.k),str(a.width),a.family) for a in anchors]==[(r["query_id"],r["query_seq"],r["k"],r["interval_width"],r["family"]) for r in h], len(anchors), len(h), sum(a.width for a in anchors))
pickle.dump(anchors,open(WORK/"g1hr_anchors.pkl","wb"))
