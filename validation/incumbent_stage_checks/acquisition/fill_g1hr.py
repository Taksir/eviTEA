from stages import *
anchors=pickle.load(open(WORK/"g1hr_anchors.pkl","rb"))
fill_intervals([a.query_seq for a in anchors]); fill_dollars([a.query_seq for a in anchors], budget=150)
