from stages import *
from evitea.incumbent.coverage_panel import _h1_family
short,h1=pickle.load(open(WORK/"g1c.pkl","rb"))
mem=sorted(set(m for l in h1 for m in _h1_family(l.source_21mer))); print(len(mem)); fill_intervals(mem, budget=150)
