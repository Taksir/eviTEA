from icache import *
from evitea.incumbent.anchor_census import *
cat=build_anchor_catalogue(parse_te_consensus(TE_FASTA), load_te_family_map(FAMMAP))
idx=CachedIndex(); census_anchor_intervals(cat, idx)
hist=ztsv("g1a","G1A_ANCHOR_INTERVALS.tsv")
bad=sum(1 for q,h in zip(cat.queries,hist) if (q.interval_lo,q.interval_hi)!=(int(h["interval_lo"]),int(h["interval_hi"])))
print("interval mismatches", bad)
