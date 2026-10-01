from icache import *
from evitea.incumbent.anchor_census import *
cat=build_anchor_catalogue(parse_te_consensus(TE_FASTA), load_te_family_map(FAMMAP))
fill_intervals([q.query_seq for q in cat.queries], budget=160)
