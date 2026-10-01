from stages import *
cat=build_anchor_catalogue(parse_te_consensus(TE_FASTA), load_te_family_map(FAMMAP))
qs=short_queries(cat); print(len(qs)); fill_intervals(qs, budget=150)
