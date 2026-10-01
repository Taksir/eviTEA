from common import *
from evitea.incumbent.anchor_census import *
t=time.time()
cat=build_anchor_catalogue(parse_te_consensus(TE_FASTA), load_te_family_map(FAMMAP)); print("built", time.time()-t, len(cat.queries), len(cat.sources))
hist=ztsv("g1a","G1A_ANCHOR_INTERVALS.tsv")
print(len(hist))
bad=0
for q,h in zip(cat.queries,hist):
    mine=(q.query_id,str(q.k),q.query_seq,",".join(q.source_te_ids),",".join(q.source_families),",".join(q.source_regions),str(int(q.has_terminal)),str(int(q.has_internal)),str(int(q.low_complexity)),str(q.source_family_count))
    his=(h["query_id"],h["k"],h["query_seq"],h["source_te_ids"],h["source_families"],h["source_regions"],h["has_terminal_provenance"],h["has_internal_provenance"],h["low_complexity"],h["source_family_count"])
    if mine!=his:
        bad+=1
        if bad<4: print(mine,his)
print("catalogue mismatches", bad, len(cat.queries)==len(hist))
prov=ztsv("g1a","G1A_ANCHOR_PROVENANCE.tsv")
bad=0
for s,h in zip(cat.sources,prov):
    mine=(s.query_id,s.te_id,s.family,str(s.k),str(s.start),s.region,s.orientation,s.query_seq,str(int(s.queryable)))
    his=(h["query_id"],h["te_id"],h["family"],h["k"],h["start"],h["region"],h["source_orientation"],h["query_seq"],h["queryable"])
    if mine!=his:
        bad+=1
        if bad<4: print(mine,his)
print("provenance mismatches", bad, len(cat.sources), len(prov))
