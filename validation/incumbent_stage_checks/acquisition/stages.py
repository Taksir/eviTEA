from icache import *
from evitea.incumbent.anchor_census import *
from evitea.incumbent.sparse_panel import *
def catalogue():
    cat=build_anchor_catalogue(parse_te_consensus(TE_FASTA), load_te_family_map(FAMMAP))
    census_anchor_intervals(cat, CachedIndex()); return cat
from evitea.incumbent.coverage_panel import *
from evitea.incumbent.coverage_panel import _short_positions, _h1_family, SHORT_K_VALUES
def short_queries(cat):
    out=set()
    for te,s in cat.records.items():
        for k in SHORT_K_VALUES:
            for st in _short_positions(len(s),k):
                src=s[st:st+k]
                if set(src)<=set("ACGT"): out.add(src); out.add(revcomp(src))
    return sorted(out)
from evitea.incumbent.sequence import revcomp
