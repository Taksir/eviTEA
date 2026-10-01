import sys, os, csv, time
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from evitea.incumbent import mate_rescue as mr
from evitea.incumbent.msbwt_queries import MsbwtPairLocator
from evitea.incumbent.sequence import load_reference_contigs
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017")); R=A/"route_b_r9e_validation_full55x"/"shadow_v2"/"reports"
def rd(p): return list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"))
def ids(s): return tuple(int(x) for x in s.split(",") if x)
ledger=[mr.MateRescueLocus(r["locus_id"],r["chr"],r["family"],float(r["candidate_start"]),float(r["candidate_end"]),float(r["candidate_midpoint"]),ids(r["left_fragment_ids"]),ids(r["right_fragment_ids"]),r["pass_two_sided_rescue"]=="True") for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_LOCUS_LEDGER.tsv")]
targets=mr.one_sided_targets(ledger)
ht=rd(R/"V2_REFERENCE_DISCRIMINATIVE_RESCUE_R6"/"V2R6_TARGET_LOCI.tsv")
print("targets equal",[(t.locus_id,t.missing_side,",".join(map(str,t.existing_left_ids))) for t in targets]==[(r["locus_id"],r["missing_side"],r["existing_left_ids"]) for r in ht],len(targets))
reserved=set()
for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_COMBINED_CALLS.tsv"): reserved|=set(ids(r["fragment_ids"]))
sub=[t for i,t in enumerate(targets) if i%130==5][:20]
ref=load_reference_contigs(Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1/Reference/dm6_TEMP2_compatible.fa")),{t.chrom for t in sub})
loc=MsbwtPairLocator(A/"msbwt_SRR9951105")
t0=time.time(); res=mr.route_one_sided_targets(targets=sub,reserved_pair_ids=reserved,reference=ref,locator=loc); print("secs",round(time.time()-t0,1))
want={t.locus_id for t in sub}
mine=sorted((p,l,res.seed_count[(p,l)]) for p,ls in res.pair_to_loci.items() for l in ls)
hist=sorted((int(r["pair_id"]),r["locus_id"],int(r["supporting_host_seed_count"])) for r in rd(R/"V2_REFERENCE_DISCRIMINATIVE_RESCUE_R6"/"V2R6_PAIR_LOCUS_ROUTING.tsv") if r["locus_id"] in want)
print("routing equal",mine==hist,len(mine),len(hist))
