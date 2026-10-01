import sys, os, csv, json, time
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from collections import defaultdict
from evitea.incumbent import mate_rescue as mr
from evitea.incumbent.te_library import load_te_records, load_record_family_map, LocalTEIndex, load_family_kmer_panels
from evitea.incumbent.sequence import load_reference_contigs
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017")); R=A/"route_b_r9e_validation_full55x"/"shadow_v2"/"reports"; D=Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1"))
out=Path(os.path.expanduser("~/scratch/r6"))
def rd(p): return list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"))
def ids(s): return tuple(int(x) for x in s.split(",") if x)
ledger=[mr.MateRescueLocus(r["locus_id"],r["chr"],r["family"],float(r["candidate_start"]),float(r["candidate_end"]),float(r["candidate_midpoint"]),ids(r["left_fragment_ids"]),ids(r["right_fragment_ids"]),r["pass_two_sided_rescue"]=="True") for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_LOCUS_LEDGER.tsv")]
targets=mr.one_sided_targets(ledger)
want=set(out.joinpath("loci.txt").read_text().split())
p2l=defaultdict(set); sc={}
for r in rd(R/"V2_REFERENCE_DISCRIMINATIVE_RESCUE_R6"/"V2R6_PAIR_LOCUS_ROUTING.tsv"):
    if r["locus_id"] in want: p2l[int(r["pair_id"])].add(r["locus_id"]); sc[(int(r["pair_id"]),r["locus_id"])]=int(r["supporting_host_seed_count"])
routing=mr.RoutingResult(targets,{p:tuple(sorted(v)) for p,v in sorted(p2l.items())},sc,frozenset())
reads={int(k):tuple(v) for k,v in json.loads((out/"reads.json").read_text()).items()}
te=D/"TE consensus library"/"transposon_sequence_set.fa"; fmap=D/"fbte_to_flybase_family.tsv"
fams={t.family for t in targets if t.locus_id in want}
panels,_=load_family_kmer_panels(te,fmap,fams,21,1.10)
index=LocalTEIndex(load_te_records(te,load_record_family_map(fmap)))
ref=load_reference_contigs(D/"Reference"/"dm6_TEMP2_compatible.fa",{t.chrom for t in targets if t.locus_id in want})
t0=time.time()
res=mr.certify_strict_pairs(routing=routing,reads=reads,reference=ref,family_kmers=panels,te_index=index,bowtie2=os.path.expanduser("~/bowtie2-2.5.5-linux-x86_64/bowtie2"),host_index_prefix=os.path.expanduser("~/scratch/bt2/idx/dm6"),work_directory=out/"work")
print("secs",round(time.time()-t0,1),dict(res.counts))
mine=[(e.locus_id,str(e.pair_id),e.side,str(e.te_mate),str(e.host_mate),e.host_view,str(e.host_start),str(e.host_end),str(e.host_midpoint),str(e.host_edits),str(e.host_candidate_distance),str(e.te_seed_hits),str(e.te_seed_fraction),str(e.te_edits),str(-e.extra_strength[0])) for e in res.evidence]
hist=[(r["locus_id"],r["pair_id"],r["side"],r["te_mate"],r["host_mate"],r["host_view"],r["host_start"],r["host_end"],r["host_midpoint"],r["host_global_best_nm"],r["host_candidate_distance"],r["te_seed_hits"],r["te_seed_fraction"],r["te_edits"],r["recruited_by_seed_count"]) for r in rd(R/"V2_REFERENCE_DISCRIMINATIVE_RESCUE_R6"/"V2R6_CERTIFIED_NEW_FRAGMENT_EVIDENCE.tsv") if r["locus_id"] in want]
print("evidence equal",sorted(mine)==sorted(hist),len(mine),len(hist))
if sorted(mine)!=sorted(hist): print(set(mine)^set(hist))
