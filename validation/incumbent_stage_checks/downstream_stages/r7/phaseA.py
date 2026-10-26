import sys, os, csv, json, time
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from collections import defaultdict
from evitea.incumbent import mate_rescue as mr, junction_groups as jg
from evitea.incumbent.msbwt_queries import MsbwtPairLocator
from evitea.incumbent.te_library import load_family_kmer_panels
from evitea.incumbent.sequence import load_reference_contigs
A=Path(os.path.expandvars("$EVITEA_DATA/AKA017")); R=A/"route_b_r9e_validation_full55x"/"shadow_v2"/"reports"; D=Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1"))
out=Path(os.path.expanduser("~/scratch/r7"))
def rd(p): return list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"))
def ids(s): return tuple(int(x) for x in s.split(",") if x)
ledger=[mr.MateRescueLocus(r["locus_id"],r["chr"],r["family"],float(r["candidate_start"]),float(r["candidate_end"]),float(r["candidate_midpoint"]),ids(r["left_fragment_ids"]),ids(r["right_fragment_ids"]),r["pass_two_sided_rescue"]=="True") for r in rd(R/"V2_ORIENTATION_CORRECTED_MATE_RESCUE_R4"/"V2R4_LOCUS_LEDGER.tsv")]
want=set(out.joinpath("loci.txt").read_text().split())
targets=[t for t in mr.one_sided_targets(ledger) if t.locus_id in want]
routed=defaultdict(list)
for r in rd(R/"V2_REFERENCE_DISCRIMINATIVE_RESCUE_R6"/"V2R6_PAIR_LOCUS_ROUTING.tsv"):
    if r["locus_id"] in want: routed[r["locus_id"]].append(int(r["pair_id"]))
reads={int(k):tuple(v) for k,v in json.loads((out/"reads.json").read_text()).items()}
te=D/"TE consensus library"/"transposon_sequence_set.fa"; fmap=D/"fbte_to_flybase_family.tsv"
panels,_=load_family_kmer_panels(te,fmap,{t.family for t in targets},21,1.10)
ref=load_reference_contigs(D/"Reference"/"dm6_TEMP2_compatible.fa",{t.chrom for t in targets})
counts={}
t0=time.time(); disc=jg.discover_proposals(targets,routed,reads,ref,panels,counts); print("discovery",len(disc),round(time.time()-t0,1))
hist=[r for r in rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_DISCOVERY_JUNCTIONS.tsv") if r["locus_id"] in want]
key=lambda p:(p.locus_id,str(p.source_pair_id),str(p.source_mate),p.read_view,p.junction_semantic,str(p.split),str(p.te_bp),str(p.te_seed_hits),str(p.te_seed_fraction),str(p.host_seed_matches_in_scope),str(p.host_seed_first_pos),p.signature,p.host_half,str(p.proposal_score))
hk=lambda r:(r["locus_id"],r["source_pair_id"],r["source_mate"],r["read_view"],r["junction_semantic"],r["split"],r["te_bp"],r["te_seed_hits"],r["te_seed_fraction"],r["host_seed_matches_in_scope"],r["host_seed_first_pos"],r["signature"],r["host_half"],r["proposal_score"])
print("discovery equal",[key(p) for p in disc]==[hk(r) for r in hist],len(hist))
# recorded qids match Bowtie2 read names
qmap={hk(r):r["qid"] for r in hist}
disc=[jg.JunctionProposal(qmap[key(p)],*[getattr(p,f) for f in jg.JunctionProposal.__slots__[1:]]) for p in disc]
acc,st=jg.screen_proposals(disc,bowtie2=os.path.expanduser("~/bowtie2-2.5.5-linux-x86_64/bowtie2"),host_index_prefix=os.path.expanduser("~/scratch/bt2/idx/dm6"),work_directory=out/"work",counts=counts)
print("status equal",[st[p.qid] for p in disc]==[r["prequery_status"] for r in hist])
loc=MsbwtPairLocator(A/"msbwt_SRR9951105")
t0=time.time(); located=jg.locate_accepted(acc,loc); print("located",len(located),round(time.time()-t0,1))
aud={r["qid"]:r for r in rd(R/"V2_MIXED_EVIDENCE_COMPLETION_R7"/"V2R7_EXACT_QUERY_AUDIT.tsv")}
ok=all(aud[p.qid]["located_pair_ids"]==",".join(map(str,pids)) and int(aud[p.qid]["raw_msbwt_rows"])==rows for p,_,_,rows,pids in located)
print("locate equal",ok, len(located)==sum(1 for p,_,_ in acc if aud[p.qid]["locate_complete"]=="True"))
json.dump(sorted({x for l in located for x in l[4]}),open(out/"hits.json","w"))
