from stages import *
from evitea.incumbent.te_library import *
from evitea.incumbent.sequence import load_reference_contigs
from evitea.incumbent import split_junction as sj
loci=pickle.load(open(WORK/"loci.pkl","rb"))
import random; random.seed(3); sub=random.sample(list(loci),60)
idx=CachedIndex()
t=time.time(); reads=idx.pair_reads(p for L in sub for p in L.pair_ids); print("reads",len(reads),time.time()-t)
te=ASSET/"TE consensus library"/"transposon_sequence_set.fa"
t=time.time(); panels,_=load_family_kmer_panels(te,FAMMAP,{L.family for L in loci},21,1.10); tei=LocalTEIndex(load_te_records(te,load_record_family_map(FAMMAP))); print("setup",time.time()-t)
t=time.time(); ref=load_reference_contigs(ASSET/"Reference"/"dm6_TEMP2_compatible.fa",{L.chrom for L in loci}); print("ref",time.time()-t)
t=time.time(); r=sj.run_split_junction_layer(sample_id="AKA017",loci=sub,reads=reads,reference=ref,family_kmers=panels,te_index=tei); dt=time.time()-t
print("fcrc1 60 loci",dt,"pairs",sum(len(L.pair_ids) for L in sub),"est total h", dt/sum(len(L.pair_ids) for L in sub)*sum(len(L.pair_ids) for L in loci)/3600)
