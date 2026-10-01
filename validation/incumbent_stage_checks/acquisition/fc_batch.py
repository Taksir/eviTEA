from stages import *
import sys
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from evitea.incumbent.te_library import *
from evitea.incumbent.sequence import load_reference_contigs
from evitea.incumbent import parallel
def main():
    layer=sys.argv[1]; T0=time.time(); BUDGET=float(sys.argv[2]) if len(sys.argv)>2 else 150
    loci=pickle.load(open(WORK/"loci.pkl","rb"))
    if layer=="mate_rescue":
        loci=pickle.load(open(WORK/"eligible.pkl","rb"))
    out=WORK/"batches"/layer; out.mkdir(parents=True,exist_ok=True)
    batches=[loci[i:i+25] for i in range(0,len(loci),25)]
    todo=[i for i in range(len(batches)) if not (out/f"{i:05d}.pkl").exists()]
    print("pending",len(todo),"of",len(batches))
    if todo:
        te=ASSET/"TE consensus library"/"transposon_sequence_set.fa"
        panels,_=load_family_kmer_panels(te,FAMMAP,{L.family for L in pickle.load(open(WORK/"loci.pkl","rb"))},21,1.10)
        tei=LocalTEIndex(load_te_records(te,load_record_family_map(FAMMAP)))
        ref=load_reference_contigs(ASSET/"Reference"/"dm6_TEMP2_compatible.fa",{L.chrom for L in loci})
        reads=pickle.load(open(WORK/"loci_reads.pkl","rb"))
        state={"reads":reads,"reference":ref,"family_kmers":panels,"te_index":tei}
        print("setup",time.time()-T0)
        with ProcessPoolExecutor(2, initializer=parallel._initialize, initargs=(state,)) as pool:
            running={}; it=iter(todo)
            def submit():
                i=next(it,None)
                if i is not None and time.time()-T0<BUDGET: running[pool.submit(parallel._run_batch,layer,batches[i])]=i
            submit(); submit()
            while running:
                done,_=wait(running,return_when=FIRST_COMPLETED)
                for f in done:
                    i=running.pop(f); pickle.dump(f.result(),open(out/f"{i:05d}.pkl","wb")); submit()
    print("remaining", sum(1 for i in range(len(batches)) if not (out/f"{i:05d}.pkl").exists()), time.time()-T0)
    
if __name__=='__main__':
    main()