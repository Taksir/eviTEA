import sys, os, json, time
sys.path.insert(0,"src"); sys.path.insert(0,".")
from pathlib import Path
from validation.legacy_adapters.fastq_random_access import SraFastqIndex
out=Path(sys.argv[1]); pairs=sorted(set(json.loads((out/"pairs.json").read_text())))
cache=out/"reads.json"; reads=json.loads(cache.read_text()) if cache.exists() else {}
F=Path(os.path.expandvars("$EVITEA_DATA/AKA017/fastq"))
a,b=SraFastqIndex(F/"SRR9951105_1.fastq"),SraFastqIndex(F/"SRR9951105_2.fastq")
t=time.time()
for p in pairs:
    if str(p) in reads: continue
    reads[str(p)]=[a.sequence(p),b.sequence(p)]
    if time.time()-t>160: break
cache.write_text(json.dumps(reads)); print(len(reads),"/",len(pairs))
