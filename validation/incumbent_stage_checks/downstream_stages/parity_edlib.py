import random, sys, time
sys.path.insert(0, "/sessions/"+__import__('os').environ.get('USER','')+"/work/EviTEA/src")
import os; sys.path.insert(0, os.path.expanduser("~/work/EviTEA/src"))
import edlib
from evitea.incumbent.alignment import hw_distance, hw_locations
rng = random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1)
def mutate(s, n, alpha):
    s=list(s)
    for _ in range(n):
        op=rng.randrange(3); i=rng.randrange(len(s)+1)
        if op==0 and i<len(s): s[i]=rng.choice(alpha)
        elif op==1: s.insert(i,rng.choice(alpha))
        elif op==2 and i<len(s): del s[i]
    return "".join(s)
bad=0; n=0; t0=time.time()
for trial in range(4000):
    alpha = rng.choice(["ACGT","AC","ACGTN","A"*3+"C"])
    tl = rng.choice([5,30,200,1601])
    target = "".join(rng.choice(alpha) for _ in range(tl))
    if rng.random()<0.5:  # tandem repeats
        unit = "".join(rng.choice(alpha) for _ in range(rng.randint(1,6)))
        pos=rng.randrange(tl); target = target[:pos] + unit*rng.randint(2,30) + target[pos:]
    ql = rng.choice([1,3,8,24,35,60,125])
    if rng.random()<0.7 and len(target)>ql:
        st=rng.randrange(len(target)-ql+1); q=mutate(target[st:st+ql], rng.randint(0,6), alpha)
    else:
        q="".join(rng.choice(alpha) for _ in range(ql))
    if not q: continue
    k = rng.choice([0,1,2,3,5,8,12])
    e = edlib.align(q, target, mode="HW", task="locations", k=k)
    ed = e["editDistance"]; exp = (None if ed<0 else ed, tuple(tuple(x) for x in e["locations"]) if ed>=0 else ())
    got = hw_locations(q, target, k)
    d = edlib.align(q, target, mode="HW", task="distance", k=k)["editDistance"]
    gd = hw_distance(q, target, k)
    n+=1
    if got!=exp or (None if d<0 else d)!=gd:
        bad+=1
        if bad<5: print("MISMATCH",repr(q),len(target),k,exp,got,d,gd)
print("trials",n,"mismatches",bad,"secs",round(time.time()-t0,1))
