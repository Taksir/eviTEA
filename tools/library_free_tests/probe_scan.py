"""Print JSON about a synthetic library-free scan: python probe_scan.py DEPTH [CAP].  Compare two source trees."""
import hashlib
import json
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import synth  # noqa: E402

depth = int(sys.argv[1])
cap = sys.argv[2] if len(sys.argv) > 2 else None
if cap is not None:
    from evitea import depth_window
    depth_window.set_depth_cap(None if cap == "none" else cap)

from evitea.library_free import _p4_lib as P  # noqa: E402
from evitea.library_free import detector as D  # noqa: E402

ref, sites, C = synth.build(depth=depth, sites=((3000, 1.0), (6000, 1.0), (9000, 1.0), (12000, 1.0), (15000, 1.0)), glen=24000)
rm = D.reference_multi({"synth": ref})
places = D.scan_contig(C, ref, "synth", rm, synth.NoFirewall())
grouped = D.group_and_call(places)
from collections import Counter  # noqa: E402

out = dict(places_sha=hashlib.sha256(D.places_tsv_bytes(places)).hexdigest(), segments_sha=hashlib.sha256(D.segments_fa_bytes(grouped)).hexdigest(),
           stages=dict(Counter(p["stage"] for p in places)), why=dict(Counter(p.get("why") for p in places if p["stage"] == "S4_abstain")),
           d_median=sorted(round(p["D"]) for p in places)[len(places) // 2] if places else None, d_max=P.D_MAX)
out["d_max"] = None if out["d_max"] == float("inf") else out["d_max"]
print(json.dumps(out, sort_keys=True))
