"""depth_cap_changes.py on two synthetic detections (depth 150 with --depth-cap 100 and with no cap)."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import synth  # noqa: E402
import depth_cap_changes as T  # noqa: E402
from evitea import depth_window  # noqa: E402
from evitea.library_free import detector as D  # noqa: E402


def places(cap):
    depth_window.set_depth_cap(cap)
    try:
        ref, sites, C = synth.build(depth=150, glen=24000, sites=((3000, 1.0), (6000, 1.0), (9000, 1.0), (12000, 1.0), (15000, 1.0)))
        rm = D.reference_multi({"synth": ref})
        return D.places_tsv_bytes(D.scan_contig(C, ref, "synth", rm, synth.NoFirewall()))
    finally:
        depth_window.set_depth_cap(None)


def main():
    with tempfile.TemporaryDirectory() as d:
        for name, cap in (("old", 100), ("new", None)):
            open(os.path.join(d, name + ".tsv"), "wb").write(places(cap))
        rows = T.main(["--old", os.path.join(d, "old.tsv"), "--new", os.path.join(d, "new.tsv"), "--out", os.path.join(d, "ch.tsv")])
    assert len(rows) == 5 and all(r["old_stage"] == "S4_abstain" and r["new_stage"] == "S4_pass" and r["old_why"] == "depth" for r in rows), rows
    print("D9 depth_cap_changes.py lists the 5 places that go from S4_abstain (depth) to S4_pass when the cap is lifted")
    print("ALL PASS")


if __name__ == "__main__":
    main()
