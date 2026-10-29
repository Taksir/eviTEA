"""Compare two library-free scans: saved pieces by value, output folders file by file.

    python tools/compare_scan.py pieces  WORK_A WORK_B  [--contigs c1,c2,...]
    python tools/compare_scan.py outputs OUT_A  OUT_B

``pieces``: for each contig, the scan piece in ``WORK/pieces/<contig>.pkl`` of both folders (or the folder itself if it holds
``<contig>.pkl``), loaded and compared by value: every place and field, the info block, dict key order, list order; NaN equals NaN.
The identity recorded with a piece must also be equal. Two pieces with equal content can differ in bytes (pickle writes a string that
appears twice once and points to it afterwards), so the file hash is not compared.

``outputs``: every file under both folders (both views), compared by sha256; ``library_free_summary.json`` is compared after the
timings, the memory use, the ``speed`` entry and the paths of the run folders are taken out.

Exit 0 when everything is equal, 1 otherwise. Prints EQUAL or DIFFERENT per item, then the first difference of each different item.
"""
from __future__ import annotations

import hashlib
import json
import math
import pickle
import sys
from pathlib import Path

# the timings, the memory use, the speed entry and the path of the run folder differ between two runs of one result
DROP = {"seconds", "rule_c_seconds", "peak_rss_gb", "speed", "scan_seconds", "source_run", "output"}


def first_difference(a, b, path="piece"):
    """The first difference of two values (None when equal): type, dict keys in order, list length, value; NaN equals NaN."""
    if type(a) is not type(b):
        return f"type at {path}: {type(a).__name__} vs {type(b).__name__}"
    if isinstance(a, dict):
        if list(a) != list(b):
            return f"keys at {path}: {list(a)[:15]} vs {list(b)[:15]}"
        for key in a:
            found = first_difference(a[key], b[key], f"{path}.{key}")
            if found:
                return found
        return None
    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return f"length at {path}: {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            found = first_difference(x, y, f"{path}[{i}]")
            if found:
                return found
        return None
    if isinstance(a, float) and math.isnan(a) and math.isnan(b):
        return None
    try:
        same = bool(a == b)
    except ValueError:  # an array: compare element by element
        import numpy as np

        same = bool(np.array_equal(a, b))
    return None if same else f"value at {path}: {a!r} vs {b!r}"


def _pieces_dir(folder: Path) -> Path:
    return folder / "pieces" if (folder / "pieces").is_dir() else folder


def compare_pieces(folder_a: Path, folder_b: Path, contigs=None) -> list[tuple[str, str, str | None]]:
    a_dir, b_dir = _pieces_dir(folder_a), _pieces_dir(folder_b)
    names = contigs or sorted({p.name[: -len(".pkl")] for p in a_dir.glob("*.pkl")} | {p.name[: -len(".pkl")] for p in b_dir.glob("*.pkl")})
    results = []
    for name in names:
        label = f"piece {name}"
        paths = [d / f"{name}.pkl" for d in (a_dir, b_dir)]
        if not all(p.is_file() for p in paths):
            results.append((label, "DIFFERENT", f"missing: {[str(p) for p in paths if not p.is_file()]}"))
            continue
        pieces = []
        for p in paths:
            with p.open("rb") as handle:
                pieces.append(pickle.load(handle))
        found = first_difference(*pieces)
        ids = []
        for d in (a_dir, b_dir):
            marker = d / f"{name}.complete.json"
            ids.append(json.loads(marker.read_text(encoding="utf-8")).get("identity") if marker.is_file() else None)
        if found is None and ids[0] != ids[1]:
            found = "identity of the piece differs: " + str({k: (ids[0] or {}).get(k) for k in (ids[0] or {}) if (ids[0] or {}).get(k) != (ids[1] or {}).get(k)})
        results.append((label, "EQUAL" if found is None else "DIFFERENT", found))
    return results


def _scrub(x):
    if isinstance(x, dict):
        return {k: _scrub(v) for k, v in x.items() if k not in DROP}
    if isinstance(x, list):
        return [_scrub(v) for v in x]
    return x


def compare_outputs(out_a: Path, out_b: Path) -> list[tuple[str, str, str | None]]:
    def files(root):
        return {str(p.relative_to(root)): p for p in sorted(root.rglob("*")) if p.is_file()}

    fa, fb = files(out_a), files(out_b)
    results = []
    for rel in sorted(set(fa) | set(fb)):
        if rel not in fa or rel not in fb:
            results.append((rel, "DIFFERENT", "only in " + ("B" if rel not in fa else "A")))
        elif Path(rel).name == "library_free_summary.json":
            found = first_difference(_scrub(json.loads(fa[rel].read_text())), _scrub(json.loads(fb[rel].read_text())), "summary")
            results.append((rel + " (without timings and speed)", "EQUAL" if found is None else "DIFFERENT", found))
        else:
            ha, hb = (hashlib.sha256(fa[rel].read_bytes()).hexdigest(), hashlib.sha256(fb[rel].read_bytes()).hexdigest())
            results.append((rel, "EQUAL" if ha == hb else "DIFFERENT", None if ha == hb else f"sha256 {ha[:12]} vs {hb[:12]}"))
    return results


def main(argv: list[str]) -> int:
    if len(argv) < 3 or argv[0] not in ("pieces", "outputs"):
        print(__doc__)
        return 2
    contigs = None
    if "--contigs" in argv:
        contigs = [c for c in argv[argv.index("--contigs") + 1].split(",") if c]
    a, b = Path(argv[1]), Path(argv[2])
    results = compare_pieces(a, b, contigs) if argv[0] == "pieces" else compare_outputs(a, b)
    for label, state, detail in results:
        print(f"{state:9s} {label}")
    bad = [(label, detail) for label, state, detail in results if state != "EQUAL"]
    for label, detail in bad:
        print(f"  first difference, {label}: {detail}")
    print(f"{len(results) - len(bad)} of {len(results)} equal")
    return 1 if bad or not results else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
