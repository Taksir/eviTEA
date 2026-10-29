"""Compare two library-free scans: saved pieces by value, output folders file by file.

    python tools/compare_scan.py pieces  WORK_A WORK_B  [--contigs c1,c2,...]
    python tools/compare_scan.py outputs OUT_A  OUT_B
    python tools/compare_scan.py commit  SUMMARY.json  COMMIT_PREFIX  [--tree CODE_DIR]
    python tools/compare_scan.py treehash CODE_DIR

``pieces``: for each contig, the scan piece in ``WORK/pieces/<contig>.pkl`` of both folders (or the folder itself if it holds
``<contig>.pkl``), loaded and compared by value: every place and field, the info block, dict key order, list order; NaN equals NaN.
The identity recorded with a piece must also be equal. Two pieces with equal content can differ in bytes (pickle writes a string that
appears twice once and points to it afterwards), so the file hash is not compared.

``outputs``: every file under both folders (both views), compared by sha256; ``library_free_summary.json`` is compared after the
timings, the memory use, the ``speed`` entry, the paths of the run folders, the ``command_line`` and the ``source_identity`` (the git commit
and dirty flag of the code that ran) are taken out. The last two are printed as NOTE lines when they differ, so that what was left out is
visible; a NOTE is not a difference.

``commit``: which code ran. The hash of the detection code is the same for the final build and the fast-scan build, so the summary's
``source_identity`` is read: ``git_commit`` must start with COMMIT_PREFIX and ``git_dirty`` be false. A container without git records no
commit; then ``--tree CODE_DIR`` (the checkout the job used) is needed: the summary's ``source_tree_sha256`` (a hash of every .py file of the
package the run imported) must equal ``treehash`` of CODE_DIR. ``treehash CODE_DIR``: that hash of ``CODE_DIR/src/evitea``; run it on a checkout
of a known commit (with git) to learn the value its runs must show.

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
DROP = {"seconds", "rule_c_seconds", "peak_rss_gb", "speed", "scan_seconds", "source_run", "output", "command_line", "source_identity"}
NOTED = ("command_line", "source_identity")


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
            a_doc, b_doc = json.loads(fa[rel].read_text()), json.loads(fb[rel].read_text())
            for key in NOTED:
                if a_doc.get(key) != b_doc.get(key):
                    results.append((f"{rel}: {key} (expected to differ, left out of the comparison)", "NOTE", f"A: {a_doc.get(key)!r} | B: {b_doc.get(key)!r}"))
        else:
            ha, hb = (hashlib.sha256(fa[rel].read_bytes()).hexdigest(), hashlib.sha256(fb[rel].read_bytes()).hexdigest())
            results.append((rel, "EQUAL" if ha == hb else "DIFFERENT", None if ha == hb else f"sha256 {ha[:12]} vs {hb[:12]}"))
    return results


def tree_hash(code_dir: Path) -> str:
    """The hash ``evitea.cli._source_tree_sha256`` records: relative path and sha256 of every .py file under CODE_DIR/src/evitea, by path."""
    package = Path(code_dir) / "src" / "evitea"
    digest = hashlib.sha256()
    for path in sorted(p for p in package.rglob("*.py") if "__pycache__" not in p.parts):
        digest.update(path.relative_to(package).as_posix().encode() + b"\0" + hashlib.sha256(path.read_bytes()).hexdigest().encode() + b"\0")
    return digest.hexdigest()


def check_commit(summary: Path, prefix: str, tree: Path | None = None) -> int:
    identity = json.loads(summary.read_text()).get("source_identity", {})
    commit, dirty, ran = identity.get("git_commit"), identity.get("git_dirty"), identity.get("source_tree_sha256")
    if commit:
        ok = commit.startswith(prefix) and dirty is False
        print(f"{'OK' if ok else 'WRONG'} git_commit {commit} (wanted {prefix}...), git_dirty {dirty}")
        return 0 if ok else 1
    if tree is None:
        print(f"WRONG no git_commit in the summary (a container without git); give --tree CODE_DIR. source_tree_sha256 {ran}")
        return 1
    want = tree_hash(tree)
    ok = ran == want
    print(f"{'OK' if ok else 'WRONG'} no git commit recorded; source_tree_sha256 of the run {ran}, of {tree} {want}")
    return 0 if ok else 1


def main(argv: list[str]) -> int:
    if len(argv) == 2 and argv[0] == "treehash":
        print(tree_hash(Path(argv[1])))
        return 0
    if len(argv) in (3, 5) and argv[0] == "commit" and (len(argv) == 3 or argv[3] == "--tree"):
        return check_commit(Path(argv[1]), argv[2], Path(argv[4]) if len(argv) == 5 else None)
    if len(argv) < 3 or argv[0] not in ("pieces", "outputs"):
        print(__doc__)
        return 2
    contigs = None
    if "--contigs" in argv:
        contigs = [c for c in argv[argv.index("--contigs") + 1].split(",") if c]
    a, b = Path(argv[1]), Path(argv[2])
    results = compare_pieces(a, b, contigs) if argv[0] == "pieces" else compare_outputs(a, b)
    notes = [r for r in results if r[1] == "NOTE"]
    results = [r for r in results if r[1] != "NOTE"]
    for label, state, detail in results:
        print(f"{state:9s} {label}")
    bad = [(label, detail) for label, state, detail in results if state != "EQUAL"]
    for label, detail in bad:
        print(f"  first difference, {label}: {detail}")
    for label, _, detail in notes:
        print(f"NOTE      {label}\n            {detail}")
    print(f"{len(results) - len(bad)} of {len(results)} equal")
    return 1 if bad or not results else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
