"""Build the reference word tables of the library-free detection once, into a cache folder that later jobs load (--reference-cache).

    python tools/build_reference_cache.py --reference-fasta REF.fa --cache-dir DIR [--no-walk-out-words]

REF.fa must be the file the library-free runs read (the prepared reference's FASTA: the key holds its sha256) and the detection code must be
the one the runs use (the key holds its hash). Needs memory for the whole reference: about 30 bytes per base for the multi-copy words (mm10:
use 96 GB). Run it in a job, not on a login node. Prints the cache folder. A second call finds the folder and loads instead of building.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--no-walk-out-words", action="store_true", help="only the multi-copy words (a run with --walk-out off)")
    args = parser.parse_args()
    from evitea.breakpoints.runtime import read_fasta
    from evitea.library_free import refcache
    from evitea.library_free.runtime import _code_sha256, _sha256

    started = time.time()
    reference_sha, code_sha = _sha256(args.reference_fasta), _code_sha256()
    print(f"reference sha256 {reference_sha}\ndetection code hash {code_sha[:8]}", flush=True)
    genome = read_fasta(args.reference_fasta)
    refcache.load_or_build(genome, reference_sha256=reference_sha, code_sha256=code_sha, cache_dir=args.cache_dir,
                           need_words=not args.no_walk_out_words, progress=lambda m: print(m, flush=True))
    folder = args.cache_dir / refcache.cache_key(reference_sha, code_sha)
    print(f"{folder}: {sum(p.stat().st_size for p in folder.iterdir()) / 1e9:.2f} GB, {time.time() - started:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
