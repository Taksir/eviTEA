"""Pairs and proposals rejected by the side tests of a finished `call` run, from its saved checkpoints.

    PYTHONPATH=src python tools/mate_side_counts.py <call work dir>      (the folder that holds engine/work/incumbent/checkpoints,
                                                                          or the checkpoints folder itself)

The same numbers are read under either rule (--mate-side-rule, --junction-side-rule), so two runs can be set side by side. With an
opt-in rule the run's manifest also holds them (incumbent counts, `side_test_rejects_*`).
"""

import json
import pickle
import sys
from pathlib import Path

from evitea.incumbent.engine import side_test_counts


def checkpoints(root: Path) -> Path:
    for candidate in (root, root / "checkpoints", root / "incumbent" / "checkpoints", root / "engine" / "work" / "incumbent" / "checkpoints"):
        if (candidate / "mate_rescue.pkl").is_file():
            return candidate
    raise SystemExit(f"no mate_rescue.pkl under {root}")


def main(argv):
    folder = checkpoints(Path(argv[0]))
    load = lambda name: pickle.loads((folder / f"{name}.pkl").read_bytes())    # noqa: E731
    mates, strict, groups = load("mate_rescue"), load("strict_pairs"), load("junction_groups")
    out = side_test_counts(mates, strict, groups)
    out["strand_rule"] = {k: v for k, v in sorted(mates.counts.items()) if k.startswith("STRAND_RULE_")}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
