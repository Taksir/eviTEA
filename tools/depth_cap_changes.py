"""List every place whose outcome differs between two library-free detections of the same sample (for example the same
reads run with --depth-cap 100 and with no cap), and every place whose local depth D is above a bound.

usage: depth_cap_changes.py --old PLACES_OLD.tsv --new PLACES_NEW.tsv [--above 100] [--out changes.tsv]
Both files are library_free_places.tsv.  Places are matched on (chrom, first, last).  Output columns: chrom, first, last, D, then
the old and the new stage, why and bL, and the kind of change (new_only, old_only, stage, why, bL).  Standard library only.
"""
import argparse
import csv
import sys


def read(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key(row):
    return row["chrom"], row["first"], row["last"]


def compare(old_rows, new_rows):
    old, new = {key(r): r for r in old_rows}, {key(r): r for r in new_rows}
    out = []
    for k in sorted(set(old) | set(new), key=lambda k: (k[0], int(k[1]), int(k[2]))):
        o, n = old.get(k), new.get(k)
        if o is None or n is None:
            kind = "new_only" if o is None else "old_only"
        else:
            kinds = [c for c in ("stage", "why", "bL") if (o.get(c) or "") != (n.get(c) or "")]
            if not kinds:
                continue
            kind = "+".join(kinds)
        ref = n or o
        out.append(dict(chrom=k[0], first=k[1], last=k[2], D=ref.get("D", ""),
                        old_stage=(o or {}).get("stage", ""), new_stage=(n or {}).get("stage", ""),
                        old_why=(o or {}).get("why", ""), new_why=(n or {}).get("why", ""),
                        old_bL=(o or {}).get("bL", ""), new_bL=(n or {}).get("bL", ""), change=kind))
    return out


def above(rows, bound):
    return [r for r in rows if r.get("D") not in ("", None) and float(r["D"]) > bound]


COLS = ["chrom", "first", "last", "D", "old_stage", "new_stage", "old_why", "new_why", "old_bL", "new_bL", "change"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--above", type=float, default=100.0)
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    old_rows, new_rows = read(a.old), read(a.new)
    changes = compare(old_rows, new_rows)
    handle = open(a.out, "w", newline="") if a.out else sys.stdout
    writer = csv.DictWriter(handle, COLS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(changes)
    if a.out:
        handle.close()
    print(f"places old {len(old_rows)}, new {len(new_rows)}; changed {len(changes)}; D above {a.above:g}: "
          f"old {len(above(old_rows, a.above))}, new {len(above(new_rows, a.above))}", file=sys.stderr)
    return changes


if __name__ == "__main__":
    main()
