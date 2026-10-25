"""Compare two tab-separated tables (or two FASTA files) and report whether they are the same.

usage: python table_diff.py A B [--ignore col,col,...] [--key col,col]

Columns named by --ignore are left out of the comparison, for example columns that exist in one table only (the recorded
one-site-check inputs). Without --key the rows are compared in file order; with --key they are matched by the key columns
first. Prints equal rows, differing rows, rows only in A, rows only in B, and the first differences.
Exit status 0 means the two are the same, 1 otherwise.
"""
import argparse
import csv
import sys


def read(path, ignore):
    if path.endswith((".fa", ".fasta")):
        recs, name = [], None
        for line in open(path):
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                recs.append({"name": name, "seq": ""})
            elif recs:
                recs[-1]["seq"] += line
        return ["name", "seq"], recs
    with open(path, encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = [c for c in rd.fieldnames if c not in ignore]
        return cols, [{c: r[c] for c in cols} for r in rd]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--ignore", default="")
    ap.add_argument("--key", default="")
    args = ap.parse_args()
    ignore = {c for c in args.ignore.split(",") if c}
    ca, ra = read(args.a, ignore)
    cb, rb = read(args.b, ignore)
    if ca != cb:
        print("columns differ:", [c for c in ca if c not in cb], [c for c in cb if c not in ca])
        return 1
    if args.key:
        keys = args.key.split(",")
        kf = lambda r: tuple(r[k] for k in keys)
        da, db = {kf(r): r for r in ra}, {kf(r): r for r in rb}
        only_a, only_b = sorted(set(da) - set(db)), sorted(set(db) - set(da))
        diff = [k for k in da if k in db and da[k] != db[k]]
        same = sum(1 for k in da if k in db and da[k] == db[k])
    else:
        n = min(len(ra), len(rb))
        diff = [i for i in range(n) if ra[i] != rb[i]]
        same = n - len(diff)
        only_a, only_b = list(range(n, len(ra))), list(range(n, len(rb)))
    print(f"equal {same} differ {len(diff)} only_in_A {len(only_a)} only_in_B {len(only_b)}")
    if diff:
        print("first differences (row number or key):", diff[:5])
    return 0 if not (diff or only_a or only_b) else 1


sys.exit(main())
