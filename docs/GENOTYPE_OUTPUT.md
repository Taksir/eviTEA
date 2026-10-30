# Genotype output (`evitea genotype`)

`evitea genotype` adds a table of hom / het genotypes to a finished library-free output. It reads the call files and writes two new
files; no call file is opened for writing and no call changes.

```
evitea genotype --library-free-dir OUT_LF --msbwt-index INDEX --reference-assets REF_ASSETS [--output DIR] [--counter plain|fast]
```

`--output` defaults to the library-free folder. The step refuses to overwrite `genotypes.tsv` or `genotype_summary.json`. Without an msBWT
index only the library-free calls are given (Step A); the library calls are left out.

## genotypes.tsv

One row per row of `calls.combined.tsv`, in its order.

| column | meaning |
|---|---|
| event_id, chr, start, end, family, source | as in `calls.combined.tsv` |
| ref_share | share of the read evidence at the break that reads straight through the reference (the empty allele); 0 = every read carries the insertion |
| n_ref, n_total | reads that continue the reference, and all reads, at the break (Step B only; `.` in Step A) |
| genotype | `hom` if ref_share <= 0.15, `het` if <= 0.85, `low` above, `.` if no count could be made |
| method | `lf_walk` (Step A) or `ref_walk` (Step B) |

The thresholds are `HOM_MAX` and `HET_MAX` of `library_free/_p4_lib.py`.

* **Step A, library-free calls.** `ref_share` is the reference-path fraction the detection recorded for the call's place: the smaller of
  `fracL` and `fracR` when the right walk broke, else `fracL` (`calling.e_value`). The table keeps these to three decimals; a call within
  0.0005 of a threshold can fall on the other side of it than the unrounded `geno` column of the places table.
* **Step B, library calls.** The library layer has insertion-side fragment counts only, so the reference-path walk of the library-free layer
  (`_p4_lib.walk`, `frac = nref / T`) is run on the msBWT at the call's break: a left walk from the first sampled word before `start`, a mirror
  walk from after `end`. `ref_share` is the smaller of the two fractions when both reach a break, else the one that does. The depth `D` is
  `_p4_lib.local_D` at the middle of the call. A call gets `.` when an anchor word is outside the contig, holds a letter other than A, C, G,
  T, or has more than twice the local depth, or when neither walk finds a divergence. No repeat map is read, so a call inside a repeat can have
  a share that reflects the other copies.

`genotype_summary.json` records the sha256 of the three input files, the genotype counts, the outcomes of the walks, the counter, the
thresholds, and the wall time of the step.
