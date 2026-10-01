# CLI and data contracts

## Preparation

All paths may be absolute or relative to the current directory. TOML-relative
paths are resolved against the directory containing the TOML file.

```powershell
evitea prepare-reference --reference dm6.fa --repeat-annotation dm6_repeatmasker_TE.bed6 `
  --bowtie2-build bowtie2-build --output prepared/reference
#   or adopt an existing Bowtie2 index of the same FASTA instead of building one:
#   --bowtie2-index path\to\bowtie2_index\dm6

evitea prepare-te --te-library transposon_sequence_set.fa `
  --family-map fbte_to_flybase_family.tsv `
  --terminal-codebook FINAL_SAFE_TERMINAL_CODEBOOK.tsv `
  --internal-codebook SAFE_INTERNAL21_INTERNAL_ONLY_CODEBOOK.tsv `
  --reference-assets prepared/reference --output prepared/te

evitea prepare-sample --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz `
  --sample-id sample --output prepared/sample
```

* `--repeat-annotation` is required by the default engine (RepeatMasker
  residency of AnchorPE clusters).
* `--family-map` accepts the two-column `te_record_id<TAB>family` form or a
  FlyBase-style table with `FBte_ID` and `FlyBase_Family` columns (the form the
  engines read; records without a family stay unassigned).
* The two codebooks are the v3 detection languages (prepared-TE schema v3,
  `languages` section).  Until their builders are EviTEA stages they are
  supplied as frozen, hash-recorded files (`source: supplied`).  A prepared
  TE asset without languages can be validated but cannot be used by `call`.
* The spec-era `--g1a-archive` / `--family-crosswalk` pair is optional and is
  not consumed by the default engine.

Prepared manifests bind source hashes, architecture identity, generated-file
hashes, and backend metadata. `validate reference`, `validate te`, and
`validate sample` recheck these bindings without accessing benchmark truth.

## Calling

Exactly one sample mode:

```powershell
evitea call --sample-assets prepared/sample `
  --reference-assets prepared/reference --te-assets prepared/te --output output/sample

evitea call --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz --sample-id sample `
  --reference-assets prepared/reference --te-assets prepared/te --output output/sample

evitea call --msbwt-index path\to\msbwt_SRR... --sample-id sample `
  --reference-assets prepared/reference --te-assets prepared/te --output output/sample
```

`--msbwt-index` uses an existing msBWT index directory (with `about.npy`,
record `2*i` / `2*i+1` = mates of pair `i`); the run manifest records the
SHA-256 of every index file.

Engineering options: `--threads` (Bowtie2), `--workers` (locus-batch
processes; default CPUs-1), `--work-dir`, `--resume`, `--bowtie2`,
`--bowtie2-build` (default: next to `--bowtie2`, else `PATH`),
`--edlib-python` (an interpreter with `edlib` for the AnchorPE fragment
decoder, used only when the running interpreter lacks it; edlib has no
CPython 3.14 Windows wheel).

### Profiles

`--profile production` (default) names AnchorPE host Bowtie2 queries
`P{pid}` and derives event-identifier labels from the sample id.
`--profile AKA017` and `--profile ISO1` are **reproduction-only**: they apply
the historical host query names (`p{pid}|h{mate}` for AKA017) and identifier
labels those samples were frozen with, and refuse any other `--sample-id`.
Bowtie2 seeds its per-read tie-breaking from the query name, so the profile is
part of the reproduction identity.  No scientific threshold is configurable.

The equivalent raw-mode TOML template is
[`../examples/evitea.raw.example.toml`](../examples/evitea.raw.example.toml)
(`[inputs] profile`, `msbwt_index`; `[runtime] bowtie2_build`,
`edlib_python`, `workers`).

## Outputs

A completed run contains:

| File | Content |
|---|---|
| `calls.native.tsv` | the callset in the frozen Route-B v3 schema (15 columns), byte-exact: R9E incumbent rows, then the rows appended by each extension layer |
| `calls.tsv` | one row per native call, genomic order: `call_origin` (INCUMBENT_SPLIT_JUNCTION / INCUMBENT_MATE_RESCUE / INCUMBENT_REPEAT_AWARE / EXTENSION_ANCHOR_PE_BILATERAL / EXTENSION_ANCHOR_PE_ASYMMETRIC / EXTENSION_FULL_MEMORY / EXTENSION_INTERNAL21), `coordinate_kind`, supports and support pair ids |
| `calls.bed` | BED view (point `p` -> `[p, p+1)`, the last base for a contig-end point; range `[start, end)`) |
| `unresolved.tsv` | fail-closed integration events (nearby different-family call, new-event multiplicity or cross-family competition) |
| `run_summary.json` | counts per origin, rows and SHA-256 after every layer, engine stage counts |
| `run.log`, `run_manifest.json` | provenance; the manifest hashes every other output |

Native coordinates are kept as the producers wrote them: `start == end` is a
breakpoint (zero-based interbase point; split-junction and all extension
events); `start < end` is a breakpoint range (mate-rescue and repeat-aware
rescue events).  An unknown native `evidence_type` fails the run closed.

Completed output directories are never overwritten. An incomplete directory
containing only EviTEA work state requires explicit `--resume`; the engine's
stage checkpoints (`<work>/engine/work`) are reused only with `--resume` and
only when their recorded input identity matches.

## Known limitation

`prepare-sample` writes a JSON read-identity map with one entry per read.  For
full-size libraries (tens of millions of reads) this file and its verification
are impractically large; use `--msbwt-index` on an existing index for such
samples until the identity map is made compact.
