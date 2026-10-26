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

* `--repeat-annotation` is required by the default engine (it flags anchor
  clusters inside TE copies).
* `--family-map` accepts a two-column `record<TAB>family` table with any record
  names (the first word of each TE FASTA header; an optional header row such as
  `te_record_id<TAB>family` is skipped), or a FlyBase-style table with `FBte_ID`
  and `FlyBase_Family` columns (records without a family stay unassigned). Every
  layer reads both forms (`evitea.family_table`).
* `--repeat-annotation` is BED6 named by TE record; see
  [`REPEAT_ANNOTATION.md`](REPEAT_ANNOTATION.md) for making it with RepeatMasker.
* The two codebooks are the terminal and internal 21-mer detection
  languages (`languages` section of the prepared TE assets).
  `--terminal-codebook` and `--internal-codebook` supply them as hash-recorded
  files (`source: supplied`); `--build-languages` builds them instead. A prepared
  TE asset without languages can be validated but not used by `call`.
* The optional `--g1a-archive` and `--family-crosswalk` pair is not used by the
  default engine.

Prepared manifests record source hashes, architecture identity, generated-file
hashes, and backend metadata. `validate reference`, `validate te`, and
`validate sample` recheck these records without reading benchmark truth.

## Calling

Exactly one sample input:

```powershell
evitea call --sample-assets prepared/sample `
  --reference-assets prepared/reference --te-assets prepared/te --output output/sample

evitea call --read1 sample_R1.fastq.gz --read2 sample_R2.fastq.gz --sample-id sample `
  --reference-assets prepared/reference --te-assets prepared/te --output output/sample

evitea call --msbwt-index path\to\msbwt_SRR... --sample-id sample `
  --reference-assets prepared/reference --te-assets prepared/te --output output/sample
```

`--msbwt-index` uses an existing msBWT index directory (with `about.npy`;
records `2*i` and `2*i+1` are the mates of pair `i`). The run manifest
records the SHA-256 of every index file.

Runtime options: `--threads` (Bowtie2), `--workers` (locus-batch
processes; default CPUs-1), `--work-dir`, `--resume`, `--bowtie2`,
`--bowtie2-build` (default: next to `--bowtie2`, else `PATH`),
`--edlib-python` (an interpreter with `edlib` for the paired-end anchor
fragment decoder, used only when the running interpreter lacks it; edlib has no
CPython 3.14 Windows wheel).

### Values taken from the sample and the reference

Each value below is written to `run_manifest.json` (`values.run_parameters`)
with where it came from.

| value | default | option |
|---|---|---|
| read length (the window of the coverage panels and of the wide-anchor step) | the prepared sample's read length, or the msBWT index's (`n` reads hold `n * (L + 1)` symbols) | `--read-length` |
| fragment-size 95th percentile (an ordinary reference pair longer than it does not explain a read pair) | estimated from the sample | `--fragment-q95` |
| fragment-size 99th percentile (mate-to-breakpoint windows; Bowtie2 `-X` of the paired reference mappings) | estimated from the sample | `--fragment-q99` |
| primary contigs (the full-memory layer looks only at these) | the reference rule below | `--primary-contigs` (comma list or file, one name per line) |

Fragment size: up to 200,000 evenly spaced read pairs are mapped to the
reference with Bowtie2 (end-to-end, no mixed or discordant pairs, `-X 5000`).
Spans of properly paired reads with both mates at MAPQ >= 20 are retained, spans
above median + 10 x MAD are dropped as outliers, and the nearest-rank 95th and
99th percentiles are taken. With fewer than 1,000 usable pairs the development
values 494 and 805 bp are used and a warning is logged.

Primary contigs, without `--primary-contigs`: every contig at least 10% as
long as the longest contig, except organelle contigs (`chrM`, `MT`, `M`,
`chrMT`, `Mito`, chloroplast names such as `ChrC` or `Pt`) and unplaced,
unlocalized or alternative contigs (`chrUn...`, `Un...`, or names with
`_random`, `_alt`, `_fix`, `_hap`, `_decoy`). The 10% scales with the genome:
it retains every yeast chromosome (chromosome I is about 15% of chromosome IV)
and every mouse or human chromosome, and leaves out unplaced scaffolds. On dm6
it gives chr2L, chr2R, chr3L, chr3R, chrX and chrY. The fly runs pass
`--primary-contigs chr2L,chr2R,chr3L,chr3R,chrX`, the list the layer was
developed with.

To reproduce a development run exactly, pass the development values:
`--read-length 125 --fragment-q95 494 --fragment-q99 805
--primary-contigs chr2L,chr2R,chr3L,chr3R,chrX`.

### Profiles

`--profile production` (default) names paired-end anchor host Bowtie2 queries
`P{pid}` and derives event-identifier labels from the sample id.
`--profile AKA017` and `--profile ISO1` are **reproduction-only**: they use the
host query names (`p{pid}|h{mate}` for AKA017) and identifier labels of the
research run, and accept only the matching `--sample-id`.
Bowtie2 seeds its per-read tie-breaking from the query name, so the profile is
part of the reproduction identity. No scientific threshold is configurable.

The equivalent raw-mode TOML template is
[`../examples/evitea.raw.example.toml`](../examples/evitea.raw.example.toml)
(`[inputs] profile`, `msbwt_index`; `[runtime] bowtie2_build`,
`edlib_python`, `workers`).

## Outputs

A completed run contains:

| File | Content |
|---|---|
| `calls.native.tsv` | the call set in the native 15-column schema, byte-exact: repeat-aware caller rows, then the rows appended by each extension layer |
| `calls.tsv` | one row per native call, genomic order: `call_origin` (INCUMBENT_SPLIT_JUNCTION / INCUMBENT_MATE_RESCUE / INCUMBENT_REPEAT_AWARE / EXTENSION_ANCHOR_PE_BILATERAL / EXTENSION_ANCHOR_PE_ASYMMETRIC / EXTENSION_FULL_MEMORY / EXTENSION_INTERNAL21), `coordinate_kind`, supports and support pair ids |
| `calls.bed` | BED view (point `p` becomes `[p, p+1)`, the last base for a contig-end point; range `[start, end)`) |
| `calls.vcf` | VCF 4.2, sites only: one record per call of `calls.tsv`, `ALT` `<INS:ME:family>`, `POS` the base before the call, `END`/`CIPOS`/`IMPRECISE` for ranges, `INFO` with `FAMILY`, `SOURCE` (library or library-free), `ORIGIN`, `EVIDENCE` and the supports (`evitea.vcf`); `library-free` writes `calls.combined.vcf` |
| `unresolved.tsv` | events not reported as calls: nearby call of another family, or competing extension events within 500 bp |
| `run_summary.json` | counts per origin, rows and SHA-256 after every layer, engine stage counts |
| `run.log`, `run_manifest.json` | provenance; the manifest hashes every other output |

Native coordinates are kept as the producers wrote them. `start == end` is a
breakpoint (zero-based interbase point; split-junction and all extension
events). `start < end` is a breakpoint range (mate-anchored search and
repeat-aware step events). An unknown native `evidence_type` stops the
run.

Completed output directories are never overwritten. An incomplete directory
containing only EviTEA work state requires explicit `--resume`. The engine's
stage checkpoints (`<work>/engine/work`) are reused only with `--resume` and
only when their recorded input identity matches.

## Known limitation

`prepare-sample` writes a JSON read-identity map with one entry per read.  For
full-size libraries (tens of millions of reads) this file and its verification
are impractically large; use `--msbwt-index` on an existing index for such
samples.
