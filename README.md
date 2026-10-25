# EviTEA

EviTEA (Evidence-integrated Transposable Element Analysis) detects non-reference transposable-element (TE) insertions
in a genome from short-read, paired-end sequencing of that genome: places where the sample carries a TE copy that the
reference genome does not have. For each insertion it reports the position and the TE family.

EviTEA reads the sample through a multi-string BWT (msBWT) of all reads, so it can look up any sequence in every read. It
maps only the reads that matter, with Bowtie2. A call requires physically independent paired-end evidence, and each
junction supports one call only. Evidence that cannot be resolved goes to `unresolved.tsv` and does not become a call.

Author: Md Taksir Hasan Majumder. License: GPL-3.0-or-later (see [LICENSE](LICENSE)).

## Two kinds of calls

- **Library calls** (`evitea call`). With a TE library (consensus sequences of the species' TE families), EviTEA calls
  insertions of those families. A repeat-aware caller runs first: split-junction reads, mate-anchored search, and a
  repeat-aware step that retains a call whose differences from the TE record recur at the locus. Four extension layers
  follow: bilateral and asymmetric paired-end anchors, full-memory anchors, and internal TE words. Output: `calls.native.tsv`,
  `calls.tsv`, `calls.bed`, `calls.vcf`.
- **Library-free calls** (`evitea library-free`, Linux or WSL). EviTEA assembles the inserted sequence at each junction in
  the reads, without using the library, and groups the insertions into families by their sequence (families named
  `LF_cluster_<n>`). Calls within 500 bp of a library call are left out. The others are added to the library calls in
  `calls.combined.tsv` and `calls.combined.vcf`.

A species without a curated TE library can be run as well. The library-free families become TE records, and the library
layer runs again on them (`evitea run --no-library`; see [docs/REPEAT_ANNOTATION.md](docs/REPEAT_ANNOTATION.md)).

## Status

Version `0.1.0.dev0`. On two *Drosophila melanogaster* development samples (dm6, 125 bp reads), the pipeline reproduces
the call sets of the research pipeline byte for byte: 838 calls with `--profile AKA017` and 194 calls with `--profile ISO1`.
That is a **method-equivalence** result, not an accuracy benchmark. See
[validation/EVITEA_V1_VALIDATION_STATUS.md](validation/EVITEA_V1_VALIDATION_STATUS.md).

EviTEA takes these values from each sample and reference: the read length, the fragment-size 95th and 99th percentiles, the
primary contigs and the mitochondrial contig. TE records may have any names (see
[docs/CLI_AND_DATA.md](docs/CLI_AND_DATA.md)). EviTEA runs on other read lengths, fragment sizes, genomes and TE libraries,
and the synthetic tests cover them. **Its accuracy on species other than *D. melanogaster* has not been measured.**

## Requirements

- **CPython 3.14** (`>=3.14,<3.15`); other Python versions are not supported.
- Linux or Windows (x86-64) for the library layer (`prepare-*`, `call`). The library-free layer, the family-record
  commands and `--refine-breakpoints` require `mappy` and run on **Linux or WSL only**.
- [Bowtie2](https://github.com/BenLangmead/bowtie2) 2.5.x on `PATH` (an external program, not installed by pip).
- [BLAST+](https://blast.ncbi.nlm.nih.gov) (`makeblastdb`, `blastn`) on `PATH`, only for `build-repeat-map` (runs without a
  TE library).

## Installation

```bash
python -m pip install evitea                        # from a release; or, from a clone: python -m pip install .
python -m pip install "evitea[anchorpe]"            # + edlib: the fragment decoder of the paired-end channel
python -m pip install "evitea[library-free]"        # + mappy and parasail: library-free layer, family records (Linux/WSL)
python -m pip install "evitea[anchorpe,library-free]"
evitea check-install
```

`pip install evitea` installs its two required Python packages: [pymsbwt](https://pypi.org/project/pymsbwt/) (version 0.5 or
later, which builds the multi-string BWT of the reads) and NumPy. The extras add optional packages:

| extra | packages | needed for |
|---|---|---|
| `anchorpe` | edlib | the whole-read fragment decoder of `call` (or point `--edlib-python` at an interpreter that has edlib; edlib has no CPython 3.14 Windows wheel) |
| `library-free` | mappy, parasail | `library-free`, `build-family-records`, `--refine-breakpoints` (Linux or WSL) |

`evitea check-install` lists every program and module EviTEA uses (Python, pymsbwt, NumPy, Bowtie2, edlib, mappy, parasail,
BLAST+) with its version, or reports it as missing. It exits with code 8 when a required one is missing.

Bowtie2, BLAST+, mappy and parasail are easiest to get from conda. `environment.yml` (repository root) creates an
environment with CPython 3.14, Bowtie2, BLAST+ and EviTEA with both extras:

```bash
conda env create -f environment.yml     # or: micromamba create -f environment.yml
conda activate evitea && evitea check-install
```

`containers/Containerfile` builds a container image from the same file (`docker build -f containers/Containerfile -t evitea .`;
for Apptainer, `apptainer build evitea.sif docker-daemon://evitea:latest`).

To run the tests: `python -m pip install pytest`, then `python -m pytest`.

## Input files

| input | what it is |
|---|---|
| reference FASTA | the reference genome |
| repeat annotation | BED6 of the TE copies in the reference (contig, start, end, TE record, score, strand); [how to make it with RepeatMasker](docs/REPEAT_ANNOTATION.md) |
| TE library | FASTA of TE consensus sequences, one record per TE (the first word of the header is the record name) |
| family map | the family of every TE record: two columns `record<TAB>family` (any names), or a FlyBase table with `FBte_ID` and `FlyBase_Family` columns |
| reads | paired-end FASTQ, read 1 and read 2 in two files (plain or `.gz`), all reads of one length |

You choose every location; no command assumes a folder. Details of each input and of every option are in
[docs/CLI_AND_DATA.md](docs/CLI_AND_DATA.md).

## A worked example

[`examples/synthetic/run_example.sh`](examples/synthetic/run_example.sh) makes a small genome (contigs `chr1`..`chr5`, `MT`),
a 120-record TE library with ordinary names and 30x paired reads of a sample with 30 non-reference TE insertions, runs
`evitea run` on them and checks the calls against the expected `calls.native.tsv` (about 10 minutes on 4 CPUs). See
[examples/README.md](examples/README.md).

## Quick start

One command, from FASTQ to calls (each step writes to its own folder under `--output`; run the command again to continue
after an interruption):

```bash
evitea run --reference data/genome.fa --repeat-annotation data/repeat_annotation.bed \
  --te-library data/te_library.fa --family-map data/family_map.tsv \
  --read1 data/s_R1.fastq.gz --read2 data/s_R2.fastq.gz --sample-id s --output runs/s --threads 8
# add --library-free for the library-free calls (Linux/WSL); use --no-library instead of the library options without one
```

The same as separate commands (the prepared reference and TE library can be reused for every sample):

```bash
evitea prepare-reference --reference data/genome.fa --repeat-annotation data/repeat_annotation.bed --output prepared/reference
#   (or adopt an existing Bowtie2 index of the same FASTA: --bowtie2-index path/to/index/genome)
evitea prepare-te --te-library data/te_library.fa --family-map data/family_map.tsv --build-languages \
  --reference-assets prepared/reference --output prepared/te
evitea prepare-sample --read1 data/s_R1.fastq.gz --read2 data/s_R2.fastq.gz --sample-id s --output prepared/s
evitea call --sample-assets prepared/s --reference-assets prepared/reference --te-assets prepared/te --output runs/s/call
evitea library-free --run-dir runs/s/call --sample-assets prepared/s --reference-assets prepared/reference \
  --output runs/s/library_free                                         # optional; Linux or WSL
```

`prepare-sample` writes a JSON map with one entry per read. The map is impractically large for full-size libraries (tens
of millions of reads). For those, build an msBWT yourself and call directly from it:

```bash
evitea call --msbwt-index path/to/msbwt_index --sample-id s \
  --reference-assets prepared/reference --te-assets prepared/te --output runs/s/call
```

Useful options for `call`: `--threads`, `--workers`, `--work-dir`, `--resume`, `--bowtie2`, `--edlib-python`, and the values
that EviTEA otherwise takes from the sample and the reference: `--read-length`, `--fragment-q95`, `--fragment-q99`,
`--primary-contigs`. Run `evitea call --help` for the full list. To reproduce a development run exactly, pass
`--read-length 125 --fragment-q95 494 --fragment-q99 805 --primary-contigs chr2L,chr2R,chr3L,chr3R,chrX`.

## Runs without a TE library

`evitea run --no-library` builds a repeat map of the reference with BLAST+ (`build-repeat-map`) and runs the library-free
layer with it. It turns the library-free families into TE records (`build-family-records`), runs the library layer on those
records, and adds the calls of that run to the library-free calls (`add-record-calls`). The result is
`record_pass/added/calls.combined.tsv`. The run requires only the reference and the reads: without `--repeat-annotation`, the
repeat map is also the reference's annotation. The same chain is in `tools/no_library.sh` (plain shell) and
`tools/slurm/no_library.sh` (Slurm). See
[docs/REPEAT_ANNOTATION.md](docs/REPEAT_ANNOTATION.md).

## Outputs

| File | Content |
|---|---|
| `calls.native.tsv` | the call set in the native 15-column schema (the scientific result) |
| `calls.tsv` | one row per call, genomic order, with `call_origin` and supporting pair IDs |
| `calls.bed` | BED view of the calls |
| `calls.vcf` | VCF 4.2 view of the calls (symbolic `<INS:ME:family>` alleles, sites only) |
| `unresolved.tsv` | events that are not reported as calls |
| `run_summary.json` | counts per origin and hashes after every layer |
| `run.log`, `run_manifest.json` | provenance: the manifest hashes every other output and records the values taken from the sample (`values.run_parameters`) |

`library-free` adds `calls.library_free.tsv`, `calls.combined.tsv` and `calls.combined.vcf`. Completed output directories
are never overwritten.

## Profiles

`--profile production` is the default. The profiles `AKA017` and `ISO1` reproduce the research call sets of the two
development samples. Each accepts only the sample ID of its own sample.

## Reproducing the validation

The scripts under `validation/` and `tools/windows/` reproduce the byte-for-byte
comparison. They require the author's private development datasets and the output of the
research pipeline, which are **not distributed**. The scripts contain **no hard-coded paths: you must
provide every path yourself**, through arguments or environment variables
(`EVITEA_DATA`, `EVITEA_ROUTE_B_V2`, `EVITEA_ROUTE_B_V3`, `EVITEA_REPO`). See
[validation/README.md](validation/README.md).

## What is validated, and limitations

- Validated: byte-for-byte method equivalence with the research pipeline on two *D. melanogaster* samples (dm6, 125 bp).
  Synthetic end-to-end runs identify every planted insertion with its family: the automated test
  (`tests/integration/test_synthetic_call.py`, 125 bp, FlyBase and plain TE names, fly and non-fly contig names) and runs at
  100 and 150 bp.
- Not validated: accuracy (recall and precision against a truth set) on any species, including on other *D. melanogaster*
  strains.
- Reads must be paired-end and of one length; trim them to one length first.
- CPython 3.14 only. The library-free layer and the family-record commands run on Linux or WSL only.
- `prepare-sample`'s read-identity map does not scale to full-size libraries; use `--msbwt-index` for those.

## Citation

See [CITATION.cff](CITATION.cff). If you use EviTEA, please cite it and retain the
copyright and license notices in any copy or derivative, as the GPL requires.

## License

EviTEA builds on msBWT (James Holt and Leonard McMillan), Bowtie2, BLAST+, mappy (minimap2), parasail, edlib and NumPy. The
licences of these programs and packages are listed in [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES).

Copyright (C) 2026 Md Taksir Hasan Majumder. Licensed under the GNU General Public
License v3.0 or later. Any distributed derivative work must also be released under the
GPL with its source and must retain the original copyright and attribution.
