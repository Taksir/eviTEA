# EviTEA

EviTEA (Evidence-integrated Transposable Element Analysis) calls transposable-element
(TE) insertions from short-read, paired-end sequencing data. It works over a
multi-string BWT (msBWT) of the reads, runs a repeat-aware incumbent caller, and then
adds four evidence-extension layers: bilateral and asymmetric AnchorPE, full-memory
AnchorPE, and INTERNAL21.

The design rules are: a junction is assigned to exactly one owner, a call needs
physically independent paired-end corroboration, and evidence that cannot be resolved is
reported as unresolved instead of being forced into a call.

Author: Md Taksir Hasan Majumder. License: GPL-3.0-or-later (see [LICENSE](LICENSE)).

## Status

Version `0.1.0.dev0`. On two *Drosophila melanogaster* development samples (AKA017 and
ISO1) the pipeline reproduces the frozen research callsets byte for byte: 838 calls for
AKA017 and 194 for ISO1. That is a **method-equivalence** result, not an accuracy
benchmark. It has been validated only on dm6 with 125 bp paired-end reads. See
[validation/EVITEA_V1_VALIDATION_STATUS.md](validation/EVITEA_V1_VALIDATION_STATUS.md).

## Requirements

- **CPython 3.14** (`>=3.14,<3.15`); other Python versions are not supported.
- Linux or Windows (x86-64) for the library layer (`prepare-*`, `call`). The library-free layer, the family-record
  commands and `--refine-breakpoints` need `mappy` and run on **Linux or WSL only**.
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

`pip install evitea` pulls in its two required Python packages, [pymsbwt](https://pypi.org/project/pymsbwt/) (version 0.5 or
later; it builds the multi-string BWT of the reads) and NumPy. The extras add optional packages:

| extra | packages | needed for |
|---|---|---|
| `anchorpe` | edlib | the whole-read fragment decoder of `call` (or point `--edlib-python` at an interpreter that has edlib; edlib has no CPython 3.14 Windows wheel) |
| `library-free` | mappy, parasail | `library-free`, `build-family-records`, `--refine-breakpoints` (Linux or WSL) |

`evitea check-install` lists every program and module EviTEA uses (Python, pymsbwt, NumPy, Bowtie2, edlib, mappy, parasail,
BLAST+) with its version, or says that it is missing; it exits with code 8 when a required one is missing.

Bowtie2, BLAST+, mappy and parasail are easiest to get from conda. `environment.yml` (repository root) creates an
environment with CPython 3.14, Bowtie2, BLAST+ and EviTEA with both extras:

```bash
conda env create -f environment.yml     # or: micromamba create -f environment.yml
conda activate evitea && evitea check-install
```

`containers/Containerfile` builds a container image from the same file (`docker build -f containers/Containerfile -t evitea .`;
for Apptainer, `apptainer build evitea.sif docker-daemon://evitea:latest`).

To run the tests: `python -m pip install pytest`, then `python -m pytest`.

## Input files, and where to put them

**You choose all file locations yourself.** Nothing in this repository assumes a
particular folder; every path is an argument you pass. A layout that works well:

```text
project/
  data/reference/   host genome FASTA (e.g. dm6) + RepeatMasker TE annotation (BED)
  data/te/          TE consensus library FASTA + family table (FBte_ID, FlyBase_Family)
  data/reads/       paired FASTQ (R1/R2), uniform read length
  prepared/         outputs of prepare-reference / prepare-te / prepare-sample (reusable)
  results/<sample>/ one `evitea call` output directory per sample
```

Details of each input are in [docs/CLI_AND_DATA.md](docs/CLI_AND_DATA.md). In short: a
host FASTA, a RepeatMasker TE annotation, a TE library FASTA, a family map
(two columns `te_record_id<TAB>family`, or a FlyBase-style table with `FBte_ID` and
`FlyBase_Family`), the two v3 language codebooks (or let EviTEA build them), and
synchronized paired FASTQ.

## Quick start

```bash
# 1. Reference (reusable across samples)
evitea prepare-reference --reference data/reference/dm6.fa \
  --repeat-annotation data/reference/dm6_repeatmasker_TE.bed6 \
  --bowtie2-build bowtie2-build --output prepared/reference
#    or adopt an existing Bowtie2 index of the same FASTA: --bowtie2-index path/to/index/dm6

# 2. TE library (reusable)
evitea prepare-te --te-library data/te/transposon_sequence_set.fa \
  --family-map data/te/fbte_to_flybase_family.tsv --build-languages \
  --reference-assets prepared/reference --output prepared/te

# 3. Sample, then call
evitea prepare-sample --read1 data/reads/s_R1.fastq.gz --read2 data/reads/s_R2.fastq.gz \
  --sample-id s --output prepared/s
evitea call --sample-assets prepared/s \
  --reference-assets prepared/reference --te-assets prepared/te --output results/s
```

`prepare-sample` writes a JSON map with one entry per read, which becomes impractically
large for full-size libraries (tens of millions of reads). For those, build an msBWT
yourself and call directly from it:

```bash
evitea call --msbwt-index path/to/msbwt_index --sample-id s \
  --reference-assets prepared/reference --te-assets prepared/te --output results/s
```

Useful options for `call`: `--threads`, `--workers`, `--work-dir`, `--resume`,
`--bowtie2`, `--edlib-python`. Run `evitea call --help` for the full list.

## Outputs

| File | Content |
|---|---|
| `calls.native.tsv` | the callset in the frozen 15-column schema (the scientific result) |
| `calls.tsv` | one row per call, genomic order, with `call_origin` and supporting pair IDs |
| `calls.bed` | BED view of the calls |
| `unresolved.tsv` | events deliberately not turned into calls |
| `run_summary.json` | counts per origin and hashes after every layer |
| `run.log`, `run_manifest.json` | provenance; the manifest hashes every other output |

Completed output directories are never overwritten.

## Profiles

`--profile production` is the default. `AKA017` and `ISO1` exist only to reproduce the
frozen development callsets and refuse any other sample ID.

## Reproducing the validation

The scripts under `validation/` and `tools/windows/` reproduce the byte-for-byte
comparison. They need the author's private development datasets and frozen research
outputs, which are **not distributed**. They also contain **no hard-coded paths: you must
supply every path yourself**, through arguments or environment variables
(`EVITEA_DATA`, `EVITEA_ROUTE_B_V2`, `EVITEA_ROUTE_B_V3`, `EVITEA_REPO`). See
[validation/README.md](validation/README.md).

## Limitations

- Host-genome specifics: dm6 and FlyBase-style TE identifiers (`FBte...`), a fixed list
  of primary chromosomes, and an insert-size constant tuned for the validation data.
- Validated only on 125 bp paired-end reads.
- CPython 3.14 only.
- Method equivalence with the research pipeline is shown; accuracy against truth sets is not claimed here.

## Citation

See [CITATION.cff](CITATION.cff). If you use EviTEA, please cite it and keep the
copyright and license notices in any copy or derivative, as the GPL requires.

## License

Copyright (C) 2026 Md Taksir Hasan Majumder. Licensed under the GNU General Public
License v3.0 or later. Any distributed derivative work must also be released under the
GPL with its source and must retain the original copyright and attribution.
