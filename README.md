<p align="center">
  <img src="assets/banner.svg" alt="EviTEA: evidence-integrated transposable element analysis" width="100%">
</p>

<p align="center">
  <a href="https://github.com/Taksir/eviTEA/releases/tag/paper-2026-11-15"><img alt="Paper release" src="https://img.shields.io/badge/paper%20release-paper--2026--11--15-ff9f1c?style=for-the-badge"></a>
  <img alt="Python 3.14" src="https://img.shields.io/badge/python-3.14-3776ab?style=for-the-badge&logo=python&logoColor=white">
  <a href="LICENSE"><img alt="License GPL-3.0-or-later" src="https://img.shields.io/badge/license-GPL--3.0--or--later-2dd4bf?style=for-the-badge"></a>
  <img alt="Linux, Windows and WSL" src="https://img.shields.io/badge/platform-Linux%20%7C%20Windows%20%7C%20WSL-8b5cf6?style=for-the-badge">
</p>

<p align="center">
  <a href="#why-evitea">Why</a> &middot;
  <a href="#how-it-works">How it works</a> &middot;
  <a href="#install">Install</a> &middot;
  <a href="#requirements">Requirements</a> &middot;
  <a href="#quick-start">Quick start</a> &middot;
  <a href="#tutorial">Tutorial</a> &middot;
  <a href="#outputs">Outputs</a> &middot;
  <a href="#settings-used-in-the-paper">Paper settings</a> &middot;
  <a href="#reproducing-the-paper">Reproduce the paper</a> &middot;
  <a href="#citation">Citation</a>
</p>

EviTEA (evidence-integrated transposable element analysis) calls non-reference transposable element (TE) insertions from the paired-end short reads of one genome. For each insertion it reports the position and the TE family. It finds insertions of the families in your TE library, and it also finds insertions of families that the library does not contain.

The method is described in *EviTEA: calling transposable element insertions from short reads when the TE library is incomplete* (Majumder, manuscript in preparation). The code used for the paper is release [`paper-2026-11-15`](https://github.com/Taksir/eviTEA/releases/tag/paper-2026-11-15) (commit `d79d205`).

<a id="why-evitea"></a>

## 🧬 Why EviTEA

A TE library holds only the families that were known when someone built it. Library-based callers cannot name an insertion of any other family. In some genomes they are a large share of all insertions. The reads still hold the evidence for them. EviTEA uses that evidence, so one run gives one call set for both kinds of family.

<a id="how-it-works"></a>

## ⚙️ How it works

<p align="center">
  <img src="assets/how-it-works.svg" alt="Reads go into one msBWT index. The library layer and the library-free layer both search it, and their calls join in one table." width="900">
</p>

EviTEA puts all reads of the sample into one index, a multi-string Burrows-Wheeler transform (msBWT; Holt and McMillan 2014). The index gives the exact count of any k-mer in the reads and the reads that contain it. EviTEA aligns with Bowtie2 only the reads it needs.

The library layer (`evitea call`) searches the index for k-mers of your TE library. It uses split reads and read pairs with one mate in unique genome sequence and the other mate in a TE. When the evidence shows a TE joined to the genome, the layer makes a call and names the family.

The library-free layer (`evitea library-free`) does not use the library. It finds sites where the reads stop following the reference, and it builds the inserted sequence from the reads on both sides. It then groups similar inserted sequences into families named `LF_cluster_<n>`. Such a family counts as a TE family when at least two independent insertions have junctions on both sides and a target-site duplication of 2 to 20 bp. Inserted sequence that the local reference, the mitochondrial genome, a segmental duplication or a tandem repeat explains is discarded. Library-free calls within 500 bp of a library call are left out, and the rest go into `calls.combined.tsv` with the library calls.

Every call needs support from independent read pairs, and each junction supports one call only. Evidence that EviTEA cannot resolve goes to `unresolved.tsv`.

<a id="install"></a>

## 📦 Install

You need CPython 3.14 (`>=3.14,<3.15`).

> [!NOTE]
> EviTEA is not on PyPI, so install it from a clone.

```bash
git clone https://github.com/Taksir/eviTEA.git
cd eviTEA
git checkout paper-2026-11-15          # the paper's code (skip this line for the newest code)
python -m pip install ".[anchorpe,library-free]"
evitea check-install
```

`pip` installs pymsbwt (0.5 or later, which builds the msBWT) and NumPy. The two extras are optional:

| Extra | Adds | Needed for |
|---|---|---|
| `anchorpe` | edlib | the paired-end step of `call`. On Windows, edlib has no CPython 3.14 wheel, so point `--edlib-python` at another interpreter that has it. |
| `library-free` | mappy, parasail | `library-free`, `build-family-records` and `--refine-breakpoints` (Linux or WSL only) |

Two programs must be on your `PATH`. pip does not install them:

- [Bowtie2](https://github.com/BenLangmead/bowtie2) 2.5.x, for every run.
- [BLAST+](https://blast.ncbi.nlm.nih.gov) (`makeblastdb`, `blastn`), only for runs without a TE library.

`evitea check-install` prints the version of every program and module that EviTEA uses. It exits with code 8 when a required one is missing.

<details>
<summary><b>Conda and container</b></summary>

With conda, `environment.yml` installs Python 3.14, Bowtie2 and BLAST+, and then installs EviTEA with both extras through pip:

```bash
conda env create -f environment.yml
conda activate evitea
evitea check-install
```

`containers/Containerfile` builds an image from the same file: `docker build -f containers/Containerfile -t evitea .` For Apptainer, run `apptainer build evitea.sif docker-daemon://evitea:latest` after that.

</details>

<a id="requirements"></a>

## 🖥️ System requirements

Use Linux (x86-64) or WSL to run everything. On native Windows only the library layer runs. The library-free layer needs Linux or WSL, and the paired-end step needs a second Python with edlib (see the table above).

On Linux, CPython 3.14 has no wheel yet for edlib and mappy, so `pip` builds them from source. You need a C compiler and `make`. The conda file installs both.

**What the paper's fly runs needed.** The fly samples are four test strains (GIM-024, MUN-009, SLA-001 and JUT-011) at 54x to 59x, with 125 bp or 151 bp read pairs. The numbers come from the supplement (Tables S2H and S9). Each run used 8 CPUs on AMD EPYC 9654 nodes. The paper did not record memory for each layer on its own. It recorded the peak memory of the whole run, so the table says "not recorded" where that is the case.

| Part of the run | Time | Peak memory | Source |
|---|---|---|---|
| Read index (once per sample; reused by every later step) | 28 to 35 min (43 min for JUT-011, on another cluster) | 34.1 GiB (MUN-009), 40.7 GiB (SLA-001), 39 GB (JUT-011). Not recorded for GIM-024 | Table S9B |
| Library layer, `call` step | 64 to 87 min | not recorded on its own | Table S9B |
| Library layer, `refine-breakpoints` | 2 to 3 min | not recorded on its own | Table S9B |
| Library-free layer, 8 workers | 29 to 37 min | not recorded on its own | Table S9B |
| Both layers, one after the other (the whole run) | 95 to 119 min, without the read index | 42.0 to 47.2 GiB | Table S9B |
| Both layers, three Kuruka fly lines | 1 h 42 min to 1 h 46 min | 44.8 to 52.4 GiB | Table S2H |

The memory figures of the full-depth runs come from the Slurm accounting record (MaxRSS). They include file pages that the job cached, so the memory that EviTEA itself holds may be lower. If you plan a fly run, ask for 8 CPUs and 64 GiB. This is our advice, not a measurement. The highest peak recorded was 52.4 GiB.

Smaller fly runs need less. Whole runs on subsampled reads of 10x to 30x took 65 to 87 min and peaked at 22.3 to 36.5 GiB (Table S9A, EviTEA only, without the read index). Preparing those samples (drawing the reads and building the index, 4 CPUs) took 9 to 19 min and peaked at 7.6 to 25.0 GiB.

**A small run with every layer measured.** The [tutorial](#tutorial) below runs each layer on its own, so it gives a memory figure per layer. We measured it on a Linux machine with 4 CPUs and 16 GB of RAM, using 1,000,000 yeast read pairs (about 25x). The peak is the largest memory of the step and its child processes, sampled while it ran.

| Step of the tutorial | Time | Peak memory |
|---|---|---|
| Read index | 41 s | 0.9 GiB |
| Reference and library preparation | about 20 s | 0.2 GiB |
| Library layer (`call` and `refine-breakpoints`) | about 3 min | 0.7 GiB |
| Library-free layer (after the library layer) | about 10 min | 2.4 GiB |
| Both layers, one after the other | about 14 min | 2.4 GiB (the largest step) |

Memory depends on the genome size and the depth (see the subsampled fly runs above), so these figures do not carry over to flies. In this run the library-free layer needed the most memory.

The library-free layer alone was recorded for the mouse. Each of the 20 jobs (one per chromosome, 8 CPUs) took 1 h 29 min to 4 h 40 min and peaked at 24.7 to 48.0 GiB, with a limit of 48 GiB per job (Table S2H).

**Disk.** The paper did not record disk use. The read index holds about one byte per sequenced base. In the yeast tutorial below, 1,000,000 read pairs (302 million bases) gave an index of 314 MB. At that ratio, a 55x fly sample (about 8 billion bases) would need an index of roughly 8 to 9 GB. This is an estimate from the tutorial, not a measurement of a fly sample.

**Versions of the paper.** The paper's runs used Python 3.14.4, Bowtie2 2.5.5, mappy 2.31, edlib 1.3.9.post1, parasail 1.3.4 and msbwt-modern3 0.3.0 ([github.com/Taksir/msbwt](https://github.com/Taksir/msbwt) commit `19eea78`), inside a container image that is not part of this repository. `pip install` takes the newest compatible release of each Python package instead, and `environment.yml` takes any Bowtie2 2.5.x. The synthetic example gives the expected output with Bowtie2 2.5.2 and pymsbwt 0.5.0.

<a id="quick-start"></a>

## 🚀 Quick start

```bash
evitea run \
  --reference genome.fa --repeat-annotation repeat_annotation.bed \
  --te-library te_library.fa --family-map family_map.tsv \
  --read1 s_R1.fastq.gz --read2 s_R2.fastq.gz \
  --sample-id s --output runs/s --threads 8
```

Add `--library-free` to run the library-free layer too (Linux or WSL). If you have no TE library, give `--no-library` in place of `--te-library` and `--family-map` (see [Runs without a TE library](#runs-without-a-te-library)).

Each step writes to its own folder under `--output`. If a run stops, run the same command again and it continues from the last finished step. The final calls are in `runs/s/call/`, or in `runs/s/library_free/calls.combined.tsv` with `--library-free`.

> [!TIP]
> To try EviTEA first on a small data set, run `bash examples/synthetic/run_example.sh`. It makes a small genome, a library of 120 TE records and about 30x of 125 bp read pairs with 30 planted insertions. Then it runs `evitea run` and compares the calls with the expected ones. It takes 6 to 10 minutes on 4 CPUs. The last line it prints is `calls.native.tsv: as expected`, and the run has 30 calls. The reads have no sequencing errors, so this checks that the programs work together, not how accurate they are. See [examples/README.md](examples/README.md).

<a id="tutorial"></a>

## 🧭 Tutorial: from a download to calls

This tutorial runs EviTEA on a small real data set: the first 1,000,000 read pairs of the yeast strain SK1 (run SRR4074258, 151 bp, about 25x). It takes the same route as the paper's runs. You build the read index yourself, and you call each layer with the options of the paper.

You need Linux or WSL, EviTEA with both extras ([Install](#install)), BLAST+ and `curl`. The run takes about 19 minutes on 4 CPUs (about 29 minutes with the optional no-library step), about 1.4 GB of disk, and the largest step needs about 2.4 GiB of memory. We measured this on a Linux machine with 4 CPUs and 16 GB of RAM. Everything happens in one new folder.

<details>
<summary><b>1. Download the reference, the TE library and the reads</b></summary>

```bash
mkdir evitea_tutorial && cd evitea_tutorial

# reference genome of S. cerevisiae (R64, called sacCer3 by UCSC)
curl -sSL -o sacCer3.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gzip -d sacCer3.fa.gz

# TE library: the six Ty sequences of the McClintock yeast set
curl -sSL -o ty_library.fa https://raw.githubusercontent.com/bergmanlab/mcclintock/master/test/sac_cer_TE_seqs.fasta

# reads: the first 1,000,000 pairs of SRR4074258 (4 lines per read, so 4,000,000 lines per file)
B=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR407/008/SRR4074258
for i in 1 2; do
  curl -sS $B/SRR4074258_$i.fastq.gz 2>/dev/null | gzip -dc 2>/dev/null | head -n 4000000 > SK1_$i.fastq || true
done
wc -l SK1_1.fastq SK1_2.fastq        # 4000000 lines each
```

`curl` stops with an error message when `head` has enough lines. That is expected, and the `|| true` hides it. This step took about 4.5 minutes on our connection. If a file has fewer than 4,000,000 lines, the download broke off early. Run the loop again.

</details>

<details>
<summary><b>2. Make the repeat annotation and the family map</b></summary>

EviTEA needs the TE copies of the reference as BED6 and the family of each library record ([Inputs](#inputs)). Here BLAST+ finds the copies (80% identity or more, 100 bp or longer). The paper used RepeatMasker for yeast, so this annotation is close to the paper's but not the same.

```bash
makeblastdb -in sacCer3.fa -dbtype nucl -out sacCer3 > /dev/null
blastn -query ty_library.fa -db sacCer3 -perc_identity 80 -evalue 1e-5 -dust no -soft_masking false \
  -outfmt "6 sseqid sstart send qseqid bitscore length" \
  | awk 'BEGIN{OFS="\t"} $6>=100 {s=$2; e=$3; st="+"; if (s>e) {t=s; s=e; e=t; st="-"} print $1, s-1, e, $4, int($5), st}' \
  | sort -k1,1 -k2,2n > repeat_annotation.bed

# TY3_1p is the long terminal repeat of Ty3, so it belongs to the family TY3
printf 'TY1\tTY1\nTY2\tTY2\nTY3\tTY3\nTY3_1p\tTY3\nTY4\tTY4\nTY5\tTY5\n' > family_map.tsv
wc -l repeat_annotation.bed          # 662 lines with BLAST+ 2.12.0
```

</details>

<details>
<summary><b>3. Build the read index</b></summary>

Every paper run reads one msBWT index of the sample. The index holds both FASTQ files. The mates of pair *i* become records 2*i* and 2*i*+1, which is the order that `--msbwt-index` requires. Save this script as `build_index.py`:

```python
"""Build the msBWT of a paired sample for `evitea call --msbwt-index`.
usage: python build_index.py R1.fastq R2.fastq OUTPUT_FOLDER THREADS"""
import logging
import sys
from MUSCython import MultiStringBWTCython

def main(read1, read2, out_dir, threads):
    MultiStringBWTCython.createMSBWTFromFastq(
        [read1, read2], out_dir, int(threads), True, logging.getLogger("msbwt"))
    MultiStringBWTCython.loadBWT(out_dir)   # also writes the count and FM-index files that EviTEA reads

if __name__ == "__main__":      # needed: the build starts worker processes
    main(*sys.argv[1:5])
```

Then run it:

```bash
mkdir index
python build_index.py SK1_1.fastq SK1_2.fastq index 4
ls index                             # about.npy  fmIndex.npy  msbwt.npy  totalCounts.npy
```

Run it from a script file, not with `python -` or `python -c`, and make the output folder first. Both mistakes stop the build. The step took 41 seconds and 0.9 GiB of memory, and the index is 314 MB. Full-size samples take longer and need more memory (see [System requirements](#requirements)). Reads must all have the same length (151 bp here), and the pairs must be in the same order in both files.

</details>

<details>
<summary><b>4. Prepare the reference and the library</b></summary>

```bash
evitea prepare-reference --reference sacCer3.fa --repeat-annotation repeat_annotation.bed \
  --output prepared/reference --threads 4
evitea prepare-te --te-library ty_library.fa --family-map family_map.tsv --build-languages \
  --reference-assets prepared/reference --output prepared/te
```

These two steps depend on the reference and the library only. For more samples of the same species, run them once and reuse the folders. They took about 20 seconds.

</details>

<details>
<summary><b>5. Call with the library layer, then refine the breakpoints</b></summary>

```bash
evitea call --msbwt-index index --sample-id SK1 \
  --reference-assets prepared/reference --te-assets prepared/te \
  --output run/call --threads 4 --mate-side-rule strand
evitea refine-breakpoints --run-dir run/call --msbwt-index index \
  --reference-assets prepared/reference --te-assets prepared/te \
  --pin-one-sided --call-annotations --output run/refine
```

`call` took about 3 minutes and `refine-breakpoints` took about 5 seconds. `call` made 50 calls (`run/call/calls.native.tsv`) and kept 48 in `run/call/calls.tsv` after the own-family copy filter.

</details>

<details>
<summary><b>6. Add the library-free layer</b></summary>

```bash
evitea library-free --run-dir run/call --msbwt-index index \
  --reference-assets prepared/reference --te-assets prepared/te \
  --contigs chrI,chrII,chrIII,chrIV,chrV,chrVI,chrVII,chrVIII,chrIX,chrX,chrXI,chrXII,chrXIII,chrXIV,chrXV,chrXVI \
  --workers 4 --weak-site-pass --weak-site-allele-depth-cut 10 --extension-stop-min 2 \
  --output run/library_free
```

This took about 10 minutes and peaked at 2.4 GiB. The contig list leaves out the mitochondrial genome (`chrM`), as in the paper. The final table is `run/library_free/calls.combined.tsv`: the library calls and the library-free calls together. It has 85 calls here, 37 of them from the library-free layer.

</details>

<details>
<summary><b>7. Optional: the library-free layer alone (no-library mode)</b></summary>

This is how the paper made its no-library results. The reference is prepared without any annotation, and there is no TE library. `library-free` needs a run folder that holds an empty set of library calls. The second command makes one:

```bash
evitea prepare-reference --reference sacCer3.fa --output prepared/reference_noannot --threads 4
python -c 'import sys; from pathlib import Path; from evitea.run_command import write_stub_run; write_stub_run(Path(sys.argv[1]), sys.argv[2])' \
  run_nolib/empty_run SK1
evitea library-free --run-dir run_nolib/empty_run --msbwt-index index \
  --reference-assets prepared/reference_noannot \
  --contigs chrI,chrII,chrIII,chrIV,chrV,chrVI,chrVII,chrVIII,chrIX,chrX,chrXI,chrXII,chrXIII,chrXIV,chrXV,chrXVI \
  --workers 4 --extension-stop-min 2 --resident-guard off --output run_nolib/library_free
```

The last step took about 10 minutes and peaked at 2.4 GiB. It made 67 calls in `run_nolib/library_free/calls.combined.tsv`.

</details>

<details>
<summary><b>8. Look at the result</b></summary>

```bash
head -3 run/library_free/calls.combined.tsv        # event_id, chr, start, end, family, source
tail -n +2 run/library_free/calls.combined.tsv | cut -f5 | sort | uniq -c | sort -k1,1nr | head
```

A call is one line with its contig, position and family. Library calls carry the family of your library (`TY1` to `TY5`). Library-free families are named `LF_cluster_<n>`. `calls.vcf` and `calls.bed` hold the library calls in standard formats ([Outputs](#outputs)). With the same inputs you should get the same counts as above. Small differences are possible if your BLAST+ version gives a different annotation.

This data set shows that the install and the commands work on real reads. It has no ground truth, so it does not measure accuracy. The ground-truth sets and the scoring scripts of the paper are described under [Reproducing the paper](#reproducing-the-paper).

</details>

<a id="inputs"></a>

## 📥 Inputs

| Input | What it is |
|---|---|
| Reference | FASTA of the reference genome |
| Repeat annotation | BED6 of the TE copies in the reference (contig, start, end, TE record, score, strand). See [docs/REPEAT_ANNOTATION.md](docs/REPEAT_ANNOTATION.md) for how to make it with RepeatMasker. |
| TE library | FASTA of TE consensus sequences. The first word of each header is the record name, and any names work. |
| Family map | the family of each TE record: two columns, `record<TAB>family`, or a FlyBase table with `FBte_ID` and `FlyBase_Family` columns |
| Reads | paired-end FASTQ in two files, plain or `.gz`. All reads must have the same length, so trim them first if they do not. |

EviTEA measures the read length and the 95th and 99th percentiles of the fragment size from the sample. It also finds the primary contigs and the mitochondrial contig of the reference. You can set each of these yourself with `--read-length`, `--fragment-q95`, `--fragment-q99` and `--primary-contigs`. [docs/CLI_AND_DATA.md](docs/CLI_AND_DATA.md) describes every input and option.

<a id="step-by-step"></a>

## 🪜 Step by step

<details>
<summary><b>The commands that <code>evitea run</code> calls</b></summary>

`evitea run` calls these commands in order. Run them yourself when you want to reuse the prepared reference and TE library for many samples:

```bash
evitea prepare-reference --reference genome.fa --repeat-annotation repeat_annotation.bed --output prepared/reference
evitea prepare-te --te-library te_library.fa --family-map family_map.tsv --build-languages \
  --reference-assets prepared/reference --output prepared/te
evitea prepare-sample --read1 s_R1.fastq.gz --read2 s_R2.fastq.gz --sample-id s --output prepared/s
evitea call --sample-assets prepared/s --reference-assets prepared/reference --te-assets prepared/te \
  --output runs/s/call
evitea library-free --run-dir runs/s/call --sample-assets prepared/s --reference-assets prepared/reference \
  --te-assets prepared/te --output runs/s/library_free
```

`prepare-reference` can adopt an existing Bowtie2 index of the same FASTA with `--bowtie2-index path/to/index/genome`.

`prepare-sample` keeps a map with one entry per read, and that map gets too large for full-size samples of tens of millions of reads. For those, build an msBWT of the reads yourself (step 3 of the [tutorial](#tutorial) gives a short script) and give it to `call`, `refine-breakpoints` and `library-free` with `--msbwt-index`. The mates of pair *i* must be records 2*i* and 2*i*+1 of the index.

```bash
evitea call --msbwt-index path/to/msbwt_index --sample-id s \
  --reference-assets prepared/reference --te-assets prepared/te --output runs/s/call
```

Run `evitea <command> --help` for all options.

</details>

<a id="outputs"></a>

## 📤 Outputs

| File | Content |
|---|---|
| `calls.tsv` | one line per call in genomic order, with the supporting read pair IDs |
| `calls.bed`, `calls.vcf` | the same calls as BED and as VCF 4.2 (symbolic `<INS:ME:family>` alleles, sites only) |
| `calls.native.tsv` | all library calls in a 15-column table, before the own-family copy filter. `library-free` and `call-tiers` read it. |
| `resident_copy_dropped.tsv` | the calls that the own-family copy filter removed, each with the reference copy that caused it |
| `unresolved.tsv` | evidence that EviTEA does not report as a call |
| `run_summary.json` | call counts per step |
| `run.log`, `run_manifest.json` | the log, and hashes of every output with the values measured from the sample |

`library-free` adds `calls.library_free.tsv`, `calls.combined.tsv` and `calls.combined.vcf`. EviTEA never overwrites a finished output folder.

The own-family copy filter removes a library call that touches or overlaps a reference copy of its own family. The reference already has that family at that spot. The copies come from your repeat annotation. The filter changes the final tables only. `calls.native.tsv` stays complete, and library-free calls are never filtered. Give `--keep-resident-copies` to turn the filter off. When you run `library-free` on a run with library calls, give either `--te-assets` or `--keep-resident-copies`.

<a id="runs-without-a-te-library"></a>

## 🧪 Runs without a TE library

`evitea run --no-library` needs only the reference and the reads:

1. `build-repeat-map` maps repeats in the reference with BLAST+.
2. The library-free layer calls insertions and groups them into families.
3. `build-family-records` turns those families into TE records.
4. The library layer runs again on those records, and `add-record-calls` adds its calls.

The final calls are in `record_pass/added/calls.combined.tsv`. The same chain is in `tools/no_library.sh` (plain shell) and `tools/slurm/no_library.sh` (Slurm). See [docs/REPEAT_ANNOTATION.md](docs/REPEAT_ANNOTATION.md).

> [!IMPORTANT]
> The no-library results in the paper come from the library-free layer alone, not from this full chain. The command is under [Settings used in the paper](#settings-used-in-the-paper).

<a id="settings-used-in-the-paper"></a>

## 📄 Settings used in the paper

All results of the paper come from commit `d79d205` (release `paper-2026-11-15`). Only the false-positive comparison of the long-read steps (Table S3E(B) of the supplement) also gives results of the earlier commit `296e725`. The supplement names its tables as follows (the numbers follow the submitted version and may change before publication):

| Supplement table | What it holds |
|---|---|
| Table S1 | every parameter of EviTEA, with its value and what it was set on |
| Table S2 | software versions (S2B), the commands (S2C), and the threads, contigs and work limit of each run (S2D) |
| Table S8, part C | recall by allele frequency |
| Table S9 | run time and peak memory |

The runs set some options that are off or different by default:

| Option | Default | In the paper |
|---|---|---|
| `--mate-side-rule` (`call`) | `midpoint` | `strand`, in all library-mode runs |
| `--pin-one-sided --call-annotations` (`refine-breakpoints`) | off | on, in all library-mode runs |
| `--weak-site-pass --weak-site-allele-depth-cut` (`library-free`) | off | on, with the value 10, in all library-mode runs |
| `--extension-stop-min` (`library-free`) | 3 | 2, in all library-mode and no-library runs |
| `--resident-guard` (`library-free`) | on | off in no-library mode |
| `--max-rescue-interval-mass` (`call`) | 1,500,000 | 2,000,000 for JUT-011 and ISO1, and 100,000,000 for the mouse (see the table below) |

The commands, with the values of the paper for the fly. `<index>` is the msBWT index of the sample ([how to build it](#tutorial)). `<ref>` and `<te>` are the folders of `prepare-reference` and `prepare-te`. `<run>` is the output folder of `call`. The paper's command lines are shortened in Table S2C, and `refine-breakpoints` also needs `--output`.

```bash
# library mode: three steps per sample
evitea call --msbwt-index <index> --sample-id <sample> --reference-assets <ref> --te-assets <te> \
    --output <run> --threads 8 --mate-side-rule strand            # add --max-rescue-interval-mass 2000000 for JUT-011 and ISO1
evitea refine-breakpoints --run-dir <run> --msbwt-index <index> --reference-assets <ref> --te-assets <te> \
    --pin-one-sided --call-annotations --output <refine>
evitea library-free --run-dir <run> --msbwt-index <index> --reference-assets <ref> --te-assets <te> \
    --contigs chr2L,chr2R,chr3L,chr3R,chrX --workers 8 \
    --weak-site-pass --weak-site-allele-depth-cut 10 --extension-stop-min 2 --output <library-free>

# no-library mode: the library-free layer alone, with no TE library and no annotation
evitea library-free --run-dir <empty run> --msbwt-index <index> --reference-assets <ref-without-annotation> \
    --contigs chr2L,chr2R,chr3L,chr3R,chrX --workers 4 --extension-stop-min 2 --resident-guard off \
    --output <library-free>
```

`<empty run>` is a folder with an empty library call set, like the `nolib_stub/` folder that `evitea run --no-library` writes. Step 7 of the [tutorial](#tutorial) shows how to make one. `<ref-without-annotation>` comes from `prepare-reference` without `--repeat-annotation`.

The other values of each run (Table S2D):

| Runs | `call` threads | `library-free` workers | Contigs of `library-free` | `--max-rescue-interval-mass` |
|---|---|---|---|---|
| Fly, four test strains at full depth and the subsampled samples | 8 | 8 | chr2L, chr2R, chr3L, chr3R, chrX | 2,000,000 for JUT-011, default for GIM-024, MUN-009 and SLA-001 |
| Fly, ISO1 and the synthetic heterozygous samples; families removed, JUT-011 | 8 | 8 | as above | 2,000,000 |
| Fly, families removed, AKA-017 | 8 | 8 | as above | default |
| Fly, no-library runs | none | 4 | as above | not used |
| Yeast, five strains | 8 | 8 | chrI to chrXVI | default |
| Yeast, no-library runs | none | 8 | chrI to chrXVI | not used |
| *Z. tritici*, three isolates | 4 | 2 | chromosomes 1 to 21 | default |
| *A. thaliana* Cvi-0 | 8 | 2 | Chr1 to Chr5 | default |
| Mouse CC055, library mode | 12 (11 workers) | 8, one job per chromosome | chr1 to chr19 and chrX (20 jobs) | 100,000,000 |
| Mouse CC055, no-library mode | none | 8, one job per chromosome | chr1 to chr19 and chrX (20 jobs) | not used |

At the default work limit of 1,500,000, the call step stopped on the JUT-011 reads (work 1,522,148). That run was repeated at 2,000,000. The ISO1 run and the synthetic heterozygous samples used 2,000,000 from the start, so that all runs with ISO1 reads share one value. The default also stopped the mouse call step, and the mouse run used 100,000,000. Breakpoint refinement stopped with an error on the synthetic heterozygous samples of JUT-011 and ISO1 at 50% and 40x and at 50% and 54x. For these two runs, `library-free` started from the calls of the `call` step. The evaluation of the paper uses the breakpoints of the call step in all runs.

The no-library results of the paper come from the library-free layer alone, with the command above. They do not come from `evitea run --no-library`, which adds a repeat map and a record pass and keeps the resident-copy test on.

<a id="reproducing-the-paper"></a>

## 🔁 Reproducing the paper

This repository holds the code and the commands. The data, the ground-truth sets and the scoring scripts of the paper are kept elsewhere. Run the code of the paper first (`git checkout paper-2026-11-15`), then follow the table below for the inputs.

> [!NOTE]
> The paper's deposit on Figshare holds the files that are not public elsewhere. **Figshare DOI: to be added at publication.** An archived copy of release `paper-2026-11-15` will be on Zenodo. **Zenodo DOI: to be added at publication.**

**Inputs.** All sequencing reads are public. The accessions are in the supplement (Table S2J).

| Species | Reads (accessions of the samples) | Reference | TE library | Repeat annotation of the reference |
|---|---|---|---|---|
| *D. melanogaster* | AKA-017 SRR9951105 (training strain), GIM-024 SRR9951106, SLA-001 SRR9951097, MUN-009 SRR12187566, JUT-011 SRR9951100, ISO1 SRR11906527; Kuruka test lines mel33 SRR23103981, mel67 SRR23103904, mel21 SRR23103920 | dm6 | the 180 transposon sequences of FlyBase | RepeatMasker annotation of dm6, BED6 |
| *S. cerevisiae* | SK1 SRR4074258, Y12 SRR4074358, DBVPG6765 SRR4074257, UWOPS03-461.4 SRR4074384, DBVPG6044 SRR4074256 | R64 (sacCer3) | six Ty sequences of the McClintock yeast set | RepeatMasker 4.1.7-p1 hits of the six Ty sequences and Tsu4 (553 intervals) |
| *Z. tritici* | ISY92 SRR4907750, OregS90 SRR3740331, Aus01 SRR3740257 | IPO323 | 42 consensus sequences of Oggenfuss et al. 2021 | blastn of the 42 sequences on IPO323 (80% identity, 100 bp or longer; 2,372 intervals) |
| *A. thaliana* | Cvi-0 ERR3624578 | TAIR10 | 348 sequences derived from the library distributed with EDTA | TAIR10 TE annotation (28,576 intervals) |
| Mouse | CC055 ERR1942452 (ENA project PRJEB14673) | mm10 | 20 sequences that cover the 11 ground-truth families | UCSC RepeatMasker track of mm10, matched by repeat name to the 20 sequences |

The reference genomes are the public assemblies named in the table. The checksum of each file that the paper used is in Table S2A. The paper's TE libraries of *A. thaliana* and the mouse were made for the study and are in the deposit. The other libraries are public (FlyBase; Nelson et al. 2017; Oggenfuss et al. 2021). The tutorial above shows how to get the yeast files.

**Ground truth and scoring.** The ground-truth sets come from long-read assemblies: Table S2I and S2J give the assembly of each strain, and Supplementary Methods (ground-truth sets) describes how each set was made. The scripts that score calls are not in this repository. They are in the deposit: the one-to-one matching at 500 bp, the paired block bootstrap, the pair-preserving subsampling (seed 20261103), the synthetic heterozygous mixing and the ISO1 exclusion mask. The folder `validation/quarantined_scoring/` of this repository is empty.

**Data that is not public.** The `validation/` folder compares this code with the research code it came from, on data that is not distributed. It is not needed to reproduce the paper.

<a id="limits"></a>

## ⚠️ Limits

- Reads must be paired-end and all of one length.
- Only CPython 3.14 is supported. The library-free layer and the family-record commands run on Linux or WSL only.
- `prepare-sample` does not scale to full-size samples. Use `--msbwt-index` for those.
- The paper measured accuracy on the species and strains above. Accuracy on other species, and on read lengths other than those tested, is not known.

<a id="tests"></a>

## ✅ Tests

Install pytest and run `python -m pytest`. The suite includes an end-to-end test on synthetic data (`tests/integration/test_synthetic_call.py`). It must find every planted insertion with its family, with both FlyBase and plain TE names and with fly and non-fly contig names. `tools/synthetic/compare_call.sh` compares the library layer of two versions on a synthetic sample, byte for byte.

The scripts under `validation/` compare this code with the research code it came from. They need data that is not distributed. See [validation/README.md](validation/README.md). The profiles `AKA017` and `ISO1` of `--profile` exist only for those checks. The default, `production`, is the one to use.

<a id="citation"></a>

## 📚 Citation

If you use EviTEA, please cite the paper:

> Majumder MTH. EviTEA: calling transposable element insertions from short reads when the TE library is incomplete. Manuscript in preparation.

For the software itself, see [CITATION.cff](CITATION.cff).

<a id="license"></a>

## ⚖️ License

Copyright (C) 2026 Md Taksir Hasan Majumder. EviTEA is free software under the GNU General Public License, version 3 or later ([LICENSE](LICENSE)). If you distribute a modified version, it must also be under the GPL, with its source and the original copyright notice.

EviTEA uses msBWT (James Holt and Leonard McMillan), Bowtie2, BLAST+, mappy (minimap2), parasail, edlib and NumPy. [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES) lists their licenses.
