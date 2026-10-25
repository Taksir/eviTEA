# Examples

## A worked example on synthetic data

`synthetic/run_example.sh` makes a small data set and calls it, start to finish:

```bash
bash examples/synthetic/run_example.sh            # about 10 minutes on 4 CPUs; Bowtie2 on PATH
```

The data set (made by `tools/synthetic/make_call_data.py`, the same files for the same arguments):

| file | content |
|---|---|
| `reference.fa` | five 200 kb contigs `chr1`..`chr5` and a 16 kb mitochondrial contig `MT` |
| `repeat_annotation.bed` | BED6 of the TE copies placed in the reference (name = TE record) |
| `te_library.fa` | 120 TE records with ordinary names (`B1_Mm_1`, `L1MdA_I_1`, `Ty1_1`, ...) in 60 families |
| `family_map.tsv` | the two-column family map, `record<TAB>family` |
| `sample_R1.fastq`, `sample_R2.fastq` | 125 bp paired reads, about 30x, fragments of 350 +- 40 bp |
| `truth_insertions.tsv` | the 30 insertions planted in the sample |

The script runs

```bash
evitea run --reference data/reference.fa --repeat-annotation data/repeat_annotation.bed \
  --te-library data/te_library.fa --family-map data/family_map.tsv \
  --read1 data/sample_R1.fastq --read2 data/sample_R2.fastq --sample-id synth --output run --threads 4
```

which is the same as the four separate commands

```bash
evitea prepare-reference --reference data/reference.fa --repeat-annotation data/repeat_annotation.bed --output run/reference
evitea prepare-te --te-library data/te_library.fa --family-map data/family_map.tsv --build-languages \
  --reference-assets run/reference --output run/te
evitea prepare-sample --read1 data/sample_R1.fastq --read2 data/sample_R2.fastq --sample-id synth --output run/sample
evitea call --sample-assets run/sample --reference-assets run/reference --te-assets run/te --output run/call
```

and checks that `run/call/calls.native.tsv` equals `synthetic/expected/calls.native.tsv`: 30 calls, one for each planted
insertion (`synthetic/expected/truth_insertions.tsv`), each with its family. `run/call/run_manifest.json` records the values
taken from the sample and the reference (`values.run_parameters`): read length 125 bp, fragment-size 95th / 99th percentiles
415 / 442 bp, primary contigs `chr1`..`chr5`.

This is a test that the programs work together, not a measure of accuracy: the reads have no sequencing errors and the
genome is small and random.

## Configuration template

[`evitea.raw.example.toml`](evitea.raw.example.toml) is a path-relative raw
paired-FASTQ configuration template. Copy it beside the data/configuration for
one run and replace the placeholder paths; do not add scientific thresholds.
