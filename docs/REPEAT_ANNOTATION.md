# The repeat annotation, and runs without a TE library

## The repeat annotation of `prepare-reference`

`prepare-reference --repeat-annotation` takes the TE copies of the reference genome as BED6:

```text
contig<TAB>start (0-based)<TAB>end<TAB>TE record<TAB>score<TAB>strand (+, - or .)
```

`call` needs it: it tells the library layer which reference places are old copies of a TE family. The contig names must be
those of the reference FASTA, and the name column must name the TE record (the first word of the TE library's FASTA header;
a FlyBase id such as `FBte0000001` anywhere in the name also works). A name such as `L1MdA_I#LINE/L1` is read as the
record `L1MdA_I`.

### Making it with RepeatMasker

Mask the reference with the same TE library that `prepare-te` gets, then convert RepeatMasker's `.out` table:

```bash
RepeatMasker -pa 8 -lib te_library.fa -nolow -norna -no_is -dir rm reference.fa
# rm/reference.fa.out: three header lines, then one line per copy. Columns used: 1 score, 5 contig, 6 begin (1-based),
# 7 end, 9 strand (+ or C for the reverse strand), 10 repeat name (the TE record).
awk 'BEGIN{OFS="\t"} NR>3 && NF>=15 {print $5, $6-1, $7, $10, $1, ($9=="C" ? "-" : "+")}' rm/reference.fa.out \
  | sort -k1,1 -k2,2n > repeat_annotation.bed
```

`-nolow -norna` leave out simple repeats, low-complexity sequence and RNA genes, which are not TE copies. If you use another
masker, write the same six columns.

The library-free layer (`library-free --repeat-annotation`) and `filter-resident-copies` read the RepeatMasker `.out` table
directly as well.

## Runs without a TE library

A species may have no curated TE library. EviTEA can still call insertions:

1. `build-repeat-map` finds the repeated part of the reference by aligning the reference to itself with BLAST+
   (`makeblastdb` and `blastn` on `PATH`). Its `repeat_map.bed` stands in for the TE annotation in the resident-copy test of
   the library-free layer.
2. `library-free` runs on an empty library call set, with the family-identity filter and `--repeat-map`.
3. `build-family-records` turns the library-free families into TE records (with their reference copies and end-word
   codebooks), and `prepare-te` prepares them.
4. `call` runs the library layer once more, on those records; `add-record-calls` adds its new calls to the library-free call
   set. The result is `record_pass/added/calls.combined.tsv`.

All four steps are one command:

```bash
evitea run --no-library --reference reference.fa --repeat-annotation repeat_annotation.bed \
  --read1 s_R1.fastq.gz --read2 s_R2.fastq.gz --sample-id s --output runs/s --threads 8
```

The same chain is `tools/no_library.sh` (a plain shell script, on prepared assets) and `tools/slurm/no_library.sh` (Slurm jobs).
The library pass of step 4 reads the reference's repeat annotation like any `call`, so `--repeat-annotation` (BED6, as above)
is needed here too. Without a curated library, make it with RepeatMasker and a library of the species' repeats from a de novo
tool (for example RepeatModeler). The repeat map cannot take its place: it is a three-column BED of repeated intervals with
no TE names.

`build-repeat-map --contigs` limits the map to some contigs (default: contigs of 1 Mb or more, the contigs the library-free
layer scans by default); `--parts` sets the number of blastn processes (the map does not depend on it).
