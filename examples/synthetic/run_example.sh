#!/bin/bash
# The public example: a small synthetic genome (five 200 kb contigs chr1..chr5 and a mitochondrial contig MT), a TE
# library of 120 records with ordinary names (B1_Mm_1, L1MdA_I_1, Ty1_1, ...) in 60 families with a two-column family map,
# and 125 bp paired reads (about 30x) of a sample that carries 30 non-reference insertions. About 10 minutes on 4 CPUs.
#   bash examples/synthetic/run_example.sh [OUTPUT_FOLDER]      (default examples/synthetic/run; PYTHON=python3)
# Requires EviTEA installed (or PYTHONPATH=src), Bowtie2 on PATH. The data are made by a script, so nothing large is stored here.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
PYTHON=${PYTHON:-python3}
OUT=${1:-$HERE/run}
DATA=$HERE/data

# 1. the input files
[ -d "$DATA" ] || "$PYTHON" "$ROOT/tools/synthetic/make_call_data.py" "$DATA" \
  --naming plain --contigs chr1,chr2,chr3,chr4,chr5 --mito MT

# 2. FASTQ to calls in one command (the same steps as prepare-reference, prepare-te, prepare-sample and call)
"$PYTHON" -m evitea run --reference "$DATA/reference.fa" --repeat-annotation "$DATA/repeat_annotation.bed" \
  --te-library "$DATA/te_library.fa" --family-map "$DATA/family_map.tsv" \
  --read1 "$DATA/sample_R1.fastq" --read2 "$DATA/sample_R2.fastq" --sample-id synth --output "$OUT" --threads 4

# 3. the result: calls.native.tsv equals the expected one; every planted insertion is identified
if cmp -s "$OUT/call/calls.native.tsv" "$HERE/expected/calls.native.tsv"; then
  echo "calls.native.tsv: as expected"
else
  echo "calls.native.tsv differs from examples/synthetic/expected/calls.native.tsv"; exit 1
fi
