#!/bin/bash
# The library layer (evitea call) of this version and of a previous one on the same synthetic sample, compared byte for byte.
#   PREVIOUS_SRC=<src folder of the previous version> bash tools/synthetic/compare_call.sh [WORK]
# Needs Bowtie2 on PATH and a python ($PYTHON, default python3) with pymsbwt and numpy. Takes about 10 minutes on 4 CPUs.
# Both versions get the same prepared sample (built by this version: every msBWT package builds the same index) and run
# prepare-reference, prepare-te --build-languages and call. This version is given the development values
# (--read-length 125 --fragment-q95 494 --fragment-q99 805 --primary-contigs chr2L,chr2R,chr3L,chr3R,chrX); the TE
# languages, their audit tables, calls.native.tsv, calls.tsv, calls.bed, unresolved.tsv, run_summary.json and run.log must
# then be identical. run_manifest.json is not compared (it records times, paths and the values used).
set -euo pipefail
: "${PREVIOUS_SRC:?set PREVIOUS_SRC to the src folder of the previous version}"
cd "$(dirname "$0")/../.."
HERE=$PWD
PYTHON=${PYTHON:-python3}
WORK=${1:-$(mktemp -d)}
mkdir -p "$WORK"
ev() { PYTHONPATH=$1 "$PYTHON" -c "import sys; from evitea.cli import main; sys.exit(main(sys.argv[1:]))" "${@:2}"; }
DATA=$WORK/data
[ -d "$DATA" ] || "$PYTHON" tools/synthetic/make_call_data.py "$DATA"
[ -d "$WORK/sample" ] || ev "$HERE/src" prepare-sample --read1 "$DATA/sample_R1.fastq" --read2 "$DATA/sample_R2.fastq" \
  --output "$WORK/sample" --sample-id synth --threads 4 > "$WORK/sample.log" 2>&1
run() {  # run <src> <out> [call options]
  local src=$1 out=$2; shift 2
  mkdir -p "$out"
  [ -d "$out/ref" ] || ev "$src" prepare-reference --reference "$DATA/reference.fa" --repeat-annotation "$DATA/repeat_annotation.bed" \
    --output "$out/ref" > "$out/ref.log" 2>&1
  [ -d "$out/te" ] || ev "$src" prepare-te --te-library "$DATA/te_library.fa" --family-map "$DATA/family_map.tsv" --build-languages \
    --reference-assets "$out/ref" --output "$out/te" > "$out/te.log" 2>&1
  rm -rf "$out/call"
  ev "$src" call --sample-assets "$WORK/sample" --reference-assets "$out/ref" --te-assets "$out/te" --output "$out/call" \
    --threads 4 --workers 2 "$@" > "$out/call.log" 2>&1
}
run "$PREVIOUS_SRC" "$WORK/previous"
run "$HERE/src" "$WORK/current" --read-length 125 --fragment-q95 494 --fragment-q99 805 \
  --primary-contigs chr2L,chr2R,chr3L,chr3R,chrX
fail=0
for f in $(cd "$WORK/previous" && ls te/languages/*.tsv te/languages/audit/*) call/calls.native.tsv call/calls.tsv call/calls.bed \
    call/unresolved.tsv call/run_summary.json call/run.log; do
  if cmp -s "$WORK/previous/$f" "$WORK/current/$f"; then echo "SAME $f"; else echo "DIFF $f"; fail=1; fi
done
[ "$fail" = 0 ] && echo "ALL SAME ($WORK)" || echo "DIFFERENCES ($WORK)"
exit $fail
