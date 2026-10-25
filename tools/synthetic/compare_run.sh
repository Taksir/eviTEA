#!/bin/bash
# `evitea run` against the separate commands on the synthetic data set, with a library (and library-free) and without one
# (against tools/no_library.sh). The outputs must be identical.
#   bash tools/synthetic/compare_run.sh [WORK]          (PYTHON=python3; Bowtie2, BLAST+, mappy and parasail needed)
set -euo pipefail
cd "$(dirname "$0")/../.."
HERE=$PWD
PYTHON=${PYTHON:-python3}
export PYTHONPATH=$HERE/src
WORK=${1:-$(mktemp -d)}
mkdir -p "$WORK"
ev() { "$PYTHON" -m evitea "$@"; }
D=$WORK/data
[ -d "$D" ] || "$PYTHON" tools/synthetic/make_call_data.py "$D"
C=chr2L,chr2R,chr3L,chr3R,chrX
COMMON=(--reference "$D/reference.fa" --repeat-annotation "$D/repeat_annotation.bed" --read1 "$D/sample_R1.fastq"
        --read2 "$D/sample_R2.fastq" --sample-id synth --threads 4 --workers 2)

# with a library
ev run "${COMMON[@]}" --te-library "$D/te_library.fa" --family-map "$D/family_map.tsv" --library-free --contigs "$C" \
  --output "$WORK/run_lib" 2> "$WORK/run_lib.log"
S=$WORK/sep_lib
if [ ! -d "$S/call" ]; then
  ev prepare-reference --reference "$D/reference.fa" --repeat-annotation "$D/repeat_annotation.bed" --output "$S/reference" --threads 4
  ev prepare-te --te-library "$D/te_library.fa" --family-map "$D/family_map.tsv" --build-languages --reference-assets "$S/reference" \
    --output "$S/te" 2> "$S.te.log"
  ev prepare-sample --read1 "$D/sample_R1.fastq" --read2 "$D/sample_R2.fastq" --sample-id synth --output "$S/sample" --threads 4
  ev call --sample-assets "$S/sample" --reference-assets "$S/reference" --te-assets "$S/te" --output "$S/call" \
    --work-dir "$S/work/call" --threads 4 --workers 2 2> "$S.call.log"
  ev library-free --run-dir "$S/call" --sample-assets "$S/sample" --reference-assets "$S/reference" --output "$S/library_free" \
    --work-dir "$S/work/library_free" --contigs "$C" 2> "$S.lf.log"
fi

# without a library
ev run "${COMMON[@]}" --no-library --contigs "$C" --output "$WORK/run_nolib" 2> "$WORK/run_nolib.log"
N=$WORK/sep_nolib
if [ ! -d "$N/record_pass" ]; then
  mkdir -p "$N"
  ev prepare-reference --reference "$D/reference.fa" --repeat-annotation "$D/repeat_annotation.bed" --output "$N/reference" --threads 4
  ev prepare-sample --read1 "$D/sample_R1.fastq" --read2 "$D/sample_R2.fastq" --sample-id synth --output "$N/sample" --threads 4
  ASSETS=$N/reference SAMPLE=$N/sample DS=synth OUT=$N CONTIGS=$C CALL_OPTS="--threads 4 --workers 2" PYTHON=$PYTHON \
    bash tools/no_library.sh 2> "$N.log"
fi

fail=0
same() {  # same <A> <B> <relative files...>
  local a=$1 b=$2; shift 2
  for f in "$@"; do
    if [ ! -e "$a/$f" ] && [ ! -e "$b/$f" ]; then continue; fi
    if cmp -s "$a/$f" "$b/$f"; then echo "SAME $f"; else echo "DIFF $f"; fail=1; fi
  done
}
CALL_FILES=(calls.native.tsv calls.tsv calls.bed calls.vcf unresolved.tsv run_summary.json run.log)
LF_FILES=(library_free_places.tsv library_free_segments.fa library_free_second_places.tsv library_free_second_segments.fa
          library_free_family_calls.tsv calls.library_free.tsv calls.combined.tsv calls.combined.vcf)
same "$WORK/run_lib/call" "$S/call" "${CALL_FILES[@]}"
same "$WORK/run_lib/library_free" "$S/library_free" "${LF_FILES[@]}"
same "$WORK/run_lib/te/languages" "$S/te/languages" $(cd "$S/te/languages" && ls *.tsv)
same "$WORK/run_nolib/repeat_map" "$N/repeat_map" repeat_map.bed
same "$WORK/run_nolib/library_free" "$N/library_free" "${LF_FILES[@]}"
same "$WORK/run_nolib/records/lib" "$N/records/lib" $(cd "$N/records/lib" && ls *.tsv *.fa *.bed 2>/dev/null)
same "$WORK/run_nolib/record_pass/call" "$N/record_pass/call" "${CALL_FILES[@]}"
same "$WORK/run_nolib/record_pass/added" "$N/record_pass/added" $(cd "$N/record_pass/added" && ls)
[ "$fail" = 0 ] && echo "ALL SAME ($WORK)" || echo "DIFFERENCES ($WORK)"
exit $fail
