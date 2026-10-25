#!/bin/bash
# Runs EviTEA without a library as one plain shell chain (no Slurm). The steps are: repeat map (once per reference), detection
# with an empty library call set, family records with the same map, and the record pass. The steps, order and arguments match
# tools/slurm/no_library.sh and `evitea run --no-library`.
#
#   ASSETS=<prepared reference assets> SAMPLE=<prepared sample assets> DS=<sample id> OUT=<new folder> \
#     [CONTIGS=chr1,chr2] [MAPDIR=<repeat map folder>] [CALL_OPTS="--threads 8"] [PYTHON=python3] bash tools/no_library.sh
# or, with no repeat annotation at all, REFERENCE=<reference FASTA> instead of ASSETS. The script builds the repeat map from the
# FASTA first. Its BED6 copy (repeat_map_bed6.bed) becomes the annotation of the reference prepared in $OUT/reference, as in
# `evitea run --no-library` without --repeat-annotation.
#
# The map is built only when $MAPDIR/repeat_map.bed (default $OUT/repeat_map) does not exist. Building it requires BLAST+ on PATH.
# The result is $OUT/record_pass/added/calls.combined.tsv.
set -euo pipefail
: "${SAMPLE:?set SAMPLE}" "${DS:?set DS}" "${OUT:?set OUT}"
[ -n "${ASSETS:-}${REFERENCE:-}" ] || { echo "set ASSETS (prepared reference) or REFERENCE (reference FASTA)"; exit 2; }
PYTHON=${PYTHON:-python3}
MAPDIR=${MAPDIR:-$OUT/repeat_map}
CONTIG_OPTS=()
[ -n "${CONTIGS:-}" ] && CONTIG_OPTS=(--contigs "$CONTIGS")
read -r -a CALL_EXTRA <<< "${CALL_OPTS:-}"
ev() { "$PYTHON" -m evitea "$@"; }
mkdir -p "$OUT"
if [ -z "${ASSETS:-}" ]; then
  if [ ! -e "$MAPDIR/repeat_map.bed" ]; then
    echo "repeat map (from the FASTA) -> $MAPDIR"
    ev build-repeat-map --reference-fasta "$REFERENCE" --output "$MAPDIR" "${CONTIG_OPTS[@]}"
  fi
  ASSETS=$OUT/reference
  if [ ! -e "$ASSETS/reference_manifest.json" ]; then
    echo "reference, annotated by the repeat map -> $ASSETS"
    ev prepare-reference --reference "$REFERENCE" --repeat-annotation "$MAPDIR/repeat_map_bed6.bed" --output "$ASSETS"
  fi
elif [ ! -e "$MAPDIR/repeat_map.bed" ]; then
  echo "repeat map -> $MAPDIR"
  ev build-repeat-map --reference-assets "$ASSETS" --output "$MAPDIR" "${CONTIG_OPTS[@]}"
fi
MAP=$MAPDIR/repeat_map.bed
echo "empty library call set -> $OUT/nolib_stub"
"$PYTHON" -c 'import sys; from pathlib import Path; from evitea.run_command import write_stub_run; write_stub_run(Path(sys.argv[1]), sys.argv[2])' \
  "$OUT/nolib_stub" "$DS"
echo "detection -> $OUT/library_free"
ev library-free --run-dir "$OUT/nolib_stub" --sample-assets "$SAMPLE" --reference-assets "$ASSETS" --output "$OUT/library_free" \
  --work-dir "$OUT/work/library_free" "${CONTIG_OPTS[@]}" --rule-c --repeat-map "$MAP"
echo "records -> $OUT/records"
SEG2=()
[ -f "$OUT/library_free/library_free_second_segments.fa" ] && SEG2=(--second-segments "$OUT/library_free/library_free_second_segments.fa")
ev build-family-records --family-calls "$OUT/library_free/library_free_family_calls.tsv" \
  --segments "$OUT/library_free/library_free_segments.fa" "${SEG2[@]}" --reference-assets "$ASSETS" --output "$OUT/records/lib" \
  --flag-rule kept --repeat-map "$MAP"
ev prepare-te --te-library "$OUT/records/lib/te_library.fa" --family-map "$OUT/records/lib/family_map_padded.tsv" \
  --terminal-codebook "$OUT/records/lib/terminal_codebook.tsv" --internal-codebook "$OUT/records/lib/internal_only_codebook.tsv" \
  --reference-assets "$ASSETS" --output "$OUT/records/assets_te"
echo "record pass -> $OUT/record_pass"
ev call --sample-assets "$SAMPLE" --reference-assets "$ASSETS" --te-assets "$OUT/records/assets_te" \
  --output "$OUT/record_pass/call" --work-dir "$OUT/work/record_pass" "${CALL_EXTRA[@]}"
ev add-record-calls --calls "$OUT/record_pass/call/calls.native.tsv" --base-calls "$OUT/library_free/calls.combined.tsv" \
  --records "$OUT/records/lib" --output "$OUT/record_pass/added"
echo "done: $OUT/record_pass/added/calls.combined.tsv"
