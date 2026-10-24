#!/bin/bash
# Family records of one sample (build-family-records, then prepare-te). The library run on the records (evitea call) is a separate
# step: copy the command line of the sample's existing call run and change only --te-assets, --output and --work-dir.
# Set before sbatch:
#   DS    sample name                     SRC   the detect.sbatch output folder whose family calls are used
#   FLAG  flag rule: kept (default), before-filters or both
#   RM    repeat annotation (RepeatMasker .out or BED) for the record guard
#   COPYREP=1  also write the copy report
#   FC    space-separated family-call tables (default: the one in SRC)      SEGFILE  segments file (default: the one in SRC)
#SBATCH --job-name=evitea_records
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
. "$EVITEA_CODE/tools/slurm/env.sh"
: "${DS:?}" "${SRC:?}"
OUT=$EVITEA_OUT/records/${DS}_${SLURM_JOB_ID}
[ -e "$OUT" ] && { echo "REFUSING $OUT exists"; exit 2; }
mkdir -p "$OUT"
EXTRA=""
[ -n "${RM:-}" ] && EXTRA="$EXTRA --repeat-annotation $RM"
[ -n "${COPYREP:-}" ] && EXTRA="$EXTRA --copy-report on"
SEG2=""
[ -z "${FC:-}" ] && [ -f "$SRC/library_free/library_free_second_segments.fa" ] && SEG2="--second-segments $SRC/library_free/library_free_second_segments.fa"
FCS=${FC:-$SRC/library_free/library_free_family_calls.tsv}
$PY -m evitea build-family-records --family-calls $FCS \
    --segments "${SEGFILE:-$SRC/library_free/library_free_segments.fa}" $SEG2 --reference-assets "$ASSETS" --output "$OUT/lib" \
    --flag-rule "${FLAG:-kept}" $EXTRA
$PY -m evitea prepare-te --te-library "$OUT/lib/te_library.fa" --family-map "$OUT/lib/family_map_padded.tsv" \
    --terminal-codebook "$OUT/lib/terminal_codebook.tsv" --internal-codebook "$OUT/lib/internal_only_codebook.tsv" \
    --reference-assets "$ASSETS" --output "$OUT/assets_te"
( cd "$OUT" && find . -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/OUTPUT_SHA256.txt" )
echo "done $(date -u +%FT%TZ)"
