#!/bin/bash
# The genotype step on an existing library-free output (evitea genotype): genotypes.tsv for every call of calls.combined.tsv. Reads the call
# files, writes only under $EVITEA_OUT/genotype/<DS>_<TAG>_<jobid>/. No scan, no change to any call file.
# Set before sbatch:
#   DS      sample name in samples.tsv (gives the msBWT index $IDX)
#   TAG     short label of the run
#   LFDIR   the finished library-free output folder (the folder that holds calls.combined.tsv)
#   COUNTER fast or plain (--counter; default fast)
#SBATCH --job-name=evitea_genotype
#SBATCH --partition=small
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
. "$EVITEA_CODE/tools/slurm/env.sh"
: "${DS:?}" "${TAG:?}" "${LFDIR:?}"
[ -f "$LFDIR/calls.combined.tsv" ] || { echo "no calls.combined.tsv in $LFDIR"; exit 2; }
OUT=$EVITEA_OUT/genotype/${DS}_${TAG}_${SLURM_JOB_ID}
[ -e "$OUT" ] && { echo "REFUSING $OUT exists"; exit 2; }
echo "host $(hostname) $(date -u +%FT%TZ) id=${SLURM_JOB_ID:-} DS=$DS TAG=$TAG LFDIR=$LFDIR"
mkdir -p "$OUT"
START=$(date +%s)
sha256sum "$IDX/about.npy" "$REFFA" "$LFDIR/calls.combined.tsv" "$LFDIR/calls.library_free.tsv" "$LFDIR/library_free_family_calls.tsv"
$PY -m evitea genotype --library-free-dir "$LFDIR" --msbwt-index "$IDX" --reference-assets "$ASSETS" --output "$OUT" \
    ${COUNTER:+--counter $COUNTER}
echo "wall_seconds $(( $(date +%s) - START ))"
( cd "$OUT" && find . -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/OUTPUT_SHA256.txt" )
echo "done $(date -u +%FT%TZ)"
