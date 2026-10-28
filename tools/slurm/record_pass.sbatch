#!/bin/bash
# The library layer once more on the family records (evitea call with the sample's own command line, only --te-assets, --output
# and --work-dir changed), then add-record-calls, which adds the record calls not yet in the existing call set.
# Set before sbatch:
#   DS    sample name (its library run's run_manifest.json gives the command line)
#   REC   the records.sbatch output folder (records/<DS>_<job>)
#   BASE  the detect.sbatch output folder whose library_free/calls.combined.tsv is the existing call set
# Resources are those of the sample's own call run plus the records; change them on the sbatch line if needed.
#SBATCH --job-name=evitea_record_pass
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=96g
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
. "$EVITEA_CODE/tools/slurm/env.sh"
: "${DS:?}" "${REC:?}" "${BASE:?}"
OUT=$EVITEA_OUT/record_pass/${DS}_${SLURM_JOB_ID}
[ -e "$OUT" ] && { echo "REFUSING $OUT exists"; exit 2; }
mkdir -p "$OUT"
mapfile -d '' ARGS < <($PY "$EVITEA_CODE/tools/slurm/call_args.py" "$RUN/run_manifest.json" "$REC/assets_te" "$OUT/call" "$OUT/work")
echo "host $(hostname) $(date -u +%FT%TZ) id=${SLURM_JOB_ID:-} DS=$DS records=$REC base=$BASE"
$PY -m evitea call "${ARGS[@]}" --keep-resident-copies   # the record names never match an annotation name: no own-family copy filter
$PY -m evitea add-record-calls --calls "$OUT/call/calls.native.tsv" --base-calls "$BASE/library_free/calls.combined.tsv" \
    --records "$REC/lib" --output "$OUT/added"
( cd "$OUT" && find . -type f -not -path './work/*' -print0 | sort -z | xargs -0 sha256sum > "$OUT/OUTPUT_SHA256.txt" )
echo "done $(date -u +%FT%TZ)"
