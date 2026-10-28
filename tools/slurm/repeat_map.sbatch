#!/bin/bash
# The repeat map of the reference (evitea build-repeat-map): BLAST+ alignment of 500 bp windows of the chosen contigs to the whole
# reference. One map per reference; every run without a library reads it. Requires BLAST+ in the image or BLAST_BIN in site.sh.
# Writes $MAPDIR/repeat_map.bed (default MAPDIR=$EVITEA_OUT/repeat_map) and refuses to replace one.
# Set before sbatch:
#   CONTIGS  comma-separated contigs to align (default: those library-free scans by default, the primary chromosomes and the
#            contigs of 1 Mb or more)
#   MAPDIR   output folder
#SBATCH --job-name=evitea_repeat_map
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
. "$EVITEA_CODE/tools/slurm/env.sh"
MAPDIR=${MAPDIR:-$EVITEA_OUT/repeat_map}
[ -e "$MAPDIR/repeat_map.bed" ] && { echo "REFUSING $MAPDIR/repeat_map.bed exists"; exit 2; }
echo "host $(hostname) $(date -u +%FT%TZ) id=${SLURM_JOB_ID:-} map=$MAPDIR contigs=${CONTIGS:-default}"
sha256sum "$REFFA"
$PY -m evitea build-repeat-map --reference-fasta "$REFFA" --output "$MAPDIR" --parts "${SLURM_CPUS_PER_TASK:-8}" \
    ${CONTIGS:+--contigs "$CONTIGS"}
sha256sum "$MAPDIR/repeat_map.bed"
echo "done $(date -u +%FT%TZ)"
