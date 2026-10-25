#!/bin/bash
# The run without a library, as one chain of jobs: repeat map (once per reference), detection with the library layer's calls left
# out, family records with the same map, and the record pass. Run it on a login node: it only submits (sbatch --dependency).
# Set EVITEA_WORK (and EVITEA_OUT, EVITEA_CODE if not at their defaults), then, from $EVITEA_OUT:
#   DS=<sample> [CONTIGS=...] [TAG=nolib] bash $EVITEA_CODE/tools/slurm/no_library.sh
# The map is built only when $MAPDIR/repeat_map.bed (default $EVITEA_OUT/repeat_map) does not exist yet. The result is
# record_pass/<DS>_<job>/added/calls.combined.tsv.
set -euo pipefail
: "${EVITEA_WORK:?set EVITEA_WORK}" "${DS:?set DS}"
EVITEA_CODE=${EVITEA_CODE:-$EVITEA_WORK/evitea_code}
EVITEA_OUT=${EVITEA_OUT:-$EVITEA_WORK/evitea_runs}
S=$EVITEA_CODE/tools/slurm
TAG=${TAG:-nolib}
MAPDIR=${MAPDIR:-$EVITEA_OUT/repeat_map}
export EVITEA_WORK EVITEA_OUT EVITEA_CODE DS CONTIGS MAPDIR
cd "$EVITEA_OUT"
DEP=""
if [ ! -e "$MAPDIR/repeat_map.bed" ]; then
  MAPJOB=$(sbatch --parsable "$S/repeat_map.sbatch")
  DEP="--dependency=afterok:$MAPJOB"
  echo "repeat map: job $MAPJOB"
fi
DETJOB=$(TAG=$TAG NOLIB=1 REPEAT_MAP=$MAPDIR/repeat_map.bed sbatch --parsable $DEP "$S/detect.sbatch")
BASE=$EVITEA_OUT/runs/${DS}_${TAG}_${DETJOB}
echo "detection: job $DETJOB -> $BASE"
RECJOB=$(SRC=$BASE FLAG=kept REPEAT_MAP=$MAPDIR/repeat_map.bed sbatch --parsable --dependency=afterok:$DETJOB "$S/records.sbatch")
REC=$EVITEA_OUT/records/${DS}_${RECJOB}
echo "records: job $RECJOB -> $REC"
PASSJOB=$(REC=$REC BASE=$BASE sbatch --parsable --dependency=afterok:$RECJOB "$S/record_pass.sbatch")
echo "record pass: job $PASSJOB -> $EVITEA_OUT/record_pass/${DS}_${PASSJOB}/added"
