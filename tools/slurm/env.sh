# Settings shared by the Slurm jobs (source this file; it changes nothing).
# Every cluster path is built from EVITEA_WORK. Site-specific locations (container, reference, samples) are read from a site file:
# copy site.example.sh to $EVITEA_OUT/site.sh and edit it.
: "${EVITEA_WORK:?set EVITEA_WORK to your working folder on the cluster}"
EVITEA_CODE=${EVITEA_CODE:-$EVITEA_WORK/evitea_code}   # a clone of this repository
EVITEA_OUT=${EVITEA_OUT:-$EVITEA_WORK/evitea_runs}     # everything the jobs write goes below this folder
case "$(hostname)" in *login*) echo "STOP: login node"; exit 9;; esac
. "${EVITEA_SITE:-$EVITEA_OUT/site.sh}"
module load apptainer; unset -f which 2>/dev/null || true
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
[ -n "${BLAST_BIN:-}" ] && export APPTAINERENV_PREPEND_PATH=$BLAST_BIN   # BLAST+ (makeblastdb, blastn) for build-repeat-map
export APPTAINERENV_PYTHONPATH=$EVITEA_CODE/src${MSBWT_PATH:+:$MSBWT_PATH}
# the python of the container that has mappy, parasail and numpy
if [ -z "${PYBIN:-}" ]; then
  for cand in $(apptainer exec "$SIF" sh -c 'ls /opt/micromamba/envs/*/bin/python 2>/dev/null'); do
    apptainer exec "$SIF" "$cand" -c 'import mappy, parasail, numpy' 2>/dev/null && { PYBIN=$cand; break; }
  done
fi
: "${PYBIN:?no python with mappy, parasail and numpy in the container}"
PY="apptainer exec --bind $EVITEA_WORK $SIF $PYBIN -B"
# sample table: DS -> MSBWT index, library run (--run-dir), default contigs
if [ -n "${DS:-}" ]; then
  IFS=$'\t' read -r _ IDX RUN CT < <(awk -F'\t' -v d="$DS" '$1==d' "$SAMPLES") || true
  [ -n "${IDX:-}" ] || { echo "unknown DS=$DS (see $SAMPLES)"; exit 2; }
  IDX=$EVITEA_WORK/$IDX; RUN=$EVITEA_WORK/$RUN
fi
CONTIGS=${CONTIGS:-${CT:-}}
# the option that drops families of old, diverged copies (family-identity filter); its present name
FAMILY_IDENTITY_FILTER=--rule-c
