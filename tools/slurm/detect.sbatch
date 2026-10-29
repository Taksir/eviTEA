#!/bin/bash
# One library-free run. Writes only under $EVITEA_OUT/runs/<DS>_<TAG>_<jobid>/ (and a work folder under $EVITEA_OUT/work/).
# Set before sbatch:
#   DS      sample name in samples.tsv          TAG    short label of the run
#   OPTS    library-free switches, quoted (default: all on)
#   CONTIGS comma-separated contigs (default: the sample's entry in samples.tsv)
#   WORK    existing work folder to reuse saved scan pieces (default a new one named <DS>_<TAG>)
#   REGROUP=1   set EVITEA_PIECES_READONLY: the run stops instead of scanning again when switches, code or inputs differ from the saved pieces
#   DET     a --detection-dir (rebuild family calls from saved tables)
#   RUNDIR  another library run for --run-dir (for example a run with some families removed from the library)
#   DEPTHCAP value of --depth-cap (default none; 100 is the previous version's window)
#   PREVIOUS=1  run the previous version's code (PREVIOUS_SRC) instead; OPTS and DEPTHCAP must then be empty
#   NOLIB=1     run without a library: the library layer's calls are an empty table (a stub of the sample's library run, kept under
#               $EVITEA_OUT/nolib/<DS>_stub); the one-site check then takes its tolerance from the layer's own two-sided calls
#   REPEAT_MAP  the repeat map (build-repeat-map output, repeat_map.bed): the resident-copy test reads it when the library layer made no call
#   WORKERS     processes that scan a contig (--workers; default the job's CPUs, SLURM_CPUS_PER_TASK). A contig is cut into tiles that the
#               processes walk; the places and every output table are the same for any number. Not given with PREVIOUS=1 (the previous
#               version has no such option)
#   COUNTER     fast or plain (--counter; default fast): how the read index is searched. fast gives the same counts with fewer searches
#   REFCACHE    folder for --reference-cache: the reference word tables, built once and loaded by later jobs (default refcache/ in the work folder)
#   TE_ASSETS   prepared TE assets of the library run (--te-assets): calls.combined.* leave out the library calls inside a reference copy of
#               their own family; required when the run has library calls, unless OPTS holds --keep-resident-copies
#SBATCH --job-name=evitea_detect
#SBATCH --partition=small
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
. "$EVITEA_CODE/tools/slurm/env.sh"
: "${DS:?}" "${TAG:?}"
if [ -n "${PREVIOUS:-}" ]; then
  [ -z "${OPTS:-}${DEPTHCAP:-}" ] || { echo "PREVIOUS=1 takes no OPTS or DEPTHCAP"; exit 2; }
  [ -z "${NOLIB:-}${REPEAT_MAP:-}" ] || { echo "PREVIOUS=1 takes no NOLIB or REPEAT_MAP"; exit 2; }
  export APPTAINERENV_PYTHONPATH=$PREVIOUS_SRC${MSBWT_PATH:+:$MSBWT_PATH}
fi
OUT=$EVITEA_OUT/runs/${DS}_${TAG}_${SLURM_JOB_ID}
[ -e "$OUT" ] && { echo "REFUSING $OUT exists"; exit 2; }
WORKDIR=${WORK:-$EVITEA_OUT/work/${DS}_${TAG}}
[ -n "${REGROUP:-}" ] && export APPTAINERENV_EVITEA_PIECES_READONLY=1
RUN=${RUNDIR:-$RUN}
if [ -n "${NOLIB:-}" ]; then
  STUB=$EVITEA_OUT/nolib/${DS}_stub
  if [ ! -e "$STUB/calls.native.tsv" ]; then
    mkdir -p "$STUB"
    head -n 1 "$RUN/calls.native.tsv" > "$STUB/calls.native.tsv"
    $PY -c 'import json,sys; json.dump({"sample_id": json.load(open(sys.argv[1]))["sample_id"]}, open(sys.argv[2], "w"))' \
        "$RUN/run_summary.json" "$STUB/run_summary.json"
  fi
  RUN=$STUB
fi
echo "host $(hostname) $(date -u +%FT%TZ) id=${SLURM_JOB_ID:-} DS=$DS TAG=$TAG OPTS=${OPTS:-} work=$WORKDIR"
mkdir -p "$OUT" "$EVITEA_OUT/work"
START=$(date +%s)
sha256sum "$IDX/about.npy" "$REFFA"
EXTRA=""
[ -n "${DET:-}" ] && EXTRA="$EXTRA --detection-dir $DET"
[ -n "${DEPTHCAP:-}" ] && EXTRA="$EXTRA --depth-cap $DEPTHCAP"
[ -n "${REPEAT_MAP:-}" ] && EXTRA="$EXTRA --repeat-map $REPEAT_MAP"
# a run with library calls: TE_ASSETS (prepared TE assets) applies the own-family reference copy filter to calls.combined.*;
# pass --keep-resident-copies in OPTS instead to keep those calls
[ -n "${TE_ASSETS:-}" ] && EXTRA="$EXTRA --te-assets $TE_ASSETS"
# scan speed: not part of the identity of a saved piece, and not given to the previous version
if [ -z "${PREVIOUS:-}" ]; then
  EXTRA="$EXTRA --workers ${WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
  [ -n "${COUNTER:-}" ] && EXTRA="$EXTRA --counter $COUNTER"
  [ -n "${REFCACHE:-}" ] && EXTRA="$EXTRA --reference-cache $REFCACHE"
fi
$PY -m evitea library-free --run-dir "$RUN" --msbwt-index "$IDX" --reference-assets "$ASSETS" --output "$OUT/library_free" \
    --work-dir "$WORKDIR" --contigs "$CONTIGS" $FAMILY_IDENTITY_FILTER ${OPTS:-} $EXTRA
echo "wall_seconds $(( $(date +%s) - START ))"
( cd "$OUT" && find . -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/OUTPUT_SHA256.txt" )
echo "done $(date -u +%FT%TZ)"
