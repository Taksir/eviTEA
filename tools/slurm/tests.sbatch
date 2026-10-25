#!/bin/bash
# The public tests, in the container. Requires EVITEA_WORK (and EVITEA_OUT, EVITEA_CODE if not at their defaults) exported.
# Submit from $EVITEA_OUT (mkdir -p logs first):  sbatch tools/slurm/tests.sbatch
#SBATCH --job-name=evitea_tests
#SBATCH --partition=small
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
. "$EVITEA_CODE/tools/slurm/env.sh"
export APPTAINERENV_PREVIOUS_SRC=$PREVIOUS_SRC   # enables the byte-for-byte comparison with the previous version
cd "$EVITEA_CODE"
for t in test_depth_window test_depth_tool test_records test_assemble test_runtime test_repeat_map test_slurm_tools; do
  echo "== $t"; $PY tools/library_free_tests/$t.py | tail -3 || echo "FAILED $t"
done
echo "done $(date -u +%FT%TZ)"
