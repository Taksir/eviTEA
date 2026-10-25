#!/bin/bash
# All public tests. Needs a python with mappy, parasail and numpy as $PYTHON.
# PREVIOUS_SRC (the src folder of the previous version) enables the byte-for-byte comparisons; they are skipped without it.
# The development comparisons with reference implementations (place, het scan, priority, check) are not part of this repository.
cd "$(dirname "$0")/../.." || exit 1
PYTHON=${PYTHON:-python3}
export PYTHONPATH=src
fail=0
for t in test_depth_window test_depth_tool test_records test_assemble test_runtime test_repeat_map test_slurm_tools; do
  echo "== $t"
  $PYTHON -B tools/library_free_tests/$t.py 2>&1 | grep -v "^PASS" | tail -5 || true
  [ "${PIPESTATUS[0]}" = 0 ] || fail=1
done
exit $fail
