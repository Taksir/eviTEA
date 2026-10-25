# Verifier for the recorded call files

This validation-only utility checks the recorded call files listed in
[`oracles.json`](oracles.json). It checks the required snapshots in order:
library layer, bilateral, asymmetric, full-memory, and final internal 21-mer.
Within each layer it checks AKA017 before ISO1. The first mismatch ends the run.

The verifier reads only the explicitly named files from the supplied research
roots. In particular, the AKA017 library-layer file is under a directory whose
name mentions truth scoring, but no other file in that directory is opened.
This utility does not read truth or other-tool data, perform scoring, or import
package code.

Run from the repository root:

```text
python validation/frozen_regression/verify_oracles.py --route-b-v2 PATH_TO_ROUTE_B_V2 --route-b-v3 PATH_TO_ROUTE_B_V3
```

Use `--output REPORT.json` to save the same deterministic JSON report that is
printed to standard output. A fully matching set returns exit code 0. A
missing, non-regular, unreadable, or mismatching file returns a nonzero exit
code and reports the earliest failed check. The report explicitly marks method
equivalence as `NOT_EVALUATED`: integrity of the recorded files is necessary
input evidence, not proof that EviTEA reproduced those calls.
