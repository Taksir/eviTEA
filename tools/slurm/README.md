# Slurm run files

Job scripts for running the tests and the library-free layer on a cluster with Slurm and Apptainer. They write only under one output folder and read every location from one variable and one site file.

## Setup

1. `export EVITEA_WORK=<your working folder>`. All cluster paths are built from it.
2. Clone this repository into a new folder (default `$EVITEA_WORK/evitea_code`; set `EVITEA_CODE` for another) and check out the branch to run. Check out the repository's checkpoint branch into a second folder: its `src` is the previous version used for the comparisons below.
3. Create the output folder `$EVITEA_WORK/evitea_runs` (`EVITEA_OUT`) with a `logs` folder in it, copy `site.example.sh` to `$EVITEA_OUT/site.sh` and fill it in (container image, reference, prepared reference assets, path of the previous version's `src`, path of the sample table).
4. Write the sample table: one line per sample, tab-separated: name, MSBWT index, folder of the library run (`evitea call` output), default contigs. The two paths are relative to `EVITEA_WORK`.
5. Submit from `$EVITEA_OUT` with `sbatch $S/<script>`, with `S=$EVITEA_CODE/tools/slurm`. Variables go in the environment.

Jobs refuse to run on a login node.

## Scripts

| Script | What it does |
|---|---|
| `tests.sbatch` | the public tests of `tools/library_free_tests`, in the container (2 CPUs, 8 GB, 2 h) |
| `detect.sbatch` | one `evitea library-free` run (2 CPUs, 16 GB, up to 12 h); variables are listed at the top of the script |
| `records.sbatch` | `build-family-records` then `prepare-te` for one sample (8 CPUs, 96 GB, 3 h) |
| `env.sh` | shared settings, sourced by the others |

## Order and checks

1. **Tests.** `sbatch $S/tests.sbatch`. Every test file must end with `ALL PASS` or `N checks, 0 failed`; the log ends with `done`.
2. **Equality with the previous version.** On one sample, with the same contigs (a few are enough), run
   `DS=<sample> TAG=previous PREVIOUS=1 sbatch $S/detect.sbatch` and
   `DS=<sample> TAG=same DEPTHCAP=100 OPTS="--previous-behaviour" sbatch $S/detect.sbatch`.
   `library_free_places.tsv`, `library_free_segments.fa`, `library_free_family_calls.tsv`, `calls.library_free.tsv` and `calls.combined.tsv` of the two runs must be identical (`python3 tools/table_diff.py A B`; exit 0 means the same). If they differ, stop and report the first differing rows. The previous version's upper depth bound is 100, which `DEPTHCAP=100` sets.
3. **Detections.** `DS=<sample> TAG=default sbatch $S/detect.sbatch` (all switches at their defaults) for each sample. In each `library_free_summary.json` read `resident_guard`: mode `family` with a `te_share` near the genome's TE share (about 0.16 for the fly reference). The mode `off` means no repeat annotation was found; pass one with `--repeat-annotation` in `OPTS`.
4. **Depth window.** Run a sample with `DEPTHCAP=100` and without, and list the places whose outcome changes with `tools/depth_cap_changes.py --old A/library_free_places.tsv --new B/library_free_places.tsv`.
5. **Single switches.** Each run with one switch off shows what that part contributes. `OPTS` examples: `--het-scan off`; `--anchored-extension off --walk-out off`; `--walk-out off`; `--check-scope none --resident-guard off`.
6. **Regrouping from saved pieces.** To change only grouping or assembly switches, reuse the scan pieces of a finished run: `WORK=$EVITEA_OUT/work/<DS>_default REGROUP=1 TAG=<label> OPTS="..." DS=<sample> sbatch $S/detect.sbatch`. With `REGROUP=1` the job stops instead of scanning again if the switches, the code or the inputs differ from the saved pieces. Switches to vary this way: `--check-scope path-dropped`, `--tol p95`, `--tol-source own-calls`, `--resident-guard per-call`, `--resident-guard off`, `--first-pass-priority off`, `--het-yield off`, `--polya-tsd on`. `RUNDIR` points a regrouped run at another library run.
7. **Family records.** `DS=<sample> SRC=$EVITEA_OUT/runs/<run folder> sbatch $S/records.sbatch`, with `RM=<RepeatMasker table>` for the record guard. Then run the library layer once more on the new assets (`evitea call` with the sample's own command line, changing only `--te-assets`, `--output` and `--work-dir`) and add its calls with `python3 -m evitea add-record-calls --calls <run>/calls.native.tsv --base-calls <run folder>/library_free/calls.combined.tsv --records <records folder>/lib --output <out>`.

## Outputs

`runs/<sample>_<label>_<job>/library_free/` holds the places tables, the segments, the family calls, `calls.library_free.tsv`, `calls.combined.tsv` and `library_free_summary.json`; `OUTPUT_SHA256.txt` lists the checksum of every file. Records go to `records/<sample>_<job>/`. A job never overwrites an existing output folder.
