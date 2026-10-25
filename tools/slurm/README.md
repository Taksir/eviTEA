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
| `tests.sbatch` | the public tests of `tools/library_free_tests`, in the container (2 CPUs, 8 GB, 2 h); the repeat-map build checks need BLAST+ (`BLAST_BIN` in `site.sh`) and are skipped without it |
| `detect.sbatch` | one `evitea library-free` run (2 CPUs, 16 GB, up to 12 h); variables are listed at the top of the script |
| `records.sbatch` | `build-family-records` then `prepare-te` for one sample (8 CPUs, 96 GB, 3 h) |
| `repeat_map.sbatch` | `build-repeat-map`: the repeat map of the reference for runs without a library (8 CPUs, 8 GB, 3 h) |
| `record_pass.sbatch` | the library layer once more on the family records (the sample's own `evitea call` command line, three options changed by `call_args.py`), then `add-record-calls` (8 CPUs, 96 GB, 12 h) |
| `no_library.sh` | run on a login node: submits repeat map (if missing), detection, records and record pass as one dependency chain |
| `env.sh` | shared settings, sourced by the others |

## Order and checks

1. **Tests.** `sbatch $S/tests.sbatch`. Every test file must end with `ALL PASS` or `N checks, 0 failed`; the log ends with `done`.
2. **Equality with the previous version.** On one sample, with the same contigs (a few are enough), run
   `DS=<sample> TAG=previous PREVIOUS=1 sbatch $S/detect.sbatch` and
   `DS=<sample> TAG=same DEPTHCAP=100 OPTS="--previous-behaviour" sbatch $S/detect.sbatch`.
   `library_free_places.tsv`, `library_free_segments.fa`, `library_free_family_calls.tsv`, `calls.library_free.tsv` and `calls.combined.tsv` of the two runs must be identical (`python3 tools/table_diff.py A B`; exit 0 means the same). If they differ, stop and report the first differing rows. The previous version's upper depth bound is 100, which `DEPTHCAP=100` sets.
   **Against saved tables.** If you hold the tables of an earlier run of the previous version on a sample, two more runs must reproduce them:
   - Same sample and the same contigs as the saved run (its job log lists them): `DS=<sample> TAG=firstpass DEPTHCAP=100 OPTS="--het-scan off --check-scope none --resident-guard off" sbatch $S/detect.sbatch`. `library_free_places.tsv` and `library_free_segments.fa` must be identical to the saved detection, and `library_free_family_calls.tsv` to the saved family calls.
   - A sample with a heterozygous mix, three contigs, with only the second scan on: `DS=<mix> CONTIGS=chr3L,chr3R,chrX TAG=hetonly DEPTHCAP=100 OPTS="--anchored-extension off --walk-out off --first-pass-priority off --second-scan-anchored off --het-yield off --check-scope path-dropped --resident-guard off" sbatch $S/detect.sbatch`. Its `library_free/library_free_second_places.tsv` must equal the saved second-scan places: `python3 tools/table_diff.py A B --ignore chkL_T,chkL_junction,chkR_T,chkR_junction`. Keep this run's `library_free_family_calls.tsv`: it is the only real-data test of the one-site check.
   If either differs, stop and report the first differing rows.
3. **Detections.** `DS=<sample> TAG=default sbatch $S/detect.sbatch` (all switches at their defaults) for each sample. In each `library_free_summary.json` read `resident_guard`: mode `family` with a `te_share` near the genome's TE share (about 0.16 for the fly reference). The mode `off` means no repeat annotation was found; pass one with `--repeat-annotation` in `OPTS`.
4. **Depth window.** Run a sample with `DEPTHCAP=100` and without, and list the places whose outcome changes with `tools/depth_cap_changes.py --old A/library_free_places.tsv --new B/library_free_places.tsv`.
5. **Single switches.** Each run with one switch off shows what that part contributes. `OPTS` examples: `--het-scan off`; `--anchored-extension off --walk-out off`; `--walk-out off`; `--check-scope none --resident-guard off`.
6. **Regrouping from saved pieces.** To change only grouping or assembly switches, reuse the scan pieces of a finished run: `WORK=$EVITEA_OUT/work/<DS>_default REGROUP=1 TAG=<label> OPTS="..." DS=<sample> sbatch $S/detect.sbatch`. With `REGROUP=1` the job stops instead of scanning again if the switches, the code or the inputs differ from the saved pieces. Switches to vary this way: `--check-scope path-dropped`, `--tol p95`, `--tol-source own-calls`, `--resident-guard per-call`, `--resident-guard off`, `--first-pass-priority off`, `--het-yield off`, `--polya-tsd on`. `RUNDIR` points a regrouped run at another library run.
7. **Family records.** `DS=<sample> SRC=$EVITEA_OUT/runs/<run folder> sbatch $S/records.sbatch`, with `RM=<RepeatMasker table>` for the record guard. Then run the library layer once more on the new assets (`evitea call` with the sample's own command line, changing only `--te-assets`, `--output` and `--work-dir`) and add its calls with `python3 -m evitea add-record-calls --calls <run>/calls.native.tsv --base-calls <run folder>/library_free/calls.combined.tsv --records <records folder>/lib --output <out>`.
   **Reproduction of the record builder.** Run it on the saved first-pass tables of a sample (`SRC=<saved run> FLAG=both FC="<the three saved family-call tables>" SEGFILE=<the saved segments file>`, no `RM`): `te_library.fa` and `family_map.tsv` must have the checksums of the earlier build. Then run it once more with `RM=<RepeatMasker table>`: `records_summary.json` lists `skipped_guard` (families that get no record) with `te_share` and `threshold`.

## Runs without a library

`DS=<sample> bash $EVITEA_CODE/tools/slurm/no_library.sh`, run from `$EVITEA_OUT` on a login node (it only submits). It needs BLAST+ in the image or `BLAST_BIN` in `site.sh`. Jobs:

1. `repeat_map.sbatch`, only if `$MAPDIR/repeat_map.bed` (default `$EVITEA_OUT/repeat_map`) does not exist. The map is per reference; it aligns the sample's default contigs (`CONTIGS` overrides; without a sample, the contigs of 1 Mb or more).
2. `detect.sbatch` with `NOLIB=1 REPEAT_MAP=<bed>`: the library layer's calls are an empty table (a stub of the sample's library run under `$EVITEA_OUT/nolib/<DS>_stub`: the header line of `calls.native.tsv` and the sample id), the family-identity filter is on, all other switches at their defaults (99th-percentile tolerance from the layer's own two-sided calls, one-site check, resident-copy test with the map).
3. `records.sbatch` with `FLAG=kept REPEAT_MAP=<bed>` on that run (consensus records of the surviving families; record guard with the map), then `prepare-te`.
4. `record_pass.sbatch`: `evitea call` with the command line recorded in the sample's library run (`run_manifest.json`; only `--te-assets`, `--output` and `--work-dir` change), then `add-record-calls` against the detection's `library_free/calls.combined.tsv`. The result is `record_pass/<DS>_<job>/added/calls.combined.tsv`.

The driver prints the job ids and output folders. To run a step alone, set the same variables by hand (`NOLIB=1 REPEAT_MAP=<bed> TAG=nolib DS=<sample> sbatch $S/detect.sbatch`, and so on).

## Outputs

`runs/<sample>_<label>_<job>/library_free/` holds the places tables, the segments, the family calls, `calls.library_free.tsv`, `calls.combined.tsv` and `library_free_summary.json`; `OUTPUT_SHA256.txt` lists the checksum of every file. Records go to `records/<sample>_<job>/`. A job never overwrites an existing output folder.
