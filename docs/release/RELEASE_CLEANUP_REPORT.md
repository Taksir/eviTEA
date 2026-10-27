# EviTEA release cleanup: report

- Base commit: `815f6c23a822c25bc0312f43ab49769f7aa48bc9` (head of `dev-merged` when the work started; the work order's
  placeholders were not filled in, so the head was used, as instructed). The relative depth floor default was taken from it
  unchanged (off); it was set to on afterwards, before the freeze.
- Branch: `release-cleanup` (the name asked for in the request; the work order says `dev-release-cleanup`).
- Last code commit: `76bcbbf` (the two fixes before the freeze, below); this report is updated in the commit after it.
- Commits are authored as `Md Taksir Hasan Majumder <Taksir@users.noreply.github.com>` with no AI attribution lines, as the work order's rule 6 asks.
- Python 3.14.6, pymsbwt 0.5.0, numpy 2.5.4, mappy 2.31, parasail 1.3.4, edlib 1.3.9, Bowtie2 2.5.4, BLAST+ 2.16.0, Linux.

## How the identity checks were done

The real fly data were not available. Every byte-for-byte check used a synthetic data set made by
`tools/synthetic/make_call_data.py` (five 200 kb fly-named contigs and `chrM`, 120 TE records in 60 families with FlyBase ids,
BED6 annotation of resident copies, 125 bp paired reads at 30x with fragments of 350 +- 40 bp, 30 planted insertions). The
fixtures under `tests/fixtures` are too small for a whole `call` run: the library layer refuses any data set whose rescue-panel
optimiser does not pick the 1,024-row panel (see "Noticed, out of scope"), so the generator adds one high-copy record.

For each check, prepare-reference, prepare-te `--build-languages`, prepare-sample and call were run with the base commit and
with the new code on the same prepared sample (the base commit's prepare-sample refuses pymsbwt, which is item A1; every
msBWT package builds the same index). Compared byte for byte: the two TE language codebooks and their 17 audit tables,
`calls.native.tsv`, `calls.tsv`, `calls.bed`, `unresolved.tsv`, `run_summary.json` and `run.log` (25 files).
`run_manifest.json` is not compared: it records times, paths and now the sample-derived values. The base commit makes 30
calls on this set, one per planted insertion. `tools/synthetic/compare_call.sh` repeats this check for any previous version.

## Part A

| item | done | commit | acceptance checks |
|---|---|---|---|
| A1 pymsbwt | yes | `1083819` | pytest passes; `EVITEA_RUN_MSBWT_LIVE=1` live test passes with pymsbwt 0.5.0 (and now checks the prepared sample names pymsbwt); run manifest `backend_versions.exact_query` = `msbwt pymsbwt 0.5.0` |
| A2 read length | yes | `7453d75` | 125 bp with `--read-length 125`: 25/25 files identical to base. 100 bp and 150 bp synthetic sets: runs complete (30 calls each), manifest `read_length` 100 / 150 (`prepared-sample`). Unit tests of the windows at 100 and 150 bp |
| A3 fragment size | yes | `e748366` | `--fragment-q95 494 --fragment-q99 805`: 25/25 identical. Unit tests: known spans give the expected nearest-rank percentiles, outliers, the fallback with its warning, option precedence, a Bowtie2 estimate on synthetic pairs (same result with 1 and 2 threads). On the synthetic set (fragments 350 +- 40) the estimate is 415 / 442 bp (theory 416 / 443) from 110,984 pairs |
| A4 primary contigs | yes | `2fac50c` | with the five arms passed: identical (combined check below). Unit tests: mouse-like, yeast I..XVI, Arabidopsis-like and dm6 contig sets |
| A5 bare names, mitochondrion | yes | `f7744ea` | dm6-named synthetic set: identical (combined check). Unit tests: `MT` found; a reference without one runs and logs it; bare `1` maps to `chr1` only when `chr1` is in the reference |
| A6 TE record names | yes | `b4f092e` | fly-named set: identical (combined check). New test runs the anchor catalogue and the TE languages on a tiny library named `B1_Mm`, `L1MdA_I`, `Ty1`, ... with a two-column map and gets the expected families. Whole `call` on a synthetic set with plain record names, a two-column map, contigs `chr1`..`chr5` and `MT`, every value from the sample: all 30 insertions called with their family |

Supporting commits: `6076bb3` (a pre-existing test failure, see below), `26ab1ed` (the synthetic generator), `495cc2e`
(`tools/synthetic/compare_call.sh`).

**Combined Part A identity check** (all six items, at `b4f092e`, with `--read-length 125 --fragment-q95 494
--fragment-q99 805 --primary-contigs chr2L,chr2R,chr3L,chr3R,chrX`): 25/25 files identical to the base commit.

**Library-free comparisons** (`tools/library_free_tests/run_all.sh` with `PREVIOUS_SRC` = the base commit's `src`, at
`b4f092e`): depth tool, records (40/40), assembly (35/35), runtime (32/32, U1 compared with the base commit and identical;
U8 skipped, `BASELINE_SRC` unset), repeat map (13/13) and Slurm tools (5/5) pass. `test_depth_window` fails at D5; it fails
in the same way at the base commit itself with `PREVIOUS_SRC` set to the base commit: D5 assumes that `PREVIOUS_SRC` is the
older checkpoint version, whose depth window had a fixed upper bound of 100. Run directly, the depth-window probes (depth 30
and 150, with and without `--depth-cap 100`) give the same places and segments with the base commit and the new code.

### What each item changed (files)

- A1: `pyproject.toml`, `src/evitea/backends/msbwt.py`, `src/evitea/cli.py`, `src/evitea/incumbent/msbwt_queries.py`,
  `README.md`, `tools/slurm/site.example.sh`, tests. pymsbwt >= 0.5 is accepted; msbwt-modern3 0.3.0 still works as a fallback.
- A2: new `src/evitea/sample_geometry.py`; `incumbent/coverage_panel.py`, `pair_aware_panel.py`, `recall_repair.py`,
  `acquisition.py`, `incumbent/engine.py`, `engine.py`, `pipeline.py`, `config.py`, `cli.py`, `callset.py`, `provenance.py`;
  `languages/internal.py`, `languages/builder.py`, `assets/te.py` (`prepare-te --read-length`, default 125, for the block
  column; the column keeps its name). Option `call --read-length`. Manifest `values.run_parameters`.
- A3: `sample_geometry.py` (estimate), `incumbent/mate_rescue.py`, `incumbent/complete_routing.py` (via acquisition),
  `anchor_pe/anchors.py`, `anchor_pe/alignments.py`, `anchor_pe/engine.py`, `full_memory/strict_a.py`, `proposals.py`,
  `full_memory/engine.py`, `internal21/acquisition.py`, `extension_engine.py`, `backends/bowtie2.py`, `constants.py`
  (q99 out of the frozen check), `engine.py`, `pipeline.py`, `config.py`, `cli.py`. Options `--fragment-q95`,
  `--fragment-q99`.
- A4: new `src/evitea/contigs.py`; `full_memory/common.py`, `proposals.py`, `full_memory/engine.py`, `anchor_pe/engine.py`,
  `engine.py`, `pipeline.py`, `config.py`, `cli.py`, `docs/CLI_AND_DATA.md`. Option `--primary-contigs`.
- A5: `contigs.py`, `resident_filter.py`, `library_free/resident.py`, `library_free/runtime.py`, `library_free/_p4_lib.py`
  (docstrings), `cli.py`. Options `library-free --mito-contig`, `filter-resident-copies --reference-assets`.
- A6: new `src/evitea/family_table.py`; `anchor_pe/clusters.py`, `anchor_pe/decoder.py`, `anchor_pe/engine.py`,
  `incumbent/anchor_census.py`, `incumbent/te_library.py`, `languages/common.py`, `inputs.py`, `breakpoints/runtime.py`,
  `tiers.py`, `resident_filter.py`, `family_records/run.py`, `families.py`, `cli.py`.

### Decisions to check

- **Bowtie2 `-X` of the paired-end channel.** The work order names `FROZEN_REFERENCE_MAX_INSERT` in `backends/bowtie2.py`.
  That command is used only by the tests; the engine's paired reference mappings use `REFERENCE_OPTIONS` in
  `anchor_pe/alignments.py`, which also hard-coded `-X 805` (the same 99th percentile). Both now take the per-sample value
  (anchor_pe, full memory, INTERNAL21). Revert `reference_options` if that was not meant.
- **Default primary contigs.** Rule: contigs at least 10% as long as the longest, without organelle, unplaced or
  alternative contigs. On dm6 it gives the five arms plus `chrY` (chr4 is 4% of chr3R), so the fly runs must pass the five
  arms, as the work order says they will.
- **Fragment estimate.** Up to 200,000 evenly spaced pairs mapped end-to-end with `-X 5000`; both mates MAPQ >= 20; spans
  above median + 10 MAD dropped; nearest-rank percentiles; fallback 494 / 805 below 1,000 pairs. These numbers are new and
  are explained in `sample_geometry.py` and `docs/CLI_AND_DATA.md`. The estimate takes one Bowtie2 run (about 2 minutes on
  the synthetic set under load). On the fly samples the frozen run will use the estimate unless `--fragment-q95/--fragment-q99`
  are passed; the work order expects calls to change where the estimate is not 494 / 805.
- **`prepare-te --read-length`.** The internal-21 block column is built at `prepare-te` time, before any sample, so it gets
  its own option (default 125). It is an audit column; no call reads it.

## Part B

All of Part B leaves the calls unchanged (checked below).

| item | done | commit | acceptance |
|---|---|---|---|
| B1 VCF | yes | `6744043` | `tests/unit/test_vcf.py` converts a small call table and library-free calls and checks every record; the TSV and BED bytes are unchanged. `call` and `refine-breakpoints` write `calls.vcf`; `library-free` writes `calls.combined.vcf` |
| B2 input errors | yes | `9874070` | `tests/negative/test_unsupported_input.py`: one negative test per check (paired input, one read length, contig names, TE ids vs map) |
| B3 run command | yes | `b83dedd`, `7a01d88` | `evitea run` (with a library, `--library-free`, or `--no-library`) and `tools/no_library.sh`; `tools/synthetic/compare_run.sh`. With a library and `--library-free`, `evitea run` and the separate commands give identical `call/` (7 files), `library_free/` (8 tables) and TE languages. Without a library, `evitea run --no-library` and `tools/no_library.sh` give identical repeat map, stub, library-free tables (8), records (8 files) and prepared record assets; **both then stop at the record-pass `call` with the same error** ("pair-aware reselection is frozen for the B=1024 rescue panel; got B=256"), because the records made from 30 single-copy synthetic insertions are all low-copy. So the no-library chain could not be run to its end on synthetic data (see "Noticed, out of scope") |
| B4 install | yes | `f50797d` | `evitea check-install`; `environment.yml`, `containers/Containerfile`; README states CPython 3.14 and Linux/WSL for the library-free layer, and the extras. The environment was created with micromamba and works (Python 3.14.9, Bowtie2 2.5.5, BLAST+ 2.17.0, all modules; `check-install` lists all). The container image was not built: no container runtime here |
| B5 CI | yes | `43f2718` | `.github/workflows/tests.yml` (Linux: check-install, pytest, live msBWT, the synthetic end-to-end calls, `run_all.sh`; Windows: pytest). `.github/` removed from `.gitignore` in the same commit. The synthetic end-to-end test (`35e6327` fixed its count) passes in both cases with `EVITEA_RUN_SYNTHETIC_CALL=1` (fly names; plain names with chr1..chr5 and MT): 2 passed in 19 minutes. The workflow itself has not run on GitHub yet |
| B6 example | yes | `702f257` | `examples/synthetic/run_example.sh` and `expected/calls.native.tsv`. Run at `3c4e295`: "calls.native.tsv: as expected" (30 calls, all 30 planted insertions with their family) |
| B7 docs | yes | `0b23c68` | README rewritten; `docs/REPEAT_ANNOTATION.md` (RepeatMasker command and `.out` to BED6, checked on a RepeatMasker-format table; `build-repeat-map` and runs without a library); `docs/README.md`, `docs/CLI_AND_DATA.md`. No accuracy claim for other species |
| B8 notices | yes | `3c4e295` | `THIRD_PARTY_NOTICES`. All licences fit GPL-3.0-or-later; parasail's licence asks for a citation (see the file) |

**Part B leaves the calls unchanged**: `tools/synthetic/compare_call.sh` at `3c4e295` against the base commit (development values passed) gives 25/25 identical files. The follow-up `ddfbd0f` changes only the library-free VCF writer.

One Part B regression was found and fixed: B1 made `library-free` write `calls.combined.vcf` through the `calls.tsv` views, which failed on the six-column library table of `tools/library_free_tests/test_runtime.py` (the failure showed in `run_all.sh` at `3c4e295`). Fixed in `ddfbd0f`; `test_runtime.py` then passes (30 checks).

## Two fixes before the freeze (after the first report)

| fix | commit | checks |
|---|---|---|
| The 1,024-row rescue panel is always kept when it fits the row limit | `68968c6` | `tests/unit/test_rescue_panel.py`: on a library whose words are all rare, the 256-, 512- and 1,024-row panels score the same and the old `max()` would keep 256 rows; the 1,024-row panel is kept and the pair-aware reselection accepts it. A 1,024-row panel over `MATERIALIZE_ROW_HARD_CAP` stops with a plain error |
| A run without a library takes no repeat annotation: the repeat map's BED6 copy stands in | `76bcbbf` | `build-repeat-map` also writes `repeat_map_bed6.bed` (name `repeat`, score 0, strand `.`; `repeat_map.bed` unchanged). `evitea run --no-library` without `--repeat-annotation` builds the map from the FASTA first and prepares the reference with that copy; `tools/no_library.sh` does the same with `REFERENCE=<FASTA>`. Unit test of the step order; `test_repeat_map.py` 15/15 (the copy holds the same intervals and prepare-reference reads it) |

Results on the synthetic data, at `76bcbbf`:

- **Identity:** `tools/synthetic/compare_call.sh` against the base commit: 25/25 files identical.
- **The no-library chain now runs to its end**, with and without an annotation. `tools/synthetic/compare_run.sh`: 69 files
  compared, all identical: `evitea run` against the separate commands with a library (and library-free); `evitea run
  --no-library --repeat-annotation ...` against `tools/no_library.sh` on prepared assets; and `evitea run --no-library` with
  no annotation against `tools/no_library.sh` with `REFERENCE=`. Each no-library chain gives 19 library-free calls, 23 record-
  pass calls and 20 calls in `record_pass/added/calls.combined.tsv` (one record-pass call added).
- **With an explicit annotation the outputs match the earlier ones:** the repeat map, the library-free tables and the records
  of the shell chain, made fresh with the fixed code, equal those `evitea run` made before the fixes (the record pass did not
  finish then, so it has no earlier output).
- **Curated annotation against the repeat map:** on this data set the two no-library chains give byte-identical library-free
  calls, record-pass calls and final calls, as expected: the same-family tests never match self-built records.
- pytest: 226 passed, 0 failed, 3 skipped (the same three opt-in tests).

The first entry of "Noticed, out of scope" (the library layer refusing genomes whose TE words are all rare) is fixed by
`68968c6`; the entry about `--repeat-annotation` being needed by `run --no-library` is fixed by `76bcbbf`.

## Tests

| | passed | failed | skipped |
|---|---|---|---|
| pytest at the base commit | 184 | 1 | 1 |
| pytest at the last commit | 223 | 0 | 3 |

- The failure at the base commit was `tests/unit/test_phase8_checkpoint_pipeline.py::test_cli_raw_build`: the test's
  argument namespace lacked `depth_cap`, which `_call` reads since the depth cap option was added. Fixed in `6076bb3`.
- Skips: `tests/integration/test_phase8_msbwt_live.py` (needs `EVITEA_RUN_MSBWT_LIVE=1`; passes when set) and the two cases
  of `tests/integration/test_synthetic_call.py` (needs `EVITEA_RUN_SYNTHETIC_CALL=1`; passes in both cases with `EVITEA_RUN_SYNTHETIC_CALL=1` (fly names; plain names with chr1..chr5 and MT): 2 passed in 19 minutes).
- `tools/library_free_tests/run_all.sh` without `PREVIOUS_SRC`: all pass at the base commit (U1 and U8 skipped).
  At `3c4e295` (no `PREVIOUS_SRC`): depth window, depth tool, records (40/40), assembly (35/35), repeat map (13/13) and Slurm tools (5/5) pass; `test_runtime` failed (the B1 regression above) and passes after `ddfbd0f` (30 checks, U1 and U8 skipped).

## Not done

- The container image was not built or run (no Docker daemon or Apptainer here); the conda environment it installs was.
- The GitHub Actions workflow has not run yet; it runs on the first push.
- The no-library chain (`evitea run --no-library`, `tools/no_library.sh`) was compared up to the record pass only: on the synthetic data its record-pass `call` stops at the frozen rescue-panel check (both forms identically).
- Nothing was run on real data, as instructed.

## Noticed, out of scope

- **The library layer refuses some genomes outright.** `incumbent/pair_aware_panel.py` stops with "pair-aware reselection
  is frozen for the B=1024 rescue panel" when the rescue-panel optimiser (`coverage_panel.PANEL_CAPS` 256 / 512 / 1024)
  picks a smaller panel. That happens when no TE word occurs more than 512 times in the reads, for example in a small genome,
  at low depth or with a library of low-copy families. It also stopped the record pass of the no-library chain on the
  synthetic data (records made from single-copy insertions). It is a frozen check and a fitted value, so it was not changed,
  but it will matter for yeast, for low-coverage samples and for runs without a library.
- `evidence/frozen_candidates.py:769` uses `start // 125` for `independent_block_id` (the `--g1a-archive` path of
  `prepare-te`, not used by the default engine). Not changed; not in A2's list.
- The frozen G1A crosswalk reader (`evidence/frozen_candidates.py` `_crosswalk`, FlyBase columns) is sha256-pinned to the
  FlyBase file and used only with `--g1a-archive`, so it was left as it is.
- `evitea run --no-library` (like the Slurm chain) needs a BED6 repeat annotation, because the record pass runs `call`,
  which reads it. A species with no TE library needs one from a de novo repeat library; the repeat map (three columns, no
  names) cannot stand in. Documented in `docs/REPEAT_ANNOTATION.md`.
- `filter-resident-copies` now maps any bare copy-table contig name (also `4`, `Y`) to its `chr` contig. If the FlyBase copy
  table that the tool lists uses bare names and has copies on chr4 or chrY, calls there could now be dropped where they were
  not before; I could not check that table. The tool is an ablation post-filter and does not change `call` output.
- `tiers.py` and `breakpoints/runtime.py` now read the family map through the shared reader, which drops rows with an empty
  family (they used to map the record to an empty name). Neither changes the default calls.
- `tools/library_free_tests/test_depth_window.py` D5 cannot pass with `PREVIOUS_SRC` set to a commit after the depth-cap
  change (see above); its comment should say which version `PREVIOUS_SRC` must be.
- `prepare-sample`'s per-read JSON identity map does not scale to full-size libraries (already in the docs).
- The rename-list items (`--rule-c`, `_p4_lib`, "frozen" in docstrings) and the two comments that name other tools
  (`incumbent/complete_routing.py:13`, `resident_filter.py:4`) were left unchanged.
- The work order is still in `docs/release/`; whether it stays is your decision.
