# Wording cleanup report

This report covers the wording rewrite of the comments, docstrings, messages, help text, README and docs. Only wording changed. The start commit is `0b485b9c50c1a85562c31b30cd6da07c8e648b7c`.

## 1. Scan counts, before and after

`wording_scan.py scan --summary`:

| Group | Before | After |
|---|---|---|
| Files scanned | 255 | 255 |
| Hits, all levels | 2246 | 1041 |
| MUST, internal label | 948 | 142 |
| MUST, history | 16 | 0 |
| MUST, claim word | 7 | 0 |
| MUST, tool name | 6 | 1 |
| MUST, person name | 1 | 3 |
| FIX, internal label | 142 | 56 |
| FIX, history | 126 | 77 |
| FIX, dev data name | 124 | 114 |
| FIX, register | 52 | 3 |
| FIX, stock phrase | 1 | 0 |
| LOOK, register | 542 | 452 |
| LOOK, field word | 279 | 192 |
| LOOK, stock word | 2 | 1 |

`wording_scan.py stats --against $START` (words in comments, docstrings, prose and messages):

| Category | Before | After |
|---|---|---|
| comment | 4029 | 3975 |
| docstring | 21257 | 21167 |
| prose | 9765 | 9684 |
| string | 8233 | 8232 |

No category has more words than before. The prose count leaves out this report, which the scanner also counts (3877 words, and its labels add scan hits); the scan counts above leave it out too. The diff touches 198 files and adds no file except this report.

## 2. Checks from section 8 of the guide

The checks ran on Python 3.14.6 with pymsbwt 0.5.0, mappy 2.31, parasail 1.3.4, edlib, numpy 2.5.4, Bowtie2 2.5.2 and BLAST+. The start numbers were taken in the same environment from a second checkout of `$START`.

| # | Check | Start | End |
|---|---|---|---|
| 1 | `guard $START` | not applicable | Exit 1, 1 failure: `evidence/frozen_candidates.py`, where Taksir asked for a code change (answer B). Without that file, exit 0. 201 REVIEW lines, explained in section 4 |
| 2 | `scan --min-level MUST` | 978 hits | 146 hits. They are names, option values, ids and `say()` strings (sections 5 and 6) |
| 3 | `stats --against $START` | not applicable | No category has more words (section 1) |
| 4 | `python -m pytest -ra` | 226 passed, 0 failed, 3 skipped | 226 passed, 0 failed, 3 skipped |
| 5 | `bash tools/library_free_tests/run_all.sh` (no `PREVIOUS_SRC`) | Exit 0; depth_window ALL PASS, depth_tool ALL PASS, records 40/0, assemble 35/35, runtime 30/0, repeat_map 15/0, slurm_tools 5/0 | Exit 0; the same results line for line, run again after the last edits |
| 6 | `compare_call.sh` (byte identity of the library layer) | not applicable | 25 of 25 identical with `CURRENT_OPTS=""`; see below |
| 7 | `run_all.sh` with `PREVIOUS_SRC` | not applicable | Same results as the start commit run against itself; see below |
| 8 | Help pages | 16 pages saved | Only wording differs (section 4) |
| 9 | `git diff --stat $START` | not applicable | 198 files changed, no file added or removed except this report |

Check 6, in detail. `PREVIOUS_SRC` was the `src` folder of the start checkout. `compare_call.sh` now takes the extra current-side options from `CURRENT_OPTS`; its default keeps the old fixed values (`--read-length 125 --fragment-q95 494 --fragment-q99 805 --primary-contigs ...`). With `CURRENT_OPTS=""` against the start commit, all 25 compared files are identical: the TE codebooks and their audit tables, `calls.native.tsv`, `calls.tsv`, `calls.bed`, `unresolved.tsv`, `run_summary.json` and `run.log`. With the default `CURRENT_OPTS` the first pass gave 24 of 25: `run_summary.json` differed in a few extension-layer counts. The same difference appeared with the start commit on both sides, so the edits did not cause it. The old fixed values are given to the current side only, so the two sides use different read-length and fragment values.

Check 7, in detail. With `PREVIOUS_SRC` set, `test_depth_window` fails at `t_scan_above_window` (assertion on `places_sha`). The guide lists this as known. It fails the same way when the start commit is compared with itself. Taksir's answer: the check expects the previous version to abstain at 150x with its default cap, so `PREVIOUS_SRC` must be a version from before the depth cap was lifted. That is a limit of the test, and the docstring and a comment now say so. All other scripts pass: `test_depth_tool` ALL PASS, `test_records` 40 checks 0 failed, `test_assemble` 35 of 35, `test_runtime` 32 checks 0 failed (32 checks; the two extra checks that need `PREVIOUS_SRC` run, and U8 is skipped because `BASELINE_SRC` is not set), `test_repeat_map` 15 checks 0 failed (BLAST+ installed), `test_slurm_tools` 5 checks 0 failed.

Environment note: before Bowtie2 was installed, the first pytest run gave 225 passed and 4 skipped. The start numbers above were taken again after the install. The scan counts in section 1 were taken before this report was added. The report itself quotes labels that the scanner flags.

## 3. Messages that tests match

Match strings changed, each in the same commit as its message:

| Test file | Old `match=` | New `match=` |
|---|---|---|
| `tests/unit/test_rescue_panel.py` | `1024-row rescue panel needs .* materialization limit` | `1024-row coverage panel needs .* materialization limit` |
| `tests/unit/test_s20_frozen_candidate_codebooks.py` | `G1A archive SHA-256 mismatch` | `archive SHA-256 does` |
| `tests/unit/test_s20_frozen_candidate_codebooks.py` | `crosswalk SHA-256 mismatch` | `crosswalk SHA-256 does` |

The new messages read "anchor archive SHA-256 does not match the expected value" and "FBte family crosswalk SHA-256 does not match the expected value". The guard needs at least three words in a match string, so the two new strings end in "does".

Messages that tests match and that keep their matched words:

- `cli.py`: "run manifest must bind every canonical non-manifest output". The test matches `every canonical`. The word "canonical" stays, because the guard rejects a change of a two-word match string.
- `preflight.py`: "same file", "do not pair up", "more than one length (100, 120 bp)". Only "needs" became "requires".
- `pipeline.py`, `provenance.py`: "cannot be overwritten", "explicit resume", "nonempty", "recursively", "together", checkpoint identity and output hash mismatch.
- `contigs.py`: "not a contig", "not in the reference". `family_table.py`: "neither a FlyBase table", "two tab-separated fields", "is empty".
- `engine.py`: "reproduces only sample", "unknown run profile". `sample_geometry.py`: "does not hold reads of one length". `constants.py`: "scientific".
- `evidence/codebook.py`: "multiple TE families". `evidence/frozen_candidates.py`: "lacks catalogue/provenance". `backends/msbwt.py`: "every physical pair", "uniform read length".
- `assets/te.py`: "--terminal-codebook", "together". `assets/reference.py`: "--repeat-annotation".
- `library_free/runtime.py`: the `EVITEA_PIECES_READONLY` message, which `tools/library_free_tests/test_runtime.py` (U3) checks.
- `library_free/assemble.py`: the value "off: no repeat annotation" is data. `test_assemble.py` (A8) checks it.

## 4. Guard REVIEW lines

`guard` exits 1 with 1 failure, and 201 REVIEW lines. The failure is `evidence/frozen_candidates.py`: the two checks `repeat_context is None` and the reference hash now run in the other order, as Taksir asked (answer B). That is a code change, in its own commit. Because the guard stops at the code change, it does not list the message edits of that file. Before the commit, the guard passed for the file with all its message edits. Every REVIEW line is one of these kinds.

- 110 changed message or help strings (78 in `raise`, 27 in `help=`, 2 in print or message calls, 3 in `match=`). Each is a wording change that keeps the value, path or option name the old message showed. Typical changes: "frozen G1A archive" to "anchor archive", "mismatch" to "does not match", "needs" to "requires", "Phase 3 Bowtie2 profile" to "the Bowtie2 profile", a dropped label such as "incumbent" in a message.
- 73 "number X appears less". Every number left with a history sentence or a project label: dates (2026, 09), stage and phase codes (Phase 0 to 9, S20 to S70, "layers 1-2", "class 1", "gate-1"), finding ids (F-031, F-035, I-165), "Section 128", "spec section 8", and the strain name GIM-024. Three are in `README.md`: 125, 194 and 838 left with the development-sample sentence that Taksir asked me to drop from the Status section (answer 14). The counts stay in the validation status file.
- 16 "code span no longer present". Six are hash prefixes whose ellipsis character became `...` (three in `tools/windows/README.md`, three in the validation status file; one of them, `0a45bc9`, is a source commit hash deleted as history). Six are the names of removed modules (`integration.py`, `proposal.py`, `exact_junction.py`, `paired_end_anchor.py`, `certification.py`, `internal_probe.py`) in the "Retired" bullet of the validation status file, which answer 12 told me to delete. `git_dirty` and `evitea call --msbwt-index ... --profile S` went with the sentences that answer 12 removed. One is an artifact of line-wrapped backticks (`prepare-reference --bowtie2-index` is still there). One is the combined command `library-free --previous-behaviour --depth-cap 100` in `docs/LIBRARY_FREE_OPTIONS.md`, which the new text writes as two separate spans.
- 2 "lines that are not comments changed". One is in `environment.yml`: a comment after a code line, with the package string identical. The other is in `tools/synthetic/compare_call.sh`: the `CURRENT_OPTS` change (answer 19), made in its own commit.

Message changes the guard shows that need a human check:

- `cli.py` `--profile` help: now "profile: 'production' (default), or AKA017 or ISO1 for the read names and event labels of those two development samples". Taken from the `engine.py` docstring.
- `cli.py` `--rule-c` help: now "the family-identity filter: also drop family calls whose family members have a median pairwise identity of <= 98.31%". The words "median pairwise" come from the style guide.
- `cli.py` `--depth-cap` help: "keeps the previous window" became "sets an upper bound of 100", read from `depth_window.py`.
- `cli.py` `--build-languages` help: "both v3 languages" became "the terminal and internal 21-mer codebooks".
- `library_free/options.py`: `--identity-min-comparisons` now says "members with plain-extension segments" (was "frozen-extension members"). `--previous-behaviour` now says "every rule switch off".
- `family_records/run.py`: "how many of" became "what fraction of", and "share" became "fraction" in the `--te-share` and `--copy-report` help.

Help pages: all 16 help outputs were saved at the start and compared after. Only wording differs. Every default, unit and value is still there (1500000, 125, 500 bp, 98.31%, depth 10 and 100, 494 and 805).

## 5. Labels I could not understand and left

Taksir's answers (section 11) explained most labels of the first pass, and the text now uses the plain names. These are still left.

- `run_command.py` docstring: "record guard" and "residency test". `shadow_layers.py`, `anchor_pe` and `full_memory`: "relaxed provisional" clusters (they match function names).
- `_p4_lib.py`: "(synthetic)" on constants and `class_strata cs_lib.rep_fracs`. The `library_free` docs and code: "walk-start word", "valid words", "haplotype scan's trusted members", "duplicate chain", "second-chance", "strides", "walked-out".
- Stage labels S1 to S7, `S2_*`, `S4_*` stay because they are written to the places table.
- Research file names in docstrings (`calling.py`: `p4_library_free/v2/code/v2_lib.py`, `precision_diag/code/precb_run.py`, `v2_diag/code/v2diag_lib.py`; `detector.py`: `p4_genome.py`, `p4_run.py`). Answer 11 says to delete such pointers where they only confuse readers. These name the source of the method, so I left them for the next pass.
- `candidates.py`: the `source` values "G1F" and "G1HR_RESCUE" are data. They wait for the rename table.
- `extension_engine.py` title: "801/820/835/838-style layers" was dropped. `contracts.py` ~634: the argument "certification rotation_id". `resident_filter.py` ~59: a data string with `dm6_TEMP2_compatible.fa.out` and "wrong-call anatomy".
- `tests/unit/test_internal21.py` line 26: "A" and "B". `tools/library_free_tests/synth_elements.py` line 75: "must stay exactly as before". "Section 128" in the validation status file ("recorded research tables").
- Strings the guard treats as data, so they keep old labels. You decided to leave the progress strings (answer 7 is marked DECIDE and is not done yet):
  - `incumbent/acquisition.py` 98, 104, 142: "coverage rescue census", "rescue panel optimization", "candidate recall repair".
  - `incumbent/engine.py` 266, 305, 325: "two-sided mate rescue", "full TE-segment family certification", "repeat-aware divergence rescue". `engine.py` 313, 337: "incumbent: sample ...", "incumbent calls: ...; v3 extension".
  - `library_free/runtime.py` 284: "family calls (v2 filter, rule c)".
  - `full_memory/engine.py`, `anchor_pe/engine.py`, `internal21/acquisition.py`, `languages/builder.py`, `terminal.py`, `internal.py`: progress lines such as "finite state:", "AnchorPE counts:", "INTERNAL21 projection", "languages: G1A anchor catalogue", "gap-fill targets".
  - `full_memory/strict_a.py`: the message "B fold drift" (suggested text: "paired-end anchor support row has a fold that does not match its pair ID").
  - `library_free/options.py` lines 101 to 102: the `--polya-cross` and `--polya-tsd` help sits in a table of tuples, which the guard treats as data. It still says "(set aside; off)". Answer 10 asks for "(off by default; experimental)". That edit fails the guard.
  - `--g1a-archive` and `--family-crosswalk` have no help text. Adding `help=` changes the syntax tree, which the guard rejects. Answer 4 suggests a text.
  - `validation/` scripts: `require("frozen ISO transfer helper")`, `say("hashing AKA FASTQs (recorded, not historically frozen)")`, the `[ISO->AKA` log prefix, `_check()` messages in `run_shadow_layer_equivalence.py`, the print "= frozen hashes" in `run_call_equivalence.py`, and the stage-check prints "T0 rows equal", "R8E equal", "CV equal", "R9E rescues", "fcrc1 roundtrip". You said to leave these.
  - `tools/library_free_tests/test_assemble.py`: check labels A3, A6, A7 (garbled: "the family's own the family-identity filter drops it") and four A10 labels with "frozen"; `test_repeat_map.py` lines 120 and 131.
  - `test_depth_window.py`: the D5 result string still says "identical to the previous version". Only the docstring and a comment now state that `PREVIOUS_SRC` must be a version from before the depth cap was lifted.

## 6. Names I left

Remaining MUST-level scan hits (142) are all names or ids that may not change, or `say()` strings from section 5. The test-case ids A1 to A10 and D1 to D9 in `test_assemble.py` and `test_depth_window.py`, the tier numbers in `tiers.py` (real tier names), the window step 250, the probe names L1, L2, R1, R2 and B1, B2 in tests, and `incumbent.MsbwtSample` module references are not labels to rewrite.

By kind:

- **Option names and values:** `--rule-c`, `--max-rescue-interval-mass`, `--g1a-archive`, `--family-crosswalk`, `--build-languages`, `--profile` (values `production`, `AKA017`, `ISO1`), `--complete-routing`, `--pin-one-sided`, `--previous-behaviour`, `--te-share`, `--depth-cap`, `--route-b-v2`, `--route-b-v3`, `--incumbent-aka017`, `--languages built|supplied`, `--host-read-names iso|aka-historical`, `pid`, `-RouteBV2`, `-RouteBV3`.
- **Constants, classes, functions, variables:** `RULE_C_MAX`, `T_E`, `rule_c_dropped`, `v2_dropped`, `feat_c`, `e_value`, `MAX_RESCUE_INTERVAL_MASS`, `RESCUE_*`, `FROZEN_RESCUE_CAP`, `FROZEN_ALIASES`, `FROZEN_FAMILY_ALIASES`, `FROZEN_LANGUAGES`, `FROZEN_INPUT_SHA`, `FROZEN_AKA`, `FAMILY_IDENTITY_FILTER`, `LEGACY_EVENT_FIELDS`, `IncumbentCall*`, `certify_*`, `CertifiedEvent`, `GroupCertification`, `MateRescueLocus`, `Shadow*`, `AnchorPE*`, `compose_frozen_shadow`, `build_frozen_candidate_codebooks`, `G1A_ARCHIVE_SHA256`, `ledger_row`, `optimize_rescue_panel`, `strict_a`, `channel_a`, `channel_b_universe`.
- **Output and data values:** statuses `DROPPED_V2`, `DROPPED_RULE_C`, `DROPPED_CHECK`, `DROPPED_RESIDENT_COPY`, `NEAR_EVITEA`, `NEAR_FIRST_PASS`, `NOT_HET`, `ADDED`; call-id prefixes `V2FCRC1_EVENT_`, `V2R4_MATE_`, `V2R9E_REPEATDIV_`, `V3FULLMEM_AP_`, `V3ANCHORPE_`, `V3INTA_AP_`; call classes `MATE_RESCUE`, `REPEAT_AWARE_RESCUE`; tier names `REPEAT_RESCUE`, `RESIDENT`; layer names `INTERNAL21`; evidence types `msbwt_r9e_repeat_aware_crossvalidated_te_divergence`; report verdict strings (`*_PASS`, `ISO_PROCEDURE_*`, `FROZEN_ORACLE_INTEGRITY_ONLY`); JSON and table keys and columns (`v2_dropped`, `rule_c_similarity_max`, `frozen_final_sha256`, `A_pairs`, `B_state`, `proposal_only_A` and similar); schemas `evitea.incumbent.checkpoint.v1`, `evitea.prepared-te.v3`, `evitea.library_free.piece.v2`; environment variables `EVITEA_ROUTE_B_V2`, `EVITEA_ROUTE_B_V3`, `EVITEA_DEPTH_CAP`, `EVITEA_PIECES_READONLY`; stage values S1 to S7.
- **File, folder and module names:** `incumbent/`, `recall_repair.py`, `mate_rescue.py`, `te_certification.py`, `shadow_layers.py`, `frozen_regression/`, `v3_layer_checks/`, `oracles.json`, `route_b_v3_*.py`, `fbte_to_flybase_family_FROZEN.tsv`, `G1A_*` files, `P4_PLACES.tsv`, `ISO1_STAGE_CHECKS_VM_20260930.json`, `g1*`/`t_g1*` scripts, `test_phase*_*.py`, `test_s20_frozen_candidate_codebooks.py`, `tools/windows/*.ps1`.
- **Competitor names in paths or code:** `dm6_TEMP2_compatible.fa`, `TEMP2_ISO1`. `tests/negative/test_production_import_boundary.py` line 15 holds `"temp2"` and `"teforest"` as strings in a forbidden-module check. The scanner's `tests/conftest.py` hit on "temp" is a false positive.
- **Person names:** the author field in `pyproject.toml`, and the dependency credit in the README (msBWT by James Holt and Leonard McMillan).
- **Strain names (115 FIX hits):** `AKA017` and `ISO1` stay where they are option values, profile values, paths or file names, in the validation scripts and records, and in the Windows launchers. `GIM-024` was removed from docstrings.
- **SHA-256 prefixes** in docstrings (such as `6432c8e1`) are flagged by the scanner as commit hashes.

## 7. TODO items and comment/code disagreements

No `TODO` or `FIXME` was found.

Places where a comment and the code disagreed were fixed by the answers of Taksir (stale hash claim in `library_free/__init__.py`, the `tests/test_place.py` and `tests/reference_*` pointers, the `MappyFirewall` docstring, the `--cache` text, the `fcrc1/check.py` order comment, the reserved-folder READMEs, and the order of the `None` check in `frozen_candidates.py`). These remain:

- `evidence/frozen_candidates.py`: the validation status file says this module is unused by the default engine, yet `assets/te.py` calls it when the anchor archive options are given.
- `validation/EVITEA_V1_VALIDATION_STATUS.md` lists the ISO1 host read names as `iso`, but `anchor_pe/alignments.py` line 52 has `pid`.
- `docs/CLI_AND_DATA.md` said the codebook builders are not yet stages. `prepare-te --build-languages` builds both. The doc now follows the help text.
- `constants.py`: the old comment said q99 "is not frozen". The code does not check q99, so the comment now says "is not checked".
- `full_memory/__init__.py` called the layer the "third" extension layer and `internal21` the "fourth". `shadow_layers.py` lists four layers. I dropped the ordinal words.
- `library_free/runtime.py` hashes eight edited source files (`_code_sha256`) for the cache identity in `.complete.json`. Comment edits change that hash, so pieces saved by the old code are scanned again. No output byte depends on it.

## 8. Claim words changed in the README and docs

- `README.md`: "finds new transposable-element (TE) insertions" became "detects non-reference transposable-element (TE) insertions". "30 new insertions" became "30 non-reference TE insertions" (also in the `examples/synthetic/run_example.sh` comment). "finds every planted insertion" became "identifies". "its new calls" became "the calls of that run".
- `README.md`: removed the promise "it will be reported separately" and the history clause "made during the release cleanup".
- `README.md` License section: added the credit line "EviTEA builds on msBWT (James Holt and Leonard McMillan), Bowtie2, BLAST+, mappy (minimap2), parasail, edlib and NumPy".
- Code docstrings: "new TE insertions" became "non-reference TE insertions" in `tools/synthetic/make_call_data.py` (its `--help` text). "new insertions" became "non-reference insertions" in `family_records/select.py` and `resident.py`. "beats ... by >= 4 edits" became "has >= 4 fewer edits than" in `split_junction.py` and `sequence.py`.
- The TEMP2 and TEforest sentence in `incumbent/complete_routing.py` became "This is the usual discordant-pair entry rule." "TEMP2 dm6 index" became "dm6 index" in the validation status file.
- No "first library-free", "only tool", "outperforms" or "superior" was found in the README or docs.

## 9. Plain names I chose

These plain names come from reading the code. Please check them.

- "rescue panel" is now "coverage panel" (matches `coverage_panel.py`). "compact screen" in the wide-anchor step docstring is now "the coverage panel".
- Validation stage codes: G1A to G1HR are the anchor census, sparse panel, coverage panel, pair-aware panel, host-mate filter with resident pair families, and the wide-anchor step. H1 is Hamming-1 loci. FCRC1 is split junction. R4 is mate-anchored search. R6 is host-seed routing. R7 is junction groups. R8C and R8D are TE-segment family confirmation and trim-variant runs. R9E is the repeat-aware library layer.
- `errors.py` "frozen architecture" is now "a rule of the program design". `config.py` messages now read "exact_query_backend accepts only msbwt" and "scientific setting is fixed and not configurable".
- The repeat-aware step is described as "a repeat-aware step that retains a call whose differences from the TE record recur at the locus" (from `repeat_rescue.py`).
- `unresolved.tsv` is described as holding competing extension events within 500 bp.
- "valid-TSD" means a family with at least two two-sided places and at least two places with a 2 to 20 letter TSD (`_p4_lib.families`).
- The tools and tests files use "k-mer" where the old text said "word" for a 21-mer or 30-mer. "end-word step" in `test_records.py` became "codebook-building step".
- The `engine.py` profile table heading "wide-anchor names" names the read names the old text called "recall-repair".
- Other wording edits: `tools/windows` README and scripts changed "frozen" to "fixed", "Route-B artifacts" to "reference artifacts", "R9E incumbent" to "library layer", "Verdict" to "Result", "oracles" to "expected stage outputs". The validation status file lost its dated session headings and the commit hash.
- Docstrings are a little shorter in places. The trim pass dropped a few asides (`sparse_panel.py` "(terminal, internal or fallback bits)", `mate_rescue.py` "Tie-free" and "locus-record quantization", and "exact junction" from the one-line docstring of `internal21/observability.py`) to keep the word counts under the start. Two sentences in `docs/CLI_AND_DATA.md` lost detail ("within 500 bp of each other" in the `unresolved.tsv` row).

## 10. CHANGELOG

All facts are kept; only the wording changed. The `--previous-behaviour` entry now says "turns off every added rule". The Unreleased section reads as a development log. Merging it into one 0.1 entry is your decision.

## 11. Answers applied, and what is still open

Taksir answered the 19 questions of the first report. This section lists what was done for each answer. Items marked DECIDE in the answers were not started.

Done:

1. "every canonical" message: kept as is.
2. `--complete-routing` help now describes the "TE word" as a TE anchor k-mer, and `build-family-records` help says "TE-end and TE-interior k-mer codebooks". Docstrings in `family_records/` and `docs/REPEAT_ANNOTATION.md` say the same. "TE word" in `coverage_panel.py` and the README became "TE anchor k-mer" and "TE k-mers".
3. "segment firewall" is now "segment filter" in the `--mito-contig` help, in `place.py`, `extension.py`, `_p4_lib.py` and `tools/library_free_tests/synth.py`. The class `MappyFirewall` and variables named `firewall` keep their names.
4. Not done (DECIDE): hiding `--g1a-archive` and `--family-crosswalk`. The help text could not be added (see section 5).
5. Confirmed: "coverage panel". No change needed.
6. Checked: no text calls GIM-024 or any other strain development data. GIM-024 is left only in two lists of forbidden names in `validation/` scripts.
7. Not done (DECIDE): the `say()` pass. The old and new lines are in the answers.
8. `library_free/__init__.py` now says "adapted from the development script `p4_lib.py`", without the hash or "byte for byte".
9. `iso1_acquisition_vm.py` now names the `ISO1_CACHE` environment variable.
10. Done: "untrimmed segment rows" (T0), "unexplained-pair cluster", the full-memory event docstring, the extra-clusters docstring, "family flag rule 'kept'", "a read-routing option of the library layer", "hostile" in `_p4_lib.py`. Not done: the `--polya-*` help (the guard blocks it).
11. Pointers to files that are not in the repository were removed in `breakpoints/` and `place.py`. The research file names in `calling.py` and `detector.py` are still there.
12. Validation status file: removed `--profile S`, the "Retired" bullet, the Windows Bowtie2 caveat and both mentions of private development notes. The regression line now says 226 passed, 3 skipped. The numbers are marked as results recorded at commit d7587b1. Not done (DECIDE): whether the file ships.
13. "checkpoint branch" is now "the version to compare against" in `tools/slurm/README.md`, `site.example.sh` and `test_runtime.py`.
14. README Status no longer gives the development-sample sentence or the counts 838 and 194 (they stay in the validation file). `docs/LIBRARY_FREE_OPTIONS.md` describes `--previous-behaviour` as turning off every added rule and giving the basic library-free scan, without "previous release". Not done (DECIDE): merging the Unreleased CHANGELOG section. The CHANGELOG still has the sentence with 838 and 194, and the README section "What is validated" still states the byte-for-byte equivalence. Please decide for both. Add `--previous-behaviour` to the rename list.
15. README credit line kept.
16. The repeat-aware step now has the longer description in the README. The `unresolved.tsv` row in `docs/CLI_AND_DATA.md` has the new text with "within 500 bp".
17. Prose says "candidate region" only where the old text had "proposal" in a sentence I touched. `full_memory` code names stay. The output-name renames (DECIDE) wait for the rename table. Candidates for the table: call-id prefixes, statuses such as `DROPPED_V2`, call classes `MATE_RESCUE` and `REPEAT_AWARE_RESCUE`, tier `REPEAT_RESCUE`, source values `G1F` and `G1HR_RESCUE` (suggested `WORKING_PAIR` and `WIDE_ANCHOR`), `--rule-c`, `--max-rescue-interval-mass`, `--previous-behaviour`.
18. The commits stay authored and committed as Md Taksir Hasan Majumder.
19. `compare_call.sh` now takes the extra current-side options from `CURRENT_OPTS`, with today's values as the default. The run with `CURRENT_OPTS=""` against the start commit is in section 2.
- `test_depth_window.py`: the docstring and a comment state that `PREVIOUS_SRC` must be a version from before the depth cap was lifted. The assertion is unchanged.
- A: the 2 unreported trim commits were not used.
- B: `evidence/frozen_candidates.py` now tests `repeat_context is None` before it reads the context. This is a code change in its own commit. The guard lists it as a failure for that file, as it should.
- C: pointers to tests that do not exist were removed from `place.py`, `breakpoints/__init__.py`, `walker.py` and `pinning.py`.
- D: the `MappyFirewall` docstring says `genome_fa` is a FASTA path.
- E: the `fcrc1/check.py` comment says the order of loci in the recorded table is not checked.
- F: `tests/integration/README.md` and `tests/parity/README.md` say what the tests in each folder check.

Still open for Taksir:

1. The DECIDE items: hide `--g1a-archive` and `--family-crosswalk`; the `say()` pass; whether the validation status file ships; the CHANGELOG merge; the output-name renames.
2. Allow an exception in the guard (or a code change) for the `--polya-*` help table and for new `help=` text on `--g1a-archive` and `--family-crosswalk`? Both are wording, but the guard cannot tell.
3. The README section "What is validated" and the CHANGELOG still state the byte-for-byte equivalence and the counts 838 and 194. Drop or keep?
4. "record guard", "residency test", "Section 128", "relaxed provisional": plain names?
5. Delete the research file names in `calling.py` and `detector.py`?
6. The commits are unsigned. GitHub may show them as Unverified.
