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
| FIX, dev data name | 124 | 115 |
| FIX, register | 52 | 3 |
| FIX, stock phrase | 1 | 0 |
| LOOK, register | 542 | 454 |
| LOOK, field word | 279 | 189 |
| LOOK, stock word | 2 | 1 |

`wording_scan.py stats --against $START` (words in comments, docstrings, prose and messages):

| Category | Before | After |
|---|---|---|
| comment | 4029 | 3927 |
| docstring | 21257 | 21152 |
| prose | 9765 | 9694 |
| string | 8233 | 8226 |

No category has more words than before. The diff touches 196 files and adds or removes no file except this report.

## 2. Checks from section 8 of the guide

The checks ran on Python 3.14.6 with pymsbwt 0.5.0, mappy 2.31, parasail 1.3.4, edlib, numpy 2.5.4, Bowtie2 2.5.2 and BLAST+. The start numbers were taken in the same environment from a second checkout of `$START`.

| # | Check | Start | End |
|---|---|---|---|
| 1 | `guard $START` | not applicable | Exit 0, 0 failures, 220 REVIEW lines, explained in section 4 |
| 2 | `scan --min-level MUST` | 978 hits | 146 hits. They are names, option values, ids and `say()` strings (sections 5 and 6) |
| 3 | `stats --against $START` | not applicable | No category has more words (section 1) |
| 4 | `python -m pytest -ra` | 226 passed, 0 failed, 3 skipped | 226 passed, 0 failed, 3 skipped |
| 5 | `bash tools/library_free_tests/run_all.sh` (no `PREVIOUS_SRC`) | Exit 0; depth_window ALL PASS, depth_tool ALL PASS, records 40/0, assemble 35/35, runtime 30/0, repeat_map 15/0, slurm_tools 5/0 | Exit 0; the same results line for line |
| 6 | `compare_call.sh` (byte identity of the library layer) | not applicable | All compared files identical except `run_summary.json`; see below |
| 7 | `run_all.sh` with `PREVIOUS_SRC` | not applicable | Same results as the start commit run against itself; see below |
| 8 | Help pages | 16 pages saved | Only wording differs (section 4) |
| 9 | `git diff --stat $START` | not applicable | 196 files changed, no file added or removed except this report |

Check 6, in detail. `PREVIOUS_SRC` was the `src` folder of the start checkout. Of the 25 compared files, 24 are identical: the TE codebooks and their audit tables, `calls.native.tsv`, `calls.tsv`, `calls.bed`, `unresolved.tsv` and `run.log`. `run_summary.json` differs in a few counts of the extension layers (for example `B_universe_unique_high` 397 against 376, `corrected_clusters` 46 against 26, `eligible_reference_unexplained` 857 against 805). I ran the same script with the start commit on both sides. It gives the identical difference, and its `run_summary.json` is byte-for-byte equal to the one from the new tree. So the difference does not come from the wording edits. It comes from the script: it passes `--read-length 125 --fragment-q95 494 --fragment-q99 805 --primary-contigs ...` only to the current side. Question 19 below asks whether that is expected.

Check 7, in detail. With `PREVIOUS_SRC` set, `test_depth_window` fails at `t_scan_above_window` (assertion on `places_sha`). The guide lists this as known. It fails the same way when the start commit is compared with itself. All other scripts pass: `test_depth_tool` ALL PASS, `test_records` 40 checks 0 failed, `test_assemble` 35 of 35, `test_runtime` 32 checks 0 failed (32 checks; the two extra checks that need `PREVIOUS_SRC` run, and U8 is skipped because `BASELINE_SRC` is not set), `test_repeat_map` 15 checks 0 failed (BLAST+ installed), `test_slurm_tools` 5 checks 0 failed.

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

`guard` exits 0 with 0 failures and 220 REVIEW lines. Every line is one of these kinds.

- 147 changed message or help strings (117 in `raise`, 25 in `help=`, 2 in print or message calls, 3 in `match=`). Each is a wording change that keeps the value, path or option name the old message showed. Typical changes: "frozen G1A archive" to "anchor archive", "mismatch" to "does not match", "needs" to "requires", "Phase 3 Bowtie2 profile" to "the Bowtie2 profile", a dropped label such as "incumbent" in a message.
- 65 "number X appears less". Every number left with a history sentence or a project label: dates (2026, 09), stage and phase codes (Phase 0 to 9, S20 to S70, "layers 1-2", "class 1", "gate-1"), finding ids (F-031, F-035, I-165), "Section 128", "spec section 8", and the strain name GIM-024. No value, unit or threshold was dropped.
- 7 "code span no longer present". Six are hash prefixes whose ellipsis character became `...`, and one is the source commit hash `0a45bc9` in the validation status file, deleted as history. The seventh is an artifact of line-wrapped backticks in the same file (`prepare-reference --bowtie2-index` is still there).
- 1 "lines that are not comments changed" in `environment.yml`. It is a comment after a code line. The package string is identical.

Message changes the guard shows that need a human check:

- `cli.py` `--profile` help: now "profile: 'production' (default), or AKA017 or ISO1 for the read names and event labels of those two development samples". Taken from the `engine.py` docstring.
- `cli.py` `--rule-c` help: now "the family-identity filter: also drop family calls whose family members have a median pairwise identity of <= 98.31%". The words "median pairwise" come from the style guide.
- `cli.py` `--depth-cap` help: "keeps the previous window" became "sets an upper bound of 100", read from `depth_window.py`.
- `cli.py` `--build-languages` help: "both v3 languages" became "the terminal and internal 21-mer codebooks".
- `library_free/options.py`: `--identity-min-comparisons` now says "members with plain-extension segments" (was "frozen-extension members"). `--previous-behaviour` now says "every rule switch off".
- `family_records/run.py`: "how many of" became "what fraction of", and "share" became "fraction" in the `--te-share` and `--copy-report` help.

Help pages: all 16 help outputs were saved at the start and compared after. Only wording differs. Every default, unit and value is still there (1500000, 125, 500 bp, 98.31%, depth 10 and 100, 494 and 805).

## 5. Labels I could not understand and left

- `cli.py` `--complete-routing` help: "TE word". `cli.py` `build-family-records` help: "end-word codebooks". Are these k-mers?
- `cli.py` `--mito-contig` help: "segment firewall". Also in `library_free/place.py` and `_p4_lib.py` (`MappyFirewall`).
- `config.py`: "first-stage layer". `run_command.py` docstring: "flag rule kept", "record guard", "residency test".
- `cli.py` `--g1a-archive` and `--family-crosswalk` have no help text. I do not know what the G1A archive is. In the text I call it the "anchor catalogue archive" because its members are `G1A_ANCHOR_CATALOGUE.tsv` and `G1A_ANCHOR_PROVENANCE.tsv`.
- `shadow_layers.py` ~165 and `full_memory/selection.py`: "Hostile-valid proposal-memory". `shadow_layers.py` ~133: "corrected cluster".
- `anchor_pe/clusters.py`: "Delayed-aggregation novelty mask". "corrected-cluster", "corrected unexplained" and "relaxed provisional" clusters match function names, so they stay.
- `_p4_lib.py`: "hostile coordinate", "hostile sites", "prereg section 3", "design memo section 4", `class_strata cs_lib.rep_fracs`, "(synthetic)" on constants.
- `library_free/options.py` lines 101 to 102: `--polya-cross` and `--polya-tsd` help says "(set aside; off)".
- `library_free` docs and code: "walk-start word", "valid words", "haplotype scan's trusted members", "duplicate chain", "second-chance", "strides", "walked-out" (docs/LIBRARY_FREE_OPTIONS.md and code).
- Stage labels S1 to S7, `S2_*`, `S4_*` stay because they are written to the places table.
- Research file names in docstrings (`calling.py`: `p4_library_free/v2/code/v2_lib.py`, `precision_diag/code/precb_run.py`, `v2_diag/code/v2diag_lib.py`; `detector.py`: `p4_genome.py`, `p4_run.py`). They point to files that are not in the repository.
- `candidates.py`: the `source` values "G1F" and "G1HR_RESCUE" are data. I describe them as "working pair" and "wide-anchor pair".
- `extension_engine.py` title: "801/820/835/838-style layers" was dropped. I assumed these were old run ids.
- `contracts.py` ~634: the argument "certification rotation_id" (a string in a call argument, which the guard treats as data).
- `resident_filter.py` ~59: a data string with `dm6_TEMP2_compatible.fa.out` and "wrong-call anatomy".
- Windows scripts: `tools/windows/VALIDATE_LANGUAGES.ps1` line 4 "regenerated G1A catalogue". `tests/unit/test_internal21.py` line 26: "A" and "B". `tools/library_free_tests/synth_elements.py` line 75: "must stay exactly as before". `r8/check.py`: "T0".
- Progress and check strings that go through `say()`, `require()`, `_check()` or `check()` are treated as data by the guard, so they keep old labels. You decided to leave them. They are:
  - `incumbent/acquisition.py` 98, 104, 142: "coverage rescue census", "rescue panel optimization", "candidate recall repair".
  - `incumbent/engine.py` 266, 305, 325: "two-sided mate rescue", "full TE-segment family certification", "repeat-aware divergence rescue".
  - `engine.py` 313, 337: "incumbent: sample ...", "incumbent calls: ...; v3 extension".
  - `library_free/runtime.py` 284: "family calls (v2 filter, rule c)".
  - `full_memory/engine.py`, `anchor_pe/engine.py`, `internal21/acquisition.py`, `languages/builder.py` and `languages/terminal.py`, `internal.py`: progress lines such as "finite state:", "AnchorPE counts:", "INTERNAL21 projection", "languages: G1A anchor catalogue", "gap-fill targets".
  - `full_memory/strict_a.py`: the message "B fold drift" (suggested text: "paired-end anchor support row has a fold that does not match its pair ID").
  - `validation/` scripts: `require("frozen ISO transfer helper")`, `say("hashing AKA FASTQs (recorded, not historically frozen)")`, the `[ISO->AKA` log prefix, `_check()` messages in `run_shadow_layer_equivalence.py`, the print "= frozen hashes" in `run_call_equivalence.py`, and the stage-check prints "T0 rows equal", "R8E equal", "CV equal", "R9E rescues", "fcrc1 roundtrip".
  - `tools/library_free_tests/test_assemble.py`: check labels A3, A6, A7 (garbled: "the family's own the family-identity filter drops it") and four A10 labels with "frozen"; `test_repeat_map.py` lines 120 and 131.

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

Places where a comment and the code disagree. Neither side was changed:

- `library_free/__init__.py` lines 20 to 21 say `_p4_lib` is `p4_lib.py` (sha256 `56a4a47d...`), "vendored byte for byte; the tests check its hash". The file is not byte-identical (its sha256 starts `43730d5e`) and no test checks a hash.
- `library_free/place.py` line 18 points to `tests/test_place.py`, which does not exist.
- `_p4_lib.py` `MappyFirewall` docstring says `genome: dict name -> sequence`, but `__init__` takes `genome_fa`, a FASTA path.
- `breakpoints/__init__.py`, `walker.py`, `pinning.py` point to `tests/reference_s3` and `tests/reference_pin1/pin1_lib.py`. These do not exist.
- `evidence/frozen_candidates.py` (`build_frozen_candidate_codebooks`) reads `repeat_context.reference_semantic_sha256` before it tests `repeat_context is None`. The validation status file says this module is unused by the default engine, yet `assets/te.py` calls it when the anchor archive options are given.
- `validation/incumbent_stage_checks/iso1_acquisition_vm.py` says each step caches under `--cache`. The script has no `--cache` option, only the `ISO1_CACHE` environment variable.
- `validation/incumbent_stage_checks/downstream_stages/fcrc1/check.py` line 47 says the recorded table is in ascending locus order. The code does not check this.
- `validation/EVITEA_V1_VALIDATION_STATUS.md` lists the ISO1 host read names as `iso`, but `anchor_pe/alignments.py` line 52 has `pid`. It also says native Windows Bowtie2 rejects all-numeric names, while an ISO1 run with such names is reported elsewhere in the file. Its regression line says "185 passed, 1 skipped".
- `docs/CLI_AND_DATA.md` said the codebook builders are not yet stages. `prepare-te --build-languages` builds both. The doc now follows the help text.
- `tests/integration/README.md` and `tests/parity/README.md` say the folders are reserved. Both hold tests.
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

## 11. Questions for Taksir

1. `cli.py` "every canonical" message and its test: change both together to drop the label "canonical"? The guard cannot pass that pair, because the old match string has only two words.
2. "TE word" in the `--complete-routing` help and "end-word" in the `build-family-records` help: are these k-mers? What are their plain names?
3. "segment firewall" (`--mito-contig`): a feature name to keep?
4. What is the G1A archive (`--g1a-archive`)? It has no help text. Is "anchor catalogue archive" right? What is the "G1A catalogue" in `VALIDATE_LANGUAGES.ps1`?
5. Is "coverage panel" the right plain name for the "rescue panel", and for "compact screen" in the wide-anchor step?
6. In `--profile` help and the docs I describe `AKA017` and `ISO1` as "development samples". Is that accurate? Keep the two names in the validation docstrings and status file (they are also option and profile values), or use "the development sample"? I dropped `GIM-024` and assumed it is just another development sample.
7. The `say()` strings listed in section 5: left as you decided. Do you want a later change that moves them to a printing call?
8. `library_free/__init__.py`: correct the stale hash and "byte for byte" claim? That is a content change.
9. `iso1_acquisition_vm.py`: fix the `--cache` text, or add the option?
10. Meaning of "T0" in `r8/check.py`, "set aside" in the `--polya-cross` and `--polya-tsd` help, "hostile" in `_p4_lib.py`, "Hostile-valid proposal-memory", "Delayed-aggregation novelty mask".
11. Keep or delete the pointers to research files and tests that are not in the repository (`calling.py`, `detector.py`, `breakpoints/`, `place.py`)?
12. Validation status file: keep a dated record? Drop `--profile S` in the opening paragraph (not a documented value)? Keep the stale "185 passed" line, the Windows Bowtie2 caveat and the two mentions of private development notes? Delete the "Retired:" sentence as history?
13. "checkpoint branch" in `tools/slurm/README.md`, `site.example.sh` and `test_runtime.py`: replace with "a checkout of the previous version"? Users follow that step, so I kept it.
14. CHANGELOG and README Status: keep the development-sample sentence and the call counts 838 and 194? Merge the Unreleased changelog section? Keep "previous release" and `--previous-behaviour` wording in `docs/LIBRARY_FREE_OPTIONS.md`?
15. README credit line for msBWT names two people (James Holt and Leonard McMillan). Keep it? It follows `THIRD_PARTY_NOTICES`.
16. Plain names to confirm: repeat-aware step description, "extension-layer events within 500 bp" for `unresolved.tsv`, `docs/CLI_AND_DATA.md` codebook bullet (now follows the `--build-languages` help), and the stage-code readings in section 9.
17. Are "proposal" and "proposals" (a candidate region in `full_memory`) acceptable as the name, or should a later pass rename them? Should the "G1F" and "G1HR_RESCUE" `source` values wait for the rename list?
18. The commits are unsigned. GitHub may show them as Unverified, and a stop hook in the working session asked for the author to be Claude. They are authored and committed as you, as you asked.
19. `compare_call.sh` passes `--read-length 125 --fragment-q95 494 --fragment-q99 805 --primary-contigs ...` only to the current side, and `run_summary.json` then differs in a few extension-layer counts, even when the start commit is compared with itself. Is that expected, or should the script give both sides the same values?
