# Wording style guide for the EviTEA code

This guide is for one job: rewrite the wording of the code and docs in plain, human language. Nothing else changes.

It is a working file. Remove `docs/release/` before the public release.

## 1. The one hard rule

**Wording only.** Change the words in comments, docstrings, error and warning messages, command-line help, and the prose of the README and docs. Change nothing else.

Do not change any of these, even when the name looks bad:

- Logic, defaults, constants, numbers, units, types, imports, order of statements.
- Names of options and switches (`--rule-c`, `--max-rescue-interval-mass`, `--profile`, and the rest).
- Names of files, folders, modules, classes, functions, variables and constants (`incumbent/`, `recall_repair.py`, `RULE_C_MAX`, `T_E`, `feat_c`).
- Text that the program writes to a file or reads back: column headers, status values (`DROPPED_V2`), event-id prefixes, manifest keys, VCF fields, the values you type for an option (`production`, `AKA017`).
- Fenced commands and code examples in the docs. Commands, flags and paths inside backticks.
- Directive comments: `# noqa`, `# type:`, `# fmt:`, `# pragma`, the shebang line, encoding lines.
- Test names, test fixtures, and the expected output files in `examples/` and `tests/fixtures/`.

The renames in `RELEASE_RENAMES` wait for Taksir's go. So the names below stay exactly as they are. Next to each one, describe the thing in plain words.

| Keep this name | Describe it as |
|---|---|
| `--rule-c`, `RULE_C_MAX`, `rule_c_dropped`, `feat_c` | the family-identity filter (drops a family whose members have median pairwise identity of 0.98309 or less) |
| `T_E`, `e_value` | the reference-path threshold (0.0405) and the reference-path fraction (the fraction of reads at the breakpoint that continue along the reference) |
| `v2_dropped`, status `DROPPED_V2` | the reference-path filter output, a status that marks a call as excluded |
| `V2FCRC1_EVENT_`, `V2R4_MATE`, `V2R9E_REPEATDIV`, `V3FULLMEM_AP_` | call-id prefixes in the combined call table |
| `--max-rescue-interval-mass` | the wide-anchor work limit |

If a label is in a name you may not change, put the name in backticks and explain it in plain words once. Do not repeat the label as if it were a word.

## 2. Before you start

1. Record the start commit: `START=$(git rev-parse HEAD)`. Every check below compares with it.
2. Run the baseline and keep the numbers:
   - `python -m pytest -ra` (the previous report gives 223 passed, 0 failed, 3 skipped; trust your own run).
   - `bash tools/library_free_tests/run_all.sh`.
3. Make a second checkout of `$START` in another folder. Its `src` folder is `PREVIOUS_SRC` for the identity checks in section 8.
4. Save every help page before you touch anything, so you can read the before and after:
   `for c in "" prepare-reference prepare-te prepare-sample call refine-breakpoints library-free build-family-records add-record-calls build-repeat-map call-tiers filter-resident-copies run check-install; do echo "== $c"; python -m evitea $c --help; done > /tmp/help_before.txt`
   (If `python -m evitea --help` lists another subcommand, add it to the loop.)
5. Run `python docs/release/wording_scan.py scan --summary > /tmp/scan_before.txt` and `python docs/release/wording_scan.py stats > /tmp/stats_before.txt`.

## 3. What to rewrite, and in what order

Do the user-facing text first. Commit after each area.

1. `src/evitea/cli.py`, `config.py`, `errors.py`, `run_command.py`, `preflight.py`, `inputs.py`: help text and error messages.
2. `README.md` and `docs/*.md` (not `docs/release/`).
3. The rest of `src/evitea/`: module docstrings first, then function docstrings, then comments.
4. `tools/`, `examples/`, `containers/`, `.github/`, `environment.yml`, `pyproject.toml`: comments and echo text only.
5. `tests/`: comments and docstrings only. Test and fixture names stay.
6. `validation/`: last. Same rules.

Skip: `LICENSE`, `THIRD_PARTY_NOTICES`, `CITATION.cff` fields, the mapping table in `docs/NAMING.md`, `docs/release/`, and any generated or fixture data.

## 4. How to write it

### Sentences

- Short, plain sentences. Aim for under 20 words. Do not go over 30 outside a table.
- One idea per sentence. Split at a semicolon. Do not nest a parenthesis inside a parenthesis.
- Use active verbs and everyday words. "Reads the file", not "performs a read of the file".
- Name each thing by what it does. If you cannot, the label is wrong; see section 5.
- Use one word for one thing in a file. Do not switch between "call", "hit" and "event" for the same object.
- Plain ASCII punctuation. No em or en dashes, no arrows (write "to" or "maps to"), no curly quotes, no emoji or check marks.
- Spell numbers and units as the code does. US spelling.

### Comments

Write a comment only when it adds something the code does not say: why a value was chosen (with its number and unit), a condition that is easy to miss, or what a block is for.

- Match the density of the neighbors. A file with few comments stays that way. Never add comments to bare code.
- Delete a comment that only repeats the next line.
- Delete a clause that tells the history: what an earlier version did, why something was changed, who decided, when, which test run it came from.
- A comment says what the code does now. Never write "now", "no longer", "fixed", "previously" or "new".
- Keep every number, unit and stated reason. If a comment and the code disagree, change neither. Write it in the report.
- Do not write a new reason you have not read in the code. If you cannot tell why a value is set, say nothing more than the old comment did.

### Docstrings

- First line: one sentence, starting with a verb ("Return", "Raise", "Drop", "Read"). Keep the Python verbs Return, Raise, Yield.
- Keep what a reader needs: units, the shape of the result, the failure cases, the reason for a fixed value.
- Do not turn a one-line docstring into a paragraph. A docstring may be shorter than before. It may not be longer without a reason.
- Some tools use their module docstring as the help text (`description=__doc__`). Read the rendered help after you edit one.

### Error and warning messages

Say what is wrong, what was expected, and what to do. At most two short sentences. Keep the offending value, path or option name that the old message showed.

| Before | After |
|---|---|
| `Phase 3 Bowtie2 profile requires max_reported=2` | `the Bowtie2 profile requires max_reported=2` |
| `incumbent fragment pair IDs must be sorted and unique` | `fragment pair IDs must be sorted and unique` |
| `frozen G1A archive SHA-256 mismatch` | name the archive by what it holds, and say the checksum does not match (read the code to find out what the archive is) |

Messages are text that other code can match.

- Before you change a message, search the whole repository for it: `grep -rn "part of the message" src tests tools validation examples .github`.
- Tests match messages with `pytest.raises(..., match="...")`. If your new message still contains the matched words, nothing else is needed. If not, change only the `match=` string in the same commit and list it in the report.
- A message that is written into an output file or compared with `==` or `in` is data, not wording. Leave it.

### Command-line help

- Follow the style of the neighboring options. In `cli.py` it is a lowercase phrase without a final period.
- Keep every default, unit and value that the old help gave.
- Say what the option does. Do not give the project history of the option.
- A choice value stays as typed. Explain it: `read-name profile: 'production' (default), or AKA017 or ISO1 for the read names of those two samples`.

### README and docs

- Keep the structure and every command. Change the prose.
- Say what the tool does and how to run it. No history, no development-sample stories, no promises.
- Credit the tools EviTEA depends on (Bowtie2, msBWT, mappy, parasail, edlib, numpy, BLAST+, RepeatMasker). Keep those names and their credits.
- Describe an input format by its name ("RepeatMasker `.out`"), not by the tool that made your test file.

## 5. Words to remove, and what to write instead

Run the scanner (section 7). Read every hit in context. A hit is a question, not a verdict. A word inside a name you may not change stays.

### 5.1 Internal labels and process words

These are project shorthand. A reader outside the project cannot know them. Write what the thing does. `docs/NAMING.md` already maps many of them; use its descriptive names.

| Label in the text | Write instead |
|---|---|
| frozen (a default value) | default |
| frozen (a list or table built into the program) | built-in |
| frozen (behavior unchanged from an earlier version) | delete the word |
| frozen (a file the run only reads, with a recorded hash) | fixed input file, or "hash-recorded file" |
| frozen (a check that refuses other values) | "accepts only this value", and give the value |
| incumbent, R9E | the library layer, or the layer's real name: exact junction evidence, paired-end anchor support, repeat-aware calls |
| strict A, exact A, channel A | exact junction evidence |
| AnchorPE, channel B, B-state | paired-end anchor support |
| INTERNAL21, TERMINAL21 | internal 21-mer, terminal 21-mer |
| shadow layers, v3 extension | the extension layers, described by what each adds |
| Route-B, v2, v3 (as an architecture or schema label) | delete, or name the format by its columns |
| rule c, detector v2 | the family-identity filter, the reference-path filter |
| recall repair, recall-repair step | the wide-anchor step |
| rescue budget, rescue limit | the wide-anchor work limit |
| mate rescue | mate-anchored search |
| repeat-aware rescue | the step that retains a call after credit for recurrent edits (read the code, then name it) |
| certified, certification | confirmed, confirmation |
| oracle, verdict, headline, ledger | delete, or say "expected value", "result", "main", "record" |
| adjudication | name what is decided: "long-read validation", "choice between candidate calls" (read the code first) |
| fails closed, fail-closed | "the call is not reported", or "the run stops with an error" |
| production (profile) | default (the option value `production` stays) |
| gate (noun, in prose) | filter, threshold, condition |
| phase 3, step 2b, part B, round 2, S100, G1A, H1 | describe the step by what it does |
| work order, handoff, dev-merged, a commit hash, a branch name | delete |

Terms that are fine and should stay as they are: library layer, library-free layer, full-memory layer, resident-copy filter, family-identity filter, reference-path filter, codebook, confidence tier, panel, one-sided completion, mate-anchored search, wide-anchor step, held-out, k-mer, breakpoint, target-site duplication (TSD).

**If you cannot tell what a label means from the code, do not guess.** Leave it, and list it in the report with the file and line. A wrong plain name is worse than a project label.

Apply this test to every term: *would a reader who has never seen this project understand it on the first read?* Write as if the development strains and TE families were unknown. State the rule first. Give a development number only where the old text did.

### 5.2 Stock wording

Delete the filler word or use the plain one.

| Do not write | Write |
|---|---|
| robust, seamless, comprehensive, sophisticated, powerful, elegant, holistic | say what it does, with a number if there is one |
| leverage, utilize | use |
| ensure, ensuring | "so that", or say which check enforces it |
| in order to | to |
| serves as, acts as, is used to, is responsible for | is, or start with the verb |
| it is worth noting, note that, importantly, notably, essentially, basically, simply, just | delete |
| gracefully, under the hood, out of the box, battle-tested | delete, or say what happens |
| delve, intricate, pivotal, meticulous, foster, bolster, realm, interplay, nuance, groundbreaking, underscore, showcase | delete, or use a plain verb |
| sanity check | check |
| a wide range of, a variety of | the number, or "several" |

### 5.3 History and log wording

A comment, message or doc is not a log. Delete the clause.

- Words: previously, originally, formerly, used to, no longer, earlier or older version, was changed, superseded, carried over, workaround, hack, bug fix, "we decided".
- Dates and clock times of work (`frozen 2026-09-20`), person names, "after the fix", "this passed", "verified on".
- Check on every "now", "new", "still", "already", "currently", "finally": remove it if the sentence reads the same without it.
- `TODO` and `FIXME`: do not delete them silently. Reword to a plain sentence and list each in the report.

Before: `"""Full-data exact-A + bilateral-B AnchorPE events (frozen 2026-09-20)."""`
After: name the events by what they are, drop the date and "frozen". Read the function first.

### 5.4 Names of other tools

No competitor tool is named in code, comments, messages, help text or docs. Name the mechanism instead.

The scanner knows these names: TEMP2, TEforest, ELITE, INSurVeyor, TranSurVeyor, GenomeDelta, MELT, xTea, McClintock, Mobster, RetroSeq, TEPID, TEFLoN, ngs_te_mapper, PoPoolationTE, TEPeak, DeviaTE, Jitterbug, STEAK, Pantera, dnaPipeTE, RepeatModeler. Other tools may exist too. Use your reading.

Before: `This is the discordant-pair entry rule that TEMP2 and TEforest use.`
After: `This is the usual discordant-pair entry rule.`

Before: `This is the mechanism by which TEMP2 avoids such calls.`
After: delete the sentence. The sentence before it already says what the filter does.

Rules for the edge cases:

- A tool name inside a file name, a folder name or a data path (`dm6_TEMP2_compatible.fa`, `TEMP2_ISO1/`) is a path. Leave it. List it in the report.
- A dependency stays by name. A credit in `README.md` or `THIRD_PARTY_NOTICES` stays (msBWT by Holt and McMillan).
- Citations belong in the paper, not in the code.

### 5.5 Development data names

AKA017, ISO1, GIM-024, SLA-001, MUN-009, JUT-011 and the other strain names are development or test data. Keep one only when it is an option value (`--profile AKA017`) or a path. Elsewhere write "the development sample" or describe the property that matters (125 bp reads, dm6).

### 5.6 Claim words

These words make a claim. Do not add them. If one already stands in the README or docs, replace it with the field's term and tell Taksir in the report. Do not decide a claim alone.

- "new insertions", "true insertions": write "non-reference TE insertions".
- "first library-free", "only tool", "beats", "outperforms", "superior": delete. Give the measured number if one is in the text.
- "rescue" as a word for what a step does: use the step's plain name (table in 5.1).

### 5.7 Field words (optional)

In docs and docstrings that describe results, the field's usual word reads better. Apply these only when the sentence stays clear. Never change a name.

| Instead of | Write |
|---|---|
| recovers, finds (a tool's result) | detects, identifies |
| matched, unmatched call | true positive, false positive (only where truth is compared) |
| settings | parameters |
| cutoff | threshold |
| share of | fraction of |
| supplied | provided |
| rejects, keeps (a filter) | excludes, retains |
| needs (software) | requires |
| wrapper | script |
| by hand | manually |

## 6. A short check for each file

After you rewrite a file, read it once more as a stranger.

1. Is every sentence short and does it say what the code does?
2. Does the comment density match the file as it was?
3. Is any label left that only the project knows?
4. Is any number, unit or reason missing compared with the old text?
5. Does the help or message still show the same default, path or value?

## 7. The scanner

`docs/release/wording_scan.py` uses only the Python standard library. Run it from the repository root. It never reads identifiers, only comments, docstrings, message text and prose.

```
python docs/release/wording_scan.py scan --summary              # counts per group and per file
python docs/release/wording_scan.py scan src/evitea/cli.py      # every hit in one file, with line numbers
python docs/release/wording_scan.py scan --min-level MUST       # only the hits that must go
python docs/release/wording_scan.py stats --against $START      # words in comments, docstrings and messages, before and after
python docs/release/wording_scan.py guard $START                # proves that only wording changed
```

Levels: **MUST** means remove it unless it sits inside a name you may not change. **FIX** means read it and change it unless it is a real option value or a path. **LOOK** means read it; most are fine.

The scanner knows words, not meaning. It cannot find a coined word it has never seen. Read each file once. A short list of other searches:

```
grep -rn -E "[A-Z][A-Z0-9]*_[A-Z0-9_]{3,}" docs README.md    # code constants in prose
grep -rn -E "\b[RSGHV][0-9][A-Z0-9]*\b" src                    # stage codes
grep -rn -i -E "frozen|incumbent|rescue|certif|ledger|oracle|verdict|headline" src docs README.md
```

### What `guard` checks

`guard REF` compares every changed file with its version at `REF`.

- A `.py` file: the syntax tree with docstrings removed and every string replaced by a placeholder must be identical. So no name, number, operator, call or statement may differ. Strings that changed are listed. A changed string must be in a `raise`, `assert`, `help=`, `description=`, `print`, logging, `warn` or `match=` position and must read like a sentence. Any other changed string fails: it may be data.
- A shell, YAML or TOML file: lines that are not comments must be identical. Differences are listed for you to judge (only echo text and help text may differ).
- A Markdown file: fenced code blocks must be identical. Backtick spans that disappeared are listed.
- Every file: a number that was in a comment, docstring or message and is no longer there is listed. Delete a number only with the history sentence that held it.
- A file added or removed fails. File names do not change. (The files in `docs/release/` named in this guide are allowed.)

Exit status 0 means no failure. Read every REVIEW line anyway and justify it in the report.

## 8. Checks that must pass at the end

All of these must match the start.

1. `python docs/release/wording_scan.py guard $START` exits 0, and every REVIEW line is explained.
2. `python docs/release/wording_scan.py scan --min-level MUST` lists only names you may not change. Each remaining line is explained in the report.
3. `python docs/release/wording_scan.py stats --against $START` shows no category with more words than before.
4. `python -m pytest -ra` gives the same passed, failed and skipped counts as at the start.
5. `bash tools/library_free_tests/run_all.sh` gives the same result as at the start. Known: `test_depth_window` fails with `PREVIOUS_SRC` set to a commit after the depth-cap change, and also fails that way at the base commit.
6. Byte identity of the library layer, once, at the end (needs Bowtie2, about 10 minutes on 4 CPUs):
   `PREVIOUS_SRC=<src folder of the $START checkout> bash tools/synthetic/compare_call.sh`
   Every compared file must be identical. The development values are passed by the script.
7. The same for the library-free tests: `PREVIOUS_SRC=<same folder> bash tools/library_free_tests/run_all.sh`.
8. Read `/tmp/help_before.txt` against the new help pages. Only the wording may differ. Every default, unit and value is still there.
9. `git diff --stat $START` shows no added or removed file except `docs/release/` files.

If a check fails, undo the last change that could cause it. Do not adjust a test, a threshold or an expected file to make a check pass.

## 9. Commits and the report

- Commit as `Md Taksir Hasan Majumder <Taksir@users.noreply.github.com>`: `git -c user.name="Md Taksir Hasan Majumder" -c user.email="Taksir@users.noreply.github.com" commit`.
- No `Co-Authored-By` line, no "generated by" line, no model name, in the message, the code or the docs.
- Plain commit messages, one area each, for example `Reword the help text and error messages in cli.py and errors.py`.
- Push to `release-cleanup`. Do not open a pull request. Do not touch `main`, tags, `LICENSE` or the version string.
- Write `docs/release/WORDING_CLEANUP_REPORT.md` in plain sentences:
  - the scan counts before and after (`scan --summary`);
  - the result of every check in section 8;
  - each message whose text a test matches, and the `match=` strings you changed;
  - each REVIEW line from `guard`, with a reason;
  - each label you could not understand and left (file and line);
  - each name you left because it is an identifier, option, path or file name, grouped by kind;
  - each `TODO`, and each place where a comment and the code disagree (do not fix either);
  - each claim word changed in the README or docs;
  - anything you were unsure of, as a short list of questions for Taksir.
