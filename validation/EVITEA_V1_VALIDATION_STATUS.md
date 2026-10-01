# EviTEA v1 validation status

**Recorded:** 2026-10-01  
**Verdict:** `EVITEA_V1_METHOD_EQUIVALENCE_PASS`  
The public pipeline (`prepare-reference` adopting the dm6 Bowtie2 index,
`prepare-te --build-languages`, `evitea call --msbwt-index ... --profile S`)
reproduced every frozen Route-B v3 layer byte for byte on both development
samples, natively on the Windows laptop, with both codebooks built by the
EviTEA language stages: `EVITEA_CALL_END_TO_END_EQUIVALENCE_PASS`
(`<work-dir>\EVITEA_CALL_END_TO_END_EQUIVALENCE.json`;
source commit `0a45bc9`, no tracked modifications — `git_dirty` reflects
untracked development notes only).  
**Historical pair-ID verdict:** `PAIR_ID_LEGACY_ADAPTER_ONLY`

This is a fail-closed implementation status, not a benchmark result. No truth,
TP/FP/FN identity, comparator callset, or score was used to reach it.

## Method-equivalence evidence (2026-10-01)

| Sample | R9E incumbent | bilateral | asymmetric | full memory | INTERNAL21 (final) | seconds |
|---|---|---|---|---|---|---:|
| AKA017 | 792 `1d1110f9` ✓ | 801 `ac352333` ✓ | 820 `8c12a561` ✓ | 835 `5efade8f` ✓ | **838 `6432c8e1`** ✓ | 11,354 |
| ISO1 | **192 `b3e48681`** ✓ | 192 `b3e48681` ✓ | 193 `31e15da3` ✓ | 194 `9e538d74` ✓ | **194 `9e538d74`** ✓ | 7,207 |

`calls.native.tsv` equals the final layer for both samples; languages bound
as `built:evitea.languages.v1` (terminal `d631f16a`, internal-only `44b98960`);
no unresolved integration events; no truth file opened.  This closes A
(default CLI), B (ISO1 R9E from the clean incumbent, with the ISO1 `{pid}`
recall-repair names under native Windows Bowtie2) and C (EviTEA-built
codebooks).

## Session 2026-09-30 (second): A–C

### A. Default `evitea call`

* `evitea.engine.run_engine` = `run_incumbent_caller` then `run_extension` on one
  msBWT; `pipeline.default_call_engine` is the CLI default (no injection).
* Canonical outputs: `calls.native.tsv` (frozen v3 schema, byte-exact; the
  equivalence target), `calls.tsv` / `calls.bed` / `unresolved.tsv` views
  (`evitea.callset`: explicit `call_origin`, unknown evidence types fail
  closed), `run_summary.json` with rows and SHA-256 after every layer.
* Profiles (`--profile`): `production` (default; AnchorPE host names
  `P{pid}`, recall-repair names `pair_{pid}`, sample-derived labels) and the
  reproduction-only `AKA017` (`p{pid}|h{mate}`) and `ISO1` (recall-repair
  `{pid}`), each bound to its own sample id.  Host profile `iso` was renamed
  `pid`.
* Prepared-TE schema v3 carries the two languages; `prepare-reference
  --bowtie2-index` adopts an existing index; `call --msbwt-index` uses an
  existing msBWT directory (index file hashes go into the run manifest).
* Retired: `integration.py`, `proposal.py`, `exact_junction.py`,
  `paired_end_anchor.py`, `certification.py`, `internal_probe.py` and their
  tests (the v3 layer invariants are covered by `test_shadow_layers.py`).
* VM: unit tests of profiles, wiring (incumbent bytes, codebooks, host names,
  labels, work layout), resume policy, language binding; a real-asset wiring
  check (prepared reference/TE from the frozen files + ISO1 msBWT) started the
  incumbent through the public CLI.

### B. ISO1 incumbent

Checked on the VM against the historical ISO1 tables
(`TE_locator_route_b/reports/G1*`, `TE_locator_route_b_v2/reports/V2_*`),
production functions, Linux Bowtie2 2.5.5, TEMP2 dm6 index
(`validation/incumbent_stage_checks/ISO1_STAGE_CHECKS_VM_20260930.json`):

| Stage | Check | Result |
|---|---|---|
| G1A | catalogue, provenance, 177,128 msBWT intervals | exact |
| G1B | panel 13,980, candidate pairs 153,405 | exact |
| G1C / G1D | short loci 678,502, H1 1,181; 13,903 actions, groups, sequences | exact |
| G1E | 994,474 pairs + read counts, groups, sequences | exact |
| G1F / G1G | 514,805 host mappings, 116,230 working pairs, pair families | exact |
| G1HR | 521 anchors, 365,738 evidence rows, 158,559 mappings, 127,075 repaired | exact **with `{pid}` query names** |
| candidate loci | 9,064 | exact |
| split junction (FCRC1) | all events 255, calls 105 | byte-identical |
| mate rescue (R4) | combined 162 calls, ledger 8,966 | byte-identical |
| R6 routing | 368,627 pair-locus routes | exact |

Finding: the historical ISO1 recall-repair Bowtie2 run named queries `{pid}`
(the AKA017 orchestrator later patched this to `pair_{pid}`); with `pair_{pid}`,
98,069 of 158,559 ISO1 mappings change their tie-broken placement and the loci
drift (9,060 vs 9,064).  This is now the `ISO1` profile's
`rescue_read_names = "numeric"`.  The downstream ISO1 Bowtie2 runs (WSL
Bowtie2 2.5.4, a separately built dm6 index) were replayed from their
historical inputs with Linux Bowtie2 2.5.5 and the TEMP2 index: identical dm6
and TE indexes; R6 global host (888,431 records), R7 host halves, R7 signatures
(all mapped records; the historical SAM lacks 153 unmapped rows), R8C and R8D
SAM bodies identical.  R6 certification through R9E (same code as AKA017) runs
on the laptop.  Launcher: `tools/windows/VALIDATE_INCUMBENT.ps1 -Sample ISO1`
(stage-by-stage diagnostics) or the end-to-end launcher above.

Native-Windows caveat: the research notes state that native Windows Bowtie2
2.5.5 rejects all-numeric FASTA names; if the ISO1 run stops at candidate
recall repair with a Bowtie2 error, run it with a WSL/Linux Bowtie2 2.5.x.

### C. Language (probe-set) builders

`evitea.languages` ports the four research builders (spec §8;
`route_b_v3_rebuild_safe_terminal21_codebook.py`,
`route_b_v3_fill_underpowered_terminal_codebook.py`,
`route_b_v3_finalize_safe_terminal_codebook.py`,
`route_b_v3_safe_internal21_codebook_gate.py`); the G1A candidate universe is
regenerated from the TE library by `incumbent.anchor_census` (no sample).
`validation/language_checks/run_language_equivalence.py`
(`evidence/LANGUAGE_EQUIVALENCE_VM_20260930.json`):

| Stage | Tables | Result |
|---|---|---|
| G1A catalogue / provenance | 177,128 / 178,870 rows, consumed columns | equal |
| rebuild terminal-21 | candidate audit 7,038, source audit 7,130, codebook 685 `64032d25`, record 180 / family 169 summaries | byte-identical |
| mixed-k gap fill | source audit 892, target summary 6, hybrid 685 `a0e08d92` | byte-identical |
| finalizer | final codebook 685 **`d631f16a`**, record status `7cd8509e`, repair log | byte-identical |
| internal-21 gate | candidate 15,193, source 15,215, full 12,280, internal-only 12,279 rows / 12,268 units **`44b98960`**, summaries | byte-identical |
| chain from regenerated G1A | all of the above | byte-identical |

`prepare-te --build-languages` binds them (`source: built:evitea.languages.v1`,
intermediates under `languages/audit/`).  Species hazards retained from the
research: FlyBase `FBte` identifier parsing, the 125-bp read block column,
k in {21, 25, 31, 35}.

## Completed evidence

- Unit, negative, fixture, output-determinism, checkpoint/resume, and package
  boundary tests cover the implemented scientific kernels.
- The installed `msbwt-modern3==0.3.0` path builds a prepared sample, recovers
  exact physical pair/mate identity from `about.npy`, and verifies frozen index
  file hashes on reload.
- Prepared-TE v2 accepts the externally supplied, exact-hash G1A archive and
  family crosswalk, preserves its candidate universe, and rebuilds terminal,
  full-internal, and internal-only languages with exhaustive TE-family and
  RepeatMasker-aware reference-context checks.
- Paired-reference S30 observations and ordered terminal/INTERNAL21 integration adapters are implemented
  with bounded tests.
- The validation-only frozen-oracle verifier checks the five recorded layers for
  AKA017 (`792 -> 801 -> 820 -> 835 -> 838`) and ISO1
  (`192 -> 192 -> 193 -> 194 -> 194`) plus their recorded final hashes.
- Oracle integrity is explicitly distinct from clean-engine method equivalence.

## Incumbent (S40) producer — resolved on AKA017

The clean incumbent branch (`evitea.incumbent.run_incumbent_caller`) derives
the repeat-aware incumbent callset from raw prepared inputs only: the sample
msBWT and `about.npy`, the TE consensus FASTA, the FBte family map, the host
reference FASTA and its Bowtie2 index. It reads no historical report file.
Every stage was run on the full AKA017 data (Linux Bowtie2 2.5.5, VM) and
compared with the frozen truth-blind artifacts:

| Stage | Frozen artifact | Result |
|---|---|---|
| Anchor census .. recall repair (G1A–G1H-R) | G1A..G1HR archives | all tables exact |
| Candidate loci | `V2FCRC1_CANDIDATE_LOCI.tsv` (5,971) | exact |
| Split-junction layer | `V2FCRC1_ALL_EVENTS.tsv` (806), `V2FCRC1_CALLS_SCORE4.tsv` (407) | byte-identical |
| Two-sided mate rescue | `V2R4_COMBINED_CALLS.tsv` (+250), locus ledger (5,585) | byte-identical / exact |
| Host-seed routing, strict pairs | `V2R6_PAIR_LOCUS_ROUTING.tsv` (350,476), 721 winners | exact |
| Junction groups | discovery (38,525, incl. query IDs), screen, 2,604 winners, ledger | exact |
| TE certification + repeat-aware rescue | `R9E_COMBINED_CALLS.tsv` (792) | byte-identical, SHA256 `1d1110f9...1757` |

In the VM the per-locus layers and the long stages were first run in
resumable chunks through the same production functions, and the assembled
engine reproduced the SHA over those checkpoints.

**Native end-to-end confirmation (2026-09-30).** A single uninterrupted run of
`tools/windows/VALIDATE_INCUMBENT.ps1` on the Windows laptop (native Bowtie2,
fresh work directory `<work-dir>`, 6,172 s) returned
`INCUMBENT_EQUIVALENCE_PASS` with every check true (`candidate_loci`,
`split_junction_all_events`, `split_junction_calls`, `mate_rescue_combined`,
`repeat_aware_final`, `frozen_final_sha256`) and final SHA256
`1d1110f939996223379edf74847d626d4bc397a5cac32ad16d5d985dca051757`.
ISO1 has not yet been run through the clean incumbent.

## v3 extension integration (S100) — reproduced exactly

`evitea.shadow_layers` implements the four frozen integration freezes
(bilateral, asymmetric, full-memory, INTERNAL21), each appending to the
previous callset as an exact byte prefix in the native R9E schema.  Fed with the
frozen research input tables,
`validation/v3_layer_checks/run_shadow_layer_equivalence.py` returns
`SHADOW_LAYER_INTEGRATION_EQUIVALENCE_PASS`: all ten layer SHA-256s (AKA017
`792 -> 801 -> 820 -> 835 -> 838`, ISO1 `192 -> 192 -> 193 -> 194 -> 194`) and
every integration decision / nearest-call attribution in the four historical
ledgers match.  This isolates integration semantics; the layer inputs are
still frozen research tables.  Source lineage and kernel discrepancies were recorded in
the author's private development notes.

## AnchorPE candidate procedure — one specification for both samples

The self-contained ISO-1 AnchorPE transfer procedure (the one reused unchanged
for GIM-024), run unchanged on AKA017 with the AKA historical Bowtie2 host
read-name profile (`p{pid}|h{host_mate}`), reproduces AKA017's frozen AnchorPE
candidate table exactly, including cluster ids
(`ISO_PROCEDURE_REPRODUCES_AKA017_CANDIDATES_EXACTLY`; 1,052 corrected
clusters, 168 targets, 91 anchors, 8,783 B pairs, 32 bilateral).  Bowtie2 query
names seed its per-read tie-breaking and are part of the reproduction profile.

## AnchorPE producer port (`evitea.anchor_pe`)

The procedure above is ported as typed production stages (codebook projection,
pair classification, fragment decoder, host / LOCAL paired-reference Bowtie2
with an explicit host read-name profile, evidence-complete ledger, clustering,
novelty/residency, Channel A, Channel B).  Reads come from the msBWT; no FASTQ
is needed.  Component equivalence on AKA017 (full data, VM, against the
laptop probe's artefacts):

| Component | Check | Result |
|---|---|---|
| msBWT codebook projection | all 175,818 pairs, unit sets and order | identical |
| read recovery from msBWT | all 128,588 decoder-universe pairs vs FASTQ | identical |
| fragment decoder (edlib 1.3.9.post1) | all 128,588 rows, every field | identical |
| Bowtie2 inputs | host (`p{pid}|h{hm}`), decoded-host, R1/R2 FASTAs | byte-identical |
| clustering .. candidates, from the Bowtie2 SAMs | primary 817, relaxed 399, corrected 1,052, targets 168, anchors 91, B 8,783 | identical; typed candidates equal the frozen table |
| bilateral / asymmetric layers from those candidates | 801 / 820 | byte-identical SHA-256 |

**Native end-to-end confirmation (2026-09-30).**
`tools/windows/VALIDATE_V3_EXTENSION.ps1` ran the engine uninterrupted from each
prepared msBWT on the Windows laptop and returned
`ANCHOR_PE_ENGINE_EQUIVALENCE_PASS`:

| Sample | host read names | candidates | bilateral | asymmetric | seconds |
|---|---|---|---|---|---:|
| AKA017 | `aka-historical` | = frozen (91; 32 bilateral, 38 one-sided) | 801 ✓ | 820 ✓ | 511 |
| ISO1 | `iso` | = frozen (17; 5 bilateral, 8 one-sided) | 192 ✓ | 193 ✓ | 568 |

Layers R9E → bilateral → asymmetric are therefore derived from raw prepared
inputs for both samples (given the frozen safe terminal codebook as an asset).

## Full-memory branch port (`evitea.full_memory`)

Typed port of the frozen full-memory chain: finite host/reference state and
pair accounting (research "Section 128"), terminal proposal memory (full data
and ten rotations), global strict exact A with the full-data event ledger, and
A-only proposal selection with full-data Channel-B support.  Stage checks
against the frozen tables (VM):

| Stage | Input | Check | AKA017 | ISO1 |
|---|---|---|---|---|
| finite state | probe's `P{pid}` host SAM + LOCAL reference SAM | host alternatives (210,014), reference alternatives (1,756,195), accounting | identical¹ | engine run |
| proposals | frozen Section-128 tables | full-data / rotation proposals | 755 / 6,979 identical | 148 / 1,353 identical |
| strict A + event ledger | frozen Section-128 tables + AnchorPE witnesses and B universe | ledger rows | 780 identical | 149 identical |
| selection + full-data B | frozen event ledger | selected events | 458 identical | 85 identical |
| full-memory layer | selected events on the asymmetric shadow | SHA-256 | 835 ✓ | 194 ✓ |

¹ Except the `family` field of the 572 `FAMILY_CONFLICT` rows (not
consumed downstream): the AKA table inherited `Dm88;invader1`-style joined
names from its historical eligibility input; EviTEA writes an empty family as
the ISO table does.  The historical A-witness keys that decide
`proposal_only_A` equal the AnchorPE engine's own Channel-A witness keys
(155/155, pair-level identical), so no historical file is needed.

The finite-state host replay uses `P{pid}` query names for both samples
(so AKA017 needs a second host mapping besides its `aka-historical` AnchorPE
one).

**Native end-to-end confirmation (2026-09-30).** The chained engines
(AnchorPE -> full memory) ran uninterrupted on the laptop from each prepared
msBWT and returned `V3_ANCHOR_PE_AND_FULL_MEMORY_ENGINE_EQUIVALENCE_PASS`:
full-data and rotation proposals, the event ledger and the selected events
equal the frozen tables for both samples, and the layers are byte-identical —
AKA017 801 / 820 / 835, ISO1 192 / 193 / 194.  The ISO1 finite state was
therefore reproduced from its own Bowtie2 output.

## INTERNAL21 port (`evitea.internal21`) and the extension engine

`internal21.acquisition` projects the frozen internal-only safe 21-mer
codebook (12,268 probes) on the msBWT, keeps one-mate family-resolved pairs,
and gives them the same finite host / LOCAL reference state as the terminal
branch — evidence completeness judged on each pair's own observed internal
probes — yielding internal strict-A witnesses.  `internal21.observability`
certifies, under the fixed terminal proposal memory, the exact events that
OLD + INTERNAL A certifies and OLD A alone does not, per 3/2 fold rotation,
and classifies their 500-bp same-family components.
`evitea.extension_engine.run_extension` chains AnchorPE -> full memory ->
INTERNAL21 and composes the four layers.

Stage checks (VM, frozen inputs):

| Stage | Check | AKA017 | ISO1 |
|---|---|---|---|
| internal projection | observed-probe sets on a 1,500-probe random subset (449,526 pairs) | — | identical |
| host alternatives | internal host SAM -> rows / per-pair summaries | — | 888,110 identical |
| reference alternatives (pair-own probes) | first 10,558 pairs of the internal reference SAM | — | 149,981 rows identical |
| internal strict-A witnesses | frozen internal accounting + reference tables | 11,550 identical (order too) | 3,758 identical |
| observability | frozen OLD / internal witnesses, proposals, B, prior ledger | 415 instances, 119 loci identical | 15 instances, 6 loci identical |
| INTERNAL21 layer | on the full-memory shadow | **838 ✓** (`6432c8e1…f701`) | **194 ✓** (`9e538d74…92b0`) |

**Native end-to-end confirmation (2026-09-30).**
`tools/windows/VALIDATE_V3_EXTENSION.ps1` ran `run_extension` uninterrupted on
the laptop from each prepared msBWT and returned
`V3_EXTENSION_ENGINE_EQUIVALENCE_PASS`: every intermediate equals the frozen
tables and every layer is byte-identical.

| Sample | bilateral | asymmetric | full-memory | INTERNAL21 | seconds |
|---|---|---|---|---|---:|
| AKA017 | 801 `ac352333` | 820 `8c12a561` | 835 `5efade8f` | **838 `6432c8e1`** | 2,383 |
| ISO1 | 192 `b3e48681` | 193 `31e15da3` | 194 `9e538d74` | **194 `9e538d74`** | 1,890 |

INTERNAL21 acquisition: AKA017 2,128,697 hit pairs / 541,722 one-mate
family-resolved / 11,550 internal witnesses; ISO1 1,395,976 / 506,538 / 3,758.

## Remaining (not method-equivalence blockers)

* `prepare-sample` writes a JSON identity map with one entry per read,
  impractical for full-size libraries; use `call --msbwt-index`.
* `evidence/frozen_candidates.py` (spec-era G1A codebooks) is unused by the
  engine.
* Species transfer hazards: `full_memory.common.PRIMARY_CHROMS`, `Q99 = 805`,
  FlyBase `FBte` parsing in the language builders, PE125 validation only.

## Regression suite

`python -m pytest` passes (185 passed, 1 live-msBWT test skipped) in the Linux
Python 3.14 environment.  Wheel/sdist creation remains unverified in the
managed Windows/Drive environment; run `tools/windows/VALIDATE_LOCAL.ps1` on an
ordinary local shell.
