# EviTEA v1 validation status

**Result:** `EVITEA_V1_METHOD_EQUIVALENCE_PASS`  
The public pipeline (`prepare-reference` adopting the dm6 Bowtie2 index,
`prepare-te --build-languages`, `evitea call --msbwt-index ... --profile S`)
reproduced every recorded layer byte for byte on both development
samples, natively on the Windows laptop, with both codebooks built by the
EviTEA codebook builders: `EVITEA_CALL_END_TO_END_EQUIVALENCE_PASS`
(`<work-dir>\EVITEA_CALL_END_TO_END_EQUIVALENCE.json`).
The source had no tracked modifications. `git_dirty` reflects untracked
development notes only.  
**Pair-ID result:** `PAIR_ID_LEGACY_ADAPTER_ONLY`

This is a status of the implementation, not a benchmark result. It uses no
truth, TP/FP/FN identity, other-tool callset, or score.
Method equivalence means that EviTEA reproduces the recorded layers byte for
byte from prepared inputs.

## Method-equivalence evidence

| Sample | library layer | bilateral | asymmetric | full memory | internal 21-mer (final) | seconds |
|---|---|---|---|---|---|---:|
| AKA017 | 792 `1d1110f9` match | 801 `ac352333` match | 820 `8c12a561` match | 835 `5efade8f` match | **838 `6432c8e1`** match | 11,354 |
| ISO1 | **192 `b3e48681`** match | 192 `b3e48681` match | 193 `31e15da3` match | 194 `9e538d74` match | **194 `9e538d74`** match | 7,207 |

`calls.native.tsv` equals the final layer for both samples. The codebooks are
bound as `built:evitea.languages.v1` (terminal `d631f16a`, internal-only
`44b98960`). No integration events are unresolved, and no truth file was
opened. This closes three items: the default CLI, the ISO1 library layer from
the clean library-layer caller (with the ISO1 `{pid}` wide-anchor query names
under native Windows Bowtie2), and the EviTEA-built codebooks.

## Default call, ISO1 library layer and codebook builders

### Default `evitea call`

* `evitea.engine.run_engine` runs `run_incumbent_caller`, then `run_extension`,
  on one msBWT. `pipeline.default_call_engine` is the CLI default (no
  injection).
* Main outputs: `calls.native.tsv` (native layer schema, byte-exact; the
  equivalence target) and the views `calls.tsv`, `calls.bed` and
  `unresolved.tsv` (`evitea.callset`: explicit `call_origin`; an unknown
  evidence type raises an error). `run_summary.json` has the rows and SHA-256
  after every layer.
* Profiles (`--profile`): `production` (default; paired-end anchor host names
  `P{pid}`, wide-anchor step names `pair_{pid}`, sample-derived labels) and the
  reproduction-only `AKA017` (`p{pid}|h{mate}`) and `ISO1` (wide-anchor step
  `{pid}`). Each is bound to its own sample id. The host profile `pid` is
  listed as `iso` in the tables below.
* The prepared-TE schema v3 contains the two codebooks.
  `prepare-reference --bowtie2-index` adopts an existing index.
  `call --msbwt-index` uses an existing msBWT directory. Index file hashes go
  into the run manifest.
* The modules `integration.py`, `proposal.py`, `exact_junction.py`,
  `paired_end_anchor.py`, `certification.py` and `internal_probe.py` and their
  tests are not part of the package. `test_shadow_layers.py` covers the
  invariants of the extension layers.
* VM: unit tests cover profiles, wiring (library-layer bytes, codebooks, host
  names, labels, work layout), resume policy and codebook binding. A wiring
  check on real assets (prepared reference and TE from the recorded files, plus
  the ISO1 msBWT) started the library layer through the public CLI.

### ISO1 library layer

Run on the VM and compared with the recorded ISO1 tables
(`TE_locator_route_b/reports/G1*`, `TE_locator_route_b_v2/reports/V2_*`),
the package functions, Linux Bowtie2 2.5.5 and the dm6 index
(`validation/incumbent_stage_checks/ISO1_STAGE_CHECKS_VM_20260930.json`):

| Stage | Check | Result |
|---|---|---|
| anchor census | catalog, provenance, 177,128 msBWT intervals | exact |
| sparse panel | panel 13,980, candidate pairs 153,405 | exact |
| coverage panel | short loci 678,502, Hamming-1 loci 1,181; 13,903 actions, groups, sequences | exact |
| pair-aware panel | 994,474 pairs + read counts, groups, sequences | exact |
| host-mate gate, resident pair families | 514,805 host mappings, 116,230 working pairs, pair families | exact |
| wide-anchor step | 521 anchors, 365,738 evidence rows, 158,559 mappings, 127,075 repaired | exact **with `{pid}` query names** |
| candidate loci | 9,064 | exact |
| split junction | all events 255, calls 105 | byte-identical |
| mate-anchored search | combined 162 calls, locus record 8,966 | byte-identical |
| host-seed routing | 368,627 pair-locus routes | exact |

Finding: the recorded ISO1 wide-anchor Bowtie2 run named its queries `{pid}`.
With `pair_{pid}`, 98,069 of 158,559 ISO1 mappings change their tie-broken
placement and the loci drift (9,060 vs 9,064). The `ISO1` profile therefore
sets `rescue_read_names = "numeric"`. The downstream ISO1 Bowtie2 runs (WSL
Bowtie2 2.5.4, a separately built dm6 index) were replayed from their recorded
inputs with Linux Bowtie2 2.5.5 and the dm6 index. The dm6 and TE indexes are
identical. These outputs are identical: the host-seed routing global host
mapping (888,431 records), the junction-group host halves, the junction-group
signatures (all mapped records; the recorded SAM does not contain 153 unmapped rows), and
the SAM bodies of the TE-segment family confirmation and trim-variant runs. The
stages after host-seed routing, up to the repeat-aware step (same code as for
AKA017), run on the laptop. Launcher:
`tools/windows/VALIDATE_INCUMBENT.ps1 -Sample ISO1` (stage-by-stage
diagnostics) or the end-to-end launcher above.

Native-Windows caveat: the research notes state that native Windows Bowtie2
2.5.5 rejects all-numeric FASTA names. If the ISO1 run stops in the wide-anchor
step with a Bowtie2 error, run it with a WSL or Linux Bowtie2 2.5.x.

### Codebook (probe-set) builders

`evitea.languages` ports the four research builders
(`route_b_v3_rebuild_safe_terminal21_codebook.py`,
`route_b_v3_fill_underpowered_terminal_codebook.py`,
`route_b_v3_finalize_safe_terminal_codebook.py`,
`route_b_v3_safe_internal21_codebook_gate.py`). `incumbent.anchor_census`
regenerates the anchor-census candidate universe from the TE library (no
sample). `validation/language_checks/run_language_equivalence.py`
(`evidence/LANGUAGE_EQUIVALENCE_VM_20260930.json`) compares the results:

| Stage | Tables | Result |
|---|---|---|
| anchor census catalog / provenance | 177,128 / 178,870 rows, consumed columns | equal |
| rebuild terminal 21-mer codebook | candidate audit 7,038, source audit 7,130, codebook 685 `64032d25`, record 180 / family 169 summaries | byte-identical |
| mixed-k gap fill | source audit 892, target summary 6, hybrid 685 `a0e08d92` | byte-identical |
| finalizer | final codebook 685 **`d631f16a`**, record status `7cd8509e`, repair log | byte-identical |
| internal 21-mer filter | candidate 15,193, source 15,215, full 12,280, internal-only 12,279 rows / 12,268 units **`44b98960`**, summaries | byte-identical |
| chain from regenerated anchor census | all of the above | byte-identical |

`prepare-te --build-languages` binds the codebooks
(`source: built:evitea.languages.v1`, intermediates under `languages/audit/`).
Species-specific assumptions: FlyBase `FBte` identifier parsing, the 125-bp
read block column, k in {21, 25, 31, 35}.

## Completed evidence

- Unit, negative, fixture, output-determinism, checkpoint/resume, and package
  boundary tests cover the implemented scientific kernels.
- The installed `msbwt-modern3==0.3.0` path builds a prepared sample, reads the
  exact physical pair/mate identity from `about.npy`, and verifies fixed index
  file hashes on reload.
- Prepared-TE accepts the externally provided, exact-hash anchor-census archive
  and family crosswalk, preserves its candidate universe, and rebuilds the
  terminal, full-internal, and internal-only codebooks with exhaustive
  TE-family and RepeatMasker-aware reference-context checks.
- Paired-reference observations and ordered terminal/internal 21-mer
  integration adapters are implemented with bounded tests.
- The validation-only verifier of the recorded call files checks the five
  recorded layers for AKA017 (`792 -> 801 -> 820 -> 835 -> 838`) and ISO1
  (`192 -> 192 -> 193 -> 194 -> 194`) plus their recorded final hashes.
- Integrity of the recorded files is distinct from method equivalence of the
  clean engine.

## Library layer producer, reproduced on AKA017

The clean library-layer branch (`evitea.incumbent.run_incumbent_caller`)
derives the repeat-aware library-layer callset from raw prepared inputs only:
the sample msBWT and `about.npy`, the TE consensus FASTA, the FBte family map,
the host reference FASTA and its Bowtie2 index. It reads no report file from the
research runs. Every stage was run on the full AKA017 data (Linux Bowtie2 2.5.5, VM) and
compared with the recorded truth-blind files:

| Stage | Recorded file | Result |
|---|---|---|
| Anchor census to wide-anchor step | the archives of these stages | all tables exact |
| Candidate loci | `V2FCRC1_CANDIDATE_LOCI.tsv` (5,971) | exact |
| Split-junction layer | `V2FCRC1_ALL_EVENTS.tsv` (806), `V2FCRC1_CALLS_SCORE4.tsv` (407) | byte-identical |
| Two-sided mate-anchored search | `V2R4_COMBINED_CALLS.tsv` (+250), locus record (5,585) | byte-identical / exact |
| Host-seed routing, strict pairs | `V2R6_PAIR_LOCUS_ROUTING.tsv` (350,476), 721 winners | exact |
| Junction groups | discovery (38,525, incl. query IDs), screen, 2,604 winners, record | exact |
| TE-family confirmation + repeat-aware step | `R9E_COMBINED_CALLS.tsv` (792) | byte-identical, SHA256 `1d1110f9...1757` |

In the VM, the per-locus layers and the long stages ran first in resumable
chunks through the same package functions. The assembled engine reproduced the
SHA over those checkpoints.

**Native end-to-end confirmation.** A single uninterrupted run of
`tools/windows/VALIDATE_INCUMBENT.ps1` on the Windows laptop (native Bowtie2,
fresh work directory `<work-dir>`, 6,172 s) returned
`INCUMBENT_EQUIVALENCE_PASS` with every check true (`candidate_loci`,
`split_junction_all_events`, `split_junction_calls`, `mate_rescue_combined`,
`repeat_aware_final`, `frozen_final_sha256`) and final SHA256
`1d1110f939996223379edf74847d626d4bc397a5cac32ad16d5d985dca051757`.
This run covered AKA017 only.

## Extension layer integration, reproduced exactly

`evitea.shadow_layers` implements the four recorded integration steps
(bilateral, asymmetric, full-memory, internal 21-mer). Each appends to the
previous callset as an exact byte prefix in the native library-layer schema.
Fed with the recorded research input tables,
`validation/v3_layer_checks/run_shadow_layer_equivalence.py` returns
`SHADOW_LAYER_INTEGRATION_EQUIVALENCE_PASS`: all ten layer SHA-256s (AKA017
`792 -> 801 -> 820 -> 835 -> 838`, ISO1 `192 -> 192 -> 193 -> 194 -> 194`) and
every integration decision and nearest-call attribution in the four recorded
decision tables match. This isolates integration semantics. The layer inputs
are the recorded research tables. Differences in source lineage and kernels are
in the author's private development notes.

## Paired-end anchor support candidates, one specification for both samples

The self-contained ISO1 transfer procedure for paired-end anchor support
candidates, run unchanged on AKA017 with the AKA017 Bowtie2 host read-name
profile (`p{pid}|h{host_mate}`), reproduces the recorded paired-end anchor
candidate table of AKA017 exactly, including cluster ids
(`ISO_PROCEDURE_REPRODUCES_AKA017_CANDIDATES_EXACTLY`; 1,052 corrected
clusters, 168 targets, 91 anchors, 8,783 paired-end support pairs, 32
bilateral). Bowtie2 query names seed its per-read tie-breaking and are part of
the reproduction profile.

## Paired-end anchor support producer (`evitea.anchor_pe`)

The procedure above is ported as typed stages of the package (codebook
projection, pair classification, fragment decoder, host / LOCAL
paired-reference Bowtie2 with an explicit host read-name profile,
evidence-complete record, clustering, novelty/residency, exact junction
evidence, paired-end support). Reads come from the msBWT; no FASTQ is needed.
Component equivalence on AKA017 (full data, VM, against the laptop probe's
artifacts):

| Component | Check | Result |
|---|---|---|
| msBWT codebook projection | all 175,818 pairs, unit sets and order | identical |
| read retrieval from msBWT | all 128,588 decoder-universe pairs vs FASTQ | identical |
| fragment decoder (edlib 1.3.9.post1) | all 128,588 rows, every field | identical |
| Bowtie2 inputs | host (`p{pid}|h{hm}`), decoded-host, R1/R2 FASTAs | byte-identical |
| clustering .. candidates, from the Bowtie2 SAMs | primary 817, relaxed 399, corrected 1,052, targets 168, anchors 91, paired-end support pairs 8,783 | identical; typed candidates equal the recorded table |
| bilateral / asymmetric layers from those candidates | 801 / 820 | byte-identical SHA-256 |

**Native end-to-end confirmation.**
`tools/windows/VALIDATE_V3_EXTENSION.ps1` ran the engine uninterrupted from each
prepared msBWT on the Windows laptop and returned
`ANCHOR_PE_ENGINE_EQUIVALENCE_PASS`:

| Sample | host read names | candidates | bilateral | asymmetric | seconds |
|---|---|---|---|---|---:|
| AKA017 | `aka-historical` | = recorded (91; 32 bilateral, 38 one-sided) | 801 match | 820 match | 511 |
| ISO1 | `iso` | = recorded (17; 5 bilateral, 8 one-sided) | 192 match | 193 match | 568 |

The library layer, the bilateral layer and the asymmetric layer are therefore
derived from raw prepared inputs for both samples, given the fixed safe
terminal codebook as an input file.

## Full-memory branch (`evitea.full_memory`)

Typed port of the recorded full-memory chain: finite host/reference state and
pair accounting, terminal proposal memory (full data and ten rotations), global
strict exact junction evidence with the full-data event record, and proposal
selection from exact junction evidence alone with full-data paired-end
support. Stage checks against the recorded tables (VM):

| Stage | Input | Check | AKA017 | ISO1 |
|---|---|---|---|---|
| finite state | probe's `P{pid}` host SAM + LOCAL reference SAM | host alternatives (210,014), reference alternatives (1,756,195), accounting | identical [1] | engine run |
| proposals | recorded research tables | full-data / rotation proposals | 755 / 6,979 identical | 148 / 1,353 identical |
| strict exact junction evidence + event record | recorded research tables + paired-end anchor witnesses and support universe | record rows | 780 identical | 149 identical |
| selection + full-data paired-end support | recorded event record | selected events | 458 identical | 85 identical |
| full-memory layer | selected events on the asymmetric layer | SHA-256 | 835 match | 194 match |

[1] Except the `family` field of the 572 `FAMILY_CONFLICT` rows (not consumed
downstream). The AKA017 table has joined names such as `Dm88;invader1` from its
recorded eligibility input, and EviTEA writes an empty family as the ISO1 table
does. The recorded exact-junction witness keys that decide `proposal_only_A`
equal the paired-end anchor engine's own exact-junction witness keys (155/155,
pair-level identical), so no recorded file is needed.

The finite-state host replay uses `P{pid}` query names for both samples. AKA017
therefore needs a second host mapping besides its `aka-historical` paired-end
anchor one.

**Native end-to-end confirmation.** The chained engines (paired-end anchor
support -> full memory) ran uninterrupted on the laptop from each prepared
msBWT and returned `V3_ANCHOR_PE_AND_FULL_MEMORY_ENGINE_EQUIVALENCE_PASS`:
full-data and rotation proposals, the event record and the selected events
equal the recorded tables for both samples, and the layers are byte-identical:
AKA017 801 / 820 / 835, ISO1 192 / 193 / 194. The ISO1 finite state was
therefore reproduced from its own Bowtie2 output.

## Internal 21-mer branch (`evitea.internal21`) and the extension engine

`internal21.acquisition` projects the fixed internal-only safe 21-mer codebook
(12,268 probes) onto the msBWT and retains pairs with one mate resolved to a
family. These pairs get the same finite host / LOCAL reference state as the
terminal branch, with evidence completeness judged on each pair's own observed
internal probes. The result is a set of internal exact junction witnesses.
`internal21.observability` confirms, under the fixed terminal proposal memory
and per 3/2 fold rotation, the exact events that the existing plus internal
exact junction evidence confirms and the existing evidence alone does not. It
classifies their 500-bp same-family components.
`evitea.extension_engine.run_extension` chains paired-end anchor support ->
full memory -> internal 21-mer and composes the four layers.

Stage checks (VM, recorded inputs):

| Stage | Check | AKA017 | ISO1 |
|---|---|---|---|
| internal projection | observed-probe sets on a 1,500-probe random subset (449,526 pairs) | - | identical |
| host alternatives | internal host SAM -> rows / per-pair summaries | - | 888,110 identical |
| reference alternatives (pair-own probes) | first 10,558 pairs of the internal reference SAM | - | 149,981 rows identical |
| internal exact junction witnesses | recorded internal accounting + reference tables | 11,550 identical (order too) | 3,758 identical |
| observability | recorded existing / internal witnesses, proposals, paired-end support, prior record | 415 instances, 119 loci identical | 15 instances, 6 loci identical |
| internal 21-mer layer | on the full-memory layer | **838 match** (`6432c8e1...f701`) | **194 match** (`9e538d74...92b0`) |

**Native end-to-end confirmation.**
`tools/windows/VALIDATE_V3_EXTENSION.ps1` ran `run_extension` uninterrupted on
the laptop from each prepared msBWT and returned
`V3_EXTENSION_ENGINE_EQUIVALENCE_PASS`: every intermediate equals the recorded
tables and every layer is byte-identical.

| Sample | bilateral | asymmetric | full-memory | internal 21-mer | seconds |
|---|---|---|---|---|---:|
| AKA017 | 801 `ac352333` | 820 `8c12a561` | 835 `5efade8f` | **838 `6432c8e1`** | 2,383 |
| ISO1 | 192 `b3e48681` | 193 `31e15da3` | 194 `9e538d74` | **194 `9e538d74`** | 1,890 |

Internal 21-mer acquisition: AKA017 2,128,697 hit pairs / 541,722 one-mate
family-resolved / 11,550 internal witnesses; ISO1 1,395,976 / 506,538 / 3,758.

## Open items (none affects method equivalence)

* `prepare-sample` writes a JSON identity map with one entry per read, which is
  impractical for full-size libraries. Use `call --msbwt-index` instead.
* `evidence/frozen_candidates.py` (anchor-census codebooks from an earlier
  specification) is not used by the engine.
* Species-specific assumptions: `full_memory.common.PRIMARY_CHROMS`,
  `Q99 = 805`, FlyBase `FBte` parsing in the codebook builders, PE125
  validation only.

## Regression suite

`python -m pytest` passes (185 passed, 1 live-msBWT test skipped) in the Linux
Python 3.14 environment. Wheel and sdist creation is untested in the managed
Windows/Drive environment. Run `tools/windows/VALIDATE_LOCAL.ps1` on an
ordinary local shell.
