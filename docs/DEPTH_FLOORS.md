# Depth-aware read floors (opt-in)

Three switches change a minimum number of reads of the pipeline for a low-depth sample. All are **off by default**. With none
given, every output is the freeze's (`296e725`) byte for byte; the saved scan pieces and library-layer checkpoints of an earlier run
are reused as before. No value is tuned here: the switches only make the minimum a setting.

The depth D is the sample's single-copy depth, the `genome_depth` of `library_free/single_copy.py` (`measure`): the median sample
count of 200,000 sampled 21-mers that occur once in the reference, over the default contigs. It is measured only when a depth cut
is given, once per work directory (`copy_profile.json`).

## 1. Library-free extension floor (`library-free`, `run`)

| Switch | Default | Effect |
|---|---|---|
| `--extension-stop-min N` | 3 | The extension of an inserted segment goes on while the best next letter has at least `max(N, 0.10 x the support at the break)` reads. 3 is the research rule. |
| `--extension-stop-min-depth-cut DCUT` | none | N applies only when D is at or under DCUT; above it the floor is 3. |

The floor used, the requested value, the cut and D are written to `library_free_summary.json` (`detection.extension_floor`).
The detection code is not edited (its hash, part of a saved piece's identity, is unchanged): the run sets
`library_free/_p4_lib.STOP_MIN` for the scan and restores it after. A floor other than 3 enters the identity of the saved pieces,
so it never reuses the default's pieces; it cannot be combined with `--detection-dir`.

## 2. Weak-site pass (`call`, `run`)

| Switch | Default | Effect |
|---|---|---|
| `--weak-site-pass` | off | After the four extension layers, call each evidence site with no row of the final table within 500 bp that has at least one pooled read on each side, or at least two split reads on one side. |
| `--weak-site-depth-cut DCUT` | required with the pass | The pass runs only when D is at or under DCUT. |

Evidence is grouped as in the desk replay of 3 Nov 2026 (`incumbent/weak_sites.py`): split events (one-sided left or right), the
left and right pair ids of the mate-rescue loci and the strict ids of the junction-group loci; items within 100 bp join one site;
per side the distinct fragment ids are pooled. The family of the new call is the family with the most supporting ids at the site;
its coordinates are the span of that family's items. The calls are named `WEAK_<n>`, have `evidence_type`
`weak_site_pooled_floor`, call origin `INCUMBENT_WEAK_SITE`, and tier 4 (`LOW_SUPPORT`) at any support. They go through the same
own-family copy filter and breakpoint refinement as every library call. The pass runs after the extension layers, the 500 bp cover is every row of the final table (library-layer and extension rows), and the rows are appended at the end (layer `weak_site` in `layer_rows`), so no extension layer reads a `WEAK_` row: a weak row cannot absorb an anchor-pair call of its family or stand beside a call of another family.
D, DCUT, the sites seen, the rows of the final table before the pass and the rows added are in `run_summary.json` (`engine_counts.depth_floors.weak_site_pass`).

## 2b. Weak-site pass in the library-free step, switched by the allele depth (`library-free`, `run`)

A second route to the same pass, off by default, separate from section 2 (give one of the two, not both). The sample depth D counts
every read, but a heterozygous or pooled sample carries the inserted allele in a share of them, so the pass is switched by the allele
depth `Da = D x a`.

| Switch | Default | Effect |
|---|---|---|
| `--weak-site-pass` | off | Run the pass of section 2 after the library-free calls and statuses are final, before `calls.combined.tsv` and `calls.combined.vcf` are written. |
| `--weak-site-allele-depth-cut DA` | required with the pass | The pass runs only when `Da` is at or under DA. |

- D is `genome_depth` of `library_free/single_copy.py` (as in section 1). `a = 1 - median(e)` over the added library-free calls of the main view: every row of
  `calls.library_free.tsv` (first pass and haplotype scan) with a measured e; `a = 1` when fewer than 10 have an e. (A pooled or
  heterozygous sample shows its allele share in the e of the added calls better than in the confirmed places of the het guard.)
- Evidence: the library layer's saved stage tables `split_junction`, `mate_rescue`, `junction_groups` from the work directory of the
  `call` run, named in the run folder's `run_manifest.json` (`normalized_config.runtime.work_dir`, else `<run>/.evitea-work`), under
  `engine/work/incumbent/checkpoints/`. A missing table stops the run with a message. `run --no-library` and
  `--weak-site-allele-depth-cut` together are refused; `run` needs `--library-free` for it, and `--weak-site-pass` with the allele cut
  goes to the library-free step, not to `call`.
- Grouping, rule, family, coordinates and `WEAK_<n>` names are those of section 2 (`incumbent/weak_sites.py`). The 500 bp cover is every
  row the combined table will hold: the library rows (after the own-family copy filter) and the added library-free calls.
- The new rows go through the same own-family copy filter (not with `--keep-resident-copies`); the rows it drops are listed in
  `resident_copy_dropped.tsv`. The kept rows are written to `calls.combined.tsv` and `calls.combined.vcf` as library rows (source
  EVITEA, after the other library rows), with `evidence_type` `weak_site_pooled_floor` and tier 4 in `call-tiers`.
  **Breakpoint refinement is not run on them: they keep the span of their evidence.**
- Main view only; `merge_view/` is unchanged. The pass reads the scan's results and the stage tables only and is not part of the
  identity of a saved scan piece: a run with the pass reuses the pieces of an earlier run with the same scan options
  (`EVITEA_PIECES_READONLY=1` works).
- `library_free_summary.json`, `weak_site_pass`: `single_copy_depth`, `e_source` ("added library-free calls"), `added_calls_used`, `median_e`, `a`, `allele_depth`, `cut`, `ran`,
  and, when it ran, `sites`, `uncovered_sites`, `passing_rule`, `dropped_by_copy_filter`, `rows_added`, `cover_rows`.

## 3. Candidate place minimum (`call`, `run`)

| Switch | Default | Effect |
|---|---|---|
| `--host-cluster-support N` | 2 | Host-mapped mates needed in a 500 bp cluster for a candidate pair (`incumbent/host_mate_gate.py`). |
| `--candidate-min-pairs N` | 2 | Candidate pairs needed for a candidate place (`incumbent/candidates.py`). |
| `--candidate-floor-depth-cut DCUT` | none | N applies only when D is at or under DCUT, else 2. |

The values used, the numbers of candidate pairs and places, and the wall time of the candidate step are in `run_summary.json`
(`engine_counts.depth_floors.candidate_floor`) and in the log. A value below 2 can raise the run time a lot at low depth.

## Checks

`tests/unit/test_depth_floors.py` and `tools/library_free_tests/test_depth_floors.py` (with `FREEZE_SRC` set to the `src` folder of
the freeze, it runs the freeze and the current code on the same sample and compares every output file).

## 4. Side of a host mate: `--mate-side-rule`, `--junction-side-rule` (opt-in, both default to today's rule)

Why. In the mate-anchored search (`incumbent/mate_rescue.py`, `mate_rescue_evidence`) a host mate is LEFT when it is a forward
read left of the middle of the place's own host mates, RIGHT when it is a reverse read right of it. At a place whose pairs all
come from one side, about half of the correct host mates lie past that middle and are dropped (`PAIR_GEOMETRY_INCONSISTENT`).
At low depth that leaves fewer than 2 pairs on the side, and the missing-side search never starts. The same middle is the test in
the two missing-side searches of the strict-pair step and in the junction groups.

`call --mate-side-rule strand`: the main loop first collects the local host anchors of a place (view and midpoint), then takes one
split point x for the place: only forward anchors, x = +infinity (all LEFT); only reverse anchors, x = -infinity (all RIGHT); both
strands, the x that maximises #(forward midpoint < x) + #(reverse midpoint > x). The score is constant between two neighbouring
distinct midpoints, so the gaps between them are the candidates; among the best ones the widest bounded gap wins and x is its
middle, then the one closest to the place's middle, then the smaller x. The two outer gaps have no middle: they are taken only when no
bounded gap reaches the best score (the smaller x first). Everything after the side is unchanged (TE-mate window, one side per
pair, two pairs per side, ledger, calls).

`call --junction-side-rule split`: the targeted missing-side tests (`mate_rescue.py` strict pairs, local and global; `junction_groups.py`
host half) compare the host midpoint with a per-locus point in place of the middle of the candidate envelope. Both strands among the
place's host anchors: the same x as `strand`. One strand: the present side's edge, the rightmost forward host midpoint (LEFT present)
or the leftmost reverse one (RIGHT present). A place without a recorded point keeps the middle. The scope window stays centred on the
middle.

The two options are independent. Under `junction-side-rule split` alone the main loop still sorts the pairs by the middle; it only
records the points (`LocusSplit` rows, taken out of the evidence list by `finalize_mate_rescue`, which then returns a
`MateRescueSplitResult` with `splits`). At the defaults the search is called as before and returns the plain result, so every output
and the checkpoint bytes are those of 18bae44.

Checkpoints. The candidate and split-junction steps do not depend on the rules; the mate-anchored search and the steps after it do.
With a rule on, those steps' checkpoints carry the rule in their identity (`StageStore.later`), so a work folder can keep the
`candidate_acquisition` and `split_junction` checkpoints and recompute the rest; a later checkpoint made under another rule is
refused. At the defaults every checkpoint has the one fingerprint it always had.

Records. With a rule on, `run_manifest.json` holds `depth_floors.side_rules` and, in the incumbent counts, the pairs given a side
that `midpoint` would have dropped (`mate_side_strand_rule_anchors_added_one_strand_loci`, `..._both_strand_loci`, and the pairs lost
the other way), plus the pairs and proposals each side test rejects (`side_test_rejects_*`: main loop, local and global missing-side
geometry and not-missing-side, junction side unresolved and not missing side). `tools/mate_side_counts.py <work dir>` prints the
`side_test_rejects_*` numbers of any finished run, so the two rules can be set side by side. `tools/mate_side_fixture.py` is a
synthetic fixture (no bowtie2): run with the src folder of 18bae44 and the current one it gives the same hashes at the defaults.

## 5. Floor-1 add-on: `--floor1-addon` (`library-free`, `run`; opt-in, default off)

A run with the weak-site pass (section 2b) uses the extension floor of the main view (`--extension-stop-min`, section 1). With
`--floor1-addon` the library-free step also runs the same scan at floor 1 in its own folder and adds the library-free rows of that run
that lie away from every row of the main combined table.

| Switch | Default | Effect |
|---|---|---|
| `--floor1-addon` | off | When the weak-site pass ran (allele depth at or under the cut), run a second library-free step at `--extension-stop-min 1`, no weak-site pass, no add-on, and append its far-away library-free rows. |

- Needs `--weak-site-pass` (with `--weak-site-allele-depth-cut` under `run`) and an msBWT index; without them the run stops at the start with a
  message. Under `run` it is refused with `--no-library`.
- Not run when the allele depth is over the cut: `library_free_summary.json` records `floor1_addon.ran = false` and every output file is the file a run
  without the option writes.
- Rule (the desk rule of `scripts/union_desk.py`, setting U1): a library-free row of the floor-1 run is added when, on its chromosome, it is more
  than 500 bp from every row of the main `calls.combined.tsv` (of either source). A gap of exactly 500 bp counts as near. Add-on rows are not
  compared to each other. The added rows are named `F1_<event id>`, source `LIBRARY_FREE`, and come after the other rows in `calls.combined.tsv`,
  `calls.library_free.tsv` and `calls.combined.vcf`. The VCF keeps its sort order and gains the `##ALT` lines of any new family.
- The second step writes to `<output>/floor1_addon/` (kept) and works in `<work>/floor1_addon/`; its scan pieces are reused when the scan options match.
  The code is in `library_free/floor1_addon.py`, outside the hashed detection files, so the detection hash is that of 18bae44.
- `library_free_summary.json`, `floor1_addon`: `requested`, `ran`, `allele_depth`, `cut`, `addon_output`, `addon_library_free_rows`, `rows_added`, `rule`.
  `output_hashes` are those of the final three files.
