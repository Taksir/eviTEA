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
| `--weak-site-pass` | off | After the library-layer calls exist, call each evidence site with no call within 500 bp that has at least one pooled read on each side, or at least two split reads on one side. |
| `--weak-site-depth-cut DCUT` | required with the pass | The pass runs only when D is at or under DCUT. |

Evidence is grouped as in the desk replay of 3 Nov 2026 (`incumbent/weak_sites.py`): split events (one-sided left or right), the
left and right pair ids of the mate-rescue loci and the strict ids of the junction-group loci; items within 100 bp join one site;
per side the distinct fragment ids are pooled. The family of the new call is the family with the most supporting ids at the site;
its coordinates are the span of that family's items. The calls are named `WEAK_<n>`, have `evidence_type`
`weak_site_pooled_floor`, call origin `INCUMBENT_WEAK_SITE`, and tier 4 (`LOW_SUPPORT`) at any support. They go through the same
own-family copy filter and breakpoint refinement as every library call, and the extension layers read them as base calls.
D, DCUT, the sites seen and the rows added are in `run_summary.json` (`engine_counts.depth_floors.weak_site_pass`).

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
