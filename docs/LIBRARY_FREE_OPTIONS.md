# Library-free layer: switches and outputs

Each added rule is a switch. `library-free --previous-behaviour` turns off every added rule and gives the basic library-free
scan (for comparison runs). With `--depth-cap 100` as well (the basic scan has the upper depth bound 100), the layer writes the
files of the basic scan byte for byte. `tools/library_free_tests/test_runtime.py` checks this on a synthetic sample.

## Switches

| Part | Switch | Default | What it does |
|---|---|---|---|
| Contigs | `--contigs a,b,c` | primary chromosomes and contigs of 1 Mb or more | the contigs scanned, in this order. The default takes every primary chromosome (at least a tenth of the longest contig, not an organelle, not an unplaced sequence) and every contig of at least 1 Mb, in reference order. It differs from the contigs of at least 1 Mb only where the longest contig is under 10 Mb: a yeast reference is scanned on all 16 chromosomes, and dm6 and mm10 on their contigs of 1 Mb or more. |
| Depth window | `--depth-cap N` | none | upper bound of the local depth D (the lower bound 10 stays). `--depth-cap 100` is the window of the basic scan. Also `EVITEA_DEPTH_CAP`. |
| Anchored extension | `--anchored-extension` | on | builds the inserted DNA from reads of one site only, where the plain extension stops short |
| Walk-out | `--walk-out` | on | moves a repeated walk-start word outward (up to 50 strides) to the nearest unique word |
| First-pass priority | `--first-pass-priority` | on | second-chance sources add calls and never change a first-pass duplicate choice, grouping or family-identity filter |
| Poly-A crossing | `--polya-cross` | off | the extension crosses a poly-A run (requires anchored extension) |
| Poly-A TSD | `--polya-tsd` | off | a poly-A tail counts toward the two TSDs of the family rule |
| Het scan | `--het-scan off\|safe-start\|second-walk\|full` | full | a second scan that follows both letters at a het SNP |
| | `--early-start`, `--right-fallback`, `--second-chance`, `--second-scan-anchored`, `--het-yield` | on | parts of the full mode; `--het-yield`: an early-start or right-side-alone place yields to the plain or second-chance walk of its duplicate chain |
| Identity evidence | `--identity-min-comparisons N` | 3 | the family-identity filter judges a family only on at least N member comparisons (each of at least 50 aligned columns): the first-pass priority members when they give N, otherwise all members. A family with fewer is not judged (its similarity is NaN, which the filter never drops). The one-site check and the resident-copy test decide. Applies to the first pass and to the haplotype scan's trusted members. 0 judges on any comparisons the members give (the previous rule). |
| Relative depth floor | `--relative-depth-floor on\|off` | on | the lower bound of the local-depth window is the smaller of 10 reads and half the median local depth of the contig (valid words); at a median of 20 or more it is 10. The upper bound is `--depth-cap`. The floor is part of the identity of a saved scan piece. |
| One-site check | `--check-scope none\|path-dropped\|unconfirmed` | unconfirmed | which calls must show one-site evidence: none (the rule of the basic scan), the calls the reference-path filter drops, or every call the library layer does not make |
| | `--tol p95\|p99` | p99 | depth tolerance of the check |
| | `--tol-source library\|own-calls` | library | tolerance from library-confirmed calls (own two-sided calls when there are none) |
| Resident-copy test | `--resident-guard off\|per-call\|family` | family | drops a call inside an annotated reference TE copy. `family`: only when the family's calls sit in copies more than twice the genome's TE share. Requires `--repeat-annotation` (RepeatMasker .out or BED) or the prepared reference's annotation; off without one. |
| | `--repeat-map FILE` | none | the repeat map of the reference (`build-repeat-map`): the resident-copy test reads it when the library layer made no call and no `--repeat-annotation` is given. With library calls it is not used, and the files are those of a run without it |
| Own-family copy filter | `--te-assets DIR` | none | prepared TE assets of the library run (the TE record to family map). With them, `calls.combined.tsv` and `calls.combined.vcf` leave out the library calls that touch a reference copy of their own family (the reference's BED6 annotation), and `resident_copy_dropped.tsv` lists them. The detection and the added calls read all library calls (`calls.native.tsv`), so they do not change. Required when the run has library calls, unless `--keep-resident-copies` is given; not needed on an empty library call set (`run --no-library`) |
| | `--keep-resident-copies` | off | keep the library calls inside a reference copy of their own family: the outputs are those of a run without the filter |
| Shared-word depth | `--shared-word-depth on\|off` | on | first pass: where the median word count around a run is below the depth floor, the run is judged with the median count of the words in the same window that the sample carries (count above zero) |
| Repeat cutoff | `--repeat-cutoff genome-share\|fixed` | genome-share | first pass: the repeat-content cutoff turns off the same share of this reference's windows as the fixed cutoff 0.5 does in the reference it was set on (never below 0.5; 0.5 on that reference) |
| Poly-A tails | `--polya-tails on\|off` | on | first pass: the extension crosses a poly-A run and a poly-A tail counts toward the two TSDs of the family rule (requires anchored extension) |
| Recurring sequence | `--recurring-sequence on\|off` | on | first pass: a place whose segment the filter rejected as a local or segmental copy is kept when the place has a TSD (or a poly-A tail) and the segment groups with a family that passes the TSD rule (requires the one-site check) |
| Merge view | `--merge-view on\|off` | on | also writes the files of the same run under the merge rules (the four rules above off) to `merge_view/` |
| Adapter test | `--adapter-test on\|off` | on | before the grouping, a place whose left or right segment holds a sequencing-adapter core is set aside (stage `S4_adapter`), so it forms no family and no family record. Cores: AGATCGGAAGAGC (TruSeq) and CTGTCTCTTATACACATCT (Nextera) with their reverse complements anywhere in a segment, and at the junction end also with the first 1 to 4 letters missing (the walk followed the reference into the adapter). A read that runs past the end of a short insert into the adapter is new sequence at single-copy depth, which the one-site check accepts. Applies to every pass and to `merge_view/`. Not part of the identity of a saved scan piece; not run with `--detection-dir` (its grouping stands). The detection summary counts the places it set aside |
| Family records | `build-family-records`, `add-record-calls` | commands | see below |
| Record guard | `--record-guard on\|off` | on | a family of old copies gets no record |
| Copy report | `--copy-report on\|off` | off | per record: fraction of its words absent from the reference, reference copies |

## Cross-species rules

The detection constants were set on the fly genome: few repeats, few copies per family, flanks close to the reference. The four
rules above keep each test and measure what it assumes, and act only where the merge rules give up on a place: a place they leave
alone is written as the merge writes it. They apply to the first pass; the haplotype scan keeps the merge rules.

- The first pass makes each place with the merge rules, then the same place under the cross-species rules. The candidate run is
  walked again only where a rule can change the place (depth below the floor, repeat content between 0.5 and the cutoff, a side
  rejected as short or satellite). A place a rule changed names the rule in the `cross_rule` column of the places table, and the
  summary counts places by rule and stage (`detection.cross_rules`).
- The merge view (`merge_view/`) holds the files of the merge rules from the same scan: its first-pass places are the merge's,
  and the haplotype scan is shared. With all four rules off (`Options.merge()`) the layer writes the merge's files byte for byte.
- `--repeat-cutoff genome-share` needs the window counts of the anchor reference (`cross.REPEAT_ANCHOR`: windows and windows at
  0.5, contigs of 1 Mb or more, measured once on EviTEA's prepared dm6). `evitea repeat-profile --reference-fasta <fa> --output
  <json>` measures them on any reference and, once the anchor is set, gives that reference's cutoff. A run keeps the profile of its
  reference in the work directory (`repeat_profile.json`) and writes the cutoff to the summary (`detection.repeat_cutoff`).

## What the one-site check does

A real junction is non-reference sequence read at single-copy depth. Each side of a call that the detection recorded is measured:
`T / D <= tolerance` (T = reads continuing the walk's context at the break, D = local depth) and the junction (30 host + first 30
inserted letters) has at most 10 of its 31-letter windows in the reference and no ungapped alignment at >= 95% identity to it.
A call passes when it has a side and every side passes. The detection records `chkL_T`, `chkL_junction`, `chkR_T`, `chkR_junction`
(extra columns of the places table, present only when `--check-scope` is not `none`).

The tolerance is the p95 or p99 of the side depth ratios of the first-pass family calls with a library call within 500 bp.

Call statuses in the family-calls table: `DROPPED_V2` (reference-path filter), `DROPPED_RULE_C` (family-identity filter), `NEAR_EVITEA`, `DROPPED_CHECK`, `DROPPED_RESIDENT_COPY`, `ADDED`.
Second-scan calls: `NOT_HET`, `DROPPED_V2`, `NEAR_FIRST_PASS`, `DROPPED_CHECK`, `DROPPED_RULE_C`, `NEAR_EVITEA`,
`DROPPED_RESIDENT_COPY`, `ADDED`. A second-scan call joins the first pass only where the first pass has no family call within 500 bp.

## Outputs of `library-free`

The files of the basic scan, unchanged columns, plus: `library_free_second_places.tsv`, `library_free_second_segments.fa` (with the het scan);
extra columns `walked`, `bL0` (walk-out), `origin` (second-scan places), the four `chk*` columns, `cross_rule` (first-pass places
with a cross-species rule on), the folder `merge_view/` (merge view), and in the family-calls table
`scan`, `origin`, `check`, `ratio_L`, `ratio_R`, `inside_copy`, `family_inside` (when a switch requires them).
Second-scan families are named `LF_cluster_h<N>` in `calls.library_free.tsv`.

Places set aside by the adapter test keep their row in the places table with stage `S4_adapter`.

Re-running with other grouping or assembly switches in the same `--work-dir` reuses the saved scan pieces. Set
`EVITEA_PIECES_READONLY=1` on such a run: it then refuses, instead of scanning again and replacing the saved pieces, if a
detection switch, the code or an input differs. `--detection-dir` rebuilds the family calls from saved tables (assembly switches only).

## Family records

A library-free family with at least 2 calls away from the library layer's calls, and more away than near, is a family the
library lacks.

```
evitea build-family-records --family-calls <library_free_family_calls.tsv> [...] --segments library_free_segments.fa \
    [--second-segments library_free_second_segments.fa] --reference-assets <ref> --output <dir> \
    [--flag-rule kept|before-filters|both] [--record-guard on|off] [--repeat-annotation <file>] [--copy-report on|off]
evitea prepare-te --te-library <dir>/te_library.fa --family-map <dir>/family_map_padded.tsv \
    --terminal-codebook <dir>/terminal_codebook.tsv --internal-codebook <dir>/internal_only_codebook.tsv --reference-assets <ref> --output <assets>
evitea call ... --te-assets <assets> ...                       # the library layer once more, on the records
evitea add-record-calls --calls <run>/calls.native.tsv --base-calls <existing>/calls.combined.tsv --records <dir> --output <out>
```

- Flag rule `kept` (default): away = calls the filters left (`ADDED`, `DROPPED_RESIDENT_COPY`); near = `NEAR_EVITEA`
  (`NEAR_FIRST_PASS` counts as neither). `before-filters`: away = calls not near a library call that the family-identity filter does not drop. `both`: either.
- Record names: first-pass families `FBte99<cluster>` / `LF<cluster>`; second-scan families `FBte97<cluster>` / `LFh<cluster>`.
  `family_map_padded.tsv` adds inert `FBte98` / `PAD` rows up to 100 entries (the anchor-pair stage requires that many); `prepare-te`
  retains only the rows of library records.
- Record guard: a family whose calls lie inside annotated copies more than `max(0.5, 2 x genome TE share)` of the time gets no
  record. The share is measured with `--repeat-annotation`, or read from the `family_inside` column with `--te-share`.
- `add-record-calls`: drops a call inside a reference copy of its own record (`annotation.bed`), adds a call only where no call of
  the existing set is within 500 bp, and retains one of calls added within 500 bp of each other (the one nearest the contig start).
  Writes `record_calls_added.tsv`, `calls.combined.tsv` (existing + added, source `RECORD`), `record_calls_report.json`.

## Runs without a library

A run whose library layer made no call has no curated TE annotation to read. Three parts replace it, all from the sequence of the
reference and the sample:

- **Repeat map.** `evitea build-repeat-map --reference-fasta <fa> --output <dir> [--contigs a,b,c] [--parts N]` aligns windows of
  500 bp (step 250) of the chosen contigs (default: the primary chromosomes and the contigs of 1 Mb or more) to the whole
  reference with BLAST+ (`makeblastdb` and `blastn` on PATH). A base is repeated when a hit of at least 100 aligned bases at 80%
  identity or more covers it, other than the window's own locus. The map is `<dir>/repeat_map.bed`; `repeat_map_summary.json`
  holds the constants, the counts and checksums. Build it once per reference. It does not depend on the sample.
- **Tolerance of the one-site check.** With no library call, the tolerance (99th percentile by default) comes from the layer's own
  two-sided calls. `--tol-source library` falls back to them on its own, so no switch is needed.
- **Resident-copy test and record guard.** `library-free --repeat-map <bed>` reads the map for the resident-copy test when the library
  layer made no call (an explicit `--repeat-annotation` wins; with library calls the map is not used and the run equals one without
  the option). `build-family-records --repeat-map <bed>` reads it for the record guard.

The sequence, with the library layer's calls left out (an `evitea call` output with the header line only, or a stub of one with the
sample id in `run_summary.json`):

```
evitea build-repeat-map --reference-fasta <fa> --output <map dir>
evitea library-free --run-dir <run without library calls> --msbwt-index <index> --reference-assets <ref> --rule-c \
    --repeat-map <map dir>/repeat_map.bed --output <lf dir>
evitea build-family-records --family-calls <lf dir>/library_free_family_calls.tsv --segments <lf dir>/library_free_segments.fa \
    [--second-segments <lf dir>/library_free_second_segments.fa] --reference-assets <ref> --flag-rule kept \
    --repeat-map <map dir>/repeat_map.bed --output <records dir>
evitea prepare-te ...                                          # as in Family records
evitea call ... --te-assets <assets> ...                       # the library layer on the family records
evitea add-record-calls --calls <run>/calls.native.tsv --base-calls <lf dir>/calls.combined.tsv --records <records dir> --output <out>
```

With no library call every family with at least 2 calls away from them is flagged, so every library-free family that survives the
filters gets a consensus record (unless the record guard drops it), and the library layer then calls the sites of those records.
`tools/slurm/no_library.sh` submits this chain on Slurm (see `tools/slurm/README.md`).

## Differences to know

- `--family-records` is not a switch of `library-free`: the second library pass is a separate `evitea call`, so the records are two commands.
- The one-site check reads the junction recorded at detection, so a call never has a missing junction.
- A run with anchored extension off cannot be rebuilt from saved tables; run a detection with it off.
- `RefWords` (walk-out) and `word_copies` (records) hold the reference in Python objects; they are sized for the fly genome. A
  mouse-size reference requires a lower-memory implementation.

## Tests (`tools/library_free_tests/`, run with `PYTHONPATH=src`)

`test_depth_window`, `test_depth_tool`, `test_assemble`, `test_records`, `test_runtime` (`PREVIOUS_SRC`: the src of the version to compare against, for U1; `BASELINE_SRC`: the src of the commit before the identity-evidence and relative-floor switches for U8; `MERGE_SRC`: the src of the merge, for U9), `test_repeat_map` (requires BLAST+ on PATH for its build checks), `test_slurm_tools`; `tools/library_free_tests/run_all.sh` runs them all.
