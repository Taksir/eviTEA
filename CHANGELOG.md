# Changelog

## [version] ([date])

- `evitea call` runs the full pipeline on one sample: the repeat-aware caller followed by four evidence-extension layers
  (bilateral, asymmetric and full-memory anchors, internal 21-mer).
- `prepare-reference`, `prepare-te` (including `--build-languages`), `prepare-sample` and `validate` commands.
- Library-free layer (`library-free`): the upper bound of the local-depth window is a switch (`--depth-cap`; default none).
  The lower bound of 10 reads stays. `--depth-cap 100` caps the window at 100.
- Library-free layer: anchored extension, walk-out, first-pass priority, a second scan for heterozygous sites (`--het-scan`),
  the one-site check and the resident-copy test, each a switch (`docs/LIBRARY_FREE_OPTIONS.md`). `--previous-behaviour` turns
  off every added rule.
- Library-free layer: `--identity-min-comparisons` (default 3; 0 judges on any comparisons). The family-identity filter judges a
  family only on at least that many member comparisons. It uses the first-pass priority members, or all members when those give
  fewer. A family with fewer comparisons goes to the one-site check and the resident-copy test.
  `--relative-depth-floor` (default on): the lower bound of the local-depth window follows the contig's own depth below 20 reads.
- `build-family-records` and `add-record-calls`: library records for the library-free families that the library lacks.
- `build-repeat-map` and the switch `--repeat-map` (`library-free`, `build-family-records`). For a run whose library layer made
  no call, the resident-copy test and the record guard read a repeat map built from the reference's alignment to itself
  (BLAST+ on PATH). Runs with library calls are unchanged. Slurm files for the chain: `tools/slurm/no_library.sh`.
- Library-free layer: one place per break. A walked-out place whose final break is the break of a place in another duplicate
  chain no longer gives a second grouped place at that break (stage `S4_same_break`); the family calls stopped on such a pair
  with a duplicate segment name.
- Tests for the library-free layer are in `tools/library_free_tests/`.
- Own-family copy filter, on by default for every run with a TE library: `call` (`calls.tsv`, `calls.bed`, `calls.vcf`) and
  `library-free` (`calls.combined.tsv`, `calls.combined.vcf`) leave out the library calls that touch a reference copy of their own
  family, read from the BED6 repeat annotation, and list them in `resident_copy_dropped.tsv`. `calls.native.tsv` stays the full
  table and `library-free` reads it, so its added calls are unchanged. `--keep-resident-copies` (`call`, `library-free`,
  `refine-breakpoints`, `run`) turns it off. `library-free` takes `--te-assets` for the family map. `run --no-library` is unchanged.
- Library-free layer: adapter test, on by default (`--adapter-test on|off`). Before the grouping, places whose inserted segment
  holds a TruSeq or Nextera adapter core (also with the junction up to 4 letters inside the adapter) are set aside with stage
  `S4_adapter` and form no family or family record. Saved scan pieces are reused with it on or off.
- Library-free layer: single-copy contigs (`--single-copy-contigs auto|all|none|LIST`, default auto). On a contig the sample
  carries in one copy (the X and Y of a male, a haploid genome), a call that would be added is not added (status
  `DROPPED_SINGLE_COPY`) when it is a haplotype-scan call or a first-pass call with the heterozygous genotype (e above 0.15).
  `auto` flags a contig at no more than 0.625 of the genome depth (sampled words that occur once in the reference, counted in
  the msBWT; `copy_profile.json` in the work directory) whose library-confirmed calls (at least 10) are at most 25%
  heterozygous. The family records count these calls as kept. Saved scan pieces are reused with any choice.
- Library-free layer: the family-identity filter is on by default in every route (`--rule-c on|off`; a bare `--rule-c` still means
  on), with a TSD skip (`--rule-c-tsd-skip on|off`, `--rule-c-tsd-share`, `--rule-c-tsd-min`; defaults on, 0.3, 5). A family that
  the filter would drop is kept when at least 5 of its trusted members (near a library call, or passing the one-site check) have a
  TSD of 2 to 20 bp at a share of 0.3 or more; it then goes on to the other statuses as if the filter were off. The filter
  candidates are sent to the one-site check as well (the tolerance still leaves them out). The summary lists every judged family
  with its TSD share. Saved scan pieces are reused with the skip on or off.
- Library-free layer: local heterozygous guard (`--het-guard on|off`, default on). After the single-copy rule, a call that would
  be added with e above 0.15 is not added (status `DROPPED_HOM_REGION`) when fewer than 25% of the 10 confirmed places nearest to
  it on the contig (the call itself counts) are heterozygous-looking. A contig with fewer than 5 confirmed places uses the share
  of the sample; a sample with fewer than 10 has no guard. The merge view takes the main view's confirmed places. Use
  `--het-guard off` for a screen for new insertions in an inbred line, where a recent private insertion is heterozygous. The
  family records count these calls as kept. Saved scan pieces are reused with it on or off. `evitea run` takes the same switch.
- Released under GPL-3.0-or-later.
- Library-free layer: faster scans that change no call. `--workers N` scans a contig on N processes (tiles of 250 kb), `--counter fast`
  (the default; `plain` is the basic search) searches the read index with fewer steps for the same counts, `--reference-cache DIR` keeps the reference word tables between jobs,
  and the read index and the minimap2 index are loaded once per job. None is part of the identity of a saved scan piece; the
  detection code is unchanged. `evitea run` gives the library-free step its CPUs (`--workers` from `--threads`).
  `tools/compare_scan.py` compares two scans.
