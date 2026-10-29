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
- Released under GPL-3.0-or-later.
