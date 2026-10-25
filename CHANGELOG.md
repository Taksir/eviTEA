# Changelog

## Unreleased

- Library-free layer: the upper bound of the local-depth window is a switch (`--depth-cap`; default none, the lower bound of 10
  reads stays; `--depth-cap 100` is the previous window).
- Library-free layer: anchored extension, walk-out, first-pass priority, a second scan for heterozygous sites (`--het-scan`),
  the one-site check and the resident-copy test, each a switch (`docs/LIBRARY_FREE_OPTIONS.md`). `--previous-behaviour` turns
  every new rule off.
- New commands `build-family-records` and `add-record-calls`: library records for the library-free families that the library lacks.
- New command `build-repeat-map` and switch `--repeat-map` (`library-free`, `build-family-records`): for a run whose library layer
  made no call, the resident-copy test and the record guard read a repeat map built from the reference's alignment to itself
  (BLAST+ on PATH). Runs with library calls are unchanged. Slurm files for the chain: `tools/slurm/no_library.sh`.
- Library-free layer: `--identity-min-comparisons` (default 3; 0 is the previous rule): the family-identity filter judges a family only on
  at least that many member comparisons, using all members when the frozen-extension members give fewer, and leaves a family with fewer
  comparisons to the one-site check and the resident-copy test. `--relative-depth-floor` (default off): the lower bound of the
  local-depth window follows the contig's own depth below 20 reads.
- Checks for the new code in `tools/library_free_tests/`.

## 0.1.0.dev0 (first public release, tagged v0.1)

- `evitea call` runs the full pipeline on one sample: the clean repeat-aware
  incumbent caller followed by four evidence-extension layers (bilateral and
  asymmetric AnchorPE, full-memory AnchorPE, INTERNAL21).
- `prepare-reference`, `prepare-te` (including `--build-languages`),
  `prepare-sample` and `validate` commands.
- Reproduces the frozen research callsets byte for byte on two *Drosophila
  melanogaster* development samples (AKA017: 838 calls, ISO1: 194 calls). This
  is a method-equivalence result, not an accuracy benchmark.
- Released under GPL-3.0-or-later.
