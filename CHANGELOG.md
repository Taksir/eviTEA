# Changelog

## Unreleased

- Library-free layer: the upper bound of the local-depth window is a switch (`--depth-cap`; default none, the lower bound of 10
  reads stays; `--depth-cap 100` is the previous window).
- Library-free layer: anchored extension, walk-out, first-pass priority, a second scan for heterozygous sites (`--het-scan`),
  the one-site check and the resident-copy test, each a switch (`docs/LIBRARY_FREE_OPTIONS.md`). `--previous-behaviour` turns
  every new rule off.
- New commands `build-family-records` and `add-record-calls`: library records for the library-free families that the library lacks.
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
