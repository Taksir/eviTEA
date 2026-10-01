# Changelog

## 0.1.0.dev0 (first public release, "EviTEA v1")

- `evitea call` runs the full pipeline on one sample: the clean repeat-aware
  incumbent caller followed by four evidence-extension layers (bilateral and
  asymmetric AnchorPE, full-memory AnchorPE, INTERNAL21).
- `prepare-reference`, `prepare-te` (including `--build-languages`),
  `prepare-sample` and `validate` commands.
- Reproduces the frozen research callsets byte for byte on two *Drosophila
  melanogaster* development samples (AKA017: 838 calls, ISO1: 194 calls). This
  is a method-equivalence result, not an accuracy benchmark.
- Released under GPL-3.0-or-later.
