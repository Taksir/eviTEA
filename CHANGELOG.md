# Changelog

## Unreleased

- Library-free layer: the upper bound of the local-depth window is a switch (`--depth-cap`; default none). The lower bound of 10
  reads always applies. `--depth-cap 100` gives the previous window.
- Library-free layer: anchored extension, walk-out, first-pass priority, a second scan for heterozygous sites (`--het-scan`),
  the one-site check and the resident-copy test are each a switch (`docs/LIBRARY_FREE_OPTIONS.md`). `--previous-behaviour` turns
  off every rule added after the previous release.
- Added the commands `build-family-records` and `add-record-calls`: library records for the library-free families that the library does not have.
- Added the command `build-repeat-map` and the switch `--repeat-map` (`library-free`, `build-family-records`). For a run whose
  library layer made no call, the resident-copy test and the record guard read a repeat map built from the reference's
  alignment to itself (BLAST+ on PATH). Runs with library calls are unchanged. Slurm files for the chain: `tools/slurm/no_library.sh`.
- Library-free layer: `--identity-min-comparisons` (default 3; 0 is the previous rule). The family-identity filter judges a family
  only on at least that many member comparisons. It uses the members that the first-pass priority picks, or all members when
  those give fewer. A family with fewer comparisons is left to the one-site check and the resident-copy test.
  `--relative-depth-floor` (default off): the lower bound of the local-depth window follows the contig's own depth below 20 reads.
- Checks for these additions in `tools/library_free_tests/`.

## 0.1.0.dev0 (first public release, tagged v0.1)

- `evitea call` runs the full pipeline on one sample: the repeat-aware caller
  followed by four evidence-extension layers (bilateral and asymmetric
  paired-end anchor support, full-memory anchor support, internal 21-mer).
- `prepare-reference`, `prepare-te` (including `--build-languages`),
  `prepare-sample` and `validate` commands.
- Reproduces the call sets of the research pipeline byte for byte on two
  *Drosophila melanogaster* development samples (838 calls with
  `--profile AKA017`, 194 calls with `--profile ISO1`). This is a
  method-equivalence result, not an accuracy benchmark.
- Released under GPL-3.0-or-later.
