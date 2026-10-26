"""Full-memory layer: exact junction evidence inside candidate regions, with paired-end anchor support.

Stages:

* ``finite_state``  finite host (-k 2) and LOCAL paired-reference (-k 20)
  alternatives with explicit reporting-cap censoring, and per-pair
  accounting;
* ``proposals``     candidate regions that confirm no event by themselves:
  q99 host-geometry regions from >= 2 pairs over >= 2 folds on the primary
  chromosomes, for the full data and for all ten 3-discovery / 2-held-out
  fold rotations;
* ``strict_a``      global exact junction evidence ("A") from finite
  reference alternatives, exact anchors, and the full-data event table;
* ``selection``     one anchor per proposal, chosen by exact junction evidence
  alone (proposal-only keys), with full-data q99 paired-end anchor support
  ("B"), giving ``shadow_layers.FullMemoryEvent``.
"""
