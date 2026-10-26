"""Library-free insertion detector and its precision filters.

``evitea library-free`` runs on a finished ``evitea call`` output. It finds
insertions without a TE library (places where the reads leave the reference,
the inserted DNA read from the msBWT, places grouped into families by shared
inserted DNA) and keeps the family calls of families with a valid
target-site duplication (TSD):

* The reference-path filter drops a call when the fraction of reads that
  still follow the reference straight across the break is ``T_E`` (0.0405) or
  more (homozygous insertions only: intended for inbred or selfing lines).
* The family-identity filter (``--rule-c``) drops a call when the median
  identity between the inserted DNA of its family's members is
  ``RULE_C_MAX`` (0.98309) or less: old, diverged resident copies rather
  than a young active family.
* A surviving call is added only when no EviTEA call lies within
  ``NEAR_EVITEA_BP`` (500) bases.

``_p4_lib`` is the research module ``p4_lib.py`` (sha256 56a4a47d...),
vendored byte for byte; the tests check its hash. Linux/WSL only (mappy and
parasail).
"""

from __future__ import annotations

__all__ = ["NEAR_EVITEA_BP", "RULE_C_MAX", "T_E"]

T_E = 0.0405
RULE_C_MAX = 0.9830917645572818
NEAR_EVITEA_BP = 500
