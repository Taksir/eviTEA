"""Library-free insertion detector and its precision filters.

``evitea library-free`` runs on a finished ``evitea call`` output. It finds
insertions without a TE library. It looks for places where the reads leave the
reference, reads the inserted DNA from the msBWT, and groups places into
families by shared inserted DNA. It keeps the family calls of families with a
valid target-site duplication (TSD) and applies these filters:

* The reference-path filter drops a call when the fraction of reads that
  still follow the reference straight across the break is ``T_E`` (0.0405) or
  more. It applies to homozygous insertions only, so it is meant for inbred
  or selfing lines.
* The family-identity filter (``--rule-c``) drops a call when the median
  identity between the inserted DNA of its family's members is
  ``RULE_C_MAX`` (0.98309) or less. Such a family is old, diverged resident
  copies rather than a young active family.
* A surviving call is added only when no EviTEA call lies within
  ``NEAR_EVITEA_BP`` (500) bases.

``_p4_lib`` is the research module ``p4_lib.py`` (sha256 56a4a47d...),
vendored byte for byte; the tests check its hash. It runs on Linux/WSL only
(mappy and parasail).
"""

from __future__ import annotations

__all__ = ["NEAR_EVITEA_BP", "RULE_C_MAX", "T_E"]

T_E = 0.0405
RULE_C_MAX = 0.9830917645572818
NEAR_EVITEA_BP = 500
