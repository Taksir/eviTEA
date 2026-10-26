"""Family records: library records built from the library-free layer's own families, and the second library pass on them.

A library-free family that the library layer does not call (two or more calls away from every library call, more than near
one) is a family the library lacks. ``build-family-records`` writes one library record per such family (the consensus of its
member junction segments, oriented and voted), the family map, the reference copies of each record, and the TE-end and TE-interior k-mer
codebooks, in the files ``prepare-te`` takes. The library layer then runs once more on those assets (``evitea call``).
``add-record-calls`` adds its calls where the existing call set has none within 500 bp and drops calls inside a reference copy
of their own record.
"""

from .build import (  # noqa: F401
    Record,
    build_record,
    join_ends,
    word_copies,
)
from .select import FlagRule, flagged_families  # noqa: F401
