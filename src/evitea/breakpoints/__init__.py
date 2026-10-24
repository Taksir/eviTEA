"""Optional breakpoint refinement of wide (range) calls.

Off by default.  When the ``refinement.breakpoints`` switch is on, every call
whose native range is wider than ``MIN_WIDTH`` bases is walked base by base
over the sample msBWT from both sides; a TE junction found inside the range
narrows the reported call to it.  The native callset is never changed; the
refined positions go to ``calls.tsv``, ``calls.bed`` and ``breakpoints.tsv``.

The method is a port of the S3a/S3 research walker (J1 junction placement);
``tests/reference_s3`` holds the frozen research code it is tested against.
"""

from __future__ import annotations

__all__ = ["MIN_WIDTH"]

MIN_WIDTH = 50
