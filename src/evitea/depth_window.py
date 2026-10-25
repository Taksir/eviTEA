"""Local-depth window shared by the breakpoint walk and the library-free layer.

Both layers judge a site only when the median 31-mer count D around it is inside the window ``[D_MIN, upper]``.
``D_MIN`` (10 reads) is the minimum evidence.  Every other depth rule of both layers is relative to D
(0.65 D, 0.15 D, the anchor test's multiple of D), so by default the window has no absolute upper bound.  Repeat copies
that collapse onto one reference place are caught by the anchor test and the multi-copy word tests, not by the depth.
``--depth-cap N`` (``set_depth_cap(N)``) sets an absolute upper bound, for example ``--depth-cap 100``.

The cap is mirrored into the environment variable ``EVITEA_DEPTH_CAP`` so that worker processes inherit it.
"""

from __future__ import annotations

import math
import os
import sys

D_MIN = 10
ENV_VAR = "EVITEA_DEPTH_CAP"
MODULES = ("evitea.breakpoints.walker", "evitea.library_free._p4_lib")

__all__ = ["D_MIN", "ENV_VAR", "parse_depth_cap", "set_depth_cap", "upper_bound"]


def parse_depth_cap(text) -> float:
    """Return the upper bound for ``text``: ``None``, ``""`` and ``"none"`` give no bound (inf); a positive number gives that bound."""

    if text is None:
        return math.inf
    if isinstance(text, (int, float)):
        value = float(text)
    else:
        text = str(text).strip().lower()
        if text in ("", "none", "off", "inf"):
            return math.inf
        value = float(text)
    if not value >= D_MIN:
        raise ValueError(f"depth cap must be at least {D_MIN} (the lower bound of the window), got {text!r}")
    return value


def upper_bound() -> float:
    """Return the upper bound of the window from the environment (inf when unset)."""

    return parse_depth_cap(os.environ.get(ENV_VAR))


def set_depth_cap(cap) -> float:
    """Set the upper bound for this process and its workers, and return it.  Also set ``D_MAX`` in the loaded modules."""

    value = parse_depth_cap(cap)
    if math.isinf(value):
        os.environ.pop(ENV_VAR, None)
    else:
        os.environ[ENV_VAR] = repr(int(value) if value == int(value) else value)
    for name in MODULES:
        module = sys.modules.get(name)
        if module is not None:
            module.D_MAX = value
    return value
