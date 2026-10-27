"""Places and segments tables of a pass. The base detector columns come first and unchanged. A switch adds its columns after them,
only when that switch is on, so a run with every switch off writes the base detector's bytes."""

from __future__ import annotations

import io
from collections.abc import Sequence

from . import detector as D
from .options import Options

__all__ = ["extra_columns", "places_tsv_bytes", "segments_fa_bytes"]

CHECK_COLUMNS = ["chkL_T", "chkL_junction", "chkR_T", "chkR_junction"]


def extra_columns(opts: Options, second: bool) -> list[str]:
    cols: list[str] = []
    if second:
        cols.append("origin")
    elif opts.walk_out:
        cols += ["walked", "bL0"]
    if opts.check_scope != "none":
        cols += CHECK_COLUMNS
    if not second and opts.any_cross:
        cols.append("cross_rule")
    return cols


def places_tsv_bytes(places: Sequence[dict], opts: Options, second: bool = False) -> bytes:
    extra = extra_columns(opts, second)
    if not extra:
        return D.places_tsv_bytes(places)
    out = io.StringIO(newline="")
    out.write("\t".join(D.COLS + ["lseg_len", "rseg_len"] + extra) + "\n")
    for p in places:
        if p["stage"] == "S2_through":
            continue
        out.write("\t".join("" if p.get(c) is None else str(p.get(c)) for c in D.COLS)
                  + f"\t{len(p.get('lseg', ''))}\t{len(p.get('rseg', ''))}\t"
                  + "\t".join("" if p.get(c) is None else str(p.get(c)) for c in extra) + "\n")
    return out.getvalue().encode("utf-8")


segments_fa_bytes = D.segments_fa_bytes
