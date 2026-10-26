"""Extension of the inserted DNA from a break: the plain extension, the anchored extension and the poly-A crossing.

``Extender`` has the signature of ``_p4_lib.extend`` (``extender(C, ctx, alt, n_alt, maxlen, fork) -> (segment, stop)``).

anchored   Where the plain extension stops short (the firewall would drop its segment by length or repeat content:
           under ``SEG_MIN`` letters or satellite), the extension is run again from the junction. At each letter, if the last
           K - 1 = 30 letters occur more than ``ANCHOR_RC`` x d0 times, where d0 is the count of the 30 host letters at the
           junction (the local depth at this junction), the next letter is chosen from the counts of the anchored word: the 30
           host letters, every inserted letter so far and the candidate letter. Only reads of this site contain that word.
           Otherwise the plain counts are used. Stop and fork rules are the plain ones. An anchored word cannot be longer
           than a read, so a re-extended segment reaches at most about one read length minus 30 letters past the junction
           while the anchored counts are in use. A re-extended segment's stop reason carries the prefix ``a_`` when the
           anchored counts were used for at least one letter; when they never were, the result equals the plain one.
cross      (requires anchored) When the extension would stop and its segment ends in a run of one letter at least ``LC_RUN``
           long, the letter after the run is counted over every run length from ``LC_RUN`` up (anchored words: host letters +
           segment before the run + run + letter). If the top letter has at least the extension's stop count and ``FORK`` of
           the summed counts, the segment takes the run at its most common length and that letter, then continues under the
           plain stop and fork rules. Stop reasons of such segments start with ``a_x_``.

Constants come from ``_p4_lib``. With both off the extender is ``_p4_lib.extend``.
"""

from __future__ import annotations

from . import _p4_lib as P

__all__ = ["Extender", "is_anchored"]

PLAIN_EXTEND = P.extend


def is_anchored(stop) -> bool:
    return str(stop or "").startswith("a_")


def _tail_run(s):
    if not s:
        return "", 0
    c = s[-1]
    return c, len(s) - len(s.rstrip(c))


def _cross(C, host, s, stop, fork):
    """(extended segment or None, record, state) for the letter after the trailing run, counted over every run length."""
    R, n = _tail_run(s)
    pre = s[:len(s) - n]
    sums = {x: 0 for x in "ACGT" if x != R}
    per = {x: {} for x in sums}
    m = P.LC_RUN
    while len(host) + len(pre) + m < P.EXT_MAX + P.K:
        base = host + pre + R * m
        if C.cnt(base) == 0:
            break
        for x in sums:
            v = C.cnt(base + x)
            if v:
                per[x][m] = v
                sums[x] += v
        m += 1
    tot = sum(sums.values())
    best = max(sorted(sums), key=lambda x: sums[x])
    rec = dict(run_letter=R, run_len=n, run_lengths_read=m - P.LC_RUN, counts=sums, letter=best)
    if sums[best] < stop or sums[best] < fork * tot:
        return None, rec, None
    mstar = max(sorted(per[best]), key=lambda k: (per[best][k], -abs(k - n), -k))
    rec.update(run_len_used=mstar, run_lens=sorted(per[best]))
    return pre + R * mstar + best, rec, dict(head=host + pre, R=R, lens=sorted(per[best]), post=best)


def _xcount(C, X, b):
    """Reads of this site that continue the crossed segment with letter b, summed over the supporting run lengths."""
    return sum(C.cnt(X["head"] + X["R"] * m + X["post"] + b) for m in X["lens"])


class Extender:
    def __init__(self, anchored: bool = False, cross: bool = False):
        if cross and not anchored:
            raise ValueError("the poly-A crossing requires the anchored extension")
        self.anchored, self.cross = bool(anchored), bool(cross)
        self.where: dict = {}
        self.stats = dict(calls=0, retried=0, anchored=0, anchored_letters=0, now_long_not_sat=0, cross_tried=0, cross_done=0)
        self.cross_log: list[dict] = []

    def __call__(self, C, ctx, alt, n_alt, maxlen=P.EXT_MAX, fork=P.FORK):
        if not self.anchored:
            return PLAIN_EXTEND(C, ctx, alt, n_alt, maxlen=maxlen, fork=fork)
        st = self.stats
        st["calls"] += 1
        seq, why = PLAIN_EXTEND(C, ctx, alt, n_alt, maxlen=maxlen, fork=fork)
        if len(seq) >= P.SEG_MIN and not P.satellite(seq):
            return seq, why
        st["retried"] += 1
        host = ctx[-(P.K - 1):]
        lim = P.ANCHOR_RC * C.cnt(host)
        stop = max(P.STOP_MIN, P.STOP_FRAC * n_alt)
        s, cur, used, reason, X, rec = alt, (ctx + alt)[-(P.K - 1):], 0, "max", None, None
        while len(s) < maxlen:
            if X is not None and len(X["post"]) >= P.K - 1:
                n = [C.cnt(cur + b) for b in "ACGT"]          # past the run: the plain counts
            elif C.cnt(cur) > lim:
                n = [C.cnt(host + s + b) if X is None else _xcount(C, X, b) for b in "ACGT"]
                used += 1
            else:
                n = [C.cnt(cur + b) for b in "ACGT"]
            bi = max(range(4), key=lambda q: n[q])
            T = sum(n)
            r = "support" if n[bi] < stop else ("fork" if n[bi] < fork * T else None)
            if r:
                if self.cross and X is None and rec is None and _tail_run(s)[1] >= P.LC_RUN:
                    st["cross_tried"] += 1
                    new, rec, state = _cross(C, host, s, stop, fork)
                    rec.update(self.where, len_before=len(s))
                    self.cross_log.append(rec)
                    if new is not None and len(new) <= maxlen:
                        st["cross_done"] += 1
                        s, X = new, state
                        cur = (ctx + s)[-(P.K - 1):]
                        continue
                reason = r
                break
            s += "ACGT"[bi]
            cur = cur[1:] + "ACGT"[bi]
            if X is not None:
                X["post"] += "ACGT"[bi]
                if _tail_run(s)[1] >= P.LC_RUN:
                    reason = "run"
                    break
        crossed = X is not None
        if rec is not None:
            rec.update(crossed=crossed, len_after=len(s), stop=reason)
        if not used and not crossed:
            return seq, why
        st["anchored"] += 1
        st["anchored_letters"] += used
        st["now_long_not_sat"] += len(s) >= P.SEG_MIN and not P.satellite(s)
        return s, ("a_x_" if crossed else "a_") + reason
