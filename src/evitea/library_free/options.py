"""Switches of the library-free layer.  Every new rule is a switch; ``Options.previous()`` has them all off, and with them all off
the layer with ``--depth-cap 100`` (the previous version's window) gives the previous version's files byte for byte.  Names say what each part does.

Rules (see LIBRARY_FREE_OPTIONS.md):
  anchored extension   build the inserted DNA from reads of one site only, where the plain extension stops short
  walk-out             move a repeated walk-start word outward to the nearest unique word
  first-pass priority  second-chance sources add calls and never change a first-pass duplicate choice, grouping or the family-identity filter
  het scan             a second scan that follows both letters at a het SNP (modes safe-start, second-walk, full)
  one-site check       a library-free call is kept only if its junction is new sequence at single-copy depth
  resident-copy test   a library-free call inside an annotated reference TE copy is dropped when its family sits mostly in such copies
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field, replace

__all__ = ["Options", "add_arguments", "from_args"]

HET_MODES = ("off", "safe-start", "second-walk", "full")
CHECK_SCOPES = ("none", "path-dropped", "unconfirmed")
TOLERANCES = ("p95", "p99")
TOL_SOURCES = ("library", "own-calls")
GUARDS = ("off", "per-call", "family")


@dataclass(frozen=True)
class Options:
    anchored_extension: bool = True
    walk_out: bool = True
    first_pass_priority: bool = True
    polya_cross: bool = False
    polya_tsd: bool = False
    het_scan: str = "full"
    early_start: bool = True
    right_fallback: bool = True
    second_scan_anchored: bool = True
    second_chance: bool = True
    het_yield: bool = True
    check_scope: str = "unconfirmed"
    tol: str = "p99"
    tol_source: str = "library"
    resident_guard: str = "family"

    def __post_init__(self):
        for name, allowed in (("het_scan", HET_MODES), ("check_scope", CHECK_SCOPES), ("tol", TOLERANCES),
                              ("tol_source", TOL_SOURCES), ("resident_guard", GUARDS)):
            if getattr(self, name) not in allowed:
                raise ValueError(f"{name} must be one of {allowed}, got {getattr(self, name)!r}")
        if self.polya_cross and not self.anchored_extension:
            raise ValueError("--polya-cross needs --anchored-extension on")

    @classmethod
    def previous(cls) -> "Options":
        """Every new rule off: the previous version's behaviour."""
        return cls(anchored_extension=False, walk_out=False, first_pass_priority=False, polya_cross=False, polya_tsd=False,
                   het_scan="off", early_start=False, right_fallback=False, second_scan_anchored=False, second_chance=False,
                   het_yield=False,
                   check_scope="none", tol="p99", tol_source="library", resident_guard="off")

    @property
    def is_previous(self) -> bool:
        return self == Options.previous()

    @property
    def first_pass_changed(self) -> bool:
        """Does the first-pass detection differ from the previous version (needs the extra place columns)?"""
        return self.anchored_extension or self.walk_out or self.polya_cross

    @property
    def detection_identity(self) -> dict:
        """The switches that change what a detection writes (part of the identity of a saved scan piece)."""
        keys = ("anchored_extension", "walk_out", "polya_cross", "het_scan", "second_chance", "early_start", "right_fallback",
                "second_scan_anchored")
        return {**{k: getattr(self, k) for k in keys}, "record_check": self.check_scope != "none"}

    def as_dict(self) -> dict:
        return asdict(self)

    def with_(self, **kw) -> "Options":
        return replace(self, **kw)


def _onoff(text: str) -> bool:
    return {"on": True, "off": False}[text]


SWITCHES = (
    # (flag, dest, default, help)
    ("--anchored-extension", "anchored_extension", "on", "build the inserted DNA from reads of one site only where the plain extension stops short"),
    ("--walk-out", "walk_out", "on", "move a repeated walk-start word outward to the nearest unique word"),
    ("--first-pass-priority", "first_pass_priority", "on", "second-chance sources add calls and never change a first-pass duplicate choice, grouping or the family-identity filter"),
    ("--polya-cross", "polya_cross", "off", "cross a poly-A run in the extension (set aside; off)"),
    ("--polya-tsd", "polya_tsd", "off", "count a poly-A tail like a target-site duplication in the family rule (set aside; off)"),
    ("--early-start", "early_start", "on", "het scan (full): walk again from 30 letters further left"),
    ("--right-fallback", "right_fallback", "on", "het scan (full): walk the right side alone when the left side shows nothing"),
    ("--second-scan-anchored", "second_scan_anchored", "on", "use the anchored extension in the het scan"),
    ("--second-chance", "second_chance", "on", "het scan (second-walk and full): walk again, following the other letter, a run that ends with no left break"),
    ("--het-yield", "het_yield", "on", "an early-start or right-only place yields to another place of its duplicate group"),
)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("library-free rules (each on its own switch)")
    for flag, dest, default, text in SWITCHES:
        g.add_argument(flag, dest=dest, choices=("on", "off"), default=None, help=f"{text} (default {default})")
    g.add_argument("--het-scan", choices=HET_MODES, default=None,
                   help="second scan that follows both letters at a het SNP: off | safe-start (safe left start) | second-walk (+ second-chance walk) | "
                        "full (+ early start and right-side fallback) (default full)")
    g.add_argument("--check-scope", choices=CHECK_SCOPES, default=None,
                   help="one-site check: none | path-dropped (a second look at the calls the reference-path filter drops) | unconfirmed "
                        "(every library-free call that the library layer does not make) (default unconfirmed)")
    g.add_argument("--tol", choices=TOLERANCES, default=None, help="depth tolerance of the one-site check (default p99)")
    g.add_argument("--tol-source", choices=TOL_SOURCES, default=None,
                   help="library: from library-confirmed calls (the layer's own two-sided calls when there are none); own-calls: always from the layer's own two-sided calls "
                        "(default library)")
    g.add_argument("--resident-guard", choices=GUARDS, default=None,
                   help="resident-copy test: family | per-call | off (default family; off when no annotation is available)")
    g.add_argument("--repeat-annotation", default=None, help="reference TE annotation for the resident-copy test (default: the prepared reference's)")
    g.add_argument("--repeat-map", default=None,
                   help="repeat map of the reference (evitea build-repeat-map) for the resident-copy test of a run without a library: used when "
                        "the library layer made no call and no --repeat-annotation is given; ignored when it made calls")
    g.add_argument("--previous-behaviour", action="store_true",
                   help="every new rule off (the previous version's behaviour); other switches given on the command line still apply")


def from_args(args: argparse.Namespace) -> Options:
    base = Options.previous() if getattr(args, "previous_behaviour", False) else Options()
    kw = {}
    for _flag, dest, _default, _text in SWITCHES:
        v = getattr(args, dest, None)
        if v is not None:
            kw[dest] = _onoff(v)
    for dest in ("het_scan", "check_scope", "tol", "tol_source", "resident_guard"):
        v = getattr(args, dest, None)
        if v is not None:
            kw[dest] = v
    return replace(base, **kw)
