#!/usr/bin/env python3
"""Find wording problems in the repository text, and check that a wording edit changed nothing else.

usage:
  python docs/release/wording_scan.py scan  [--min-level LEVEL] [--summary] [PATH ...]
  python docs/release/wording_scan.py stats [--against REF] [PATH ...]
  python docs/release/wording_scan.py guard REF

scan   Lists comments, docstrings, messages, help text and prose lines that contain a word or phrase from
       docs/release/WORDING_STYLE_GUIDE.md. Identifiers are never scanned. Levels: MUST, FIX, LOOK.
stats  Counts the words in comments, docstrings and messages per file (with --against, next to the count at REF).
guard  Compares every file that differs from REF with its version at REF. It fails when anything other than
       comments, docstrings and message text changed, when a number disappeared from a comment, or when a file
       was added or removed. Exit status 0 means the edit is wording only (still read the REVIEW lines).

Only the Python standard library is used. Run it from the repository root.
"""

from __future__ import annotations

import argparse
import ast
import collections
import difflib
import io
import json
import re
import shlex
import subprocess
import sys
import tokenize
from pathlib import Path

PY = {".py"}
PROSE = {".md", ".txt", ".rst"}
SCRIPT = {".sh", ".bash"}
CONFIG = {".yml", ".yaml", ".toml", ".cfg", ".ini", ".cff"}
CONFIG_NAMES = {"Containerfile", "Dockerfile", "environment.yml"}
SKIP_PREFIXES = (
    "tests/fixtures/",
    "examples/synthetic/expected/",
    "docs/release/RELEASE_CLEANUP_",
    "docs/release/WORDING_STYLE_GUIDE.md",
    "docs/release/wording_scan.py",
)
SKIP_NAMES = {"LICENSE", "THIRD_PARTY_NOTICES"}
DIRECTIVE = re.compile(r"^#\s*(!|-\*-|noqa|type:|fmt:|pragma|pylint|ruff|isort|mypy|flake8|yapf|nosec|shellcheck|syntax=)", re.I)
LEVEL_ORDER = {"MUST": 0, "FIX": 1, "LOOK": 2}

# (level, group, regex, advice). Case-insensitive; every regex is applied to comment, docstring, message and prose text only.
P = []


def add(level, group, pattern, advice):
    P.append((level, group, re.compile(pattern, re.I), advice))


# 1. AI and stock words
add("MUST", "stock word", r"\b(delv\w*|intricat\w*|pivotal|meticulous\w*|foster\w*|bolster\w*|realm|garner\w*|interplay|nuanc\w*|groundbreaking|underscor\w*|showcas\w*|surpass\w*|boasts?)\b", "delete or use a plain verb")
add("MUST", "stock word", r"\b(robust\w*|seamless\w*|leverag\w*|comprehensive\w*|crucial\w*|streamlin\w*|utili[sz]\w*|holistic\w*|cutting-edge|state-of-the-art|sophisticated|powerful|elegant\w*|gracefully|under the hood|battle-tested|rock-solid|game-changer|paradigm|landscape|tapestry|testament|myriad|plethora)\b", "say what it does, with a number if there is one")
add("FIX", "stock phrase", r"\b(it(?:'s| is) worth noting|worth noting|note that|importantly|notably|essentially|fundamentally|in order to|serves as|acts as|plays a (?:crucial|key|vital) role|a (?:wide )?(?:variety|range) of|sanity[- ]check|honest(?:ly)?|simply|just works|out of the box)\b", "delete the filler, or use: to, is, a check")
add("LOOK", "stock word", r"\b(ensur(?:e|es|ed|ing)|vital|critical|key insight|handle[sd]?|handling|robustly|appropriately|carefully)\b", "name the condition or the check instead")
# 2. Project labels and process words
add("MUST", "internal label", r"\bfail(?:s|ed|ing)?[- ](?:closed|open)\b|\bfail-closed\b|\bfails? the run closed\b", "say what happens: the run stops with an error; the call is not reported")
add("MUST", "claim word", r"\b(?:new|true|novel) (?:transposable[- ]element |TE )?insertions?\b|\bfirst library-free\b|\bonly tool\b|\b(?:outperform\w*|superior|beats?)\b|\blevel with\b", "non-reference TE insertions; never a ranking word; a changed claim is reported to Taksir, not decided in the edit")
add("MUST", "internal label", r"\b(frozen|unfrozen|unfreeze|freeze|thaw\w*)\b", "say what is unchanged or fixed and why: default, built-in, unchanged, fixed value")
add("MUST", "internal label", r"\b(incumbent|oracles?|verdicts?|headline|ledger|adjudicat\w*|certif\w*|rescue[sd]?|rescuing|recall[- ]repair|mate rescue|rule[- ]c|detector v2|route[- ]b|strict[- ]a|exact a|channel [ab]|shadow(?: layers?)?)\b", "use the plain name in docs/NAMING.md and the guide's table; if none fits, describe what it does")
add("MUST", "internal label", r"\b(?:phase|step|round|stage|part|addendum|amendment|arm|look|case|rule|item|gate|panel|tier|batch)[ -](?:[0-9]+[a-z]?(?:-[ivx]+)?(?![\w])|(?-i:[A-Z])(?![\w]))", "numbered or lettered project step: describe what the step does")
add("MUST", "internal label", r"\b(dev-merged|dev-release-cleanup|release-cleanup|work order|handoff|hand-off|checkpoint[- ]?[0-9]+|(?=[0-9a-f]*[0-9])(?=[0-9a-f]*[a-f])[0-9a-f]{7,40})\b", "remove branch names, commit hashes and work-order references")
add("MUST", "internal label", r"\b(?:R[0-9][A-Z0-9]*|S[0-9]{2,3}|G[0-9][A-Z]|H[0-9]|V[0-9]+[A-Z]*|v[0-9]|[ABCEL][0-9]{1,2}[a-z]?)\b(?![-_.a-z])", "stage or case code: describe the stage by what it does (allowed: Q95, Q99, F1, TE family names)")
add("FIX", "internal label", r"\b(production|legacy|historical|research[- ]stage|dev(?:elopment)? (?:data|strain|set)|reproduction[- ]only|incumbent path|decision|strategy|memo|thread|session|prompt)\b", "production -> default; describe the profile, drop the project history")
add("FIX", "internal label", r"\b(canonical)\b", "keep only for a strand-canonical k-mer or probe; otherwise main, standard, usual")
add("MUST", "AI or process word", r"\b(claude|anthropic|openai|chat-?gpt|gpt-?[0-9]\w*|codex|copilot|opus|sonnet|fable|haiku|llm|ai[- ]generated|as an ai|co-authored-by|generated (?:by|with))\b", "remove; wording about the writer does not belong in code")
add("MUST", "person name", r"\b(taksir|fernando|martin|mcmillan|holt|ferris)\b", "remove person names from comments; credits belong in CITATION.cff and THIRD_PARTY_NOTICES")
add("FIX", "dev data name", r"\b(AKA-?017|ISO1|GIM-?024|SLA-?001|MUN-?009|JUT-?011|TEN-?015|KIE-?094|LUN-?007|Canton-?S|CC0[0-9]{2}|A3|A4)\b", "keep only where it is an option value or a file name; elsewhere say 'the development sample'")
# 3. History and log narration
add("MUST", "history", r"\b(previously|originally|formerly|used to|no longer|earlier version|older version|old version|new version|was changed|were changed|has been changed|changed from|renamed from|moved from|superseded|carried over|backward[- ]compat\w*|turned out|we decided|decided to|after the fix|bug ?fix|regressed|hot-?fix|workaround|work-around|hack)\b", "describe what the code does now; delete what it did before")
add("FIX", "history", r"\b(now|new|newly|currently|already|still|finally|recently|at the time|for now|temporar\w*|TODO|FIXME|XXX)\b", "delete 'now' and 'currently' if the sentence reads the same without them; list TODO items in the report, do not delete them silently")
add("FIX", "history", r"\b(passed|failed|pass(?:es|ing)? (?:criteri\w*|line)|reached|acceptance criteri\w*|met the criteri\w*|verified (?:on|against|with)|checked (?:on|against|with)|confirmed (?:on|with)|measured on|tested on)\b", "a comment states what the code does, not whether a check passed")
add("MUST", "history", r"\b20[0-9]{2}-[0-9]{2}-[0-9]{2}\b|\b[0-9]{1,2}:[0-9]{2}(?::[0-9]{2})?Z\b|\bUTC\b|\b[0-9]{1,2} (?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\b", "no dates or clock times of work")
# 4. Other tools (competitors); dependencies such as Bowtie2, msBWT, mappy, parasail, RepeatMasker, BLAST are fine
add("MUST", "tool name", r"\b(temp2?|te-?forest|elite|insurveyor|transurveyor|genome-?delta|melt|xtea|tea|mcclintock|mobster|retroseq|tepid|teflon|ngs[_-]te[_-]mapper|popoolation\w*|tepeak|te-?peak|deviate|jitterbug|steak|pantera|dna-?pipete|repeatmodeler|tangram|tracker|t-lex\w*)\b(?<!\bTE-tracker)", "name the mechanism instead of the tool; dependencies and input formats stay (guide section 4)")
# 5. Register
add("FIX", "register", r"\b(simply|just|obviously|clearly|of course|basically|actually|really|very|quite|rather|pretty|fairly|a bit|kind of|sort of|in a sense)\b", "delete; add nothing in its place")
add("FIX", "register", r"[—–→⇒✓✗•…“”‘’]|\bi\.e\.|\be\.g\.|\betc\.", "plain ASCII: use a comma, a colon, 'to', 'for example', or end the sentence")
add("FIX", "register", r"\b(this (?:function|method|class|module|section) (?:is responsible for|is used to|provides|serves)|responsible for|is used to|allows? (?:us|you|the user) to|enables? (?:us|you|the user) to|aims to|seeks to|strives to|designed to|intended to)\b", "start with the verb: 'Reads', 'Returns', 'Drops'")
add("LOOK", "register", r";\s", "one idea per sentence; split at the semicolon")
add("LOOK", "register", r"\((?:[^()]*\([^()]*\)[^()]*)\)", "parenthesis inside a parenthesis: split the sentence")
# 6. Field words from FINAL_WRITING_RULE.md 3.1 that matter in code text
add("LOOK", "field word", r"\b(recovers?|recovered|recovery|finds)\b", "detects, reports, identifies")
add("LOOK", "field word", r"\b(unmatched|settings|cutoffs?|share of|comparator|supplied|confirmable|screen|rejects?|rebuilt|wrapper|proposals|by hand|whatever|favour\w*)\b", "false positives; parameters; threshold; fraction of; tool; provided; validated by long reads; filter; excludes; constructed; script; candidates; manually; regardless of; favor")
add("LOOK", "field word", r"\b(keeps?|carries|carrier|removes|lacks|stays|sits|lies|enters|needs|makes no call)\b", "retains, contains, excludes, is absent from, remains, is located in, is added to, requires, is not reported (docstring verbs Return, Raise, Yield stay)")


def git(*args, check=True):
    r = subprocess.run(["git", *args], capture_output=True, text=True, encoding="utf-8")
    if check and r.returncode:
        sys.exit(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return r.stdout


def tracked_files(paths):
    files = git("ls-files", "--", *(paths or ["."])).splitlines()
    out = []
    for f in files:
        name = Path(f).name
        if name in SKIP_NAMES or f.startswith(SKIP_PREFIXES):
            continue
        if kind_of(f):
            out.append(f)
    return out


def kind_of(path):
    p = Path(path)
    if p.suffix in PY:
        return "py"
    if p.suffix in PROSE:
        return "prose"
    if p.suffix in SCRIPT:
        return "script"
    if p.suffix in CONFIG or p.name in CONFIG_NAMES:
        return "config"
    return None


# ---------- text units: (line, kind, text) ----------

def is_message(text):
    words = re.findall(r"[A-Za-z]{2,}", text)
    return len(words) >= 3 and " " in text.strip()


def py_units(source):
    units = []
    try:
        tree = ast.parse(source)
    except SyntaxError as e:
        return [(e.lineno or 1, "error", f"cannot parse: {e.msg}")]
    doc_ids = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
                doc_ids.add(id(body[0].value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            kind = "docstring" if id(node) in doc_ids else "string"
            if kind == "string" and not is_message(node.value):
                continue
            for i, line in enumerate(node.value.splitlines() or [""]):
                units.append((node.lineno + i, kind, line))
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        if tok.type == tokenize.COMMENT and not DIRECTIVE.match(tok.string):
            units.append((tok.start[0], "comment", tok.string.lstrip("# ").rstrip()))
    return sorted(units)


def prose_units(text):
    units, fenced = [], False
    for n, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if s.startswith("```"):
            fenced = not fenced
            continue
        if fenced:
            if s.startswith("#") and not s.startswith("#!"):
                units.append((n, "comment", s.lstrip("# ")))
            continue
        if s.startswith("|") and set(s) <= set("|-: "):
            continue
        units.append((n, "prose", re.sub(r"`[^`]*`", "``", line)))
    return units


def script_units(text, config=False):
    units = []
    for n, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if s.startswith("#") and not DIRECTIVE.match(s):
            units.append((n, "comment", s.lstrip("# ")))
        elif re.search(r"\s#\s", line) and not DIRECTIVE.match(s):
            units.append((n, "comment", line.split(" # ", 1)[-1].strip()))
        if re.search(r"\b(echo|printf|die|usage|description|summary|name)\b", s) and not s.startswith("#"):
            units.append((n, "string", re.sub(r"`[^`]*`", "``", line)))
    return units


def units_of(path, text):
    k = kind_of(path)
    if k == "py":
        return py_units(text)
    if k == "prose":
        return prose_units(text)
    return script_units(text, config=(k == "config"))


# ---------- scan ----------

def scan(args):
    files = tracked_files(args.paths)
    floor = LEVEL_ORDER[args.min_level]
    hits, per_group, per_file = [], collections.Counter(), collections.Counter()
    for f in files:
        text = Path(f).read_text(encoding="utf-8", errors="replace")
        for line, kind, unit in units_of(f, text):
            if kind == "error":
                hits.append((f, line, "MUST", "parse", unit, unit, ""))
                continue
            for level, group, rx, advice in P:
                if LEVEL_ORDER[level] > floor:
                    continue
                for m in rx.finditer(unit):
                    if f == "docs/NAMING.md" and group in ("internal label", "AI or process word"):
                        continue
                    if group == "person name" and re.search(r"copyright|author:|licen[sc]e", unit, re.I):
                        continue
                    hits.append((f, line, level, group, m.group(0), unit.strip()[:110], advice))
    hits.sort(key=lambda h: (h[0], h[1], LEVEL_ORDER[h[2]]))
    for h in hits:
        per_group[(h[2], h[3])] += 1
        per_file[h[0]] += 1
    if args.json:
        print(json.dumps([dict(zip(("file", "line", "level", "group", "match", "text", "advice"), h)) for h in hits], indent=1))
        return 1 if any(h[2] == "MUST" for h in hits) else 0
    if not args.summary:
        for f, line, level, group, match, text, advice in hits:
            print(f"{f}:{line}: [{level} {group}] '{match}'  | {text}")
        print()
    print(f"{len(files)} files scanned, {len(hits)} hits")
    for (level, group), n in sorted(per_group.items(), key=lambda x: (LEVEL_ORDER[x[0][0]], -x[1])):
        print(f"  {level:4s} {group:20s} {n}")
    if args.summary:
        print("\nfiles with most hits:")
        for f, n in per_file.most_common(25):
            print(f"  {n:5d}  {f}")
    return 1 if any(h[2] == "MUST" for h in hits) else 0


# ---------- stats ----------

def word_counts(path, text):
    c = collections.Counter()
    for line, kind, unit in units_of(path, text):
        if kind != "error":
            c[kind] += len(unit.split())
    return c


def stats(args):
    rows = []
    for f in tracked_files(args.paths):
        now = word_counts(f, Path(f).read_text(encoding="utf-8", errors="replace"))
        before = None
        if args.against:
            old = git("show", f"{args.against}:{f}", check=False)
            before = word_counts(f, old) if old else collections.Counter()
        rows.append((f, now, before))
    tot_now, tot_before = collections.Counter(), collections.Counter()
    for f, now, before in rows:
        tot_now.update(now)
        if before is not None:
            tot_before.update(before)
        if args.against:
            if now != before:
                print(f"{f}: " + ", ".join(f"{k} {before[k]}->{now[k]}" for k in sorted(set(now) | set(before))))
        else:
            print(f"{f}: " + ", ".join(f"{k} {v}" for k, v in sorted(now.items())))
    print("\ntotal now   : " + ", ".join(f"{k} {v}" for k, v in sorted(tot_now.items())))
    if args.against:
        print("total before: " + ", ".join(f"{k} {v}" for k, v in sorted(tot_before.items())))
        grew = [k for k in tot_now if tot_now[k] > tot_before[k]]
        if grew:
            print("WARNING: more words than before in: " + ", ".join(grew) + " (the edit should not add words)")
    return 0


# ---------- guard ----------

MESSAGE_CALLS = {"print", "warn", "warning", "error", "info", "debug", "exception", "critical", "exit", "die", "fail", "skip", "echo", "write", "log"}
MESSAGE_KEYWORDS = {"help", "description", "epilog", "usage", "reason", "metavar", "message", "msg"}


class Normalise(ast.NodeTransformer):
    """Remove docstrings, and replace every string by a placeholder while recording it with its context."""

    def __init__(self):
        self.units = []
        self.ctx = []

    def visit(self, node):
        self.ctx.append(node)
        try:
            return super().visit(node)
        finally:
            self.ctx.pop()

    def _body(self, node):
        if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:] or [ast.Pass()]
        return node

    def visit_Module(self, node):
        self._body(node)
        return self.generic_visit(node)

    visit_ClassDef = visit_FunctionDef = visit_AsyncFunctionDef = lambda self, node: self.generic_visit(self._body(node))

    def context(self):
        for anc in reversed(self.ctx[:-1]):
            if isinstance(anc, (ast.Raise, ast.Assert)):
                return "raise"
            if isinstance(anc, ast.keyword) and anc.arg in MESSAGE_KEYWORDS:
                return "help"
            if isinstance(anc, ast.keyword) and anc.arg == "match":
                return "test-match"
            if isinstance(anc, ast.Call):
                f = anc.func
                name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
                if name in MESSAGE_CALLS or name == "error":
                    return "message"
            if isinstance(anc, ast.stmt):
                break
        return "other"

    def visit_Constant(self, node):
        if isinstance(node.value, str):
            self.units.append(("S", node.value, self.context()))
            return ast.copy_location(ast.Constant("§"), node)
        return node

    def visit_JoinedStr(self, node):
        template = "".join(v.value if isinstance(v, ast.Constant) else "{}" for v in node.values)
        self.units.append(("F", template, self.context()))
        keep = [self.visit(v) for v in node.values if isinstance(v, ast.FormattedValue)]
        return ast.copy_location(ast.JoinedStr(keep), node)


def normalise(source):
    tree = ast.parse(source)
    n = Normalise()
    tree = n.visit(tree)
    return ast.dump(tree, indent=1), n.units


NUMBER = re.compile(r"(?<![\w.])[0-9][0-9,._]*[0-9]%?|(?<![\w.])[0-9]%?")


def number_bag(path, text):
    bag = collections.Counter()
    for _, kind, unit in units_of(path, text):
        if kind != "error":
            bag.update(NUMBER.findall(unit))
    return bag


def code_lines(text, kind):
    """Lines that are not blank and not a whole-line comment; for scripts, comments after a command are cut off too."""
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s or (s.startswith("#") and not s.startswith("#!")):
            continue
        if kind == "script":
            try:
                lex = shlex.shlex(line, posix=True)
                lex.whitespace_split = True
                lex.commenters = "#"
                s = " ".join(lex)
            except ValueError:
                pass
        out.append(s)
    return out


def fenced_code(text):
    blocks, cur, fenced = [], [], False
    for line in text.splitlines():
        if line.strip().startswith("```"):
            if fenced:
                blocks.append([x for x in cur if x.strip() and not x.strip().startswith("#")])
                cur = []
            fenced = not fenced
        elif fenced:
            cur.append(line)
    return blocks


def guard(args):
    base = args.ref
    git("rev-parse", "--verify", base)
    status = git("diff", "--name-status", "--no-renames", base, "--").splitlines()
    status += ["A\t" + f for f in git("ls-files", "--others", "--exclude-standard").splitlines()]
    fail, review = [], []
    for row in status:
        flag, _, f = row.partition("\t")
        if f.startswith("docs/release/") and Path(f).name in ("WORDING_STYLE_GUIDE.md", "wording_scan.py", "WORDING_CLEANUP_REPORT.md"):
            continue
        if flag != "M":
            fail.append(f"{f}: file {'added' if flag == 'A' else 'removed' if flag == 'D' else flag}; file names must not change")
            continue
        k = kind_of(f)
        old = git("show", f"{base}:{f}")
        new = Path(f).read_text(encoding="utf-8")
        if k == "py":
            try:
                d_old, u_old = normalise(old)
            except SyntaxError as e:
                fail.append(f"{f}: cannot parse the old version: {e}")
                continue
            try:
                d_new, u_new = normalise(new)
            except SyntaxError as e:
                fail.append(f"{f}: new version does not parse: {e}")
                continue
            if d_old != d_new:
                diff = list(difflib.unified_diff(d_old.splitlines(), d_new.splitlines(), "old", "new", n=2, lineterm=""))
                fail.append(f"{f}: the code changed (not only comments, docstrings and messages):\n    " + "\n    ".join(diff[:24]))
            else:
                for (ko, to, co), (kn, tn, cn) in zip(u_old, u_new):
                    if to == tn:
                        continue
                    where = f"{f}: string in {co} context"
                    if co not in ("raise", "help", "message", "test-match") or not (is_message(to) and is_message(tn)):
                        fail.append(f"{where} changed but does not look like message text: {to!r} -> {tn!r}")
                    else:
                        review.append(f"{where}: {to[:80]!r} -> {tn[:80]!r}")
        else:
            if k == "prose":
                if fenced_code(old) != fenced_code(new):
                    fail.append(f"{f}: a fenced code block changed (commands and examples must stay as they are)")
                for span in set(re.findall(r"`([^`\n]+)`", old)) - set(re.findall(r"`([^`\n]+)`", new)):
                    review.append(f"{f}: code span no longer present: `{span}`")
            elif code_lines(old, k) != code_lines(new, k):
                a, b = code_lines(old, k), code_lines(new, k)
                diff = [x for x in difflib.unified_diff(a, b, "old", "new", n=0, lineterm="") if x[:1] in "+-" and x[:3] not in ("+++", "---")]
                review.append(f"{f}: lines that are not comments changed (allowed only for echo, printf and help text):\n    " + "\n    ".join(diff[:20]))
        lost = number_bag(f, old) - number_bag(f, new)
        for num, cnt in sorted(lost.items()):
            review.append(f"{f}: number {num} appears {cnt} time(s) less in comments, docstrings and messages (keep every value and unit; delete only with the sentence that was history)")
    print("GUARD against", base)
    for x in fail:
        print("FAIL  ", x)
    for x in review:
        print("REVIEW", x)
    print(f"\n{len(status)} changed files, {len(fail)} failures, {len(review)} items to review")
    return 1 if fail else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scan")
    s.add_argument("--min-level", choices=list(LEVEL_ORDER), default="LOOK")
    s.add_argument("--summary", action="store_true")
    s.add_argument("--json", action="store_true")
    s.add_argument("paths", nargs="*")
    t = sub.add_parser("stats")
    t.add_argument("--against")
    t.add_argument("paths", nargs="*")
    g = sub.add_parser("guard")
    g.add_argument("ref")
    args = ap.parse_args(argv)
    return {"scan": scan, "stats": stats, "guard": guard}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
