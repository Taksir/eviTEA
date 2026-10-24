"""``build-family-records`` and ``add-record-calls``."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from ..errors import InputContractError
from ..library_free import detector as D
from ..library_free.resident import CopyIndex
from . import build as B
from . import select as S

SCAN_ID = {"first": ("FBte99%05d", "LF%s"), "second": ("FBte97%05d", "LFh%s")}
PAD_ID, PAD_NAME = "FBte98%05d", "PAD%03d"
NEAR_BP = 500
MAP_HEADER = "FBte_ID\tFlyBase_Family\tStructural_Class\tSubtype\tSource_Accession\tFeature_Type\n"


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _segments(path):
    if path is None:
        return {}
    return D.parse_segments(Path(path).read_bytes())


def _members(rows, seg, by_pos):
    out = []
    for r in rows:
        d = seg.get((r["chrom"], str(r["bL"]), r["cluster"])) or by_pos.get((r["chrom"], str(r["bL"])))
        if d is not None:
            out.append((d.get("lseg", ""), d.get("rseg", "")))
    return out


def _by_pos(seg):
    out = {}
    for (c, b, _cl), d in seg.items():
        out.setdefault((c, b), d)
    return out


def _fasta_total(path):
    return sum(len(s) for _, s in B._reference_chunks(path))


def _ids(scan, cluster):
    fmt, name = SCAN_ID[scan]
    return fmt % int(cluster), name % cluster


def padded_map(records):
    """The family map followed by inert rows up to PAD_TO entries (the anchor-pair stage needs that many entries; prepare-te keeps
    only rows of library records, so the assets still hold the records alone)."""
    body = "".join("%s\t%s\t\t\t\t\n" % (te, fam) for te, fam, _ in records)
    n = len(records)
    pad = "".join("%s\t%s\t\t\t\t\n" % (PAD_ID % k, PAD_NAME % k) for k in range(1, max(0, B.PAD_TO - n) + 1))
    return MAP_HEADER + body + pad


def build_family_records(*, family_calls, segments, second_segments=None, reference_fasta, output, flag_rule="kept",
                         record_guard=True, repeat_annotation=None, te_share=None, copy_report=False, progress=print):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    seg = {"first": _segments(segments), "second": _segments(second_segments)}
    pos = {k: _by_pos(v) for k, v in seg.items()}
    tables = [read_tsv(p) for p in family_calls]
    flagged: dict = {}
    for rows in tables:
        for key, v in S.flagged_families(rows, flag_rule).items():
            flagged.setdefault(key, v)
    progress(f"flagged families: {len(flagged)}")
    copies = None
    genome_bp = 0
    if repeat_annotation is not None:
        names = {h.split()[0] for h, _ in _fasta_headers(reference_fasta)}
        copies = CopyIndex.from_file(repeat_annotation, names)
        genome_bp = _fasta_total(reference_fasta)
        te_share = copies.share(genome_bp)
    threshold = S.guard_threshold(te_share)
    info, built = {}, {}
    for (scan, cl), v in flagged.items():
        te, name = _ids(scan, cl)
        inside = S.inside_share(v["rows"], copies)
        row = dict(record_id=te, name=name, scan=scan, cluster=cl, flagged_by=",".join(v["flagged_by"]),
                   away_kept=v["kept"][0], near_kept=v["kept"][1], away_before=v["before-filters"][0],
                   near_before=v["before-filters"][1], calls=len(v["rows"]),
                   family_inside="" if inside is None else f"{inside:.4f}", status="")
        if record_guard and inside is not None and inside > threshold:
            row["status"] = "skipped_guard"
            info[(scan, cl)] = row
            continue
        if scan == "second" and not seg["second"]:
            raise InputContractError("a haplotype-scan family is flagged: give --second-segments")
        if not seg[scan]:
            raise InputContractError("flagged families need --segments")
        mem = _members(v["rows"], seg[scan], pos[scan])
        rec = B.build_record(mem)
        if rec is None:
            row["status"] = "no_record"
        else:
            row.update(status="built", usable=rec.members, placed=rec.placed, a_len=rec.a_len, b_len=rec.b_len,
                       joined=rec.joined, length=len(rec.seq))
            built[(scan, cl)] = (te, name, rec)
        info[(scan, cl)] = row
    recs = [(te, name, rec.seq) for te, name, rec in built.values()]          # in the order the families were flagged
    ref = str(reference_fasta)
    (out / "te_library.fa").write_text("".join(">%s\n%s\n" % (te, s) for te, _, s in recs))
    (out / "family_map.tsv").write_text(MAP_HEADER + "".join("%s\t%s\t\t\t\t\n" % (te, fam) for te, fam, _ in recs))
    (out / "family_map_padded.tsv").write_text(padded_map(recs))
    progress(f"records: {len(recs)}; finding reference copies")
    seqs = {te: s for te, _, s in recs}
    copy_rows = B.word_copies(seqs, ref) if recs else []
    (out / "annotation.bed").write_text("".join("%s\t%d\t%d\t%s\n" % row for row in sorted(copy_rows)))
    if recs:
        tb, ib = B.languages(out / "te_library.fa", out / "family_map.tsv", out / "annotation.bed", ref)
    else:
        tb = ib = b""
    (out / "terminal_codebook.tsv").write_bytes(tb)
    (out / "internal_only_codebook.tsv").write_bytes(ib)
    term = Counter(r["te_id"] for r in csv.DictReader(tb.decode().splitlines(), delimiter="\t")) if tb else Counter()
    intl = Counter(r["te_id"] for r in csv.DictReader(ib.decode().splitlines(), delimiter="\t")) if ib else Counter()
    for row in info.values():
        row["terminal_words"], row["internal_words"] = term.get(row["record_id"], ""), intl.get(row["record_id"], "")
    cols = ["record_id", "name", "scan", "cluster", "flagged_by", "away_kept", "near_kept", "away_before", "near_before", "calls",
            "family_inside", "status", "usable", "placed", "a_len", "b_len", "joined", "length", "terminal_words", "internal_words"]
    order = list(info.values())
    (out / "records_report.tsv").write_text(
        "\t".join(cols) + "\n" + "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in order))
    if copy_report:
        cnt = B.count_record_words(seqs, ref) if recs else {}
        own = Counter(row[3] for row in copy_rows)
        ccols = ["record_id", "name", "length", "reference_absent_share", "reference_copies", "family_inside"]
        lines = ["\t".join(ccols)]
        for te, name, s in recs:
            r = next(x for x in info.values() if x["record_id"] == te)
            lines.append("\t".join(map(str, (te, name, len(s), B.reference_absent_share(s, cnt), own.get(te, 0), r["family_inside"]))))
        (out / "copy_report.tsv").write_text("\n".join(lines) + "\n")
    summary = dict(flag_rule=flag_rule, flagged=len(flagged), built=len(recs),
                   skipped_guard=sorted(r["name"] for r in info.values() if r["status"] == "skipped_guard"),
                   no_record=sorted(r["name"] for r in info.values() if r["status"] == "no_record"),
                   record_guard=dict(on=bool(record_guard), te_share=None if te_share is None else round(te_share, 4),
                                     threshold=round(threshold, 4),
                                     evaluated=sum(r["family_inside"] != "" for r in info.values())))
    (out / "records_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    lines = [f"records {len(recs)}"]
    names = ["te_library.fa", "family_map.tsv", "family_map_padded.tsv", "annotation.bed", "terminal_codebook.tsv",
             "internal_only_codebook.tsv"]
    lines += ["%s %s" % (f, B.sha12(out / f)) for f in names]
    lines += ["%s %s %s/%s" % (te, fam, term.get(te, 0), intl.get(te, 0)) for te, fam, _ in recs]
    (out / "sha12.txt").write_text("\n".join(lines) + "\n")
    progress("\n".join(lines))
    return summary


def _fasta_headers(path):
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].strip(), None


def _pos(r):
    c = r.get("chrom") or r.get("chr")
    return c, int(r["start"]), int(r["end"])


def _near(a, b, tol=NEAR_BP):
    return a[0] == b[0] and a[1] - tol <= b[2] and b[1] - tol <= a[2]


class _Index:
    def __init__(self, items=()):
        self.d = defaultdict(list)
        for x in items:
            self.add(x)

    def add(self, x):
        self.d[(x[0], x[1] // 10000)].append(x)

    def hit(self, x, tol=NEAR_BP):
        for k in range((x[1] - tol) // 10000, (x[2] + tol) // 10000 + 1):
            if any(_near(x, y, tol) for y in self.d.get((x[0], k), ())):
                return True
        return False


def add_record_calls(*, calls, base_calls, annotation, family_map, output):
    """Calls of the second library pass (on the record assets) that are added to the existing call set: a call inside a reference
    copy of its own record is dropped; a call is added only where no call of the existing set lies within 500 bp; calls added
    within 500 bp of each other count once (the one nearest the start of the contig is kept)."""
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    fam_of = {r["FBte_ID"]: r["FlyBase_Family"] for r in read_tsv(family_map)}
    names = {f for f in fam_of.values() if not f.startswith("PAD")}
    copies = defaultdict(list)
    for line in Path(annotation).read_text().splitlines():
        if line.strip():
            c, s, e, te = line.split("\t")[:4]
            if te in fam_of:
                copies[(c, fam_of[te])].append((int(s), int(e)))
    rows = read_tsv(calls)
    pad = [r for r in rows if r.get("family", "").startswith("PAD")]
    mine = [r for r in rows if r.get("family") in names]

    def resident(r):
        c, s, e = _pos(r)
        return any(s <= b and e >= a for a, b in copies.get((c, r["family"]), ()))

    kept = [r for r in mine if not resident(r)]
    base_rows = read_tsv(base_calls)
    bi = _Index(_pos(r) for r in base_rows)
    away = [r for r in kept if not bi.hit(_pos(r))]

    def support(r):
        for k in ("score", "total_support"):
            if r.get(k) not in (None, ""):
                return float(r[k])
        return 0.0

    added, ai = [], _Index()
    for r in sorted(away, key=lambda r: (_pos(r)[0], _pos(r)[1], -support(r))):
        if not ai.hit(_pos(r)):
            added.append(r)
            ai.add(_pos(r))
    report = dict(record_calls=len(mine), inside_own_copy=len(mine) - len(kept), near_existing=len(kept) - len(away),
                  same_place=len(away) - len(added), added=len(added), padding_calls=len(pad),
                  added_by_record=dict(sorted(Counter(r["family"] for r in added).items())))
    cols = ["event_id", "chr", "start", "end", "family", "source"]
    body = [["RECORD_EVENT_%06d" % (i + 1), _pos(r)[0], r["start"], r["end"], r["family"], "RECORD"] for i, r in enumerate(added)]
    (out / "record_calls_added.tsv").write_text("\t".join(cols) + "\n" + "".join("\t".join(map(str, b)) + "\n" for b in body))
    head = list(base_rows[0]) if base_rows else cols
    comb = ["\t".join(head)] + ["\t".join(r.get(h, "") for h in head) for r in base_rows]
    comb += ["\t".join({"event_id": b[0], "chr": b[1], "start": b[2], "end": b[3], "family": b[4], "source": b[5]}.get(h, "") for h in head)
             for b in body]
    (out / "calls.combined.tsv").write_text("\n".join(comb) + "\n")
    (out / "record_calls_report.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    return report


def add_arguments_build(p):
    p.add_argument("--family-calls", type=Path, nargs="+", required=True,
                   help="library_free_family_calls.tsv of one or more runs (families flagged in any are built; the first run "
                        "that flags a family names it)")
    p.add_argument("--segments", type=Path, help="library_free_segments.fa of the first-pass families")
    p.add_argument("--second-segments", type=Path, help="library_free_second_segments.fa (haplotype-scan families)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--reference-assets", type=Path)
    src.add_argument("--reference-fasta", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--flag-rule", choices=S.FLAG_RULES, default="kept")
    p.add_argument("--record-guard", choices=("on", "off"), default="on")
    p.add_argument("--repeat-annotation", type=Path,
                   help="RepeatMasker .out or BED of the reference: measures how many of a family's calls lie inside old copies")
    p.add_argument("--te-share", type=float, help="the genome's TE share, when the annotation is not given")
    p.add_argument("--copy-report", choices=("on", "off"), default="off",
                   help="also write copy_report.tsv: for each record the share of its words absent from the reference and its "
                        "reference copies")


def add_arguments_add(p):
    p.add_argument("--calls", type=Path, required=True, help="calls.native.tsv of the library run on the record assets")
    p.add_argument("--base-calls", type=Path, required=True, help="calls.combined.tsv of the existing call set (library + library-free)")
    p.add_argument("--records", type=Path, required=True, help="the build output folder (annotation.bed, family_map.tsv)")
    p.add_argument("--output", type=Path, required=True)
