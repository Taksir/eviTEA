"""``evitea run``: one command from FASTQ to calls, with a TE library or without one.

It runs the separate commands in their usual order, with the same arguments, and writes each step's output to its own
folder under ``--output``.  A step whose folder holds a completion marker (``steps/<step>.done``) is not run again, so an
interrupted run continues with ``evitea run`` and the same arguments.  ``run_steps.json`` lists the command line of every
step.

With a library (``--te-library`` and ``--family-map``)::

    reference/  prepare-reference          te/      prepare-te --build-languages
    sample/     prepare-sample             call/    call
    library_free/   library-free (with --library-free)

Without a library (``--no-library``; the chain of ``tools/no_library.sh`` and the Slurm ``tools/slurm/no_library.sh``)::

    reference/  prepare-reference          sample/  prepare-sample
    repeat_map/     build-repeat-map (BLAST+ on PATH)

``--repeat-annotation`` is optional here.  Without it, ``build-repeat-map`` runs first, on the reference FASTA, and its BED6
copy (``repeat_map/repeat_map_bed6.bed``: every repeated interval named ``repeat``) is the reference's annotation for
``prepare-reference`` and so for the record pass's ``call``.  The run then takes only the reference and the reads.  A
self-built record never shares a name with an annotation entry, so the same-family tests that read the annotation (the
paired-end channel's residency test, the reference check of the record codebooks) behave the same with the map as with a
curated annotation; what the map changes is the near-any-TE field of ``call-tiers`` and the segment classes of
``--refine-breakpoints``, which read any annotated interval as a TE.  An explicit ``--repeat-annotation`` wins.
    nolib_stub/     an empty library call set (the header of calls.native.tsv and the sample id)
    library_free/   library-free on the stub, family-identity filter on, resident-copy test with the repeat map
    records/lib/, records/assets_te/   build-family-records (flag rule kept, record guard with the map), prepare-te
    record_pass/call/   call on the record assets       record_pass/added/   add-record-calls

The final calls are ``call/`` (with a library), ``library_free/calls.combined.tsv`` (with ``--library-free``) and
``record_pass/added/calls.combined.tsv`` (without a library).
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Callable, Sequence
from pathlib import Path

from .errors import InputContractError
from .incumbent.events import LEGACY_EVENT_FIELDS
from .repeat_map import BED6

__all__ = ["add_arguments", "library_steps", "no_library_steps", "run", "write_stub_run"]

# the switch of the family-identity filter, as the Slurm runs pass it
FAMILY_IDENTITY_FILTER = "--rule-c"


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reference", type=Path, required=True, help="reference FASTA")
    parser.add_argument("--repeat-annotation", type=Path,
                        help="TE annotation of the reference (BED6; see docs/REPEAT_ANNOTATION.md); required with a library; "
                        "with --no-library the repeat map stands in for it when it is not given")
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--sample-id")
    parser.add_argument("--output", type=Path, required=True, help="run folder (created; a step already done is kept)")
    library = parser.add_mutually_exclusive_group(required=True)
    library.add_argument("--te-library", type=Path, help="TE library FASTA (with --family-map)")
    library.add_argument("--no-library", action="store_true",
                         help="run without a TE library: library-free calls turned into records, then a library pass on them")
    parser.add_argument("--family-map", type=Path, help="TE record to family map (with --te-library)")
    parser.add_argument("--library-free", action="store_true", help="with a library: also run the library-free layer")
    parser.add_argument("--contigs", help="contigs of the library-free layer and the repeat map (default: contigs >= 1 Mb)")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--workers", type=int)
    for name in ("--read-length", "--fragment-q95", "--fragment-q99", "--primary-contigs"):
        parser.add_argument(name, help=f"passed to call ({name})")


def _call_options(args: argparse.Namespace) -> list[str]:
    out = ["--threads", str(args.threads)]
    if args.workers is not None:
        out += ["--workers", str(args.workers)]
    for name in ("read_length", "fragment_q95", "fragment_q99", "primary_contigs"):
        value = getattr(args, name)
        if value is not None:
            out += ["--" + name.replace("_", "-"), str(value)]
    return out


def _common_steps(args: argparse.Namespace, out: Path, sample_id: str) -> list[tuple[str, list[str]]]:
    annotation = args.repeat_annotation if args.repeat_annotation is not None else out / "repeat_map" / BED6
    return [
        ("reference", ["prepare-reference", "--reference", str(args.reference), "--repeat-annotation",
                       str(annotation), "--output", str(out / "reference"), "--threads", str(args.threads)]),
        ("sample", ["prepare-sample", "--read1", str(args.read1), "--read2", str(args.read2), "--sample-id", sample_id,
                    "--output", str(out / "sample"), "--threads", str(args.threads)]),
    ]


def _contigs(args: argparse.Namespace) -> list[str]:
    return ["--contigs", args.contigs] if args.contigs else []


def library_steps(args: argparse.Namespace, out: Path, sample_id: str) -> list[tuple[str, list[str]]]:
    reference, sample = _common_steps(args, out, sample_id)
    steps = [
        reference,
        ("te", ["prepare-te", "--te-library", str(args.te_library), "--family-map", str(args.family_map),
                "--build-languages", "--reference-assets", str(out / "reference"), "--output", str(out / "te")]),
        sample,
        ("call", ["call", "--sample-assets", str(out / "sample"), "--reference-assets", str(out / "reference"),
                  "--te-assets", str(out / "te"), "--output", str(out / "call"), "--work-dir", str(out / "work" / "call"),
                  *_call_options(args)]),
    ]
    if args.library_free:
        steps.append(("library_free", [
            "library-free", "--run-dir", str(out / "call"), "--sample-assets", str(out / "sample"), "--reference-assets",
            str(out / "reference"), "--output", str(out / "library_free"), "--work-dir", str(out / "work" / "library_free"),
            *_contigs(args)]))
    return steps


def no_library_steps(args: argparse.Namespace, out: Path, sample_id: str) -> list[tuple[str, list[str]]]:
    reference, sample = _common_steps(args, out, sample_id)
    repeat_map = out / "repeat_map" / "repeat_map.bed"
    records = out / "records"
    second = out / "library_free" / "library_free_second_segments.fa"
    if args.repeat_annotation is None:
        # the map first, from the FASTA: its BED6 copy is the reference's annotation
        first = [("repeat_map", ["build-repeat-map", "--reference-fasta", str(args.reference), "--output",
                                 str(out / "repeat_map"), *_contigs(args)]), reference, sample]
    else:
        first = [reference, sample, ("repeat_map", ["build-repeat-map", "--reference-assets", str(out / "reference"),
                                                    "--output", str(out / "repeat_map"), *_contigs(args)])]
    return [
        *first,
        ("nolib_stub", ["@stub", str(out / "nolib_stub"), sample_id]),
        ("library_free", ["library-free", "--run-dir", str(out / "nolib_stub"), "--sample-assets", str(out / "sample"),
                          "--reference-assets", str(out / "reference"), "--output", str(out / "library_free"),
                          "--work-dir", str(out / "work" / "library_free"), *_contigs(args), FAMILY_IDENTITY_FILTER,
                          "--repeat-map", str(repeat_map)]),
        ("records", ["build-family-records", "--family-calls", str(out / "library_free" / "library_free_family_calls.tsv"),
                     "--segments", str(out / "library_free" / "library_free_segments.fa"), "@second-segments", str(second),
                     "--reference-assets", str(out / "reference"), "--output", str(records / "lib"), "--flag-rule", "kept",
                     "--repeat-map", str(repeat_map)]),
        ("records_te", ["prepare-te", "--te-library", str(records / "lib" / "te_library.fa"), "--family-map",
                        str(records / "lib" / "family_map_padded.tsv"), "--terminal-codebook",
                        str(records / "lib" / "terminal_codebook.tsv"), "--internal-codebook",
                        str(records / "lib" / "internal_only_codebook.tsv"), "--reference-assets", str(out / "reference"),
                        "--output", str(records / "assets_te")]),
        ("record_pass_call", ["call", "--sample-assets", str(out / "sample"), "--reference-assets", str(out / "reference"),
                              "--te-assets", str(records / "assets_te"), "--output", str(out / "record_pass" / "call"),
                              "--work-dir", str(out / "work" / "record_pass"), *_call_options(args)]),
        ("record_pass_added", ["add-record-calls", "--calls", str(out / "record_pass" / "call" / "calls.native.tsv"),
                               "--base-calls", str(out / "library_free" / "calls.combined.tsv"), "--records",
                               str(records / "lib"), "--output", str(out / "record_pass" / "added")]),
    ]


def write_stub_run(folder: Path, sample_id: str) -> None:
    """An empty library call set: the header of calls.native.tsv and the sample id (what library-free reads)."""

    folder.mkdir(parents=True, exist_ok=True)
    (folder / "calls.native.tsv").write_text("\t".join(LEGACY_EVENT_FIELDS) + "\n", encoding="utf-8", newline="\n")
    (folder / "run_summary.json").write_text(json.dumps({"sample_id": sample_id}), encoding="utf-8")


def _resolve(words: list[str]) -> list[str]:
    """Drop ``@second-segments <path>`` when the file is absent, else pass it as ``--second-segments``."""

    out, i = [], 0
    while i < len(words):
        if words[i] == "@second-segments":
            if Path(words[i + 1]).is_file():
                out += ["--second-segments", words[i + 1]]
            i += 2
        else:
            out.append(words[i])
            i += 1
    return out


def run(args: argparse.Namespace, main: Callable[[Sequence[str]], int], say: Callable[[str], None]) -> int:
    from .pipeline import derive_sample_id

    if args.no_library:
        if args.family_map is not None or args.library_free:
            raise InputContractError("--no-library takes no --family-map or --library-free")
    elif args.family_map is None:
        raise InputContractError("--te-library needs --family-map")
    elif args.repeat_annotation is None:
        raise InputContractError("a run with a TE library needs --repeat-annotation (see docs/REPEAT_ANNOTATION.md)")
    out = Path(args.output).expanduser().resolve()
    sample_id = derive_sample_id(args.read1, args.sample_id)
    steps = no_library_steps(args, out, sample_id) if args.no_library else library_steps(args, out, sample_id)
    marks = out / "steps"
    marks.mkdir(parents=True, exist_ok=True)
    (out / "run_steps.json").write_text(json.dumps([{"step": name, "command": ["evitea", *words]} for name, words in steps],
                                                   indent=2), encoding="utf-8")
    for name, words in steps:
        done = marks / f"{name}.done"
        if done.is_file():
            say(f"run: {name} already done")
            continue
        say(f"run: {name}")
        if words[0] == "@stub":
            write_stub_run(Path(words[1]), words[2])
        else:
            code = main(_resolve(words))
            if code != 0:
                say(f"run: step {name} failed (exit code {code}); fix the cause and run the same command again")
                return code
        done.write_text("done\n", encoding="utf-8")
    final = (out / "record_pass" / "added" / "calls.combined.tsv" if args.no_library
             else out / "library_free" / "calls.combined.tsv" if args.library_free else out / "call" / "calls.tsv")
    say(f"run: done; calls in {final}")
    return 0
