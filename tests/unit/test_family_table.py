from __future__ import annotations

import csv
import gzip
import io
import random
from pathlib import Path

import pytest

from evitea.errors import InputContractError
from evitea.family_table import FLYBASE, OTHER, TWO_COLUMN, family_pairs, family_table_kind

FLYBASE_TABLE = (
    "FBte_ID\tFlyBase_Family\tStructural_Class\n"
    "FBte0000001\troo\tLTR\n"
    "FBte0000002\tjockey\tLINE\n"
    "# a comment\n"
    "\n"
    "FBte0000003\tjockey\tLINE\n"
)
PLAIN = "B1_Mm\tB1\nL1MdA_I\tL1\nTy1\tTy1_family\n"


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_both_formats_are_read(tmp_path: Path) -> None:
    flybase = _write(tmp_path, "fb.tsv", FLYBASE_TABLE)
    assert family_table_kind(flybase) == FLYBASE
    assert family_pairs(flybase) == [("FBte0000001", "roo"), ("FBte0000002", "jockey"), ("FBte0000003", "jockey")]
    plain = _write(tmp_path, "plain.tsv", PLAIN)
    assert family_table_kind(plain) == TWO_COLUMN
    assert family_pairs(plain) == [("B1_Mm", "B1"), ("L1MdA_I", "L1"), ("Ty1", "Ty1_family")]
    headed = _write(tmp_path, "headed.tsv", "te_record_id\tfamily\n" + PLAIN)
    assert family_pairs(headed) == family_pairs(plain)
    packed = tmp_path / "plain.tsv.gz"
    with gzip.open(packed, "wt", encoding="utf-8") as handle:
        handle.write(PLAIN)
    assert family_pairs(packed) == family_pairs(plain)


def test_other_layouts_and_empty_maps_are_refused(tmp_path: Path) -> None:
    other = _write(tmp_path, "other.tsv", "a\tb\tc\n1\t2\t3\n")
    assert family_table_kind(other) == OTHER
    with pytest.raises(InputContractError, match="neither a FlyBase table"):
        family_pairs(other)
    with pytest.raises(InputContractError, match="empty"):
        family_pairs(_write(tmp_path, "empty.tsv", "te_record_id\tfamily\n"))
    with pytest.raises(InputContractError, match="two tab-separated fields"):
        family_pairs(_write(tmp_path, "ragged.tsv", "B1_Mm\tB1\nTy1\n"))


def test_every_reader_accepts_a_two_column_map_with_plain_names(tmp_path: Path) -> None:
    from evitea.anchor_pe.clusters import family_map_rows, load_family_map, unmapped_records
    from evitea.breakpoints.runtime import _family_map as breakpoint_map
    from evitea.incumbent.anchor_census import load_te_family_map
    from evitea.incumbent.te_library import load_heuristic_family_map, load_record_family_map
    from evitea.inputs import load_family_map as input_map
    from evitea.languages.common import parse_family_map_gapfill, parse_family_map_rebuild
    from evitea.resident_filter import load_family_ids
    from evitea.tiers import _family_map as tier_map

    plain = _write(tmp_path, "plain.tsv", PLAIN)
    expected = {"B1_Mm": "B1", "L1MdA_I": "L1", "Ty1": "Ty1_family"}
    for reader in (load_te_family_map, load_record_family_map, load_heuristic_family_map, parse_family_map_rebuild,
                   parse_family_map_gapfill, breakpoint_map, tier_map):
        assert reader(plain) == expected, reader.__module__ + "." + reader.__name__
    assert set(input_map(plain).assignments) == set(expected)
    # three entries are enough now (the paired-end channel used to refuse maps of fewer than 100 entries)
    mapping, families = load_family_map(family_map_rows(plain))
    assert mapping == expected and families == {"B1", "L1", "Ty1_family"}
    assert unmapped_records(["B1_Mm", "Ty1", "IAPEz"], mapping) == 1
    ids = load_family_ids(plain)
    assert ids.of("L1") == {"L1MdA_I"} and ids.of("Ty1") == {"Ty1"} and ids.of("unknown") is None


def test_flybase_tables_read_as_before(tmp_path: Path) -> None:
    from evitea.anchor_pe.clusters import family_map_rows, load_family_map
    from evitea.incumbent.anchor_census import load_te_family_map
    from evitea.incumbent.te_library import load_record_family_map
    from evitea.languages.common import parse_family_map_rebuild

    flybase = _write(tmp_path, "fb.tsv", FLYBASE_TABLE)
    expected = {"FBte0000001": "roo", "FBte0000002": "jockey", "FBte0000003": "jockey"}
    assert load_te_family_map(flybase) == expected
    assert load_record_family_map(flybase) == expected
    assert parse_family_map_rebuild(flybase) == expected
    rows = family_map_rows(flybase)
    assert rows == list(csv.DictReader(io.StringIO(FLYBASE_TABLE), delimiter="\t"))
    assert load_family_map(rows)[0] == expected


def test_annotation_names_find_plain_records() -> None:
    from evitea.languages.common import annotation_record

    families = {"L1MdA_I": "L1", "FBte0000001": "roo"}
    assert annotation_record("L1MdA_I", families) == "L1MdA_I"
    assert annotation_record("L1MdA_I#LINE/L1", families) == "L1MdA_I"
    assert annotation_record("roo", families) is None


def _rnd(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_library_layer_on_a_plain_named_library(tmp_path: Path) -> None:
    """The anchor catalogue and the TE languages of a tiny library with ordinary record names and a two-column map."""

    from evitea.incumbent.anchor_census import build_anchor_catalogue, load_te_family_map, parse_te_consensus
    from evitea.languages import build_languages

    rng = random.Random(5)
    shared = _rnd(rng, 300)
    records = {
        "B1_Mm": ("B1", _rnd(rng, 700)),
        "B1_Mm2": ("B1", None),
        "L1MdA_I": ("L1", _rnd(rng, 900)),
        "Ty1": ("Ty1", _rnd(rng, 800)),
        "IAPEz": ("IAP", shared[:120] + _rnd(rng, 500) + shared[-120:]),
        "MERVL": ("ERVL", shared[:120] + _rnd(rng, 600) + shared[-120:]),
    }
    base = records["B1_Mm"][1]
    records["B1_Mm2"] = ("B1", "".join(rng.choice("ACGT") if rng.random() < 0.03 else b for b in base))
    library = tmp_path / "library.fa"
    library.write_text("".join(f">{name} some description\n{seq}\n" for name, (_, seq) in records.items()), encoding="ascii")
    family_map = tmp_path / "families.tsv"
    family_map.write_text("".join(f"{name}\t{family}\n" for name, (family, _) in records.items()), encoding="utf-8")
    genome = _rnd(rng, 20_000)
    copy = records["L1MdA_I"][1]
    genome = genome[:5_000] + copy + genome[5_000 + len(copy):]
    reference = tmp_path / "reference.fa"
    reference.write_text(f">chr1\n{genome}\n", encoding="ascii")
    annotation = tmp_path / "repeats.bed"
    annotation.write_text(f"chr1\t5000\t{5_000 + len(copy)}\tL1MdA_I\t0\t+\n", encoding="utf-8")

    catalogue = build_anchor_catalogue(parse_te_consensus(library), load_te_family_map(family_map))
    assert catalogue.family_by_te() == {name: family for name, (family, _) in records.items()}
    built = build_languages(te_fasta=library, family_map=family_map, reference_fasta=reference, repeatmasker_bed=annotation)
    rows = list(csv.DictReader(io.StringIO(built.terminal_codebook.decode()), delimiter="\t"))
    assert rows
    by_record = {}
    for row in rows:
        by_record.setdefault(row["te_id"], set()).add(row["family"])
    assert set(by_record) <= set(records)
    assert all(families == {records[name][0]} for name, families in by_record.items())
    # the records whose ends are shared by two families get no end word of their own from the 21-mer pass
    assert {"B1_Mm", "L1MdA_I", "Ty1"} <= set(by_record)
