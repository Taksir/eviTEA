"""One place per break: duplicate groups whose chosen place shares its break with a chosen place ranked before it are set aside
before the grouping."""

from __future__ import annotations

import copy
import random

from evitea.library_free import detector as D
from evitea.library_free import priority as PR
from evitea.library_free import runtime as R
from evitea.library_free import same_break as SB
from evitea.library_free.options import Options

OPTS = Options().with_(recurring_sequence=False)
RNG = random.Random(11)
SEG = "".join(RNG.choice("ACGT") for _ in range(300))


def place(bL, stage="S4_pass", walked_from=None, two_sided=1, seg=SEG):
    p = dict(chrom="chr2", first=bL - 40, last=bL - 10, n_low=3, D=30.0, stage=stage, bL=bL)
    if stage in ("S4_pass", "S4_reject"):
        p.update(lseg=seg, rseg=seg[::-1], lstop="end", rstop="end", seg_keep="L,R", seg_reject="", two_sided=two_sided,
                 tsd=8 if two_sided else None)
    if walked_from is not None:
        p.update(walked=1, bL0=walked_from)
    return p


def test_walked_place_on_the_break_of_a_plain_place_is_set_aside_with_its_group():
    plain = place(98668328)
    walked = place(98668328, walked_from=98667650)
    near_walked = place(98667660, stage="S4_abstain")  # chained with the walked place's plain break
    places = [near_walked, walked, plain]
    assert SB.set_aside_same_break(places, OPTS) == 2
    assert [p["stage"] for p in places] == [SB.STAGE, SB.STAGE, "S4_pass"]
    grouped = PR.group_and_call(places, OPTS)
    assert grouped == [plain]
    assert D.parse_segments(D.segments_fa_bytes(grouped))  # one name per segment


def test_without_the_step_the_pair_carries_one_name():
    places = [place(98668328, walked_from=98667650), place(98668328)]
    grouped = PR.group_and_call(places, OPTS)
    assert len(grouped) == 2
    try:
        D.parse_segments(D.segments_fa_bytes(grouped))
    except ValueError as exc:
        assert "duplicate segment" in str(exc)
    else:
        raise AssertionError("expected a duplicate segment")


def test_walked_pass_stays_where_the_plain_place_abstained():
    plain = place(5000, stage="S4_abstain")
    walked = place(5000, walked_from=4000)
    places = [walked, plain]
    assert SB.set_aside_same_break(places, OPTS) == 1
    assert plain["stage"] == SB.STAGE and walked["stage"] == "S4_pass"


def test_plain_reject_keeps_its_break_against_a_walked_pass():
    plain = place(5000, stage="S4_reject")
    walked = place(5000, walked_from=4000)
    places = [plain, walked]
    assert SB.set_aside_same_break(places, OPTS) == 1
    assert walked["stage"] == SB.STAGE and plain["stage"] == "S4_reject"


def test_two_walked_places_on_one_break_keep_the_two_sided_one():
    one = place(5000, walked_from=3000, two_sided=0)
    two = place(5000, walked_from=4400)
    places = [one, two]
    assert SB.set_aside_same_break(places, OPTS) == 1
    assert one["stage"] == SB.STAGE and two["stage"] == "S4_pass"


def test_no_shared_break_changes_nothing():
    places = []
    for k in range(60):
        base = 10_000 * k
        places.append(place(base + RNG.randrange(0, 300), stage=RNG.choice(["S4_pass", "S4_reject", "S4_abstain"])))
        if RNG.random() < 0.5:  # a walked-out place chained with it: dedupe resolves the shared break itself
            places.append(place(base + RNG.randrange(0, 300), walked_from=base + RNG.randrange(0, 300)))
    before = copy.deepcopy(places)
    assert SB.set_aside_same_break(places, OPTS) == 0
    a = PR.group_and_call(places, OPTS)
    b = PR.group_and_call(before, OPTS)
    strip = lambda ps: [{k: v for k, v in p.items() if k != "_rep"} for p in ps]  # noqa: E731
    assert strip(places) == strip(before) and strip(a) == strip(b)


def test_off_without_walk_out_or_priority():
    for o in (OPTS.with_(walk_out=False), OPTS.with_(first_pass_priority=False)):
        places = [place(5000, walked_from=4000), place(5000)]
        assert SB.set_aside_same_break(places, o) == 0
        assert [p["stage"] for p in places] == ["S4_pass", "S4_pass"]


def test_set_aside_places_enter_no_duplicate_group_and_saved_pieces_do_not_depend_on_the_step():
    assert SB.STAGE not in PR.RANK
    assert "same_break.py" not in R.CODE_FILES
