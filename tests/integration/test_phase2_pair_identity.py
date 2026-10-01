"""End-to-end Phase 2 contract from paired FASTQ through physical identity."""

import shutil
from pathlib import Path

from evitea.contracts import PairIdentityMode
from evitea.identity import assign_pair_identity, fold_id
from evitea.inputs import iter_paired_fastq

FIXTURES = Path(__file__).parents[1] / "fixtures"


def test_pair_ids_are_zero_based_ordinals_and_folds_wrap_modulo_five(tmp_path):
    read1 = tmp_path / "phase2_pairs_R1.fastq"
    read2 = tmp_path / "phase2_pairs_R2.fastq"
    shutil.copyfile(FIXTURES / "phase2_pairs_R1.fastq.txt", read1)
    shutil.copyfile(FIXTURES / "phase2_pairs_R2.fastq.txt", read2)
    raw_pairs = iter_paired_fastq(
        read1,
        read2,
    )
    physical_pairs = [
        assign_pair_identity(pair, ordinal0, sample_id="fixture-sample")
        for ordinal0, pair in enumerate(raw_pairs)
    ]

    assert len(physical_pairs) == 6
    assert [pair.raw_pair_ordinal0 for pair in physical_pairs] == list(range(6))
    assert [pair.pair_id for pair in physical_pairs] == list(range(6))
    assert [pair.fold_id for pair in physical_pairs] == [0, 1, 2, 3, 4, 0]
    assert [fold_id(pair.pair_id) for pair in physical_pairs] == [0, 1, 2, 3, 4, 0]
    assert all(pair.identity_mode is PairIdentityMode.EVITEA_ORDINAL0 for pair in physical_pairs)
    assert physical_pairs[1].pair_name == physical_pairs[5].pair_name == "repeat"
    assert physical_pairs[1].pair_id != physical_pairs[5].pair_id
