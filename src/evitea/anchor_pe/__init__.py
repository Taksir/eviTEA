"""Paired-end anchor support evidence producer (``AnchorPE`` in names).

Feeds the bilateral and asymmetric layers.  Stages:

1. Exact msBWT projection of the final safe terminal codebook.
2. Physical-pair classification.
3. Whole-read TE-record fragment decoder (input of the relaxed novelty mask).
4. Opposite-host Bowtie2 ``-k 2``.
5. LOCAL paired-reference Bowtie2 ``-k 20``, giving the evidence-complete
   reference table.
6. 500-bp host-midpoint clusters (primary, relaxed provisional, corrected
   unexplained).
7. RepeatMasker residency, giving the novel, non-resident targets.
8. Exact soft-clip anchors ("A", ``channel_a``).
9. NO_LOCAL_PROPER + UNIQUE_HIGH q99 geometry ("B", ``channel_b_universe``).
10. Candidates.
"""

__all__ = ["AnchorPEInputs", "AnchorPEResult", "HOST_READ_NAME_PROFILES", "derive_candidates", "run_anchor_pe"]


def __getattr__(name: str):
    # Lazy: the decoder worker must import without the msBWT/NumPy stack.
    if name in __all__:
        from . import engine

        return getattr(engine, name)
    raise AttributeError(name)
