"""AnchorPE evidence producer (S20-S70 for the bilateral/asymmetric layers).

A typed port of the frozen, self-contained Route-B v3 AnchorPE transfer
procedure (``route_b_v3_anchor_pe_iso1_transfer_freeze.py``; also reused
unchanged for GIM-024).  Run unchanged it reproduces the frozen AnchorPE
candidate tables of both development samples, ISO-1 and AKA-017 (the latter
with the AKA historical Bowtie2 host read-name profile).

Stages::

    final safe terminal codebook -> exact msBWT projection -> physical-pair
    classification -> whole-read TE-record fragment decoder (relaxed novelty)
    -> opposite-host Bowtie2 -k 2 -> LOCAL paired reference Bowtie2 -k 20
    -> evidence-complete reference ledger -> 500-bp host-midpoint clusters
    (old primary / relaxed provisional / corrected unexplained) -> RepeatMasker
    residency -> novel+nonresident targets -> Channel-A exact soft-clip anchors
    -> Channel-B NO_LOCAL_PROPER + UNIQUE_HIGH q99 geometry -> candidates.
"""

__all__ = ["AnchorPEInputs", "AnchorPEResult", "HOST_READ_NAME_PROFILES", "derive_candidates", "run_anchor_pe"]


def __getattr__(name: str):
    # Lazy: the decoder worker must import without the msBWT/NumPy stack.
    if name in __all__:
        from . import engine

        return getattr(engine, name)
    raise AttributeError(name)
