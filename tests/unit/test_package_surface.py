"""Package metadata smoke tests."""

import evitea


def test_public_metadata_surface() -> None:
    assert evitea.__version__ == "0.1.0.dev0"
    assert (
        evitea.ARCHITECTURE_ID
        == "EVITEA_V1_REPEAT_AWARE_EVIDENCE_ARCHITECTURE_20260928_REVISION_1"
    )
