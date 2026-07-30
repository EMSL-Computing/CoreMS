"""Unit tests for corems.molecular_networking (fixture-based, no private data)."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from corems.molecular_id.search.database_interfaces import MSPInterface
from corems.molecular_networking import MolecularNetwork
from corems.molecular_networking.network_builder import MolecularNetwork as MNClass


MSP_PATH = Path.cwd() / "tests/tests_data/lcms/test_db.msp"

FE_KWARGS = {
    "normalize_intensity": True,
    "min_ms2_difference_in_da": 0.02,
    "max_ms2_tolerance_in_da": 0.01,
    "max_indexed_mz": 3000,
    "precursor_ions_removal_da": None,
    "noise_threshold": 0,
}


class MockSpectrum:
    """Minimal spectrum compatible with SimilarityEngine."""

    def __init__(self, mz_exp, abundance, name=None):
        self.mz_exp = np.asarray(mz_exp, dtype=float)
        self.abundance = np.asarray(abundance, dtype=float)
        self.name = name


def _mock_spectra_pair():
    """Two similar spectra and one dissimilar spectrum."""
    base_mz = np.array([100.0, 150.0, 200.0, 250.0, 280.0], dtype=float)
    base_ab = np.array([1.0, 0.8, 0.6, 0.4, 0.2], dtype=float)
    similar_mz = base_mz + 0.0001
    other_mz = np.array([120.0, 180.0, 220.0], dtype=float)
    other_ab = np.array([1.0, 0.5, 0.3], dtype=float)
    spectra = [
        MockSpectrum(base_mz, base_ab, name="a"),
        MockSpectrum(similar_mz, base_ab, name="b"),
        MockSpectrum(other_mz, other_ab, name="c"),
    ]
    ids = ["q0", "q1", "q2"]
    precursor_mzs = [300.0, 300.0, 250.0]
    return spectra, ids, precursor_mzs


@pytest.fixture
def msp_fe_lib():
    if not MSP_PATH.is_file():
        pytest.skip(f"MSP fixture missing: {MSP_PATH}")
    msp = MSPInterface(file_path=str(MSP_PATH))
    fe_lib, _meta = msp.get_metabolomics_spectra_library(
        polarity="negative",
        format="flashentropy",
        normalize=True,
        fe_kwargs=FE_KWARGS,
    )
    return fe_lib, msp


def test_query_vs_query_open_produces_edges():
    spectra, ids, precursor_mzs = _mock_spectra_pair()
    mn = MolecularNetwork(
        fe_lib=None,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
        use_parallel=False,
    )
    mn.run_query_vs_query_only(
        spectra, ids, query_precursor_mzs=precursor_mzs
    )
    edges = mn.get_network_edges(metric="entropy_similarity")
    assert isinstance(edges, list)
    assert len(edges) >= 1
    # similar pair q0–q1 should be present with high score
    pairs = {(min(a, b), max(a, b)): score for a, b, score in edges}
    assert ("q0", "q1") in pairs
    assert pairs[("q0", "q1")] > 0.5

    stats = mn.get_network_stats(metric="entropy_similarity")
    assert stats["n_nodes"] == 3
    assert stats["n_edges"] == len(edges)
    assert stats["threshold"] == 0.1


def test_query_vs_query_neutral_loss_smoke():
    spectra, ids, precursor_mzs = _mock_spectra_pair()
    mn = MolecularNetwork(
        fe_lib=None,
        search_type="neutral_loss",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.05, "cosine": 0.05},
        use_parallel=False,
    )
    mn.run_query_vs_query_only(
        spectra, ids, query_precursor_mzs=precursor_mzs
    )
    edges = mn.get_network_edges()
    assert isinstance(edges, list)
    stats = mn.get_network_stats()
    assert stats["n_nodes"] == 3


def test_query_vs_library_with_test_msp(msp_fe_lib):
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    assert len(df) >= 1

    # Build a noisy query from the first library row
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    assert peaks.ndim == 2 and peaks.shape[1] == 2
    rng = np.random.default_rng(0)
    mz = peaks[:, 0] + rng.normal(0, 0.0002, size=peaks.shape[0])
    ab = np.clip(peaks[:, 1] * (1 + rng.normal(0, 0.01, size=peaks.shape[0])), 0, None)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(mz, ab, name="mock_from_lib0")

    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.2, "cosine": 0.2},
        use_parallel=False,
    )
    mn.query_vs_library(
        [q],
        ["mock0"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=False,
    )
    edges = mn.get_network_edges(metric="entropy_similarity")
    assert isinstance(edges, list)
    # Expect at least one edge involving the query
    assert any("mock0" in (a, b) for a, b, _ in edges)
    stats = mn.get_network_stats()
    assert stats["n_nodes"] >= 2


def test_prepare_query_spectra_from_lcms_object():
    class FakeMS2:
        def __init__(self):
            self.mz_exp = np.array([100.0, 200.0])
            self.abundance = np.array([1.0, 0.5])

    class FakeMF:
        def __init__(self, mid, has_ms2=True):
            self.id = mid
            self.mz = 300.0 + mid
            self.best_ms2 = FakeMS2() if has_ms2 else None

    class FakeLCMS:
        def __init__(self):
            self.mass_features = {
                0: FakeMF(0, True),
                1: FakeMF(1, False),
                2: FakeMF(2, True),
            }

    lcms = FakeLCMS()
    spectra, ids, pmzs = MNClass.prepare_query_spectra_from_lcms_object(lcms)
    assert ids == ["0", "2"]
    assert len(spectra) == 2
    assert len(pmzs) == 2
    assert pmzs[0] == pytest.approx(300.0)

    spectra2, ids2, _ = MNClass.prepare_query_spectra_from_lcms_object(
        lcms, mf_ids={2}
    )
    assert ids2 == ["2"]
    assert len(spectra2) == 1


def test_fe_build_attrs_stored_on_index(msp_fe_lib):
    fe_lib, _ = msp_fe_lib
    # database_interfaces stores selected build_index kwargs for later retrieval
    assert hasattr(fe_lib, "_build_min_ms2_difference_in_da")
    assert fe_lib._build_min_ms2_difference_in_da == 0.02
