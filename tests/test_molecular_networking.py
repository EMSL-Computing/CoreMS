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
    # Library nodes use lib:<idx> IDs so they never collide with query mf_ids
    assert any(
        str(a).startswith("lib:") or str(b).startswith("lib:") for a, b, _ in edges
    ), "expected at least one query–library edge with lib: node IDs"
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


def test_prepare_query_spectra_from_lcms_collection():
    class FakeMS2:
        def __init__(self, seed):
            self.mz_exp = np.array([100.0 + seed, 200.0])
            self.abundance = np.array([1.0, 0.5])

    class FakeMF:
        def __init__(self, mid, mz, has_ms2=True):
            self.id = mid
            self.mz = mz
            self.best_ms2 = FakeMS2(mid) if has_ms2 else None

    class FakeSample:
        def __init__(self, mass_features):
            self.mass_features = mass_features

    class FakeCollection:
        def __init__(self):
            # sample 0: mf 10 (has MS2), sample 1: mf 20 (no MS2), sample 0: mf 11 (has MS2)
            self._samples = {
                0: FakeSample(
                    {
                        10: FakeMF(10, 301.0, True),
                        11: FakeMF(11, 302.0, True),
                    }
                ),
                1: FakeSample({20: FakeMF(20, 400.0, False)}),
            }
            self._reps = __import__("pandas").DataFrame(
                [
                    {
                        "cluster": 0,
                        "sample_id": 0,
                        "mf_id": 10,
                        "coll_mf_id": "0_10",
                        "has_ms2": True,
                        "intensity": 100.0,
                    },
                    {
                        "cluster": 1,
                        "sample_id": 1,
                        "mf_id": 20,
                        "coll_mf_id": "1_20",
                        "has_ms2": False,
                        "intensity": 50.0,
                    },
                    {
                        "cluster": 2,
                        "sample_id": 0,
                        "mf_id": 11,
                        "coll_mf_id": "0_11",
                        "has_ms2": True,
                        "intensity": 80.0,
                    },
                ]
            )

        def __getitem__(self, index):
            return self._samples[index]

        def get_representative_mass_features_for_all_clusters(
            self, representative_metric=None
        ):
            return self._reps.copy()

    coll = FakeCollection()
    spectra, ids, pmzs = MNClass.prepare_query_spectra_from_lcms_collection(coll)
    # cluster 1 skipped (no MS2)
    assert ids == ["0_10", "0_11"]
    assert len(spectra) == 2
    assert pmzs == pytest.approx([301.0, 302.0])

    spectra2, ids2, _ = MNClass.prepare_query_spectra_from_lcms_collection(
        coll, cluster_ids={2}
    )
    assert ids2 == ["0_11"]
    assert len(spectra2) == 1

    with pytest.raises(AttributeError):
        MNClass.prepare_query_spectra_from_lcms_collection(object())


def test_fe_build_attrs_stored_on_index(msp_fe_lib):
    fe_lib, _ = msp_fe_lib
    # database_interfaces stores selected build_index kwargs for later retrieval
    assert hasattr(fe_lib, "_build_min_ms2_difference_in_da")
    assert fe_lib._build_min_ms2_difference_in_da == 0.02
