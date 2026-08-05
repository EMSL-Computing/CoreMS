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


def test_plot_network_static_smoke(msp_fe_lib, tmp_path):
    """Static plot_network returns a Figure and can save a PNG."""
    pytest.importorskip("networkx")
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q0")

    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
    )
    mn.query_vs_library(
        [q],
        ["q0"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=False,
    )

    png = tmp_path / "network.png"
    fig = mn.plot_network(
        metric="entropy_similarity",
        path=str(png),
        return_fig=True,
        bypass_clustering=True,
        show_labels=True,
    )
    assert fig is not None
    assert hasattr(fig, "savefig")
    assert png.is_file() and png.stat().st_size > 0


def test_plot_interactive_network_smoke(msp_fe_lib, tmp_path):
    """Interactive HTML plot writes a file when ipysigma is installed."""
    pytest.importorskip("networkx")
    pytest.importorskip("ipysigma")
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q0")

    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
    )
    mn.query_vs_library(
        [q],
        ["q0"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=False,
    )
    html = tmp_path / "network.html"
    out = mn.plot_interactive_network(
        metric="entropy_similarity",
        out_path=str(html),
        bypass_clustering=True,
    )
    assert Path(out).is_file()
    assert Path(out).stat().st_size > 0


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


def test_library_vs_library_filtered_skips_empty_preserves_ids():
    """Empty mid-list library entries must not shift IDs/precursors of survivors."""
    from corems.molecular_networking.similarity_engine import SimilarityEngine

    class FakeFELib:
        def __init__(self, entries):
            self._entries = entries
            # SimilarityEngine may read ms2 tolerance from fe_lib.entropy_search
            self.entropy_search = SimpleNamespace(max_ms2_tolerance_in_da=0.01)

        def __getitem__(self, idx):
            return self._entries[idx]

        def __len__(self):
            return len(self._entries)

    peaks_a = np.array([[100.0, 1.0], [150.0, 0.8], [200.0, 0.5]], dtype=float)
    peaks_c = peaks_a + np.array([[0.0001, 0.0], [0.0001, 0.0], [0.0001, 0.0]])
    # Index 1 is empty → skipped; survivors should keep ids id0 and id2 (not id0, id1)
    entries = [
        {"peaks": peaks_a, "precursor_mz": 300.0, "spectra_id": "A"},
        {"peaks": np.empty((0, 2), dtype=float), "precursor_mz": 999.0, "spectra_id": "EMPTY"},
        {"peaks": peaks_c, "precursor_mz": 301.0, "spectra_id": "C"},
    ]
    fe = FakeFELib(entries)
    engine = SimilarityEngine(
        fe_lib=fe,
        search_type="open",
        additional_similarities=[],
        ms2_tolerance_da=0.01,
    )
    result = engine.compute_library_vs_library_filtered(
        library_indices=[0, 1, 2],
        spectrum_ids=["id0", "id1", "id2"],
        precursor_mzs=[300.0, 999.0, 301.0],
    )
    entropy_pairs = result.get("entropy_similarity", {})
    # Only the two valid spectra participate; pair keys must use id0/id2
    assert entropy_pairs, "expected at least one pair among non-empty library spectra"
    for (a, b), score in entropy_pairs.items():
        assert {a, b} == {"id0", "id2"}, f"unexpected pair ids {(a, b)} (id1 is empty)"
        assert score > 0.5


def test_drop_queries_resets_stage_flags_and_allows_rerun():
    spectra, ids, precursor_mzs = _mock_spectra_pair()
    mn = MolecularNetwork(
        fe_lib=None,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
    )
    mn.run_query_vs_query_only(spectra, ids, query_precursor_mzs=precursor_mzs)
    assert mn.stage_query_query_done is True
    assert mn._has_queries_run is True
    assert len(mn.get_network_edges()) >= 1

    mn.drop_queries()
    assert mn.stage_query_query_done is False
    assert mn.stage_query_library_done is False
    assert mn.stage_library_library_done is False
    assert mn._has_queries_run is False
    assert mn._stage2_entropy_pairs is None
    assert mn.get_network_edges() == []
    assert mn.similarity_matrices["entropy_similarity"].n_spectra == 0

    # Lifecycle: can run again after clear
    mn.run_query_vs_query_only(spectra[:2], ids[:2], query_precursor_mzs=precursor_mzs[:2])
    assert mn.stage_query_query_done is True
    assert mn.similarity_matrices["entropy_similarity"].n_spectra == 2


def test_hydrate_library_similarities_stage3(msp_fe_lib):
    """Stage 3 (library–library among matched entries) runs without error."""
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    if len(df) < 2:
        pytest.skip("need at least 2 library spectra for stage-3 edges")

    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q0")

    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
    )
    mn.query_vs_library(
        [q],
        ["q0"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=True,
        library_similarity_threshold=0.1,
    )
    assert mn.stage_query_query_done is True
    assert mn.stage_query_library_done is True
    assert mn.stage_library_library_done is True
    edges = mn.get_network_edges(metric="entropy_similarity")
    assert any("q0" in (a, b) for a, b, _ in edges)


def test_identity_search_type_smoke(msp_fe_lib):
    """Default-ish identity mode with precursor filtering smoke-tests cleanly."""
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    if pmz <= 0:
        pytest.skip("library row lacks precursor m/z for identity search")

    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q_id")
    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="identity",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
        ms1_tolerance_da=0.5,
    )
    mn.query_vs_library(
        [q],
        ["q_id"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=False,
    )
    stats = mn.get_network_stats()
    assert stats["n_nodes"] >= 1
    assert mn.stage_query_library_done is True


def test_save_edge_list_export_maps_lib_to_spectra_id(msp_fe_lib, tmp_path):
    """CSV export should map internal lib:<idx> nodes to spectra_id when present."""
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q0")

    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=[],
        similarity_thresholds={"entropy_similarity": 0.1},
    )
    mn.query_vs_library(
        [q],
        ["q0"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=False,
    )
    edges = mn.get_network_edges(metric="entropy_similarity")
    lib_nodes = []
    for a, b, _ in edges:
        for node in (a, b):
            if str(node).startswith("lib:"):
                lib_nodes.append(str(node))
    assert lib_nodes, "need internal lib: edges before export"

    # Resolve expected export IDs for every library endpoint in the edge list
    expected_export_ids = set()
    for lib_node in lib_nodes:
        lib_idx = MolecularNetwork.library_index_from_node_id(lib_node)
        assert lib_idx is not None
        entry = fe_lib[lib_idx]
        assert isinstance(entry, dict)
        export_id = str(entry.get("spectra_id") or entry.get("id") or "")
        assert export_id, f"library entry {lib_idx} has no spectra_id/id"
        expected_export_ids.add(export_id)

    out = tmp_path / "edges.csv"
    mn.save_edge_list(str(out), metric="entropy_similarity")
    text = out.read_text()
    assert "q0" in text
    for export_id in expected_export_ids:
        assert export_id in text
    # No internal lib: ids should remain in the CSV
    assert "lib:" not in text


def test_staged_api_stage2_has_no_threshold_kwarg():
    """Stage 2 no longer accepts library_similarity_threshold (lives on stage 3)."""
    import inspect

    sig = inspect.signature(MolecularNetwork.run_query_vs_library_stage)
    assert "library_similarity_threshold" not in sig.parameters
    sig3 = inspect.signature(MolecularNetwork.run_library_vs_library_stage)
    assert "library_similarity_threshold" in sig3.parameters


def test_defaults_search_type_open():
    """Constructor default search_type is open (DDA-friendly)."""
    from corems.molecular_networking.similarity_engine import SimilarityEngine

    mn = MolecularNetwork(fe_lib=None)
    assert mn.search_type == "open"
    assert not hasattr(mn._engine, "use_parallel")
    assert not hasattr(mn._engine, "n_jobs")

    eng = SimilarityEngine(fe_lib=None)
    assert eng.search_type == "open"
    assert not hasattr(eng, "use_parallel")
    assert not hasattr(eng, "n_jobs")


def test_dual_path_compute_all_vs_all_removed():
    """Slow O(n²) dual path was deleted; only FE vectorised path remains."""
    from corems.molecular_networking.similarity_engine import SimilarityEngine

    assert not hasattr(SimilarityEngine, "compute_all_vs_all")
    assert not hasattr(SimilarityEngine, "_pairwise_entropy")
    assert not hasattr(SimilarityEngine, "_entropy_score_pair")
    assert hasattr(SimilarityEngine, "compute_all_vs_all_with_lib")
    assert hasattr(SimilarityEngine, "search_queries_against_library")


def test_search_queries_against_library_public_api(msp_fe_lib):
    """Public engine Q–L API returns scores + library_size without private hooks."""
    from corems.molecular_networking.similarity_engine import SimilarityEngine

    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q0")

    engine = SimilarityEngine(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=["cosine"],
    )
    scores, lib_size = engine.search_queries_against_library(
        [q],
        ["q0"],
        query_precursor_mzs=[pmz],
        format_library_id=MolecularNetwork.library_node_id,
    )
    assert lib_size >= 1
    entropy = scores["entropy_similarity"]
    assert entropy
    assert any(qid == "q0" and str(lid).startswith("lib:") for qid, lid in entropy)


def _network_with_edges(msp_fe_lib):
    """Build a MolecularNetwork that has at least one query–library edge."""
    fe_lib, msp = msp_fe_lib
    df = msp._data_frame
    row = df.iloc[0]
    peaks = np.asarray(row.peaks, dtype=float)
    pmz = float(getattr(row, "precursormz", 0.0) or 0.0)
    q = MockSpectrum(peaks[:, 0], peaks[:, 1], name="q0")
    mn = MolecularNetwork(
        fe_lib=fe_lib,
        search_type="open",
        additional_similarities=["cosine"],
        similarity_thresholds={"entropy_similarity": 0.1, "cosine": 0.1},
    )
    mn.query_vs_library(
        [q],
        ["q0"],
        query_precursor_mzs=[pmz],
        fe_kwargs=FE_KWARGS,
        hydrate_library_similarities=False,
    )
    return mn


def test_compute_network_clusters_smoke(msp_fe_lib):
    pytest.importorskip("networkx")
    mn = _network_with_edges(msp_fe_lib)
    summary = mn.compute_network_clusters(
        metric="entropy_similarity",
        max_edges=100,
        compute_layout=True,
    )
    assert summary["metric"] == "entropy_similarity"
    assert summary["n_nodes"] >= 1
    assert "entropy_similarity" in mn._network_clusters
    artifact = mn._network_clusters["entropy_similarity"]
    assert not artifact["node_table"].empty or summary["n_nodes"] == 0
    assert "schema_version" in artifact
    assert artifact["schema_version"] == 1


def test_save_and_load_network_clusters_roundtrip(msp_fe_lib, tmp_path):
    """save_network_clusters writes CSVs; load_network_clusters restores cache."""
    pytest.importorskip("networkx")
    mn = _network_with_edges(msp_fe_lib)
    summary = mn.compute_network_clusters(
        metric="entropy_similarity",
        max_edges=100,
        compute_layout=True,
    )
    assert summary["n_nodes"] >= 1

    out_dir = tmp_path / "clusters"
    paths = mn.save_network_clusters(
        str(out_dir), metric="entropy_similarity", run_id="runA"
    )
    for key in ("nodes", "communities", "edges", "layout", "manifest"):
        assert key in paths
        assert Path(paths[key]).is_file()
        assert Path(paths[key]).stat().st_size > 0

    # Filenames include metric + run_id suffix
    assert "entropy_similarity_clusters_runA_nodes.csv" in paths["nodes"]
    assert "entropy_similarity_clusters_runA_manifest.csv" in paths["manifest"]

    # Clear cache and reload from disk
    mn.drop_network_clusters(metric="entropy_similarity")
    assert "entropy_similarity" not in mn._network_clusters

    loaded = mn.load_network_clusters(
        str(out_dir), metric="entropy_similarity", run_id="runA"
    )
    assert loaded["metric"] == "entropy_similarity"
    assert loaded["n_nodes"] == summary["n_nodes"]
    assert loaded["n_edges"] == summary["n_edges"]
    assert loaded["n_clusters"] == summary["n_clusters"]

    restored = mn._network_clusters["entropy_similarity"]
    assert restored["schema_version"] == 1
    assert list(restored["node_table"].columns)
    assert "params" in restored
    # Layout may be empty for tiny graphs but file/table should exist
    assert "layout_table" in restored


def test_save_network_clusters_requires_compute(msp_fe_lib, tmp_path):
    mn = _network_with_edges(msp_fe_lib)
    with pytest.raises(RuntimeError, match="No clusters available"):
        mn.save_network_clusters(str(tmp_path), metric="entropy_similarity")


def test_load_network_clusters_missing_files(tmp_path):
    mn = MolecularNetwork(fe_lib=None)
    with pytest.raises(FileNotFoundError, match="Missing cluster artifact"):
        mn.load_network_clusters(str(tmp_path), metric="entropy_similarity")


def test_drop_network_clusters_clears_cache(msp_fe_lib):
    pytest.importorskip("networkx")
    mn = _network_with_edges(msp_fe_lib)
    mn.compute_network_clusters(metric="entropy_similarity", max_edges=50)
    assert "entropy_similarity" in mn._network_clusters

    mn.drop_network_clusters(metric="entropy_similarity")
    assert "entropy_similarity" not in mn._network_clusters

    mn.compute_network_clusters(metric="entropy_similarity", max_edges=50)
    mn.drop_network_clusters()  # clear all
    assert mn._network_clusters == {}
