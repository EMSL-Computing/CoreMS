"""
MolecularNetwork
================

Main user-facing interface for building and querying molecular networks
from a collection of mass spectra.

Each similarity metric (entropy_similarity, cosine, …) has its own
SimilarityMatrix and can be queried independently.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from corems.molecular_networking.similarity_matrix import SimilarityMatrix
from corems.molecular_networking.similarity_engine import SimilarityEngine

# Metrics that require precursor_mzs
_PRECURSOR_REQUIRED = {"identity", "neutral_loss"}


class MolecularNetwork:
    """Build and query a molecular network from a collection of mass spectra.

    Parameters
    ----------
    fe_lib : ms_entropy.FlashEntropySearch
        Pre-built FlashEntropy search instance.
    search_type : str
        FlashEntropy search mode: ``"identity"``, ``"open"``, or
        ``"neutral_loss"``.  Default ``"identity"``.
    additional_similarities : list of str, optional
        Extra similarity metrics to compute alongside entropy similarity.
        Currently supported: ``["cosine"]``.  Default ``["cosine"]``.
    similarity_thresholds : dict, optional
        Mapping of metric name → edge threshold.  Edges are created in each
        metric's network when score >= threshold.  Metrics not listed use
        a default of 0.5.
        Example: ``{"entropy_similarity": 0.5, "cosine": 0.6}``
    peak_sep_da : float
        Minimum m/z separation between peaks (Da).  Default 0.01.
    ms1_tolerance_da : float
        Precursor m/z tolerance (Da) for identity/neutral_loss.  Default 0.01.
    ms2_tolerance_da : float
        Fragment m/z tolerance (Da) for FlashEntropy search.  Default 0.005.
    entropy_threshold_low : float
        Minimum entropy similarity score required to trigger additional metric
        computation.  Default 0.1.
    use_parallel : bool
        Enable multiprocessing for additional metric computation.  Default True.
    n_jobs : int
        Number of worker processes.  -1 uses all available cores.  Default -1.

    Attributes
    ----------
    similarity_matrices : dict of str → SimilarityMatrix
        One SimilarityMatrix per metric (``"entropy_similarity"``, ``"cosine"``, …).
    _spectra : list
        All registered spectrum objects (in registration order).
    _spectrum_ids : list of str
        All registered spectrum IDs (in registration order).
    _precursor_mzs : list of float or None
        Precursor m/z for each registered spectrum.
    _lib_indices : list of int or None
        FlashEntropy library index for each registered spectrum.
    """

    def __init__(
        self,
        fe_lib,
        search_type: str = "identity",
        additional_similarities: list[str] | None = None,
        similarity_thresholds: dict[str, float] | None = None,
        peak_sep_da: float = 0.01,
        ms1_tolerance_da: float = 0.01,
        ms2_tolerance_da: float = 0.005,
        entropy_threshold_low: float = 0.1,
        use_parallel: bool = True,
        n_jobs: int = -1,
    ):
        if additional_similarities is None:
            additional_similarities = ["cosine"]
        if similarity_thresholds is None:
            similarity_thresholds = {}

        self.fe_lib = fe_lib
        self.search_type = search_type
        self.additional_similarities = list(additional_similarities)
        self.similarity_thresholds = similarity_thresholds

        # Build the engine
        self._engine = SimilarityEngine(
            fe_lib=fe_lib,
            search_type=search_type,
            additional_similarities=additional_similarities,
            peak_sep_da=peak_sep_da,
            ms1_tolerance_da=ms1_tolerance_da,
            ms2_tolerance_da=ms2_tolerance_da,
            entropy_threshold_low=entropy_threshold_low,
            use_parallel=use_parallel,
            n_jobs=n_jobs,
        )

        # One SimilarityMatrix per metric
        all_metrics = ["entropy_similarity"] + list(additional_similarities)
        self.similarity_matrices: dict[str, SimilarityMatrix] = {
            m: SimilarityMatrix(metric_name=m) for m in all_metrics
        }

        # Internal spectrum registry
        self._spectra: list = []
        self._spectrum_ids: list[str] = []
        self._precursor_mzs: list[float | None] = []
        self._lib_indices: list[int | None] = []

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _threshold_for(self, metric: str) -> float:
        """Return the edge threshold for *metric*."""
        return self.similarity_thresholds.get(metric, 0.5)

    def _update_matrices(
        self, new_scores: dict[str, dict[tuple[str, str], float]]
    ):
        """Write new similarity scores into the appropriate SimilarityMatrix."""
        for metric, pairs in new_scores.items():
            mat = self.similarity_matrices.get(metric)
            if mat is None:
                continue
            for (id1, id2), score in pairs.items():
                mat.set_similarity(id1, id2, score)
        # Finalise all matrices for efficient querying
        for mat in self.similarity_matrices.values():
            mat.finalise()

    # ── Public API ────────────────────────────────────────────────────────────

    def add_spectra(
        self,
        spectra: list,
        spectrum_ids: list[str],
        precursor_mzs: list[float | None] | None = None,
        lib_indices: list[int | None] | None = None,
    ):
        """Add spectra to the network and compute new pairwise similarities.

        For ``"identity"`` and ``"neutral_loss"`` search types, *precursor_mzs*
        is **required** (one value per spectrum).  For ``"open"`` search type,
        *precursor_mzs* is ignored.

        Only NEW pairs are computed:
          - new-vs-existing
          - new-vs-new
        Existing-vs-existing pairs are never recomputed.

        Parameters
        ----------
        spectra : list
            Spectrum objects with ``.mz_exp`` and ``.abundance`` attributes.
        spectrum_ids : list of str
            User-provided IDs, one per spectrum.  Must be unique across all
            calls to ``add_spectra``.
        precursor_mzs : list of float or None, optional
            Precursor m/z for each spectrum.  Required for ``"identity"`` and
            ``"neutral_loss"`` search types.
        lib_indices : list of int or None, optional
            Index of each spectrum in the FlashEntropy library.  Providing
            these enables exact pairwise entropy scores; otherwise the engine
            searches each spectrum against the library and reads off the score
            at the partner's library position.

        Raises
        ------
        ValueError
            If *precursor_mzs* is None for identity/neutral_loss search types.
        ValueError
            If any spectrum ID in *spectrum_ids* is already registered.
        """
        n_new = len(spectra)
        if n_new == 0:
            return

        # Validate precursor_mzs requirement
        if self.search_type in _PRECURSOR_REQUIRED:
            if precursor_mzs is None:
                raise ValueError(
                    f"precursor_mzs is required for search_type='{self.search_type}'. "
                    "Provide a list of precursor m/z values, one per spectrum."
                )
            if len(precursor_mzs) != n_new:
                raise ValueError(
                    f"precursor_mzs length ({len(precursor_mzs)}) must match "
                    f"spectra length ({n_new})."
                )
        else:
            precursor_mzs = [None] * n_new

        if lib_indices is None:
            lib_indices = [None] * n_new

        # Check for duplicate IDs
        duplicates = set(spectrum_ids) & set(self._spectrum_ids)
        if duplicates:
            raise ValueError(
                f"The following spectrum IDs are already registered: {duplicates}"
            )

        # Register new IDs in all matrices
        for mat in self.similarity_matrices.values():
            mat.register_spectra(spectrum_ids)

        existing_spectra = list(self._spectra)
        existing_ids = list(self._spectrum_ids)
        existing_pmzs = list(self._precursor_mzs)
        existing_lib_idx = list(self._lib_indices)

        # Append to internal registry
        self._spectra.extend(spectra)
        self._spectrum_ids.extend(spectrum_ids)
        self._precursor_mzs.extend(precursor_mzs)
        self._lib_indices.extend(lib_indices)

        # ── Compute similarities ──────────────────────────────────────────────
        if not existing_spectra:
            # First batch: all-vs-all within the new batch
            new_scores = self._engine.compute_all_vs_all(
                spectra=spectra,
                spectrum_ids=spectrum_ids,
                precursor_mzs=precursor_mzs,
                lib_indices=lib_indices,
            )
        else:
            # Subsequent batches: new-vs-existing + new-vs-new
            new_scores = self._engine.compute_new_vs_existing(
                new_spectra=spectra,
                new_ids=spectrum_ids,
                existing_spectra=existing_spectra,
                existing_ids=existing_ids,
                new_precursor_mzs=precursor_mzs,
                existing_precursor_mzs=existing_pmzs,
                new_lib_indices=lib_indices,
                existing_lib_indices=existing_lib_idx,
            )

        self._update_matrices(new_scores)

    def get_network_edges(
        self, metric: str = "entropy_similarity"
    ) -> list[tuple[str, str, float]]:
        """Return all edges above the similarity threshold for *metric*.

        Parameters
        ----------
        metric : str
            Similarity metric to query.  Default ``"entropy_similarity"``.

        Returns
        -------
        list of (id1, id2, score)
            Upper-triangle pairs with score >= threshold.
        """
        mat = self._get_matrix(metric)
        threshold = self._threshold_for(metric)
        return mat.get_pairs_above_threshold(threshold)

    def get_spectrum_neighbors(
        self,
        spectrum_id: str,
        metric: str = "entropy_similarity",
    ) -> list[tuple[str, float]]:
        """Return all neighbors of *spectrum_id* above the threshold.

        Parameters
        ----------
        spectrum_id : str
            Query spectrum ID.
        metric : str
            Similarity metric to query.  Default ``"entropy_similarity"``.

        Returns
        -------
        list of (neighbor_id, score)
            Sorted by score descending.
        """
        mat = self._get_matrix(metric)
        threshold = self._threshold_for(metric)

        if spectrum_id not in mat._id_to_idx:
            return []

        idx = mat._id_to_idx[spectrum_id]
        if mat._matrix is None:
            return []

        mat.finalise()
        row = mat._matrix.getrow(idx)
        cx = row.tocoo()
        neighbors = []
        for j, v in zip(cx.col, cx.data):
            if v >= threshold and j != idx:
                neighbors.append((mat._idx_to_id[j], float(v)))
        neighbors.sort(key=lambda x: x[1], reverse=True)
        return neighbors

    def get_network_stats(
        self, metric: str = "entropy_similarity"
    ) -> dict[str, Any]:
        """Return basic network statistics for *metric*.

        Parameters
        ----------
        metric : str
            Similarity metric to query.  Default ``"entropy_similarity"``.

        Returns
        -------
        dict with keys:
            - ``n_nodes`` : number of registered spectra
            - ``n_edges`` : number of edges above threshold
            - ``avg_degree`` : average node degree
            - ``density`` : fraction of possible edges that exist
            - ``threshold`` : the threshold used
            - ``metric`` : the metric name
        """
        mat = self._get_matrix(metric)
        threshold = self._threshold_for(metric)
        edges = mat.get_pairs_above_threshold(threshold)
        n_nodes = mat.n_spectra
        n_edges = len(edges)
        max_edges = n_nodes * (n_nodes - 1) / 2 if n_nodes > 1 else 1
        degrees = {sid: 0 for sid in mat.spectrum_ids}
        for id1, id2, _ in edges:
            degrees[id1] += 1
            degrees[id2] += 1
        avg_degree = np.mean(list(degrees.values())) if degrees else 0.0
        return {
            "metric": metric,
            "threshold": threshold,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "avg_degree": float(avg_degree),
            "density": n_edges / max_edges if max_edges > 0 else 0.0,
        }

    # ── Graph export ──────────────────────────────────────────────────────────

    def to_networkx(self, metric: str = "entropy_similarity"):
        """Return a networkx Graph for *metric*.

        Requires ``networkx`` to be installed.

        Parameters
        ----------
        metric : str
            Similarity metric.  Default ``"entropy_similarity"``.

        Returns
        -------
        networkx.Graph
            Nodes are spectrum IDs; edges carry a ``"score"`` attribute.
        """
        try:
            import networkx as nx
        except ImportError:
            raise ImportError(
                "networkx is required for to_networkx(). "
                "Install it with: pip install networkx"
            )
        mat = self._get_matrix(metric)
        threshold = self._threshold_for(metric)
        G = nx.Graph()
        G.add_nodes_from(mat.spectrum_ids)
        for id1, id2, score in mat.get_pairs_above_threshold(threshold):
            G.add_edge(id1, id2, score=score)
        return G

    # ── File output ───────────────────────────────────────────────────────────

    def save_graphml(
        self,
        path: str,
        metric: str = "entropy_similarity",
    ):
        """Save the network as a GraphML file (compatible with Cytoscape).

        Parameters
        ----------
        path : str
            Output file path.
        metric : str
            Similarity metric.  Default ``"entropy_similarity"``.
        """
        G = self.to_networkx(metric=metric)
        import networkx as nx
        nx.write_graphml(G, path)
        print(f"Saved GraphML to {path}")

    def save_edge_list(
        self,
        path: str,
        metric: str = "entropy_similarity",
    ):
        """Save the edge list as a CSV file.

        Parameters
        ----------
        path : str
            Output file path.
        metric : str
            Similarity metric.  Default ``"entropy_similarity"``.
        """
        edges = self.get_network_edges(metric=metric)
        df = pd.DataFrame(edges, columns=["id1", "id2", "score"])
        df.to_csv(path, index=False)
        print(f"Saved edge list ({len(edges)} edges) to {path}")

    def save_similarity_matrix(
        self,
        path: str,
        metric: str = "entropy_similarity",
        threshold: float = 0.0,
    ):
        """Save the similarity matrix as a CSV file.

        Parameters
        ----------
        path : str
            Output file path.
        metric : str
            Similarity metric.  Default ``"entropy_similarity"``.
        threshold : float
            Minimum score to include.  Default 0.0 (all stored pairs).
        """
        mat = self._get_matrix(metric)
        df = mat.to_dataframe(threshold=threshold)
        df.to_csv(path, index=False)
        print(f"Saved similarity matrix ({len(df)} pairs) to {path}")

    # ── Visualisation ─────────────────────────────────────────────────────────

    def plot_network(
        self,
        metric: str = "entropy_similarity",
        node_label: str | None = None,
        layout: str = "spring",
        output_file: str | None = None,
        figsize: tuple[int, int] = (12, 10),
        dpi: int = 150,
    ):
        """Plot the molecular network using matplotlib + networkx.

        Parameters
        ----------
        metric : str
            Similarity metric.  Default ``"entropy_similarity"``.
        node_label : str or None
            Node attribute to use as label.  If None, uses spectrum IDs.
        layout : str
            networkx layout algorithm (``"spring"``, ``"kamada_kawai"``,
            ``"circular"``, etc.).  Default ``"spring"``.
        output_file : str or None
            If provided, save the figure to this path.
        figsize : tuple
            Figure size in inches.  Default (12, 10).
        dpi : int
            Figure resolution.  Default 150.
        """
        try:
            import matplotlib.pyplot as plt
            import networkx as nx
        except ImportError:
            raise ImportError(
                "matplotlib and networkx are required for plot_network(). "
                "Install with: pip install matplotlib networkx"
            )

        G = self.to_networkx(metric=metric)
        threshold = self._threshold_for(metric)

        layout_funcs = {
            "spring": nx.spring_layout,
            "kamada_kawai": nx.kamada_kawai_layout,
            "circular": nx.circular_layout,
            "spectral": nx.spectral_layout,
            "random": nx.random_layout,
        }
        layout_func = layout_funcs.get(layout, nx.spring_layout)
        pos = layout_func(G, seed=42)

        edge_weights = [G[u][v]["score"] for u, v in G.edges()]
        edge_widths = [w * 3 for w in edge_weights]

        labels = {n: n for n in G.nodes()}

        fig, ax = plt.subplots(figsize=figsize)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=300, node_color="steelblue", alpha=0.8)
        nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, alpha=0.5, edge_color="gray")
        nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=8)
        ax.set_title(
            f"Molecular Network – {metric} (threshold={threshold:.2f})\n"
            f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges"
        )
        ax.axis("off")
        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
            print(f"Saved network plot to {output_file}")
        else:
            plt.show()
        plt.close(fig)

    def plot_similarity_heatmap(
        self,
        metric: str = "entropy_similarity",
        output_file: str | None = None,
        figsize: tuple[int, int] = (10, 8),
        dpi: int = 150,
    ):
        """Plot a heatmap of the similarity matrix.

        Parameters
        ----------
        metric : str
            Similarity metric.  Default ``"entropy_similarity"``.
        output_file : str or None
            If provided, save the figure to this path.
        figsize : tuple
            Figure size in inches.  Default (10, 8).
        dpi : int
            Figure resolution.  Default 150.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is required for plot_similarity_heatmap(). "
                "Install with: pip install matplotlib"
            )

        mat = self._get_matrix(metric)
        dense = mat.to_dense()
        ids = mat.spectrum_ids

        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(dense, cmap="viridis", vmin=0, vmax=1, aspect="auto")
        plt.colorbar(im, ax=ax, label="Similarity score")
        ax.set_xticks(range(len(ids)))
        ax.set_yticks(range(len(ids)))
        ax.set_xticklabels(ids, rotation=90, fontsize=7)
        ax.set_yticklabels(ids, fontsize=7)
        ax.set_title(f"Similarity Heatmap – {metric}")
        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
            print(f"Saved heatmap to {output_file}")
        else:
            plt.show()
        plt.close(fig)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _get_matrix(self, metric: str) -> SimilarityMatrix:
        if metric not in self.similarity_matrices:
            raise KeyError(
                f"Unknown metric '{metric}'. "
                f"Available: {list(self.similarity_matrices.keys())}"
            )
        return self.similarity_matrices[metric]

    def __repr__(self) -> str:
        n = len(self._spectrum_ids)
        metrics = list(self.similarity_matrices.keys())
        return (
            f"MolecularNetwork(search_type='{self.search_type}', "
            f"n_spectra={n}, metrics={metrics})"
        )
