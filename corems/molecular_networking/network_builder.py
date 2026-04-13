"""
MolecularNetwork
================

Main user-facing interface for building and querying molecular networks
from a collection of mass spectra.

Usage pattern (lazy / explicit):
  1. Create a MolecularNetwork (no computation happens).
  2. Call ``query_vs_library()`` to compute similarities.
  3. Query edges, neighbors, and statistics.
  4. Save outputs.

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

    Initialization is **lazy** – no similarity computation happens until you
    explicitly call :meth:`query_vs_library`.

    Parameters
    ----------
    fe_lib : ms_entropy.FlashEntropySearch
        Pre-built FlashEntropy search instance (built from the reference
        library).  Used for query-vs-library similarity computation.
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

        # Build the engine (wraps fe_lib + search parameters)
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

        # One SimilarityMatrix per metric – empty until query_vs_library() is called
        all_metrics = ["entropy_similarity"] + list(additional_similarities)
        self.similarity_matrices: dict[str, SimilarityMatrix] = {
            m: SimilarityMatrix(metric_name=m) for m in all_metrics
        }

        # Accumulate query spectra (kept separate from the library FE index)
        self._all_query_spectra: list = []
        self._all_query_ids: list[str] = []
        self._all_query_precursor_mzs: list[float | None] = []
        self._all_query_lib_indices: list[int | None] = []
        # Whether query_vs_library has been run (disabled repeated runs until dropped)
        self._has_queries_run: bool = False

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

    def query_vs_library(
        self,
        query_spectra: list,
        query_ids: list[str],
        query_precursor_mzs: list[float | None] | None = None,
        fe_kwargs: dict | None = None,
        *,
        query_lib_indices: list[int | None] | None = None,
    ):
        """Compute query-vs-query and query-vs-library similarities.

        This is the primary method for building a molecular network.  It:

        1. Adds to (or initiates) a FlashEntropy index from the **query spectra**.
        2. Uses the pre-built library FlashEntropy index (``self.fe_lib``) to
           compute query-vs-library similarities and also computes the query-vs-query 
           similarities from the query FlashEntropy index.
        3. Stores all results in ``self.similarity_matrices``.

        For ``"identity"`` and ``"neutral_loss"`` search types,
        *query_precursor_mzs* is **required**.

        Parameters
        ----------
        query_spectra : list
            Experimental spectrum objects with ``.mz_exp`` and ``.abundance``.
        query_ids : list of str
            Unique IDs for each query spectrum.
        query_precursor_mzs : list of float or None, optional
            Precursor m/z for each query spectrum.  Required for
            ``"identity"`` and ``"neutral_loss"`` search types.
        fe_kwargs : dict, optional
            Extra keyword arguments forwarded to
            ``FlashEntropySearch`` when building the query-vs-query index.
            Defaults to the same settings used for the library index.

        Raises
        ------
        ValueError
            If *query_precursor_mzs* is None for identity/neutral_loss search.
        """
        n_query = len(query_spectra)

        if n_query == 0:
            return

        # Validate precursor_mzs requirement
        if self.search_type in _PRECURSOR_REQUIRED:
            if query_precursor_mzs is None:
                raise ValueError(
                    f"query_precursor_mzs is required for "
                    f"search_type='{self.search_type}'. "
                    "Provide a list of precursor m/z values, one per query spectrum."
                )
            if len(query_precursor_mzs) != n_query:
                raise ValueError(
                    f"query_precursor_mzs length ({len(query_precursor_mzs)}) "
                    f"must match query_spectra length ({n_query})."
                )
        else:
            query_precursor_mzs = [None] * n_query if query_precursor_mzs is None else query_precursor_mzs

        if query_lib_indices is None:
            query_lib_indices = [None] * n_query

        # Disallow repeated incremental additions — require explicit drop to run again
        if self._has_queries_run:
            raise RuntimeError(
                "query_vs_library has already been run — call drop_queries() to clear queries and results before running again."
            )

        # Store provided queries (single-run behavior)
        self._all_query_spectra = list(query_spectra)
        self._all_query_ids = list(query_ids)
        self._all_query_precursor_mzs = list(query_precursor_mzs)
        self._all_query_lib_indices = list(query_lib_indices)

        # Register query IDs now
        for mat in self.similarity_matrices.values():
            mat.register_spectra(self._all_query_ids)

        # ── Stage 1: Query-vs-Query ───────────────────────────────────────────
        # Build a temporary FE index from the query spectra, then search
        # each query against it to get all-vs-all query-vs-query scores.
        print(f"  [query_vs_library] Computing query-vs-query "
              f"({n_query} × {n_query} = {n_query * (n_query - 1) // 2} pairs) …")

        # Build a temporary FE index from the query spectra and compute all-vs-all
        query_fe_lib = self._engine.build_fe_index_from_spectra(
            spectra=self._all_query_spectra,
            precursor_mzs=self._all_query_precursor_mzs,
            fe_kwargs=fe_kwargs,
        )
        qq_scores = self._engine.compute_all_vs_all_with_lib(
            spectra=self._all_query_spectra,
            spectrum_ids=self._all_query_ids,
            precursor_mzs=self._all_query_precursor_mzs,
            fe_lib_override=query_fe_lib,
        )
        self._update_matrices(qq_scores)

        # Mark as run
        self._has_queries_run = True

        # ── Stage 2: Query-vs-Library ─────────────────────────────────────────
        # Use the pre-built library FE index (self.fe_lib) to search each
        # query spectrum against the full library.
        print(f"  [query_vs_library] Computing query-vs-library against internal FE library …")

        entropy_pairs: dict[tuple[str, str], float] = {}
        cosine_pairs: dict[tuple[str, str], float] = {}
        cross_pairs_for_cosine: list[tuple[int, int]] = []

        lib_size = 0
        # For each query, search against the FE library and collect scores
        for qi, (spec, pmz) in enumerate(zip(query_spectra, query_precursor_mzs)):
            peaks = self._engine._peaks_array(spec)
            if peaks.shape[0] == 0:
                continue
            result_vec = self._engine._clean_and_search(peaks, pmz)
            if result_vec is None:
                continue
            lib_size = max(lib_size, len(result_vec))
            # If a specific library index for this query is provided, only record that one
            if query_lib_indices and query_lib_indices[qi] is not None:
                li = query_lib_indices[qi]
                if 0 <= li < len(result_vec):
                    score = float(result_vec[li])
                    if score > 0.0:
                        entropy_pairs[(query_ids[qi], str(li))] = score
                        if score >= self._engine.entropy_threshold_low:
                            cross_pairs_for_cosine.append((qi, li))
            else:
                for lib_idx, score in enumerate(result_vec):
                    if score > 0.0:
                        entropy_pairs[(query_ids[qi], str(lib_idx))] = float(score)
                        if score >= self._engine.entropy_threshold_low:
                            cross_pairs_for_cosine.append((qi, lib_idx))

        # Register synthesized library IDs now that we know lib_size
        lib_ids = [str(i) for i in range(lib_size)]
        for mat in self.similarity_matrices.values():
            mat.register_spectra(lib_ids)

        # Compute cosine for cross pairs if requested and if we can access library spectra
        if "cosine" in self._engine.additional_similarities and cross_pairs_for_cosine:
            try:
                lib_specs = getattr(self.fe_lib, "spectra", None) or getattr(self.fe_lib, "library", None)
            except Exception:
                lib_specs = None

            if lib_specs is not None:
                cross_cos = self._engine._compute_cosine_for_pairs(
                    cross_pairs_for_cosine, query_spectra, lib_specs
                )
                for (i, j), score in cross_cos.items():
                    cosine_pairs[(query_ids[i], str(j))] = score

        combined: dict[str, dict[tuple[str, str], float]] = {"entropy_similarity": entropy_pairs}
        if cosine_pairs:
            combined["cosine"] = cosine_pairs
        self._update_matrices(combined)

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

    def drop_queries(self):
        """Clear all stored queries and computed results.

        Resets internal query storage and reinitialises the similarity
        matrices so `query_vs_library` can be called again.
        """
        # Clear stored queries
        self._all_query_spectra = []
        self._all_query_ids = []
        self._all_query_precursor_mzs = []
        self._all_query_lib_indices = []
        self._has_queries_run = False

        # Recreate empty similarity matrices for each metric
        all_metrics = list(self.similarity_matrices.keys())
        self.similarity_matrices = {m: SimilarityMatrix(metric_name=m) for m in all_metrics}

    def __repr__(self) -> str:
        n_nodes = sum(mat.n_spectra for mat in self.similarity_matrices.values()) // max(
            1, len(self.similarity_matrices)
        )
        metrics = list(self.similarity_matrices.keys())
        return (
            f"MolecularNetwork(search_type='{self.search_type}', "
            f"n_spectra={n_nodes}, metrics={metrics})"
        )
