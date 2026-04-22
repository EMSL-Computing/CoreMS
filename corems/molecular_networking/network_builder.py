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

from typing import Any

import numpy as np
import pandas as pd

from corems.molecular_networking.similarity_matrix import SimilarityMatrix
from corems.molecular_networking.similarity_engine import SimilarityEngine

# Metrics that require precursor_mzs
_PRECURSOR_REQUIRED = {"identity", "neutral_loss"}


class MolecularNetwork:
    """Build a molecular network around a reference library.

    Initialization is **lazy** - no similarity computation happens until you
    explicitly call :meth:`query_vs_library`.

    Parameters
    ----------
    fe_lib : ms_entropy.FlashEntropySearch
        Pre-built FlashEntropy search instance (built from the reference
        library).  Used for query-vs-library similarity computation and for
        building properly configured cleaning methods for input spectra.
        Tolerance parameters are automatically extracted from this library.
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
    ms1_tolerance_da : float, optional
        Precursor m/z tolerance (Da) for identity search only.
        If None (default) and identity search is selected, uses 0.01 Da.
    entropy_threshold_low : float | None
        Minimum entropy similarity score required to trigger additional metric
        computation.  Default half of lowest non-entropy similarity threshold.
    use_parallel : bool
        Enable multiprocessing for additional metric computation.  Default True.
    n_jobs : int
        Number of worker processes.  -1 uses all available cores.  Default -1.

    Attributes
    ----------
    similarity_matrices : dict of str → SimilarityMatrix
        One SimilarityMatrix per metric (``"entropy_similarity"``, ``"cosine"``, …).
    
    Notes
    -----
    The following parameters are automatically extracted from the FE library:
    
    - ``ms2_tolerance_da`` : From ``fe_lib.entropy_search.max_ms2_tolerance_in_da``
    - ``peak_sep_da`` : Computed as ``2 * ms2_tolerance_da`` (FE requirement)
    
    This ensures all similarity calculations use the same cleaning parameters
    as the FE library index.
    """

    def __init__(
        self,
        fe_lib,
        search_type: str = "identity",
        additional_similarities: list[str] | None = None,
        similarity_thresholds: dict[str, float] | None = None,
        ms1_tolerance_da: float | None = None,
        entropy_threshold_low: float | None = None,
        use_parallel: bool = True,
        n_jobs: int = -1,
    ):
        if additional_similarities is None:
            additional_similarities = ["cosine"]
        if similarity_thresholds is None:
            similarity_thresholds = {}
        
        if entropy_threshold_low is None:
            non_entropy_thresholds = [
                v for k, v in similarity_thresholds.items() if k != "entropy_similarity"
            ]
            if non_entropy_thresholds:
                entropy_threshold_low = min(non_entropy_thresholds) / 2
            else:
                entropy_threshold_low = 0.25 # Placeholder that won't be used if no additional similarities are computed

        self.fe_lib = fe_lib
        self.search_type = search_type
        self.additional_similarities = list(additional_similarities)
        self.similarity_thresholds = similarity_thresholds

        # Build the engine (wraps fe_lib + search parameters)
        # Tolerance parameters are extracted from fe_lib automatically
        self._engine = SimilarityEngine(
            fe_lib=fe_lib,
            search_type=search_type,
            additional_similarities=additional_similarities,
            ms1_tolerance_da=ms1_tolerance_da,
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
        hydrate_library_similarities: bool = False,
        library_similarity_threshold: float = 0.3,
    ):
        """Compute query-vs-query and query-vs-library similarities.

        This is the primary method for building a molecular network.  It:

        1. Builds a temporary FlashEntropy index from the **query spectra** and
           computes query-vs-query similarities.
        2. Uses the pre-built library FlashEntropy index (``self.fe_lib``) to
           compute query-vs-library similarities.
        3. Optionally (when *hydrate_library_similarities* is True), computes
           library-vs-library similarities for the subset of library spectra
           that had similarity ≥ *library_similarity_threshold* to any query.
        4. Stores all results in ``self.similarity_matrices``.

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
        hydrate_library_similarities : bool, optional
            When True, perform a third stage that computes library-vs-library
            similarities for the subset of library spectra that matched any
            query above *library_similarity_threshold*.  Default False.
        library_similarity_threshold : float, optional
            Minimum query-vs-library entropy similarity score required for a
            library spectrum to be included in the library-vs-library stage.
            Only used when *hydrate_library_similarities* is True.
            Default 0.5.

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
            if any(pmz is None for pmz in query_precursor_mzs):
                raise ValueError(
                    f"query_precursor_mzs contains None values for search_type='{self.search_type}'. "
                    "Provide a valid precursor m/z for every query spectrum."
                )
        else:
            query_precursor_mzs = [None] * n_query if query_precursor_mzs is None else query_precursor_mzs

        # Disallow repeated incremental additions — require explicit drop to run again
        if self._has_queries_run:
            raise RuntimeError(
                "query_vs_library has already been run — call drop_queries() to clear queries and results before running again."
            )

        # Store provided queries (single-run behavior)
        self._all_query_spectra = list(query_spectra)
        self._all_query_ids = list(query_ids)
        self._all_query_precursor_mzs = list(query_precursor_mzs)

        # Register query IDs now
        for mat in self.similarity_matrices.values():
            mat.register_spectra(self._all_query_ids)

        # ── Stage 1: Query-vs-Query ───────────────────────────────────────────
        # Build a temporary FE index from the query spectra, then search
        # each query against it to get all-vs-all query-vs-query scores.
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

        entropy_pairs: dict[tuple[str, str], float] = {}
        cosine_pairs: dict[tuple[str, str], float] = {}
        
        # Group library indices by query for efficient cosine computation
        query_to_lib_indices: dict[int, list[int]] = {}

        lib_size = 0
        
        # Search each query against library
        for qi, (spec, pmz) in enumerate(zip(query_spectra, query_precursor_mzs)):
            peaks = self._engine._peaks_array(spec)
            if peaks.shape[0] == 0:
                continue
            result_vec = self._engine._clean_and_search(peaks, pmz)
            if result_vec is None:
                continue
            lib_size = max(lib_size, len(result_vec))
            
            # Collect library indices that pass entropy threshold for this query
            lib_indices_for_query = []
            
            # Search full library for each query
            for lib_idx, score in enumerate(result_vec):
                if score > 0.0:
                    entropy_pairs[(query_ids[qi], str(lib_idx))] = float(score)
                    if score >= self._engine.entropy_threshold_low:
                        lib_indices_for_query.append(lib_idx)
            
            if lib_indices_for_query:
                query_to_lib_indices[qi] = lib_indices_for_query
        
        total_pairs = sum(len(v) for v in query_to_lib_indices.values())
        print(f"  [query_vs_library] Found {total_pairs} query-vs-library pairs above entropy threshold ({self._engine.entropy_threshold_low})")

        # Register synthesized library IDs now that we know lib_size
        lib_ids = [str(i) for i in range(lib_size)]
        for mat in self.similarity_matrices.values():
            mat.register_spectra(lib_ids)

        # Compute cosine for query-vs-library pairs if requested
        if "cosine" in self._engine.additional_similarities and query_to_lib_indices:
            print(f"  [query_vs_library] Computing cosine for {len(query_to_lib_indices)} queries against library...")
            if self.search_type == "neutral_loss":
                print(
                    "  [query_vs_library] neutral_loss mode: cosine computed in neutral-loss mass space "
                    "(precursor_mz - fragment_mz)."
                )
            
            # Compute cosine efficiently: one query at a time against all its matching library spectra
            # Cleaned peaks are extracted directly from the FE library
            for qi, lib_indices in query_to_lib_indices.items():
                cosine_scores = self._engine._compute_cosine_for_query_vs_library(
                    query_spectrum=query_spectra[qi],
                    query_precursor_mz=query_precursor_mzs[qi],
                    library_indices=lib_indices,
                )
                
                # Store results
                for lib_idx, score in cosine_scores.items():
                    cosine_pairs[(query_ids[qi], str(lib_idx))] = score
            
            print(f"  [query_vs_library] Stored {len(cosine_pairs)} cosine pairs")

        combined: dict[str, dict[tuple[str, str], float]] = {"entropy_similarity": entropy_pairs}
        if cosine_pairs:
            combined["cosine"] = cosine_pairs
        self._update_matrices(combined)

        # ── Stage 3: Library-vs-Library (optional, filtered) ──────────────────
        # Only executed when hydrate_library_similarities=True.
        # Identifies library spectra that matched any query above
        # library_similarity_threshold, then computes all-vs-all similarities
        # within that filtered subset.
        if hydrate_library_similarities:
            # Collect unique library indices that passed the threshold
            matched_lib_indices: list[int] = sorted(
                {
                    int(lib_id)
                    for (q_id, lib_id), score in entropy_pairs.items()
                    if q_id in set(query_ids) and score >= library_similarity_threshold
                }
            )

            n_matched = len(matched_lib_indices)
            print(
                f"  [query_vs_library] Stage 3 – library-vs-library for "
                f"{n_matched} matched library spectra "
                f"(threshold={library_similarity_threshold}) …"
            )

            if n_matched > 1:
                matched_lib_ids = [str(i) for i in matched_lib_indices]
                ll_scores = self._engine.compute_library_vs_library_filtered(
                    library_indices=matched_lib_indices,
                    spectrum_ids=matched_lib_ids,
                )
                self._update_matrices(ll_scores)
                n_ll_pairs = sum(
                    len(v) for v in ll_scores.values()
                )
                print(
                    f"  [query_vs_library] Stage 3 complete – "
                    f"{n_ll_pairs} library-vs-library pairs stored."
                )
            else:
                print(
                    "  [query_vs_library] Stage 3 skipped – "
                    "fewer than 2 library spectra matched the threshold."
                )

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

    # ── File output ───────────────────────────────────────────────────────────

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
