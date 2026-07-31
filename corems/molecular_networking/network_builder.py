"""
MolecularNetwork
================

Main user-facing interface for building and querying molecular networks
from a collection of mass spectra.

Two usage patterns are supported:

**All-in-one** (requires a reference library of a FlashEntropySearch instance):

.. code-block:: python

    mn = MolecularNetwork(fe_lib=my_fe_lib, search_type="open")
    mn.query_vs_library(spectra, ids, precursor_mzs)
    edges = mn.get_network_edges()

**Staged / query-only** (no reference library required for stage 1):

.. code-block:: python

    mn = MolecularNetwork(fe_lib=None, search_type="open")
    mn.run_query_vs_query_only(spectra, ids)          # stage 1
    mn.run_query_vs_library_stage()                   # stage 2 (needs fe_lib)
    mn.run_library_vs_library_stage()                 # stage 3 (optional)

Each similarity metric (``"entropy_similarity"``, ``"cosine"``, …) has its own
:class:`~corems.molecular_networking.similarity_matrix.SimilarityMatrix` and
can be queried independently.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from corems.molecular_networking.similarity_matrix import SimilarityMatrix
from corems.molecular_networking.similarity_engine import SimilarityEngine
from corems.molecular_networking.network_visualize import NetworkVisualizeMixin

# Metrics that require precursor_mzs
_PRECURSOR_REQUIRED = {"identity", "neutral_loss"}

# Internal node IDs for library spectra (must not collide with query mf IDs)
_LIB_NODE_PREFIX = "lib:"


class MolecularNetwork(NetworkVisualizeMixin):
    """Build a molecular network around a reference library.

    Initialization is **lazy** — no similarity computation happens until you
    call one of the computation methods.

    Parameters
    ----------
    fe_lib : ms_entropy.FlashEntropySearch or None, optional
        Pre-built FlashEntropy search instance (built from the reference
        library).  Required for query-vs-library and library-vs-library
        stages.  Pass ``None`` (default) to build a query-only network
        (stage 1 only).  When provided, fragment-tolerance parameters are
        automatically extracted from this library unless overridden by
        *ms2_tolerance_da*.
    search_type : str
        FlashEntropy search mode: ``"identity"``, ``"open"``, or
        ``"neutral_loss"``.  Default ``"open"`` (typical DDA networking;
        use ``"identity"`` when precursor filtering is required).
    additional_similarities : list of str, optional
        Extra similarity metrics to compute alongside entropy similarity.
        Currently supported: ``["cosine"]``.  Default ``["cosine"]``.
        Cosine is gated by *entropy_threshold_low* (best-effort secondary score).
    similarity_thresholds : dict, optional
        Mapping of metric name → edge threshold.  Edges are created in each
        metric's network when score >= threshold.  Metrics not listed use
        a default of 0.5.
        Example: ``{"entropy_similarity": 0.5, "cosine": 0.6}``
    ms1_tolerance_da : float, optional
        Precursor m/z tolerance (Da) used for ``"identity"`` search to filter
        library candidates by precursor m/z.  Ignored for ``"open"`` and
        ``"neutral_loss"`` (precursor filter is disabled for those modes).
        If None (default), uses 0.01 Da.
    ms2_tolerance_da : float, optional
        Fragment m/z tolerance (Da) used for spectrum cleaning and similarity
        scoring.  Priority: explicit kwarg > extracted from *fe_lib* >
        0.01 Da fallback.  Override this when *fe_lib* is ``None`` or when
        you need a tolerance different from the library's built-in value.
    entropy_threshold_low : float or None, optional
        Minimum entropy similarity score required to trigger additional metric
        computation (e.g. cosine).  If None (default), set to half the lowest
        non-entropy similarity threshold, or 0.25 when no thresholds are given.
    use_parallel : bool
        Forwarded to :class:`~corems.molecular_networking.similarity_engine.SimilarityEngine`.
        When ``True``, cosine scores for **all-vs-all pair batches** (query–query
        and optional library–library stages) may use a process pool.
        Does **not** parallelize FlashEntropy scoring or query–library cosine.
        Default ``False`` (safer for notebooks / macOS spawn / CI).  See Notes.
    n_jobs : int
        Worker process count when *use_parallel* is ``True`` (``-1`` = all
        cores).  Default ``-1``.

    Attributes
    ----------
    similarity_matrices : dict of str → SimilarityMatrix
        One :class:`~corems.molecular_networking.similarity_matrix.SimilarityMatrix`
        per metric (``"entropy_similarity"``, ``"cosine"``, …).
    stage_query_query_done : bool
        True after query-vs-query similarities have been computed (stage 1).
    stage_query_library_done : bool
        True after query-vs-library similarities have been computed (stage 2).
    stage_library_library_done : bool
        True after matched library-vs-library similarities have been computed
        (stage 3).

    Notes
    -----
    **Tolerance resolution order** (highest priority first):

    1. Explicit *ms2_tolerance_da* kwarg passed to this constructor.
    2. Value extracted from ``fe_lib.entropy_search.max_ms2_tolerance_in_da``.
    3. Hard-coded fallback of 0.01 Da (used when *fe_lib* is ``None`` and no
       explicit value is given).

    **Precursor filter behaviour by search type:**

    - ``"identity"`` : library candidates are filtered by precursor m/z within
      *ms1_tolerance_da* before scoring.
    - ``"open"`` : precursor filter is disabled; all library entries are scored.
    - ``"neutral_loss"`` : precursor filter is disabled; matching occurs in
      neutral-loss mass space (``precursor_mz − fragment_mz``), not on the
      precursor itself.

    **Parallelism (*use_parallel* / *n_jobs*):**

    These kwargs are stored only on the internal
    :class:`~corems.molecular_networking.similarity_engine.SimilarityEngine`
    (not as attributes of this class).  Scope:

    - **Uses a process pool when** ``use_parallel=True``, cosine is in
      *additional_similarities*, and there are multiple entropy-gated pairs
      to score in an all-vs-all batch (stage 1 query–query, and stage 3
      library–library when enabled).
    - **Never parallelizes** FlashEntropy / entropy similarity itself.
    - **Does not** parallelize query–library cosine (stage 2); that path is
      always serial per query.
    - **No-op** if cosine is disabled, only one pair is scored, or
      ``n_jobs <= 1``.
    """

    def __init__(
        self,
        fe_lib=None,
        search_type: str = "open",
        additional_similarities: list[str] | None = None,
        similarity_thresholds: dict[str, float] | None = None,
        ms1_tolerance_da: float | None = None,
        ms2_tolerance_da: float | None = None,
        entropy_threshold_low: float | None = None,
        use_parallel: bool = False,
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
        # When fe_lib is provided, tolerance parameters are extracted from it.
        # When fe_lib is None (queries-only use), SimilarityEngine will fall back
        # to its internal defaults and only query-vs-query stage is allowed.
        self._engine = SimilarityEngine(
            fe_lib=fe_lib,
            search_type=search_type,
            additional_similarities=additional_similarities,
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

        # Stage flags – enforce order:
        #   1) query-vs-query
        #   2) query-vs-library
        #   3) matched library-vs-library
        self.stage_query_query_done: bool = False
        self.stage_query_library_done: bool = False
        self.stage_library_library_done: bool = False

        # Whether query_vs_library has been run (disabled repeated runs until dropped)
        self._has_queries_run: bool = False

        # Cache for stage 2 entropy pairs (used by stage 3)
        self._stage2_entropy_pairs: dict[tuple[str, str], float] | None = None

        # Cache of computed network clustering artifacts keyed by metric.
        self._network_clusters: dict[str, dict] = {}

    @staticmethod
    def prepare_query_spectra_from_lcms_object(
        lcms_obj,
        mf_ids: set | None = None,
    ) -> tuple[list, list[str], list[float]]:
        """Build query lists from an ``LCMSBase`` object's mass features.

        Iterates ``lcms_obj.mass_features`` and collects each feature's
        ``best_ms2`` spectrum, skipping features with no MS2 or an empty
        ``mz_exp``.  If *mf_ids* is provided, only those mass features are
        included.

        Parameters
        ----------
        lcms_obj : LCMSBase
            LCMS object whose mass features have been populated with MS2
            spectra (e.g. via :meth:`add_associated_ms2_dda`).
        mf_ids : set, optional
            If provided, only mass features whose id is in this set are
            included.

        Returns
        -------
        (query_spectra, query_ids, query_precursor_mzs)
            Three parallel lists ready to pass to :meth:`query_vs_library`
            or :meth:`run_query_vs_query_only`.
        """
        query_spectra: list = []
        query_ids: list[str] = []
        query_precursor_mzs: list[float] = []
        for mf_id, mf in lcms_obj.mass_features.items():
            if mf_ids is not None and mf_id not in mf_ids:
                continue
            ms2 = mf.best_ms2
            if ms2 is None:
                continue
            if not hasattr(ms2, "mz_exp") or len(ms2.mz_exp) == 0:
                continue
            query_spectra.append(ms2)
            query_ids.append(str(mf_id))
            query_precursor_mzs.append(float(mf.mz))
        return query_spectra, query_ids, query_precursor_mzs

    @staticmethod
    def prepare_query_spectra_from_lcms_collection(
        lcms_collection,
        cluster_ids: set | None = None,
        representative_metric: str | None = None,
    ) -> tuple[list, list[str], list[float]]:
        """Build query lists from representative mass features in an ``LCMSCollection``.

        Uses :meth:`~corems.mass_spectra.calc.lc_calc.LCMSCollectionCalculations.get_representative_mass_features_for_all_clusters`
        to select one mass feature per consensus cluster, then collects each
        feature's ``best_ms2`` from the corresponding sample object.

        Representative mass features must already be loaded with MS2 spectra
        (e.g. via ``process_consensus_features(load_representatives=True,
        add_ms2=True, ...)`` or ``reload_representative_mass_features(add_ms2=True)``).
        Clusters whose representative has no usable MS2 are skipped.

        Parameters
        ----------
        lcms_collection : LCMSCollection
            Collection with consensus clusters and loaded representative
            mass features (with MS2 where available).
        cluster_ids : set, optional
            If provided, only representatives for these cluster IDs are
            included.
        representative_metric : str, optional
            Metric passed to
            :meth:`~corems.mass_spectra.calc.lc_calc.LCMSCollectionCalculations.get_representative_mass_features_for_all_clusters`.
            If ``None``, uses the collection parameter
            ``consensus_representative_metric``.

        Returns
        -------
        (query_spectra, query_ids, query_precursor_mzs)
            Three parallel lists ready to pass to :meth:`query_vs_library`
            or :meth:`run_query_vs_query_only`.  Query IDs are collection-level
            ``coll_mf_id`` strings (``"{sample_id}_{mf_id}"``).

        Raises
        ------
        AttributeError
            If *lcms_collection* does not expose representative lookup or
            sample indexing (not an ``LCMSCollection``-like object).
        """
        if not hasattr(lcms_collection, "get_representative_mass_features_for_all_clusters"):
            raise AttributeError(
                "lcms_collection must provide get_representative_mass_features_for_all_clusters "
                "(expected an LCMSCollection)."
            )

        representatives = lcms_collection.get_representative_mass_features_for_all_clusters(
            representative_metric=representative_metric
        )
        if cluster_ids is not None:
            representatives = representatives[
                representatives["cluster"].isin(cluster_ids)
            ]

        query_spectra: list = []
        query_ids: list[str] = []
        query_precursor_mzs: list[float] = []

        for _, row in representatives.iterrows():
            sample_id = int(row["sample_id"])
            mf_id = row["mf_id"]
            coll_mf_id = str(row["coll_mf_id"])

            try:
                sample = lcms_collection[sample_id]
            except (IndexError, KeyError, TypeError):
                continue

            mass_features = getattr(sample, "mass_features", None) or {}
            mf = mass_features.get(mf_id)
            if mf is None:
                # mf_id from dataframes may not match dict key type (int/str/np)
                try:
                    mf = mass_features.get(int(mf_id))
                except (TypeError, ValueError):
                    mf = None
            if mf is None:
                try:
                    mf = mass_features.get(str(mf_id))
                except (TypeError, ValueError):
                    mf = None
            if mf is None:
                continue

            ms2 = getattr(mf, "best_ms2", None)
            if ms2 is None:
                continue
            if not hasattr(ms2, "mz_exp") or len(ms2.mz_exp) == 0:
                continue

            query_spectra.append(ms2)
            query_ids.append(coll_mf_id)
            query_precursor_mzs.append(float(mf.mz))

        return query_spectra, query_ids, query_precursor_mzs

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _threshold_for(self, metric: str) -> float:
        """Return the edge threshold for *metric*."""
        return self.similarity_thresholds.get(metric, 0.5)

    @staticmethod
    def library_node_id(lib_idx: int) -> str:
        """Return the internal network node ID for a library spectrum index."""
        return f"{_LIB_NODE_PREFIX}{int(lib_idx)}"

    @staticmethod
    def library_index_from_node_id(node_id: str) -> int | None:
        """Parse a library matrix index from an internal node ID, if any.

        Accepts the canonical ``lib:<idx>`` form. Bare integer strings are
        also accepted for backward compatibility with older matrices.
        """
        if node_id is None:
            return None
        s = str(node_id)
        if s.startswith(_LIB_NODE_PREFIX):
            try:
                return int(s[len(_LIB_NODE_PREFIX) :])
            except ValueError:
                return None
        try:
            return int(s)
        except (TypeError, ValueError):
            return None

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

    def _prepare_queries(
        self,
        query_spectra: list,
        query_ids: list[str],
        query_precursor_mzs: list[float | None] | None,
    ) -> None:
        """Validate and store query spectra, IDs, and precursor m/z values.

        Resets all stage flags and the stage-2 entropy-pair cache so that a
        fresh computation can proceed.  Called internally by
        :meth:`run_query_vs_query_only` and :meth:`query_vs_library`.

        Parameters
        ----------
        query_spectra : list
            Spectrum objects with ``.mz_exp`` and ``.abundance`` attributes.
        query_ids : list of str
            Unique string ID for each spectrum.  Must be the same length as
            *query_spectra*.
        query_precursor_mzs : list of float or None, or None
            Precursor m/z for each spectrum.  Required (non-None list) for
            ``"identity"`` and ``"neutral_loss"`` search types.  Pass ``None``
            to auto-fill with ``None`` values for ``"open"`` search.

        Raises
        ------
        ValueError
            If *query_spectra* is empty, lengths are mismatched, or precursor
            m/z values are missing/None for a precursor-required search type.
        """
        n_query = len(query_spectra)
        if n_query == 0:
            raise ValueError("No query spectra provided.")

        if len(query_ids) != n_query:
            raise ValueError(
                f"query_ids length ({len(query_ids)}) must match query_spectra length ({n_query})."
            )

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
            if query_precursor_mzs is None:
                query_precursor_mzs = [None] * n_query
            elif len(query_precursor_mzs) != n_query:
                raise ValueError(
                    f"query_precursor_mzs length ({len(query_precursor_mzs)}) "
                    f"must match query_spectra length ({n_query})."
                )

        # Store provided queries
        self._all_query_spectra = list(query_spectra)
        self._all_query_ids = list(query_ids)
        self._all_query_precursor_mzs = list(query_precursor_mzs)

        # Register query IDs now
        for mat in self.similarity_matrices.values():
            mat.register_spectra(self._all_query_ids)

        # Reset stage flags and caches
        self.stage_query_query_done = False
        self.stage_query_library_done = False
        self.stage_library_library_done = False
        self._stage2_entropy_pairs = None

    def _run_query_vs_query(
        self,
        fe_kwargs: dict | None = None,
    ) -> None:
        """Stage 1 – compute query-vs-query similarities.

        Builds a temporary FlashEntropy index from the stored query spectra
        and computes all-vs-all pairwise similarities within that set.
        Results are written into :attr:`similarity_matrices` and
        :attr:`stage_query_query_done` is set to ``True``.

        Parameters
        ----------
        fe_kwargs : dict, optional
            Extra keyword arguments forwarded to
            :class:`ms_entropy.FlashEntropySearch` when building the
            temporary query index.

        Raises
        ------
        RuntimeError
            If no queries have been prepared (i.e. :meth:`_prepare_queries`
            has not been called).
        """
        if not self._all_query_spectra:
            raise RuntimeError("No queries prepared. Call run_query_vs_query_only() or query_vs_library() first.")

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
        self.stage_query_query_done = True

    def _run_query_vs_library(
        self,
    ) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], float]]:
        """Stage 2 – compute query-vs-library similarities.

        Searches each stored query spectrum against the reference FE library
        (``self.fe_lib``).  Entropy similarity scores are stored for all
        library entries with score > 0; cosine scores are computed for pairs
        above :attr:`~SimilarityEngine.entropy_threshold_low`.  Results are
        written into :attr:`similarity_matrices`,
        :attr:`stage_query_library_done` is set to ``True``, and the entropy
        pairs are cached in :attr:`_stage2_entropy_pairs` for use by stage 3.

        Returns
        -------
        entropy_pairs : dict mapping (query_id, lib_idx_str) → float
            All query-vs-library entropy similarity scores > 0.
        cosine_pairs : dict mapping (query_id, lib_idx_str) → float
            Cosine scores for pairs above the entropy threshold.

        Raises
        ------
        RuntimeError
            If stage 1 has not been completed or if ``self.fe_lib`` is ``None``.
        """
        if not self.stage_query_query_done:
            raise RuntimeError("query-vs-query stage must be run before query-vs-library stage.")
        if self.fe_lib is None:
            raise RuntimeError(
                "query-vs-library stage requires a reference FE library (fe_lib). "
                "Construct MolecularNetwork with fe_lib or use run_query_vs_query_only for query-only networks."
            )

        combined, lib_size = self._engine.search_queries_against_library(
            query_spectra=self._all_query_spectra,
            query_ids=self._all_query_ids,
            query_precursor_mzs=self._all_query_precursor_mzs,
            format_library_id=self.library_node_id,
        )

        lib_ids = [self.library_node_id(i) for i in range(lib_size)]
        for mat in self.similarity_matrices.values():
            mat.register_spectra(lib_ids)

        self._update_matrices(combined)

        entropy_pairs = combined.get("entropy_similarity", {})
        cosine_pairs = combined.get("cosine", {})

        self.stage_query_library_done = True
        self._stage2_entropy_pairs = entropy_pairs

        return entropy_pairs, cosine_pairs

    def _run_library_vs_library(
        self,
        entropy_pairs: dict[tuple[str, str], float],
        library_similarity_threshold: float,
    ) -> None:
        """Stage 3 – compute library-vs-library similarities for matched entries.

        Identifies the subset of library spectra that scored ≥
        *library_similarity_threshold* against any query in stage 2, then
        computes all-vs-all pairwise similarities within that subset.
        Results are written into :attr:`similarity_matrices` and
        :attr:`stage_library_library_done` is set to ``True``.

        Parameters
        ----------
        entropy_pairs : dict mapping (query_id, lib_idx_str) → float
            The entropy similarity pairs returned by :meth:`_run_query_vs_library`.
            Used to identify which library indices exceeded the threshold.
        library_similarity_threshold : float
            Minimum entropy similarity score a library entry must have against
            any query to be included in this stage.

        Raises
        ------
        RuntimeError
            If stage 2 has not been completed.
        """
        if not self.stage_query_library_done:
            raise RuntimeError("query-vs-library stage must be run before library-vs-library stage.")

        matched_lib_indices: list[int] = sorted(
            {
                lib_idx
                for (q_id, lib_id), score in entropy_pairs.items()
                if q_id in set(self._all_query_ids)
                and score >= library_similarity_threshold
                and (lib_idx := self.library_index_from_node_id(lib_id)) is not None
            }
        )

        n_matched = len(matched_lib_indices)

        if n_matched > 1:
            matched_lib_ids = [self.library_node_id(i) for i in matched_lib_indices]
            ll_scores = self._engine.compute_library_vs_library_filtered(
                library_indices=matched_lib_indices,
                spectrum_ids=matched_lib_ids,
            )
            self._update_matrices(ll_scores)

        self.stage_library_library_done = True

    def _export_node_id(self, node_id: str, query_id_set: set[str]) -> str:
        """Return node ID for file export.

        Query IDs are preserved as-is. Library nodes stored internally as
        ``lib:<idx>`` (or bare index strings) are mapped to FlashEntropy entry
        ``spectra_id`` when available (fallback to ``id``, then the internal
        node ID).
        """
        if node_id in query_id_set:
            return node_id

        lib_idx = self.library_index_from_node_id(node_id)
        if lib_idx is None or lib_idx < 0:
            return node_id

        if self.fe_lib is None:
            return node_id

        try:
            lib_entry = self.fe_lib[lib_idx]
        except (IndexError, KeyError, TypeError):
            return node_id

        if not isinstance(lib_entry, dict):
            return node_id

        lib_id = lib_entry.get("spectra_id")
        if lib_id is None:
            lib_id = lib_entry.get("id")
        if lib_id is None:
            return node_id

        lib_id_str = str(lib_id)
        return lib_id_str if lib_id_str else node_id

    # ── Public API ────────────────────────────────────────────────────────────

    def run_query_vs_query_only(
        self,
        query_spectra: list,
        query_ids: list[str],
        query_precursor_mzs: list[float | None] | None = None,
        fe_kwargs: dict | None = None,
    ) -> None:
        """Compute query-vs-query similarities (stage 1 only).

        Builds a temporary FlashEntropy index from *query_spectra* and
        computes all-vs-all pairwise similarities within that set.  Results
        are stored in :attr:`similarity_matrices`.

        After this call, :attr:`stage_query_query_done` is ``True``.  If a
        reference library was provided at construction time you may continue
        with :meth:`run_query_vs_library_stage`.

        Parameters
        ----------
        query_spectra : list
            Spectrum objects with ``.mz_exp`` and ``.abundance`` attributes.
        query_ids : list of str
            Unique string ID for each spectrum.
        query_precursor_mzs : list of float or None, optional
            Precursor m/z for each spectrum.  Required (non-None list) for
            ``"identity"`` and ``"neutral_loss"`` search types.
        fe_kwargs : dict, optional
            Extra keyword arguments forwarded to
            :class:`ms_entropy.FlashEntropySearch` when building the
            temporary query index.

        Raises
        ------
        RuntimeError
            If queries have already been run (call :meth:`drop_queries` first).
        ValueError
            If *query_spectra* is empty, lengths are mismatched, or precursor
            m/z values are missing for a precursor-required search type.
        """
        if self._has_queries_run:
            raise RuntimeError(
                "query_vs_library has already been run — call drop_queries() to clear queries and results before running again."
            )

        self._prepare_queries(
            query_spectra=query_spectra,
            query_ids=query_ids,
            query_precursor_mzs=query_precursor_mzs,
        )
        self._run_query_vs_query(fe_kwargs=fe_kwargs)
        self._has_queries_run = True

    def run_query_vs_library_stage(self) -> None:
        """Continue from stage 1 and compute query-vs-library similarities (stage 2).

        Searches each stored query spectrum against the reference FE library
        and stores entropy + cosine scores in :attr:`similarity_matrices`.
        After this call, :attr:`stage_query_library_done` is ``True``.

        Stage-3 library subset selection uses
        *library_similarity_threshold* on :meth:`run_library_vs_library_stage`
        (or the same keyword on :meth:`query_vs_library`), not this method.

        Raises
        ------
        RuntimeError
            If stage 1 has not been completed (call
            :meth:`run_query_vs_query_only` first) or if ``fe_lib`` is
            ``None``.
        """
        if not self._has_queries_run:
            raise RuntimeError(
                "No queries have been run yet — call run_query_vs_query_only() or query_vs_library() first."
            )

        entropy_pairs, _ = self._run_query_vs_library()
        self._stage2_entropy_pairs = entropy_pairs

    def run_library_vs_library_stage(
        self,
        library_similarity_threshold: float = 0.3,
    ) -> None:
        """Continue from stage 2 and compute library-vs-library similarities (stage 3).

        Identifies library spectra that scored ≥ *library_similarity_threshold*
        against any query in stage 2, then computes all-vs-all pairwise
        similarities within that subset.  Results are stored in
        :attr:`similarity_matrices` and :attr:`stage_library_library_done`
        is set to ``True``.

        Parameters
        ----------
        library_similarity_threshold : float, optional
            Minimum entropy similarity score (from stage 2) required for a
            library spectrum to be included in this stage.  Default 0.3.

        Raises
        ------
        RuntimeError
            If stage 2 has not been completed (call
            :meth:`run_query_vs_library_stage` first).
        """
        if self._stage2_entropy_pairs is None:
            raise RuntimeError(
                "query-vs-library stage must be run before library-vs-library stage."
            )

        self._run_library_vs_library(
            entropy_pairs=self._stage2_entropy_pairs,
            library_similarity_threshold=library_similarity_threshold,
        )

    def query_vs_library(
        self,
        query_spectra: list,
        query_ids: list[str],
        query_precursor_mzs: list[float | None] | None = None,
        fe_kwargs: dict | None = None,
        *,
        hydrate_library_similarities: bool = False,
        library_similarity_threshold: float = 0.3,
    ) -> None:
        """Run all three stages in one call (convenience wrapper).

        Executes stages 1–3 sequentially:

        1. **Query-vs-query** – builds a temporary FlashEntropy index from
           *query_spectra* and computes all-vs-all pairwise similarities.
        2. **Query-vs-library** – searches each query against the reference
           FE library (``self.fe_lib``).
        3. **Library-vs-library** *(optional)* – when
           *hydrate_library_similarities* is ``True``, computes pairwise
           similarities for the subset of library spectra that matched any
           query above *library_similarity_threshold*.

        All results are stored in :attr:`similarity_matrices`.

        Parameters
        ----------
        query_spectra : list
            Spectrum objects with ``.mz_exp`` and ``.abundance`` attributes.
        query_ids : list of str
            Unique string ID for each spectrum.
        query_precursor_mzs : list of float or None, optional
            Precursor m/z for each spectrum.  Required (non-None list) for
            ``"identity"`` and ``"neutral_loss"`` search types.
        fe_kwargs : dict, optional
            Extra keyword arguments forwarded to
            :class:`ms_entropy.FlashEntropySearch` when building the
            temporary query index (stage 1).
        hydrate_library_similarities : bool, optional
            When ``True``, run stage 3 (library-vs-library).  Default ``False``.
        library_similarity_threshold : float, optional
            Minimum entropy similarity score (stage 2) required for a library
            spectrum to be included in stage 3.  Default 0.3.

        Raises
        ------
        RuntimeError
            If queries have already been run (call :meth:`drop_queries` first)
            or if ``fe_lib`` is ``None`` (required for stages 2 and 3).
        ValueError
            If *query_spectra* is empty, lengths are mismatched, or precursor
            m/z values are missing for a precursor-required search type.
        """
        if self._has_queries_run:
            raise RuntimeError(
                "query_vs_library has already been run — call drop_queries() to clear queries and results before running again."
            )

        self._prepare_queries(
            query_spectra=query_spectra,
            query_ids=query_ids,
            query_precursor_mzs=query_precursor_mzs,
        )
        self._run_query_vs_query(fe_kwargs=fe_kwargs)
        entropy_pairs, _ = self._run_query_vs_library()

        if hydrate_library_similarities:
            self._run_library_vs_library(
                entropy_pairs=entropy_pairs,
                library_similarity_threshold=library_similarity_threshold,
            )

        self._has_queries_run = True

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
        query_id_set = set(self._all_query_ids)
        export_edges = [
            (
                self._export_node_id(id1, query_id_set),
                self._export_node_id(id2, query_id_set),
                score,
            )
            for id1, id2, score in edges
        ]
        df = pd.DataFrame(export_edges, columns=["id1", "id2", "score"])
        df.to_csv(path, index=False)

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

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _get_matrix(self, metric: str) -> SimilarityMatrix:
        if metric not in self.similarity_matrices:
            raise KeyError(
                f"Unknown metric '{metric}'. "
                f"Available: {list(self.similarity_matrices.keys())}"
            )
        return self.similarity_matrices[metric]

    def drop_queries(self):
        """Clear all stored queries and computed similarity results.

        Resets internal query storage, stage flags, the stage-2 entropy-pair
        cache, clustering artifacts, and reinitialises all
        :attr:`similarity_matrices` to empty so that
        :meth:`run_query_vs_query_only` or :meth:`query_vs_library` can be
        called again with a consistent cleared state.
        """
        # Clear stored queries
        self._all_query_spectra = []
        self._all_query_ids = []
        self._all_query_precursor_mzs = []
        self._has_queries_run = False
        self._stage2_entropy_pairs = None

        # Stage flags must match empty matrices (not "stale done")
        self.stage_query_query_done = False
        self.stage_query_library_done = False
        self.stage_library_library_done = False

        # Recreate empty similarity matrices for each metric
        all_metrics = list(self.similarity_matrices.keys())
        self.similarity_matrices = {m: SimilarityMatrix(metric_name=m) for m in all_metrics}

        # Drop any cached clustering artifacts tied to previous query results.
        self._network_clusters = {}

    def __repr__(self) -> str:
        n_nodes = sum(mat.n_spectra for mat in self.similarity_matrices.values()) // max(
            1, len(self.similarity_matrices)
        )
        metrics = list(self.similarity_matrices.keys())
        return (
            f"MolecularNetwork(search_type='{self.search_type}', "
            f"n_spectra={n_nodes}, metrics={metrics})"
        )
