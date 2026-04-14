"""
SimilarityEngine
================

Computes pairwise spectral similarities between query spectra using:
  1. FlashEntropy (fast, vectorized) – always computed first.
  2. Additional metrics (e.g., cosine) – computed only for pairs that
     pass a low entropy-similarity threshold.

Supports three FlashEntropy search modes:
  - "identity"     : precursor-matched (requires precursor_mzs)
  - "open"         : no precursor matching
  - "neutral_loss" : neutral-loss matched (requires precursor_mzs)

Parallel processing uses multiprocessing.Pool to match CoreMS conventions.
"""

from __future__ import annotations

import multiprocessing
from itertools import combinations
from typing import Any

import numpy as np

from corems.mass_spectra.calc.lc_calc import find_closest

# Supported additional similarity metrics
_SUPPORTED_ADDITIONAL = {"cosine"}

# Map search_type → FlashEntropy method name and result key
_FE_METHOD_MAP = {
    "identity": ("identity", "identity_search"),
    "open": ("open", "open_search"),
    "neutral_loss": ("neutral_loss", "neutral_loss_search"),
}


def _compute_cosine_pair(args):
    """Worker function for multiprocessing: compute cosine similarity for one pair.

    Uses tolerance-based peak matching to align spectra, then computes cosine on
    aligned vectors including both matched and unmatched peaks (unmatched peaks
    get 0 abundance in the other spectrum).

    Parameters
    ----------
    args : tuple
        (mz1, abun1, mz2, abun2, tolerance_da)
        where tolerance_da is the m/z matching tolerance in Da.

    Returns
    -------
    float
        Cosine similarity score in [0, 1].
    """
    if len(args) == 5:
        mz1, abun1, mz2, abun2, tolerance_da = args
    else:
        mz1, abun1, mz2, abun2 = args
        tolerance_da = 0.01

    try:
        # Convert to numpy arrays
        mz1 = np.asarray(mz1, dtype=float)
        abun1 = np.asarray(abun1, dtype=float)
        mz2 = np.asarray(mz2, dtype=float)
        abun2 = np.asarray(abun2, dtype=float)
        
        if len(mz1) == 0 or len(mz2) == 0:
            return 0.0
        
        # Sort both spectra by m/z
        idx1 = np.argsort(mz1)
        mz1_sorted = mz1[idx1]
        abun1_sorted = abun1[idx1]
        
        idx2 = np.argsort(mz2)
        mz2_sorted = mz2[idx2]
        abun2_sorted = abun2[idx2]
        
        # Build aligned vectors including all peaks
        vec1 = []
        vec2 = []
        used_spec2 = np.zeros(len(mz2_sorted), dtype=bool)
        
        # For each peak in spec1, find match in spec2 or add as unmatched
        for i in range(len(mz1_sorted)):
            # Find closest peak in spec2
            closest_idx = find_closest(mz2_sorted, np.array([mz1_sorted[i]]))[0]
            diff = abs(mz2_sorted[closest_idx] - mz1_sorted[i])
            
            if diff <= tolerance_da:
                # Matched peak
                vec1.append(abun1_sorted[i])
                vec2.append(abun2_sorted[closest_idx])
                used_spec2[closest_idx] = True
            else:
                # Unmatched peak in spec1
                vec1.append(abun1_sorted[i])
                vec2.append(0.0)
        
        # Add unmatched peaks from spec2
        for j in range(len(mz2_sorted)):
            if not used_spec2[j]:
                vec1.append(0.0)
                vec2.append(abun2_sorted[j])
        
        # Convert to numpy arrays
        vec1 = np.array(vec1, dtype=float)
        vec2 = np.array(vec2, dtype=float)
        
        # Compute cosine similarity
        norm1 = np.linalg.norm(vec1)
        norm2 = np.linalg.norm(vec2)
        
        if norm1 == 0 or norm2 == 0:
            return 0.0
        
        cosine = np.dot(vec1, vec2) / (norm1 * norm2)
        return float(np.clip(cosine, 0.0, 1.0))
        
    except Exception:
        return 0.0


class SimilarityEngine:
    """Compute pairwise spectral similarities using FlashEntropy + optional extras.

    Parameters
    ----------
    fe_lib : ms_entropy.FlashEntropySearch
        Pre-built FlashEntropy search instance (from MSPInterface._to_flashentropy).
    search_type : str
        FlashEntropy search mode: ``"identity"``, ``"open"``, or ``"neutral_loss"``.
        Default ``"identity"``.
    additional_similarities : list of str, optional
        Extra similarity metrics to compute for pairs passing the entropy threshold.
        Currently supported: ``["cosine"]``.  Default ``["cosine"]``.
    peak_sep_da : float
        Minimum m/z separation between peaks (Da).  Default 0.01.
    ms1_tolerance_da : float
        Precursor m/z tolerance (Da) for identity/neutral_loss search.  Default 0.01.
    ms2_tolerance_da : float
        Fragment m/z tolerance (Da) for FlashEntropy search.  Default 0.005.
    entropy_threshold_low : float
        Minimum entropy similarity score required to trigger additional metric
        computation.  Default 0.1.
    use_parallel : bool
        Enable multiprocessing for additional metric computation.  Default True.
    n_jobs : int
        Number of worker processes.  -1 uses all available cores.  Default -1.
    """

    def __init__(
        self,
        fe_lib,
        search_type: str = "identity",
        additional_similarities: list[str] | None = None,
        peak_sep_da: float = 0.01,
        ms1_tolerance_da: float = 0.01,
        ms2_tolerance_da: float = 0.005,
        entropy_threshold_low: float = 0.1,
        use_parallel: bool = True,
        n_jobs: int = -1,
    ):
        if search_type not in _FE_METHOD_MAP:
            raise ValueError(
                f"search_type must be one of {list(_FE_METHOD_MAP.keys())}, "
                f"got '{search_type}'."
            )
        if additional_similarities is None:
            additional_similarities = ["cosine"]
        unsupported = set(additional_similarities) - _SUPPORTED_ADDITIONAL
        if unsupported:
            raise ValueError(
                f"Unsupported additional_similarities: {unsupported}. "
                f"Supported: {_SUPPORTED_ADDITIONAL}"
            )

        self.fe_lib = fe_lib
        self.search_type = search_type
        #TODO KRH: Check that the additional similarity searches are done in the same manner (open vs neutral loss)
        self.additional_similarities = list(additional_similarities)
        self.peak_sep_da = peak_sep_da
        self.ms1_tolerance_da = ms1_tolerance_da
        self.ms2_tolerance_da = ms2_tolerance_da
        self.entropy_threshold_low = entropy_threshold_low
        self.use_parallel = use_parallel
        self.n_jobs = (
            multiprocessing.cpu_count() if n_jobs == -1 else max(1, n_jobs)
        )

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _peaks_array(self, spectrum) -> np.ndarray:
        """Return (N, 2) peaks array from a spectrum object."""
        mz = np.asarray(spectrum.mz_exp, dtype=float)
        ab = np.asarray(spectrum.abundance, dtype=float)
        return np.column_stack((mz, ab))

    def _clean_and_search(
        self,
        peaks: np.ndarray,
        precursor_mz: float | None,
        fe_lib_override=None,
    ) -> np.ndarray:
        """Run FlashEntropy search for one query spectrum.

        Parameters
        ----------
        peaks : np.ndarray of shape (N, 2)
            [[mz, abundance], ...]
        precursor_mz : float or None
            Required for identity/neutral_loss; ignored for open.
        fe_lib_override : optional
            If provided, use this FlashEntropy library instead of ``self.fe_lib``.

        Returns
        -------
        np.ndarray
            1-D array of entropy similarity scores, one per library entry.
        """
        fe_method, fe_key = _FE_METHOD_MAP[self.search_type]
        fe = fe_lib_override if fe_lib_override is not None else self.fe_lib

        # Use a dummy precursor_mz for open search
        pmz = precursor_mz if precursor_mz is not None else 0.0

        cleaned = fe.clean_spectrum_for_search(
            precursor_mz=pmz,
            peaks=peaks,
            precursor_ions_removal_da=None,
            noise_threshold=0.0,
            min_ms2_difference_in_da=self.peak_sep_da,
        )

        search_kwargs: dict[str, Any] = dict(
            peaks=cleaned,
            ms2_tolerance_in_da=self.ms2_tolerance_da,
            method={fe_method},
            precursor_ions_removal_da=None,
            noise_threshold=0.0,
            target="cpu",
        )
        if self.search_type in ("identity", "neutral_loss"):
            search_kwargs["precursor_mz"] = pmz
            search_kwargs["ms1_tolerance_in_da"] = self.ms1_tolerance_da
        else:
            # open search – precursor_mz still required by ms_entropy API
            search_kwargs["precursor_mz"] = pmz
            search_kwargs["ms1_tolerance_in_da"] = 1e9  # effectively no filter

        results = fe.search(**search_kwargs)
        return results[fe_key]

    def _compute_cosine_for_query_vs_library(
        self,
        query_spectrum,
        query_precursor_mz: float | None,
        library_indices: list[int],
        lib_specs: list,
    ) -> dict[int, float]:
        """Compute cosine similarity for one query against multiple library spectra.
        
        This is more efficient than computing cosine for each pair individually because:
        1. We clean the query spectrum once
        2. We use find_closest to match library peaks to query peaks
        3. We compute cosine including all peaks (matched and unmatched)
        
        Parameters
        ----------
        query_spectrum : spectrum object
            Query spectrum with .mz_exp and .abundance attributes
        query_precursor_mz : float or None
            Precursor m/z for the query
        library_indices : list of int
            Library indices to compute cosine against
        lib_specs : list
            Library spectra (dicts with "peaks" key)
            
        Returns
        -------
        dict mapping library_idx → cosine score
        """
        from corems.mass_spectra.calc.lc_calc import find_closest
        
        # Get and clean query peaks (same as FlashEntropy does)
        query_peaks = self._peaks_array(query_spectrum)
        if query_peaks.shape[0] == 0:
            return {}
        
        pmz = query_precursor_mz if query_precursor_mz is not None else 0.0
        cleaned_query = self.fe_lib.clean_spectrum_for_search(
            precursor_mz=pmz,
            peaks=query_peaks,
            precursor_ions_removal_da=None,
            noise_threshold=0.0,
            min_ms2_difference_in_da=self.peak_sep_da,
        )
        
        if cleaned_query.shape[0] == 0:
            return {}
        
        # Sort cleaned query peaks by m/z
        query_mz = cleaned_query[:, 0]
        query_abun = cleaned_query[:, 1]
        sort_idx = np.argsort(query_mz)
        query_mz_sorted = query_mz[sort_idx]
        query_abun_sorted = query_abun[sort_idx]
        
        # Compute cosine for each library spectrum
        cosine_scores = {}
        for lib_idx in library_indices:
            if lib_idx >= len(lib_specs):
                continue
            
            lib_entry = lib_specs[lib_idx]
            lib_peaks = np.asarray(lib_entry.get("peaks", []), dtype=float)
            
            if lib_peaks.ndim != 2 or lib_peaks.shape[1] < 2 or lib_peaks.shape[0] == 0:
                continue
            
            lib_mz = lib_peaks[:, 0]
            lib_abun = lib_peaks[:, 1]
            
            # Sort library peaks by m/z
            lib_sort_idx = np.argsort(lib_mz)
            lib_mz_sorted = lib_mz[lib_sort_idx]
            lib_abun_sorted = lib_abun[lib_sort_idx]
            
            # Build aligned vectors including all peaks
            vec_query = []
            vec_lib = []
            used_lib = np.zeros(len(lib_mz_sorted), dtype=bool)
            
            # For each query peak, find match in library or add as unmatched
            for i in range(len(query_mz_sorted)):
                closest_idx = find_closest(lib_mz_sorted, np.array([query_mz_sorted[i]]))[0]
                diff = abs(lib_mz_sorted[closest_idx] - query_mz_sorted[i])
                
                if diff <= self.ms2_tolerance_da:
                    # Matched peak
                    vec_query.append(query_abun_sorted[i])
                    vec_lib.append(lib_abun_sorted[closest_idx])
                    used_lib[closest_idx] = True
                else:
                    # Unmatched query peak
                    vec_query.append(query_abun_sorted[i])
                    vec_lib.append(0.0)
            
            # Add unmatched library peaks
            for j in range(len(lib_mz_sorted)):
                if not used_lib[j]:
                    vec_query.append(0.0)
                    vec_lib.append(lib_abun_sorted[j])
            
            # Convert to numpy arrays
            vec_query = np.array(vec_query, dtype=float)
            vec_lib = np.array(vec_lib, dtype=float)
            
            # Compute cosine similarity
            norm_query = np.linalg.norm(vec_query)
            norm_lib = np.linalg.norm(vec_lib)
            
            if norm_query == 0 or norm_lib == 0:
                continue
            
            cosine = np.dot(vec_query, vec_lib) / (norm_query * norm_lib)
            cosine_scores[lib_idx] = float(np.clip(cosine, 0.0, 1.0))
        
        return cosine_scores

    def _entropy_score_pair(
        self,
        spec_a,
        precursor_a: float | None,
        spec_b,
        precursor_b: float | None,
    ) -> float:
        """Compute entropy similarity between two spectra.

        We search spec_a against the library, then look up the library index
        that corresponds to spec_b.  Because the library is built from the
        same spectra, we can use the index directly.

        NOTE: This method is used for all-vs-all and new-vs-existing.
        The caller is responsible for mapping library indices to spectrum IDs.
        """
        raise NotImplementedError(
            "_entropy_score_pair is not used directly; "
            "use compute_all_vs_all or compute_new_vs_existing."
        )

    def _compute_entropy_matrix(
        self,
        spectra: list,
        precursor_mzs: list[float | None],
        lib_indices: list[int],
    ) -> np.ndarray:
        """Compute entropy similarity scores for all spectra against the library.

        Parameters
        ----------
        spectra : list
            Query spectra (objects with .mz_exp and .abundance).
        precursor_mzs : list of float or None
            Precursor m/z for each spectrum.
        lib_indices : list of int
            Library indices corresponding to each spectrum in *spectra*.
            Used to extract pairwise scores from the full search result vector.

        Returns
        -------
        np.ndarray of shape (len(spectra), len(spectra))
            Pairwise entropy similarity matrix.
        """
        n = len(spectra)
        scores = np.zeros((n, n), dtype=np.float32)
        for i, (spec, pmz) in enumerate(zip(spectra, precursor_mzs)):
            peaks = self._peaks_array(spec)
            if peaks.shape[0] == 0:
                continue
            result_vec = self._clean_and_search(peaks, pmz)
            # Extract scores for all spectra in our set
            for j, lib_idx in enumerate(lib_indices):
                if lib_idx is not None and lib_idx < len(result_vec):
                    scores[i, j] = result_vec[lib_idx]
        return scores

    def _compute_cosine_for_pairs(
        self,
        pairs: list[tuple[int, int]],
        spectra_a: list,
        spectra_b: list,
    ) -> dict[tuple[int, int], float]:
        """Compute cosine similarity for a list of (i, j) index pairs.

        Parameters
        ----------
        pairs : list of (i, j)
            Index pairs into spectra_a and spectra_b respectively.
        spectra_a, spectra_b : list
            Spectrum objects.

        Returns
        -------
        dict mapping (i, j) → cosine score
        """
        if not pairs:
            return {}

        # TODO: Add spectrum cleaning step before cosine calculation.
        # Currently, spectra are used directly without cleaning (e.g., via
        # clean_spectrum_for_search()). This differs from query-vs-library
        # cosine computation which cleans spectra first. Consider adding
        # cleaning here for consistency and to improve cosine accuracy.

        args = [
            (
                np.asarray(spectra_a[i].mz_exp, dtype=float),
                np.asarray(spectra_a[i].abundance, dtype=float),
                np.asarray(spectra_b[j].mz_exp, dtype=float),
                np.asarray(spectra_b[j].abundance, dtype=float),
                self.ms2_tolerance_da,
            )
            for i, j in pairs
        ]

        if self.use_parallel and len(args) > 1 and self.n_jobs > 1:
            with multiprocessing.Pool(min(self.n_jobs, len(args))) as pool:
                cosine_scores = pool.map(_compute_cosine_pair, args)
        else:
            cosine_scores = [_compute_cosine_pair(a) for a in args]

        return {pair: score for pair, score in zip(pairs, cosine_scores)}

    # ── Public API ────────────────────────────────────────────────────────────

    def compute_all_vs_all(
        self,
        spectra: list,
        spectrum_ids: list[str],
        precursor_mzs: list[float | None] | None = None,
        lib_indices: list[int | None] | None = None,
    ) -> dict[str, dict[tuple[str, str], float]]:
        """Compute all-vs-all pairwise similarities.

        Parameters
        ----------
        spectra : list
            Spectrum objects (must have .mz_exp and .abundance).
        spectrum_ids : list of str
            User-provided IDs, one per spectrum.
        precursor_mzs : list of float or None, optional
            Precursor m/z for each spectrum.  Required for ``"identity"`` and
            ``"neutral_loss"`` search types.  Ignored for ``"open"``.
        lib_indices : list of int or None, optional
            Index of each spectrum in the FlashEntropy library.  If None,
            entropy similarity is computed by searching each spectrum against
            the library and using the best match.  Providing correct indices
            gives exact pairwise scores.

        Returns
        -------
        dict
            ``{metric_name: {(id1, id2): score}}``
            where ``metric_name`` is ``"entropy_similarity"`` plus any
            additional metrics.  Only upper-triangle pairs are returned.
        """
        n = len(spectra)
        if n == 0:
            return {}

        if precursor_mzs is None:
            precursor_mzs = [None] * n
        if len(precursor_mzs) != n:
            raise ValueError("precursor_mzs must have the same length as spectra.")

        # ── Stage 1: entropy similarity ───────────────────────────────────────
        entropy_pairs: dict[tuple[str, str], float] = {}
        pairs_for_additional: list[tuple[int, int]] = []

        for i, j in combinations(range(n), 2):
            score = self._pairwise_entropy(
                spectra[i], precursor_mzs[i],
                spectra[j], precursor_mzs[j],
                lib_indices[i] if lib_indices else None,
                lib_indices[j] if lib_indices else None,
            )
            if score > 0.0:
                entropy_pairs[(spectrum_ids[i], spectrum_ids[j])] = score
                if score >= self.entropy_threshold_low:
                    pairs_for_additional.append((i, j))

        result: dict[str, dict[tuple[str, str], float]] = {
            "entropy_similarity": entropy_pairs
        }

        # ── Stage 2: additional metrics ───────────────────────────────────────
        for metric in self.additional_similarities:
            if metric == "cosine":
                cosine_idx_scores = self._compute_cosine_for_pairs(
                    pairs_for_additional, spectra, spectra
                )
                result["cosine"] = {
                    (spectrum_ids[i], spectrum_ids[j]): score
                    for (i, j), score in cosine_idx_scores.items()
                }

        return result

    def compute_new_vs_existing(
        self,
        new_spectra: list,
        new_ids: list[str],
        existing_spectra: list,
        existing_ids: list[str],
        new_precursor_mzs: list[float | None] | None = None,
        existing_precursor_mzs: list[float | None] | None = None,
        new_lib_indices: list[int | None] | None = None,
        existing_lib_indices: list[int | None] | None = None,
    ) -> dict[str, dict[tuple[str, str], float]]:
        """Compute similarities between new spectra and existing spectra.

        Also computes within-new-batch similarities.
        Does NOT recompute existing-vs-existing pairs.

        Parameters
        ----------
        new_spectra : list
            New spectrum objects.
        new_ids : list of str
            IDs for new spectra.
        existing_spectra : list
            Already-processed spectrum objects.
        existing_ids : list of str
            IDs for existing spectra.
        new_precursor_mzs : list of float or None, optional
            Precursor m/z for new spectra.
        existing_precursor_mzs : list of float or None, optional
            Precursor m/z for existing spectra.
        new_lib_indices, existing_lib_indices : list of int or None, optional
            Library indices for new/existing spectra.

        Returns
        -------
        dict
            ``{metric_name: {(id1, id2): score}}``
            Only new pairs (new-vs-existing and new-vs-new).
        """
        n_new = len(new_spectra)
        n_exist = len(existing_spectra)

        if new_precursor_mzs is None:
            new_precursor_mzs = [None] * n_new
        if existing_precursor_mzs is None:
            existing_precursor_mzs = [None] * n_exist

        entropy_pairs: dict[tuple[str, str], float] = {}
        pairs_for_additional_cross: list[tuple[int, int]] = []  # (new_i, exist_j)
        pairs_for_additional_new: list[tuple[int, int]] = []    # (new_i, new_j)

        # ── new vs existing ───────────────────────────────────────────────────
        for i in range(n_new):
            for j in range(n_exist):
                score = self._pairwise_entropy(
                    new_spectra[i], new_precursor_mzs[i],
                    existing_spectra[j], existing_precursor_mzs[j],
                    new_lib_indices[i] if new_lib_indices else None,
                    existing_lib_indices[j] if existing_lib_indices else None,
                )
                if score > 0.0:
                    entropy_pairs[(new_ids[i], existing_ids[j])] = score
                    if score >= self.entropy_threshold_low:
                        pairs_for_additional_cross.append((i, j))

        # ── new vs new ────────────────────────────────────────────────────────
        for i, j in combinations(range(n_new), 2):
            score = self._pairwise_entropy(
                new_spectra[i], new_precursor_mzs[i],
                new_spectra[j], new_precursor_mzs[j],
                new_lib_indices[i] if new_lib_indices else None,
                new_lib_indices[j] if new_lib_indices else None,
            )
            if score > 0.0:
                entropy_pairs[(new_ids[i], new_ids[j])] = score
                if score >= self.entropy_threshold_low:
                    pairs_for_additional_new.append((i, j))

        result: dict[str, dict[tuple[str, str], float]] = {
            "entropy_similarity": entropy_pairs
        }

        # ── Stage 2: additional metrics ───────────────────────────────────────
        for metric in self.additional_similarities:
            if metric == "cosine":
                metric_scores: dict[tuple[str, str], float] = {}

                # cross-batch cosine
                cross_cosine = self._compute_cosine_for_pairs(
                    pairs_for_additional_cross, new_spectra, existing_spectra
                )
                for (i, j), score in cross_cosine.items():
                    metric_scores[(new_ids[i], existing_ids[j])] = score

                # within-new cosine
                new_cosine = self._compute_cosine_for_pairs(
                    pairs_for_additional_new, new_spectra, new_spectra
                )
                for (i, j), score in new_cosine.items():
                    metric_scores[(new_ids[i], new_ids[j])] = score

                result["cosine"] = metric_scores

        return result

    def _pairwise_entropy(
        self,
        spec_a,
        pmz_a: float | None,
        spec_b,
        pmz_b: float | None,
        lib_idx_a: int | None,
        lib_idx_b: int | None,
        fe_lib_override=None,
    ) -> float:
        """Compute entropy similarity between two spectra.

        Strategy: search spec_a against the library and read off the score
        at lib_idx_b (the library position of spec_b).  If lib_idx_b is None,
        fall back to searching spec_b and reading lib_idx_a.

        If neither library index is known, returns 0.0 (cannot compute without
        library indices for the query-vs-query case).

        Parameters
        ----------
        fe_lib_override : optional
            If provided, use this FlashEntropy library instead of ``self.fe_lib``.
            Used for query-vs-query searches where a temporary index is built
            from the query spectra themselves.
        """
        peaks_a = self._peaks_array(spec_a)
        if peaks_a.shape[0] == 0:
            return 0.0

        if lib_idx_b is not None:
            result_vec = self._clean_and_search(peaks_a, pmz_a, fe_lib_override=fe_lib_override)
            if lib_idx_b < len(result_vec):
                return float(result_vec[lib_idx_b])
            return 0.0

        if lib_idx_a is not None:
            peaks_b = self._peaks_array(spec_b)
            if peaks_b.shape[0] == 0:
                return 0.0
            result_vec = self._clean_and_search(peaks_b, pmz_b, fe_lib_override=fe_lib_override)
            if lib_idx_a < len(result_vec):
                return float(result_vec[lib_idx_a])
            return 0.0

        # No library indices available – cannot compute entropy similarity
        return 0.0
    
    def build_fe_index_from_spectra(
        self,
        spectra: list,
        precursor_mzs: list[float | None] | None = None,
        fe_kwargs: dict | None = None,
    ):
        """Build a FlashEntropy search index from a list of spectrum objects.

        This is used to create a temporary index from query spectra so that
        query-vs-query similarities can be computed without requiring the
        spectra to already be in the main library.

        Parameters
        ----------
        spectra : list
            Spectrum objects with ``.mz_exp`` and ``.abundance`` attributes.
        precursor_mzs : list of float or None, optional
            Precursor m/z for each spectrum.  Used as the ``precursor_mz``
            field in the FlashEntropy library.  Defaults to 0.0 for each.
        fe_kwargs : dict, optional
            Extra keyword arguments forwarded to ``FlashEntropySearch``.
            Defaults to the same settings used for the main library.

        Returns
        -------
        ms_entropy.FlashEntropySearch
            A new FlashEntropy search instance indexed on *spectra*.
        """
        try:
            from ms_entropy import FlashEntropySearch
        except ImportError:
            raise ImportError(
                "ms_entropy is required for build_fe_index_from_spectra(). "
                "Install with: pip install ms_entropy"
            )

        if precursor_mzs is None:
            precursor_mzs = [None] * len(spectra)

        # Default FE kwargs mirror the main library settings
        default_fe_kwargs = {
            "normalize_intensity": True,
            "min_ms2_difference_in_da": self.peak_sep_da * 2,
            "max_ms2_tolerance_in_da": self.ms2_tolerance_da,
            "max_indexed_mz": 3000,
            "precursor_ions_removal_da": None,
            "noise_threshold": 0,
        }
        if fe_kwargs:
            default_fe_kwargs.update(fe_kwargs)

        # Build the spectral library list
        spectral_library = []
        for i, (spec, pmz) in enumerate(zip(spectra, precursor_mzs)):
            peaks = self._peaks_array(spec)
            if peaks.shape[0] == 0:
                continue
            spectral_library.append({
                "id": i,
                "precursor_mz": float(pmz) if pmz is not None else 0.0,
                "peaks": peaks.tolist(),
            })

        # Match the FlashEntropy constructor/init and build_index signature
        fe_init_kws = [
            "max_ms2_tolerance_in_da",
            "mz_index_step",
            "low_memory",
            "path_data",
        ]
        fe_init_kws = {k: v for k, v in default_fe_kwargs.items() if k in fe_init_kws}

        fe = FlashEntropySearch(**fe_init_kws)

        fe_index_kws = [
            "max_indexed_mz",
            "precursor_ions_removal_da",
            "noise_threshold",
            "min_ms2_difference_in_da",
            "max_peak_num",
        ]
        fe_index_kws = {k: v for k, v in default_fe_kwargs.items() if k in fe_index_kws}

        fe.build_index(spectral_library, **fe_index_kws, clean_spectra=True)
        return fe

    def compute_library_vs_library_filtered(
        self,
        library_indices: list[int],
        spectrum_ids: list[str],
        precursor_mzs: list[float | None] | None = None,
    ) -> dict[str, dict[tuple[str, str], float]]:
        """Compute library-vs-library similarities for a filtered subset of library spectra.

        Extracts the spectra at *library_indices* from ``self.fe_lib``, builds a
        temporary FlashEntropy index from them, and computes all-vs-all pairwise
        similarities within that subset.

        Parameters
        ----------
        library_indices : list of int
            Indices into ``self.fe_lib`` of the library spectra to include.
        spectrum_ids : list of str
            IDs to assign to each library spectrum (must be same length as
            *library_indices*).
        precursor_mzs : list of float or None, optional
            Precursor m/z for each library spectrum.  If None, attempts to read
            from ``self.fe_lib`` (attribute ``precursor_mz`` on each entry), or
            defaults to 0.0.

        Returns
        -------
        dict
            ``{metric_name: {(id1, id2): score}}``
            All-vs-all pairs within the filtered library subset.
        """
        n = len(library_indices)
        if n == 0:
            return {}

        # Extract spectra from the FE library
        # ms_entropy stores spectra as a list of dicts with 'peaks' and 'precursor_mz'
        lib_spectra_raw = getattr(self.fe_lib, "spectra", None) or getattr(self.fe_lib, "library", None)
        if lib_spectra_raw is None:
            return {}

        # Build lightweight spectrum objects from the raw library entries
        class _LibSpec:
            __slots__ = ("mz_exp", "abundance")
            def __init__(self, peaks_arr):
                self.mz_exp = peaks_arr[:, 0]
                self.abundance = peaks_arr[:, 1]

        lib_spectra: list = []
        lib_precursor_mzs: list[float | None] = []

        for li in library_indices:
            if li >= len(lib_spectra_raw):
                continue
            entry = lib_spectra_raw[li]
            peaks = np.asarray(entry.get("peaks", []), dtype=float)
            if peaks.ndim != 2 or peaks.shape[1] < 2 or peaks.shape[0] == 0:
                continue
            lib_spectra.append(_LibSpec(peaks))
            if precursor_mzs is not None:
                lib_precursor_mzs.append(precursor_mzs[len(lib_spectra) - 1])
            else:
                lib_precursor_mzs.append(float(entry.get("precursor_mz", 0.0) or 0.0))

        if not lib_spectra:
            return {}

        # Trim spectrum_ids to match successfully extracted spectra
        valid_ids = spectrum_ids[: len(lib_spectra)]

        # Build a temporary FE index from the filtered library spectra
        temp_fe = self.build_fe_index_from_spectra(
            spectra=lib_spectra,
            precursor_mzs=lib_precursor_mzs,
        )

        # Compute all-vs-all using the temporary index
        return self.compute_all_vs_all_with_lib(
            spectra=lib_spectra,
            spectrum_ids=valid_ids,
            precursor_mzs=lib_precursor_mzs,
            fe_lib_override=temp_fe,
        )

    def compute_all_vs_all_with_lib(
        self,
        spectra: list,
        spectrum_ids: list[str],
        precursor_mzs: list[float | None] | None = None,
        fe_lib_override=None,
    ) -> dict[str, dict[tuple[str, str], float]]:
        """Compute all-vs-all pairwise similarities using a provided FE library.

        This is the same as :meth:`compute_all_vs_all` but uses
        *fe_lib_override* instead of ``self.fe_lib``.  The library indices
        are derived from the order of *spectra* (index 0, 1, 2, …).

        Parameters
        ----------
        spectra : list
            Spectrum objects.
        spectrum_ids : list of str
            User-provided IDs, one per spectrum.
        precursor_mzs : list of float or None, optional
            Precursor m/z for each spectrum.
        fe_lib_override : ms_entropy.FlashEntropySearch, optional
            FlashEntropy library to search against.  If None, falls back to
            ``self.fe_lib``.

        Returns
        -------
        dict
            ``{metric_name: {(id1, id2): score}}``
        """
        n = len(spectra)
        if n == 0:
            return {}

        if precursor_mzs is None:
            precursor_mzs = [None] * n

        # Library indices are simply 0..n-1 (positions in fe_lib_override)
        lib_indices = list(range(n))

        # Temporarily swap the FE library if an override is provided
        original_fe_lib = self.fe_lib
        if fe_lib_override is not None:
            self.fe_lib = fe_lib_override

        try:
            entropy_pairs: dict[tuple[str, str], float] = {}
            pairs_for_additional: list[tuple[int, int]] = []

            for i, j in combinations(range(n), 2):
                score = self._pairwise_entropy(
                    spectra[i], precursor_mzs[i],
                    spectra[j], precursor_mzs[j],
                    lib_indices[i],
                    lib_indices[j],
                )
                if score > 0.0:
                    entropy_pairs[(spectrum_ids[i], spectrum_ids[j])] = score
                    if score >= self.entropy_threshold_low:
                        pairs_for_additional.append((i, j))

            result: dict[str, dict[tuple[str, str], float]] = {
                "entropy_similarity": entropy_pairs
            }

            for metric in self.additional_similarities:
                if metric == "cosine":
                    cosine_idx_scores = self._compute_cosine_for_pairs(
                        pairs_for_additional, spectra, spectra
                    )
                    result["cosine"] = {
                        (spectrum_ids[i], spectrum_ids[j]): score
                        for (i, j), score in cosine_idx_scores.items()
                    }
        finally:
            # Always restore the original FE library
            self.fe_lib = original_fe_lib

        return result
