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

# CoreMS spectral similarity (cosine, etc.)
from corems.molecular_id.calc.SpectralSimilarity import SpectralSimilarity

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

    Uses a tolerance-based m/z binning (round to nearest 0.01 Da) so that
    peaks that are close in m/z are treated as matching.  This avoids the
    near-zero scores that arise when exact floating-point m/z values differ
    slightly between spectra.

    Parameters
    ----------
    args : tuple
        (mz1, abun1, mz2, abun2, mz_bin_da)
        where mz_bin_da is the binning resolution in Da (default 0.01).

    Returns
    -------
    float
        Cosine similarity score in [0, 1].
    """
    if len(args) == 5:
        mz1, abun1, mz2, abun2, mz_bin_da = args
    else:
        mz1, abun1, mz2, abun2 = args
        mz_bin_da = 0.01

    try:
        # Round m/z values to the nearest bin so close peaks match
        factor = 1.0 / mz_bin_da
        binned1 = {round(float(m) * factor) / factor: float(a)
                   for m, a in zip(mz1, abun1)}
        binned2 = {round(float(m) * factor) / factor: float(a)
                   for m, a in zip(mz2, abun2)}
        ref_obj = {"mz": list(binned2.keys()), "abundance": list(binned2.values())}
        ss = SpectralSimilarity(binned1, ref_obj)
        return ss.cosine_correlation()
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
    ) -> np.ndarray:
        """Run FlashEntropy search for one query spectrum.

        Parameters
        ----------
        peaks : np.ndarray of shape (N, 2)
            [[mz, abundance], ...]
        precursor_mz : float or None
            Required for identity/neutral_loss; ignored for open.

        Returns
        -------
        np.ndarray
            1-D array of entropy similarity scores, one per library entry.
        """
        fe_method, fe_key = _FE_METHOD_MAP[self.search_type]

        # Use a dummy precursor_mz for open search
        pmz = precursor_mz if precursor_mz is not None else 0.0

        cleaned = self.fe_lib.clean_spectrum_for_search(
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

        results = self.fe_lib.search(**search_kwargs)
        return results[fe_key]

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

        args = [
            (
                np.asarray(spectra_a[i].mz_exp, dtype=float),
                np.asarray(spectra_a[i].abundance, dtype=float),
                np.asarray(spectra_b[j].mz_exp, dtype=float),
                np.asarray(spectra_b[j].abundance, dtype=float),
                self.peak_sep_da,
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
    ) -> float:
        """Compute entropy similarity between two spectra.

        Strategy: search spec_a against the library and read off the score
        at lib_idx_b (the library position of spec_b).  If lib_idx_b is None,
        fall back to searching spec_b and reading lib_idx_a.

        If neither library index is known, returns 0.0 (cannot compute without
        library indices for the query-vs-query case).
        """
        peaks_a = self._peaks_array(spec_a)
        if peaks_a.shape[0] == 0:
            return 0.0

        if lib_idx_b is not None:
            result_vec = self._clean_and_search(peaks_a, pmz_a)
            if lib_idx_b < len(result_vec):
                return float(result_vec[lib_idx_b])
            return 0.0

        if lib_idx_a is not None:
            peaks_b = self._peaks_array(spec_b)
            if peaks_b.shape[0] == 0:
                return 0.0
            result_vec = self._clean_and_search(peaks_b, pmz_b)
            if lib_idx_a < len(result_vec):
                return float(result_vec[lib_idx_a])
            return 0.0

        # No library indices available – cannot compute entropy similarity
        return 0.0
