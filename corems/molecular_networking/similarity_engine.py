"""
SimilarityEngine
================

Computes pairwise spectral similarities between query spectra using:

1. **FlashEntropy** (fast, vectorised) — always computed first.
2. **Additional metrics** (e.g. cosine) — computed only for pairs that
   pass a low entropy-similarity gate (*entropy_threshold_low*).

Supported FlashEntropy search modes
------------------------------------
``"identity"``
    Precursor-matched search.  Library candidates are filtered by precursor
    m/z within *ms1_tolerance_da* before scoring.  Requires precursor m/z
    for every query spectrum.

``"open"``
    No precursor filtering.  All library entries are scored against every
    query.  Precursor m/z is accepted but ignored for filtering.

``"neutral_loss"``
    Matching occurs in neutral-loss mass space
    (``precursor_mz − fragment_mz``).  Requires precursor m/z for every
    query spectrum.  The precursor filter is **disabled** (set to 1e9 Da)
    because neutral-loss matching does not operate on the precursor itself.

"""

from __future__ import annotations

from collections.abc import Callable
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


def _sort_peaks_by_mz(mz, abun):
    """Return float m/z and abundance arrays sorted by m/z (same argsort as cosine)."""
    mz = np.asarray(mz, dtype=float)
    abun = np.asarray(abun, dtype=float)
    if mz.size == 0:
        return mz, abun
    idx = np.argsort(mz)
    return mz[idx], abun[idx]


def _align_and_compute_cosine_presorted(mz1, abun1, mz2, abun2, tolerance_da):
    """Cosine similarity assuming both spectra are already sorted by m/z.

    Internal hot path for batch cosine after unique spectra have been
    pre-sorted once.  Callers must pass m/z-sorted peak lists (e.g. from
    :func:`_sort_peaks_by_mz`).  Behaviour matches sorting then aligning.

    Parameters
    ----------
    mz1, abun1, mz2, abun2 : array-like
        Peak arrays **sorted by m/z ascending** within each spectrum.
    tolerance_da : float
        m/z matching tolerance in Da.

    Returns
    -------
    float
        Cosine similarity in [0, 1].  Empty spectra return 0.0.
    """
    mz1 = np.asarray(mz1, dtype=float)
    abun1 = np.asarray(abun1, dtype=float)
    mz2 = np.asarray(mz2, dtype=float)
    abun2 = np.asarray(abun2, dtype=float)

    n1 = mz1.size
    n2 = mz2.size
    if n1 == 0 or n2 == 0:
        return 0.0

    # Vectorised closest-peak search (both arrays sorted, as required by
    # find_closest).  Same per-peak result as calling find_closest once per
    # spectrum-1 peak.
    closest_idx = find_closest(mz2, mz1)
    diffs = np.abs(mz2[closest_idx] - mz1)
    matched = diffs <= float(tolerance_da)

    # Aligned vectors: first n1 entries follow spectrum-1 peak order (matched
    # or unmatched).  Same layout as the original list-append loop.
    vec1 = np.empty(n1, dtype=float)
    vec2 = np.empty(n1, dtype=float)
    vec1[:] = abun1
    vec2[:] = 0.0
    if np.any(matched):
        vec2[matched] = abun2[closest_idx[matched]]

    used_spec2 = np.zeros(n2, dtype=bool)
    if np.any(matched):
        used_spec2[closest_idx[matched]] = True

    unmatched2 = ~used_spec2
    n_unmatched2 = int(np.count_nonzero(unmatched2))
    if n_unmatched2:
        vec1 = np.concatenate((vec1, np.zeros(n_unmatched2, dtype=float)))
        vec2 = np.concatenate((vec2, abun2[unmatched2]))

    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    if norm1 == 0.0 or norm2 == 0.0:
        return 0.0

    cosine = float(np.dot(vec1, vec2) / (norm1 * norm2))
    return float(np.clip(cosine, 0.0, 1.0))


def _align_and_compute_cosine(mz1, abun1, mz2, abun2, tolerance_da):
    """Align two spectra by m/z tolerance and compute cosine similarity.

    Sorts both spectra, then delegates to
    :func:`_align_and_compute_cosine_presorted`.  Numerically identical to
    the original sequential alignment:

    1. Sort both spectra by m/z.
    2. For each peak in spectrum 1, take the closest m/z in spectrum 2
       (:func:`~corems.mass_spectra.calc.lc_calc.find_closest`).
    3. If within *tolerance_da*, pair abundances; otherwise pair spectrum-1
       abundance with 0.  The same spectrum-2 peak may be paired more than
       once (matching does not skip already-used peaks).
    4. Append unused spectrum-2 peaks as (0, abundance).
    5. Cosine of the two aligned vectors, clipped to [0, 1].

    For batch work over unique spectra, pre-sort once with
    :func:`_sort_peaks_by_mz` and call
    :func:`_align_and_compute_cosine_presorted` to avoid repeated sorts.

    Parameters
    ----------
    mz1, abun1 : array-like
        m/z and abundance arrays for spectrum 1 (will be sorted internally)
    mz2, abun2 : array-like
        m/z and abundance arrays for spectrum 2 (will be sorted internally)
    tolerance_da : float
        m/z matching tolerance in Da

    Returns
    -------
    float
        Cosine similarity score in [0, 1].  Empty spectra return 0.0.
        Unexpected input errors are not swallowed (they propagate).
    """
    mz1_s, abun1_s = _sort_peaks_by_mz(mz1, abun1)
    mz2_s, abun2_s = _sort_peaks_by_mz(mz2, abun2)
    return _align_and_compute_cosine_presorted(
        mz1_s, abun1_s, mz2_s, abun2_s, tolerance_da
    )


class SimilarityEngine:
    """Compute pairwise spectral similarities using FlashEntropy + optional extras.

    Parameters
    ----------
    fe_lib : ms_entropy.FlashEntropySearch or None
        Pre-built FlashEntropy search instance.  Pass ``None`` for
        query-only mode (only :meth:`compute_all_vs_all_with_lib` with a
        temporary index is used; query-vs-library methods are unavailable).
        When provided, fragment-tolerance parameters are extracted from this
        library unless overridden by *ms2_tolerance_da*.
    search_type : str
        FlashEntropy search mode: ``"identity"``, ``"open"``, or
        ``"neutral_loss"``.  Default ``"open"`` (matches typical DDA
        molecular-networking use; use ``"identity"`` when precursor
        filtering is required).
    additional_similarities : list of str, optional
        Extra similarity metrics to compute for pairs that pass the entropy
        threshold.  Currently supported: ``["cosine"]``.
        Default ``["cosine"]``.  Cosine is a best-effort peak-aligned
        score gated by *entropy_threshold_low* (not a second primary engine).
    ms1_tolerance_da : float, optional
        Precursor m/z tolerance (Da) used **only** for ``"identity"`` search
        to filter library candidates by precursor m/z.  Ignored for
        ``"open"`` and ``"neutral_loss"`` (precursor filter is disabled for
        those modes).  If ``None`` (default), uses 0.01 Da.
    ms2_tolerance_da : float, optional
        Fragment m/z tolerance (Da) for spectrum cleaning and similarity
        scoring.  Resolution priority: explicit kwarg > extracted from
        *fe_lib* > 0.01 Da fallback.
    entropy_threshold_low : float
        Minimum entropy similarity score required to trigger additional
        metric computation.  Default 0.1.

    Attributes
    ----------
    ms2_tolerance_da : float
        Fragment m/z tolerance resolved at construction time.
    peak_sep_da : float
        Minimum peak separation used during spectrum cleaning
        (``2 * ms2_tolerance_da``).
    ms1_tolerance_da : float
        Precursor m/z tolerance for identity search.
    entropy_threshold_low : float
        Low-entropy gate for triggering additional metric computation.

    Notes
    -----
    **Tolerance resolution order** (highest priority first):

    1. Explicit *ms2_tolerance_da* kwarg.
    2. ``fe_lib.entropy_search.max_ms2_tolerance_in_da`` (when *fe_lib* is
       not ``None``).
    3. Hard-coded fallback of 0.01 Da.

    **Precursor filter behaviour by search type:**

    - ``"identity"`` : library candidates are filtered by precursor m/z
      within *ms1_tolerance_da* before scoring.
    - ``"open"`` : precursor filter is disabled (``ms1_tolerance_in_da``
      passed to FlashEntropy is set to ``1e9``).
    - ``"neutral_loss"`` : precursor filter is disabled for the same reason
      as ``"open"``; matching occurs in neutral-loss mass space
      (``precursor_mz − fragment_mz``), not on the precursor itself.

    Cosine (when requested) runs **single-threaded**: unique spectra are
    pre-sorted once, then pairs are scored sequentially via
    :func:`_align_and_compute_cosine_presorted`.
    """

    def __init__(
        self,
        fe_lib,
        search_type: str = "open",
        additional_similarities: list[str] | None = None,
        ms1_tolerance_da: float | None = None,
        ms2_tolerance_da: float | None = None,
        entropy_threshold_low: float = 0.1,
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

        # Tolerance priority: explicit kwarg > fe_lib extraction > 0.01 fallback.
        # When fe_lib is None (queries-only mode), a temporary FE index is built
        # from the query spectra themselves; ms2_tolerance_da controls that index.
        if ms2_tolerance_da is not None:
            self.ms2_tolerance_da = ms2_tolerance_da
        elif fe_lib is not None:
            self.ms2_tolerance_da = fe_lib.entropy_search.max_ms2_tolerance_in_da
        else:
            self.ms2_tolerance_da = 0.01
        self.peak_sep_da = 2 * self.ms2_tolerance_da
        self.ms1_tolerance_da = ms1_tolerance_da if ms1_tolerance_da is not None else 0.01

        self.entropy_threshold_low = entropy_threshold_low

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _peaks_array(self, spectrum) -> np.ndarray:
        """Return (N, 2) peaks array from a spectrum object."""
        mz = np.asarray(spectrum.mz_exp, dtype=float)
        ab = np.asarray(spectrum.abundance, dtype=float)
        return np.column_stack((mz, ab))

    def _require_precursor_mz(
        self,
        precursor_mz: float | None,
        context: str,
    ) -> float:
        """Validate precursor m/z and return it as float."""
        if precursor_mz is None:
            raise ValueError(
                f"{context} is missing precursor_mz, which is required for "
                f"search_type='{self.search_type}'."
            )

        precursor = float(precursor_mz)
        if not np.isfinite(precursor):
            raise ValueError(
                f"{context} has non-finite precursor_mz={precursor_mz!r}."
            )
        return precursor

    def _transform_to_neutral_loss_space(
        self,
        peaks: np.ndarray,
        precursor_mz: float | None,
        context: str,
    ) -> np.ndarray:
        """Transform peaks from fragment m/z to neutral-loss mass space.

        Neutral-loss mass is computed as:
            neutral_loss_mz = precursor_mz - fragment_mz
        """
        precursor = self._require_precursor_mz(precursor_mz, context)
        peaks = np.asarray(peaks, dtype=float)

        if peaks.ndim != 2 or peaks.shape[1] < 2:
            raise ValueError(
                f"{context} peaks must be a 2D array with at least 2 columns; "
                f"got shape={peaks.shape}."
            )

        if peaks.shape[0] == 0:
            return peaks.copy()

        # Neutral-loss masses are only valid for fragment_mz <= precursor_mz.
        # Drop invalid fragment peaks to mirror the neutral-loss search domain.
        valid_mask = peaks[:, 0] <= precursor
        if not np.any(valid_mask):
            return np.empty((0, peaks.shape[1]), dtype=float)

        nl_peaks = peaks[valid_mask].copy()
        nl_peaks[:, 0] = precursor - nl_peaks[:, 0]

        # Keep deterministic ordering before downstream alignment.
        sort_idx = np.argsort(nl_peaks[:, 0])
        return nl_peaks[sort_idx]

    def _clean_and_search(
        self,
        peaks: np.ndarray,
        precursor_mz: float | None,
        fe_lib_override=None,
    ) -> np.ndarray:
        """Run FlashEntropy search for one query spectrum.

        Passes raw peaks directly to FlashEntropy's ``search()`` method,
        which cleans the spectrum internally.  This avoids double-cleaning
        and ensures the query is processed identically to library entries.

        **Precursor filter behaviour:**

        - ``"identity"`` : ``ms1_tolerance_in_da`` is set to
          :attr:`ms1_tolerance_da` so only library entries within that
          precursor window are scored.
        - ``"open"`` and ``"neutral_loss"`` : ``ms1_tolerance_in_da`` is
          set to ``1e9`` (effectively no filter).  For ``"neutral_loss"``,
          matching occurs in neutral-loss mass space internally within
          FlashEntropy, so filtering on precursor m/z would incorrectly
          exclude valid matches.

        For ``"neutral_loss"``, fragment peaks with m/z > precursor m/z are
        stripped before the search call (they have no valid neutral-loss
        representation).

        Parameters
        ----------
        peaks : np.ndarray of shape (N, 2)
            Raw ``[[mz, abundance], ...]`` array.  FlashEntropy cleans this
            internally (normalisation, denoising, peak separation).
        precursor_mz : float or None
            Required for ``"identity"`` and ``"neutral_loss"``; a dummy
            value of 0.0 is used for ``"open"`` when ``None`` is passed.
        fe_lib_override : ms_entropy.FlashEntropySearch, optional
            If provided, search against this library instead of
            ``self.fe_lib``.  Used for query-vs-query searches where a
            temporary index is built from the query spectra themselves.

        Returns
        -------
        np.ndarray
            1-D array of entropy similarity scores, one per library entry.
            Returns a zero array if the query has no valid peaks after
            pre-filtering.
        """
        fe_method, fe_key = _FE_METHOD_MAP[self.search_type]
        fe = fe_lib_override if fe_lib_override is not None else self.fe_lib

        if self.search_type in ("identity", "neutral_loss"):
            pmz = self._require_precursor_mz(precursor_mz, "Query spectrum")
        else:
            # Use a dummy precursor_mz for open search
            pmz = precursor_mz if precursor_mz is not None else 0.0

        peaks_for_search = peaks
        if self.search_type == "neutral_loss":
            peaks_for_search = np.asarray(peaks, dtype=float)
            if peaks_for_search.ndim != 2 or peaks_for_search.shape[1] < 2:
                raise ValueError(
                    f"Query spectrum peaks must be 2D with at least 2 columns for neutral_loss search; "
                    f"got shape={peaks_for_search.shape}."
                )
            peaks_for_search = peaks_for_search[peaks_for_search[:, 0] <= pmz]
            if peaks_for_search.shape[0] == 0:
                try:
                    lib_size = len(fe)
                except TypeError:
                    lib_size = len(getattr(fe, "precursor_mz_array", []))
                return np.zeros(lib_size, dtype=float)

        # Pass raw peaks to search() - it will clean them internally
        # This ensures consistent cleaning between library and query spectra
        search_kwargs: dict[str, Any] = dict(
            peaks=peaks_for_search,  # Raw peaks - search() will clean them
            ms2_tolerance_in_da=self.ms2_tolerance_da,
            method={fe_method},
            precursor_ions_removal_da=None,
            noise_threshold=0.0, #TODO: get this as an attribute from FE as well
            min_ms2_difference_in_da=self.peak_sep_da,
            target="cpu",
        )
        if self.search_type == "identity":
            search_kwargs["precursor_mz"] = pmz
            search_kwargs["ms1_tolerance_in_da"] = self.ms1_tolerance_da
        else:
            # open search and neutral_loss: precursor_mz still required by ms_entropy
            # API but precursor filtering is disabled (1e9 = no filter).
            # For neutral_loss, matching occurs on neutral-loss fragment space, not
            # on precursor m/z, so the precursor filter must be disabled.
            search_kwargs["precursor_mz"] = pmz
            search_kwargs["ms1_tolerance_in_da"] = 1e9  # effectively no filter

        results = fe.search(**search_kwargs)
        return results[fe_key]

    def _compute_cosine_for_query_vs_library(
        self,
        query_spectrum,
        query_precursor_mz: float | None,
        library_indices: list[int],
    ) -> dict[int, float]:
        """Compute cosine similarity for one query against multiple library spectra.

        Extracts cleaned peaks directly from the FlashEntropy library to ensure
        consistency with entropy similarity calculations.  For
        ``"neutral_loss"`` search, both query and library peaks are
        transformed to neutral-loss mass space before alignment.

        Parameters
        ----------
        query_spectrum : spectrum object
            Query spectrum with ``.mz_exp`` and ``.abundance`` attributes.
        query_precursor_mz : float or None
            Precursor m/z for the query.  Required for ``"neutral_loss"``
            search; used as a dummy 0.0 for ``"open"`` when ``None``.
        library_indices : list of int
            Indices into ``self.fe_lib`` of the library spectra to score.

        Returns
        -------
        dict mapping int → float
            ``{library_idx: cosine_score}`` for all pairs with score > 0.
        """
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
        cleaned_query = np.asarray(cleaned_query, dtype=float)
        
        if cleaned_query.shape[0] == 0:
            return {}

        if self.search_type == "neutral_loss":
            cleaned_query = self._transform_to_neutral_loss_space(
                cleaned_query,
                query_precursor_mz,
                "Query spectrum",
            )
        
        # Sort query once; each library spectrum is sorted once before cosine
        query_mz_sorted, query_abun_sorted = _sort_peaks_by_mz(
            cleaned_query[:, 0], cleaned_query[:, 1]
        )

        cosine_scores = {}
        for lib_idx in library_indices:
            # Extract cleaned peaks from FE library (same as entropy similarity uses)
            try:
                lib_entry = self.fe_lib[lib_idx]
                lib_peaks = np.asarray(lib_entry["peaks"], dtype=float)
            except (IndexError, KeyError, TypeError):
                continue

            if self.search_type == "neutral_loss":
                lib_peaks = self._transform_to_neutral_loss_space(
                    lib_peaks,
                    lib_entry.get("precursor_mz"),
                    f"Library spectrum at index {lib_idx}",
                )

            if lib_peaks.ndim != 2 or lib_peaks.shape[1] < 2 or lib_peaks.shape[0] == 0:
                continue

            lib_mz_sorted, lib_abun_sorted = _sort_peaks_by_mz(
                lib_peaks[:, 0], lib_peaks[:, 1]
            )

            cosine = _align_and_compute_cosine_presorted(
                query_mz_sorted,
                query_abun_sorted,
                lib_mz_sorted,
                lib_abun_sorted,
                self.ms2_tolerance_da,
            )
            if cosine > 0.0:
                cosine_scores[lib_idx] = cosine

        return cosine_scores

    def _compute_entropy_matrix_and_pairs(
        self,
        spectra: list,
        spectrum_ids: list[str],
        precursor_mzs: list[float | None],
        fe_lib_override=None,
    ) -> tuple[dict[tuple[str, str], float], list[tuple[int, int]]]:
        """Compute entropy similarity and extract upper-triangle pairs.

        Uses FlashEntropy's built-in search to efficiently compute pairwise
        similarities. Each spectrum is searched once against the library,
        and pairwise scores are extracted from the result vectors.
        
        Streams results directly to sparse storage without building dense matrix.

        Parameters
        ----------
        spectra : list
            Spectrum objects with .mz_exp and .abundance attributes.
        spectrum_ids : list of str
            User-provided IDs, one per spectrum.
        precursor_mzs : list of float or None
            Precursor m/z for each spectrum.
        fe_lib_override : FlashEntropySearch, optional
            FlashEntropy library to search against. If None, uses self.fe_lib.

        Returns
        -------
        entropy_pairs : dict mapping (id1, id2) → score
            Upper-triangle pairs with score > 0.0
        pairs_for_additional : list of (i, j)
            Index pairs with score >= entropy_threshold_low (for additional metrics)
        """
        n = len(spectra)
        fe = fe_lib_override if fe_lib_override is not None else self.fe_lib
        
        # Stream results directly to sparse storage (no dense intermediate)
        entropy_pairs: dict[tuple[str, str], float] = {}
        pairs_for_additional: list[tuple[int, int]] = []
        
        for i, (spec, pmz) in enumerate(zip(spectra, precursor_mzs)):
            peaks = self._peaks_array(spec)
            if peaks.shape[0] == 0:
                continue
            # Search against the FE library
            result_vec = self._clean_and_search(peaks, pmz, fe_lib_override=fe)
            
            if result_vec is None or len(result_vec) == 0:
                continue
            
            # Extract all pairs from result_vec, but only store upper-triangle (j > i)
            # result_vec[j] contains the similarity score between spectrum i and library entry j
            # When searching within same library (all-vs-all), result_vec has n entries
            for j in range(min(len(result_vec), n)):
                if j <= i:
                    # Skip lower triangle and diagonal to avoid duplicates
                    continue
                score = float(result_vec[j])
                if score > 0.0:
                    entropy_pairs[(spectrum_ids[i], spectrum_ids[j])] = score
                    if score >= self.entropy_threshold_low:
                        pairs_for_additional.append((i, j))

        return entropy_pairs, pairs_for_additional

    def _compute_cosine_for_pairs(
        self,
        pairs: list[tuple[int, int]],
        lib_indices: list[int | None],
        fe_lib_override = None,
    ) -> dict[tuple[int, int], float]:
        """Compute cosine similarity for a list of (i, j) index pairs.

        Extracts cleaned peaks directly from the FlashEntropy library for all
        unique spectrum indices referenced in the pairs list.

        Parameters
        ----------
        pairs : list of (i, j)
            Index pairs. Both i and j are indices into the lib_indices list.
        lib_indices : list of int or None
            FlashEntropy library indices. For each index in pairs, lib_indices[index]
            gives the FE library position. Cleaned peaks are extracted via
            fe_lib[lib_indices[idx]]["peaks"]. All spectra must be indexed.
        fe_lib_override : FlashEntropySearch, optional
            FlashEntropy library to use. If None, uses self.fe_lib.

        Returns
        -------
        dict mapping (i, j) → cosine score
        """
        if not pairs:
            return {}

        # Get FE library
        fe = fe_lib_override if fe_lib_override is not None else self.fe_lib

        # Extract unique indices from all pairs
        unique_indices = {i for i, _ in pairs} | {j for _, j in pairs}

        # Extract cleaned peaks once per unique spectrum and sort by m/z once
        # so pair scoring can use _align_and_compute_cosine_presorted.
        sorted_peaks: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        for idx in unique_indices:
            if lib_indices[idx] is not None:
                spec_dict = fe[lib_indices[idx]]
                peaks = np.asarray(spec_dict["peaks"], dtype=float)
                if peaks.ndim != 2 or peaks.shape[1] < 2:
                    raise ValueError(
                        f"Spectrum at library index {lib_indices[idx]} has invalid peaks shape={peaks.shape}."
                    )

                if self.search_type == "neutral_loss":
                    peaks = self._transform_to_neutral_loss_space(
                        peaks,
                        spec_dict.get("precursor_mz"),
                        f"Library spectrum at index {lib_indices[idx]}",
                    )

                mz_s, ab_s = _sort_peaks_by_mz(peaks[:, 0], peaks[:, 1])
                sorted_peaks[idx] = (mz_s, ab_s)
            else:
                raise ValueError(
                    f"Spectrum index {idx} has no library index. "
                    "All spectra must be indexed in FlashEntropy library."
                )

        # Score pairs sequentially with pre-sorted peak arrays
        cosine_scores: dict[tuple[int, int], float] = {}
        tol = self.ms2_tolerance_da
        for i, j in pairs:
            mz_i, ab_i = sorted_peaks[i]
            mz_j, ab_j = sorted_peaks[j]
            cosine_scores[(i, j)] = _align_and_compute_cosine_presorted(
                mz_i, ab_i, mz_j, ab_j, tol
            )

        return cosine_scores

    # ── Public API ────────────────────────────────────────────────────────────

    def search_queries_against_library(
        self,
        query_spectra: list,
        query_ids: list[str],
        query_precursor_mzs: list[float | None] | None = None,
        *,
        format_library_id: Callable[[int], str] | None = None,
    ) -> tuple[dict[str, dict[tuple[str, str], float]], int]:
        """Search each query spectrum against ``self.fe_lib`` (query–library).

        Computes entropy similarity for all library hits with score > 0 and,
        when configured, cosine for pairs at or above
        :attr:`entropy_threshold_low`.

        Parameters
        ----------
        query_spectra : list
            Spectrum objects with ``.mz_exp`` and ``.abundance``.
        query_ids : list of str
            Unique ID for each query (same length as *query_spectra*).
        query_precursor_mzs : list of float or None, optional
            Precursor m/z per query.  Required for ``"identity"`` and
            ``"neutral_loss"``; optional for ``"open"``.
        format_library_id : callable, optional
            Maps library index → node ID string used in pair keys.
            Default ``str`` (bare index).  Callers that need non-colliding
            IDs (e.g. ``MolecularNetwork.library_node_id``) should pass a
            formatter.

        Returns
        -------
        scores : dict
            ``{metric_name: {(query_id, library_id): score}}`` with
            ``"entropy_similarity"`` and optionally ``"cosine"``.
        library_size : int
            Length of the FlashEntropy result vector (library size used for
            registration / indexing).  ``0`` if no successful searches ran.

        Raises
        ------
        RuntimeError
            If ``self.fe_lib`` is ``None``.
        ValueError
            If lengths of *query_ids* / *query_precursor_mzs* do not match
            *query_spectra*.
        """
        if self.fe_lib is None:
            raise RuntimeError(
                "search_queries_against_library requires a reference FE library "
                "(fe_lib)."
            )

        n_query = len(query_spectra)
        if len(query_ids) != n_query:
            raise ValueError(
                f"query_ids length ({len(query_ids)}) must match "
                f"query_spectra length ({n_query})."
            )
        if query_precursor_mzs is None:
            query_precursor_mzs = [None] * n_query
        elif len(query_precursor_mzs) != n_query:
            raise ValueError(
                f"query_precursor_mzs length ({len(query_precursor_mzs)}) must match "
                f"query_spectra length ({n_query})."
            )

        id_fmt = format_library_id if format_library_id is not None else str

        entropy_pairs: dict[tuple[str, str], float] = {}
        cosine_pairs: dict[tuple[str, str], float] = {}
        query_to_lib_indices: dict[int, list[int]] = {}
        library_size = 0

        for qi, (spec, pmz) in enumerate(zip(query_spectra, query_precursor_mzs)):
            peaks = self._peaks_array(spec)
            if peaks.shape[0] == 0:
                continue
            result_vec = self._clean_and_search(peaks, pmz)
            if result_vec is None:
                continue
            library_size = max(library_size, len(result_vec))

            lib_indices_for_query: list[int] = []
            for lib_idx, score in enumerate(result_vec):
                if score > 0.0:
                    entropy_pairs[
                        (query_ids[qi], id_fmt(lib_idx))
                    ] = float(score)
                    if score >= self.entropy_threshold_low:
                        lib_indices_for_query.append(lib_idx)

            if lib_indices_for_query:
                query_to_lib_indices[qi] = lib_indices_for_query

        if "cosine" in self.additional_similarities and query_to_lib_indices:
            for qi, lib_indices in query_to_lib_indices.items():
                cosine_scores = self._compute_cosine_for_query_vs_library(
                    query_spectrum=query_spectra[qi],
                    query_precursor_mz=query_precursor_mzs[qi],
                    library_indices=lib_indices,
                )
                for lib_idx, score in cosine_scores.items():
                    cosine_pairs[(query_ids[qi], id_fmt(lib_idx))] = score

        scores: dict[str, dict[tuple[str, str], float]] = {
            "entropy_similarity": entropy_pairs
        }
        if cosine_pairs:
            scores["cosine"] = cosine_pairs

        return scores, library_size

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
        # Extract build_index parameters from the main FE library if available
        # (these are stored as custom attributes by _build_flash_entropy_index)
        default_fe_kwargs = {
            "normalize_intensity": True,
            "min_ms2_difference_in_da": self.peak_sep_da,  # Already 2x tolerance from extraction
            "max_ms2_tolerance_in_da": self.ms2_tolerance_da,
            "max_indexed_mz": getattr(self.fe_lib, "_build_max_indexed_mz", 3000),
            "precursor_ions_removal_da": getattr(self.fe_lib, "_build_precursor_ions_removal_da", None),
            "noise_threshold": getattr(self.fe_lib, "_build_noise_threshold", 0),
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

        if self.fe_lib is None:
            raise RuntimeError(
                "compute_library_vs_library_filtered requires a reference FE library "
                "(fe_lib)."
            )

        if len(spectrum_ids) != n:
            raise ValueError(
                f"spectrum_ids length ({len(spectrum_ids)}) must match "
                f"library_indices length ({n})."
            )
        if precursor_mzs is not None and len(precursor_mzs) != n:
            raise ValueError(
                f"precursor_mzs length ({len(precursor_mzs)}) must match "
                f"library_indices length ({n})."
            )

        # Build lightweight spectrum objects from the raw library entries
        class _LibSpec:
            __slots__ = ("mz_exp", "abundance")

            def __init__(self, peaks_arr):
                self.mz_exp = peaks_arr[:, 0]
                self.abundance = peaks_arr[:, 1]

        lib_spectra: list = []
        lib_precursor_mzs: list[float | None] = []
        # IDs for successfully extracted entries only (same order as lib_spectra)
        valid_ids: list[str] = []

        for i, idx in enumerate(library_indices):
            try:
                entry = self.fe_lib[idx]
            except (IndexError, KeyError, TypeError) as exc:
                raise ValueError(
                    f"Cannot read library entry at index {idx}."
                ) from exc

            if not isinstance(entry, dict):
                raise TypeError(
                    f"Library entry at index {idx} must be a dict; "
                    f"got {type(entry).__name__}."
                )

            peaks = np.asarray(entry.get("peaks", []), dtype=float)
            # Skip empty/invalid peak tables but keep spectrum_ids/precursors
            # aligned to the original index *i* for survivors (not first-N trim).
            if peaks.ndim != 2 or peaks.shape[1] < 2 or peaks.shape[0] == 0:
                continue

            lib_spectra.append(_LibSpec(peaks))
            valid_ids.append(spectrum_ids[i])
            if precursor_mzs is not None:
                lib_precursor_mzs.append(precursor_mzs[i])
            else:
                lib_precursor_mzs.append(
                    float(entry.get("precursor_mz", 0.0) or 0.0)
                )

        if not lib_spectra:
            return {}

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

        Canonical all-vs-all path: each spectrum is searched once against the
        FlashEntropy index (*fe_lib_override* or ``self.fe_lib``), and
        upper-triangle pairwise scores are extracted from the result vectors.

        Parameters
        ----------
        spectra : list
            Spectrum objects with .mz_exp and .abundance attributes.
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

        # ── Stage 1: Compute entropy similarity matrix using FE search ────────
        entropy_pairs, pairs_for_additional = self._compute_entropy_matrix_and_pairs(
            spectra, spectrum_ids, precursor_mzs, fe_lib_override
        )

        result: dict[str, dict[tuple[str, str], float]] = {
            "entropy_similarity": entropy_pairs
        }

        # ── Stage 2: Compute additional metrics for pairs above threshold ────
        fe = fe_lib_override if fe_lib_override is not None else self.fe_lib
        for metric in self.additional_similarities:
            if metric == "cosine":
                cosine_idx_scores = self._compute_cosine_for_pairs(
                    pairs_for_additional,
                    lib_indices,
                    fe_lib_override=fe,
                )
                result["cosine"] = {
                    (spectrum_ids[i], spectrum_ids[j]): score
                    for (i, j), score in cosine_idx_scores.items()
                }

        return result
