"""Equivalence and smoke tests for peak-aligned cosine similarity.

The reference implementation freezes the original sequential alignment
(per-peak ``find_closest`` + list appends).  The production utility must
match it **exactly** on every case (bit-identical float results).
"""

from __future__ import annotations

import numpy as np
import pytest

from corems.mass_spectra.calc.lc_calc import find_closest
from corems.molecular_networking.similarity_engine import (
    _align_and_compute_cosine,
    _align_and_compute_cosine_presorted,
    _sort_peaks_by_mz,
)


def _cosine_reference_sequential(mz1, abun1, mz2, abun2, tolerance_da):
    """Original sequential algorithm (do not 'improve' — golden reference)."""
    mz1 = np.asarray(mz1, dtype=float)
    abun1 = np.asarray(abun1, dtype=float)
    mz2 = np.asarray(mz2, dtype=float)
    abun2 = np.asarray(abun2, dtype=float)

    if len(mz1) == 0 or len(mz2) == 0:
        return 0.0

    idx1 = np.argsort(mz1)
    mz1_sorted = mz1[idx1]
    abun1_sorted = abun1[idx1]

    idx2 = np.argsort(mz2)
    mz2_sorted = mz2[idx2]
    abun2_sorted = abun2[idx2]

    vec1 = []
    vec2 = []
    used_spec2 = np.zeros(len(mz2_sorted), dtype=bool)

    for i in range(len(mz1_sorted)):
        closest_idx = find_closest(mz2_sorted, np.array([mz1_sorted[i]]))[0]
        diff = abs(mz2_sorted[closest_idx] - mz1_sorted[i])

        if diff <= tolerance_da:
            vec1.append(abun1_sorted[i])
            vec2.append(abun2_sorted[closest_idx])
            used_spec2[closest_idx] = True
        else:
            vec1.append(abun1_sorted[i])
            vec2.append(0.0)

    for j in range(len(mz2_sorted)):
        if not used_spec2[j]:
            vec1.append(0.0)
            vec2.append(abun2_sorted[j])

    vec1 = np.array(vec1, dtype=float)
    vec2 = np.array(vec2, dtype=float)

    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    if norm1 == 0 or norm2 == 0:
        return 0.0

    cosine = np.dot(vec1, vec2) / (norm1 * norm2)
    return float(np.clip(cosine, 0.0, 1.0))


def _assert_exact(a, b):
    assert a == b, f"not bit-identical: {a!r} vs {b!r} (delta={a - b!r})"


@pytest.mark.parametrize(
    "name,mz1,ab1,mz2,ab2,tol",
    [
        (
            "identical",
            np.array([100.0, 150.0, 200.0, 250.0]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            np.array([100.0, 150.0, 200.0, 250.0]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            0.01,
        ),
        (
            "shift_in",
            np.array([100.0, 150.0, 200.0, 250.0]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            np.array([100.005, 150.005, 200.005, 250.005]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            0.01,
        ),
        (
            "shift_out",
            np.array([100.0, 150.0, 200.0, 250.0]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            np.array([100.02, 150.02, 200.02, 250.02]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            0.01,
        ),
        (
            "diff",
            np.array([100.0, 150.0, 200.0, 250.0]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            np.array([110.0, 180.0, 220.0]),
            np.array([1.0, 0.5, 0.3]),
            0.01,
        ),
        (
            "empty",
            np.array([100.0, 150.0]),
            np.array([1.0, 0.5]),
            np.array([]),
            np.array([]),
            0.01,
        ),
        (
            "multi_closest",
            np.array([100.0, 150.0]),
            np.array([1.0, 0.8]),
            np.array([100.0, 100.001, 150.0]),
            np.array([0.5, 0.7, 1.0]),
            0.01,
        ),
        (
            "unsorted",
            np.array([250.0, 200.0, 150.0, 100.0]),
            np.array([0.2, 0.5, 0.8, 1.0]),
            np.array([100.0, 150.0, 200.0, 250.0]),
            np.array([1.0, 0.8, 0.5, 0.2]),
            0.01,
        ),
        (
            "same_peak2_twice",
            # Two spectrum-1 peaks closer to the same spectrum-2 peak than
            # tolerance allows for any other peak → reuse of that peak2 index.
            np.array([100.0, 100.002]),
            np.array([1.0, 0.5]),
            np.array([100.001, 200.0]),
            np.array([0.9, 0.1]),
            0.01,
        ),
    ],
)
def test_cosine_matches_sequential_reference(name, mz1, ab1, mz2, ab2, tol):
    ref = _cosine_reference_sequential(mz1, ab1, mz2, ab2, tol)
    fast = _align_and_compute_cosine(mz1, ab1, mz2, ab2, tol)
    _assert_exact(fast, ref)


def test_cosine_random_spectra_exact_match():
    rng = np.random.default_rng(42)
    for n1, n2, tol in [
        (5, 5, 0.01),
        (20, 15, 0.02),
        (50, 80, 0.005),
        (100, 100, 0.01),
        (1, 30, 0.01),
        (40, 1, 0.01),
    ]:
        for _ in range(25):
            mz1 = np.sort(rng.uniform(50, 500, size=n1))
            mz2 = np.sort(rng.uniform(50, 500, size=n2))
            ab1 = rng.random(n1)
            ab2 = rng.random(n2)
            # occasionally shuffle inputs
            if rng.random() < 0.5:
                p = rng.permutation(n1)
                mz1, ab1 = mz1[p], ab1[p]
            ref = _cosine_reference_sequential(mz1, ab1, mz2, ab2, tol)
            fast = _align_and_compute_cosine(mz1, ab1, mz2, ab2, tol)
            _assert_exact(fast, ref)


def test_cosine_zero_abundance_spectra():
    mz = np.array([100.0, 200.0])
    ab0 = np.array([0.0, 0.0])
    ab1 = np.array([1.0, 0.5])
    ref = _cosine_reference_sequential(mz, ab0, mz, ab1, 0.01)
    fast = _align_and_compute_cosine(mz, ab0, mz, ab1, 0.01)
    _assert_exact(fast, ref)
    assert fast == 0.0


def test_presorted_matches_sorting_wrapper():
    """Pre-sort once then presorted path must match full _align_and_compute_cosine."""
    rng = np.random.default_rng(7)
    for _ in range(40):
        n1, n2 = int(rng.integers(1, 60)), int(rng.integers(1, 60))
        mz1 = rng.uniform(50, 500, size=n1)
        mz2 = rng.uniform(50, 500, size=n2)
        ab1 = rng.random(n1)
        ab2 = rng.random(n2)
        tol = 0.01
        full = _align_and_compute_cosine(mz1, ab1, mz2, ab2, tol)
        m1, a1 = _sort_peaks_by_mz(mz1, ab1)
        m2, a2 = _sort_peaks_by_mz(mz2, ab2)
        pre = _align_and_compute_cosine_presorted(m1, a1, m2, a2, tol)
        _assert_exact(pre, full)
        # Reference sequential as well
        _assert_exact(full, _cosine_reference_sequential(mz1, ab1, mz2, ab2, tol))
