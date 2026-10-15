"""Apex chosen by search_for_targeted_mass_features_batch."""
import numpy as np
import pandas as pd

from corems.mass_spectra.factory.lc_class import LCMSBase


class _Parent:
    """Stand-in LCMS parent. The search only uses it to build features."""


def _search(ms1df, mz_mins, mz_maxs, st_mins, st_maxs, st_aligned=False):
    set_ids = [f"w{i}" for i in range(len(mz_mins))]
    found = LCMSBase.search_for_targeted_mass_features_batch(
        _Parent(),
        ms1df,
        np.asarray(mz_mins, dtype=float),
        np.asarray(mz_maxs, dtype=float),
        np.asarray(st_mins, dtype=float),
        np.asarray(st_maxs, dtype=float),
        set_ids,
        st_aligned=st_aligned,
    )
    return [found[set_id] for set_id in set_ids]


def test_targeted_batch_search_picks_highest_intensity_and_earliest_tie():
    """Highest intensity wins. An intensity tie keeps the earliest row.

    A NaN m/z is not a candidate. A window with no row returns apex scan -99.
    """
    ms1df = pd.DataFrame(
        {
            "mz": [100.0, 100.0, 100.0, np.nan, 200.0],
            "scan_time": [1.0, 1.0, 3.0, 1.0, 1.0],
            "scan": [10, 11, 12, 14, 13],
            "intensity": [5.0, 9.0, 9.0, 100.0, 4.0],
        }
    )

    peaks = _search(
        ms1df,
        mz_mins=[100.0, 100.0, 150.0, 200.0],
        mz_maxs=[100.0, 100.0, 160.0, 200.0],
        st_mins=[0.0, 0.0, 0.0, 0.0],
        st_maxs=[2.0, 4.0, 2.0, 2.0],
    )

    # Rows 0 and 1 are in the window. Intensity 9 beats 5.
    assert peaks[0].apex_scan == 11
    assert peaks[0].mz == 100.0
    assert peaks[0].intensity == 9.0

    # Rows 1 and 2 both have intensity 9. Row 1 is earlier.
    assert peaks[1].apex_scan == 11

    assert peaks[2].apex_scan == -99
    assert np.isnan(peaks[2].mz)

    assert peaks[3].apex_scan == 13
    assert peaks[3].mz == 200.0


def test_targeted_batch_search_empty_ms1_is_a_miss():
    ms1df = pd.DataFrame(
        {"mz": [], "scan_time": [], "scan": [], "intensity": []}
    )
    peaks = _search(
        ms1df,
        mz_mins=[100.0],
        mz_maxs=[101.0],
        st_mins=[0.0],
        st_maxs=[1.0],
    )
    assert peaks[0].apex_scan == -99
    assert np.isnan(peaks[0].mz)


def test_targeted_batch_search_uses_aligned_time_column():
    ms1df = pd.DataFrame(
        {
            "mz": [100.0, 100.0],
            "scan_time": [1.0, 1.0],
            "scan_time_aligned": [0.2, 5.0],
            "scan": [10, 11],
            "intensity": [1.0, 50.0],
        }
    )
    peaks = _search(
        ms1df,
        mz_mins=[100.0],
        mz_maxs=[100.0],
        st_mins=[0.0],
        st_maxs=[1.0],
        st_aligned=True,
    )
    # The intense row is outside the aligned-time window.
    assert peaks[0].apex_scan == 10
    assert peaks[0].intensity == 1.0
