from types import SimpleNamespace

import numpy as np
import pytest

from corems.mass_spectrum.calc.NoiseCalc import NoiseThresholdCalc


class DummyNoiseCalc(NoiseThresholdCalc):
    pass


def _build_calc(mz_profile, abundance_profile, min_mz, max_mz):
    calc = DummyNoiseCalc()
    calc.is_centroid = False
    calc.settings = SimpleNamespace(
        noise_threshold_method="signal_noise",
        noise_min_mz=min_mz,
        noise_max_mz=max_mz,
    )
    calc.mz_exp_profile = np.array(mz_profile)
    calc.abundance_profile = np.array(abundance_profile)
    return calc


def test_cut_mz_domain_noise_inclusive_ascending():
    calc = _build_calc(
        mz_profile=[1.0, 2.0, 3.0, 4.0],
        abundance_profile=[10.0, 20.0, 30.0, 40.0],
        min_mz=2.0,
        max_mz=3.0,
    )

    mz_cut, abundance_cut = calc.cut_mz_domain_noise()

    assert np.array_equal(mz_cut, np.array([2.0, 3.0]))
    assert np.array_equal(abundance_cut, np.array([20.0, 30.0]))


def test_cut_mz_domain_noise_inclusive_descending():
    calc = _build_calc(
        mz_profile=[4.0, 3.0, 2.0, 1.0],
        abundance_profile=[40.0, 30.0, 20.0, 10.0],
        min_mz=2.0,
        max_mz=3.0,
    )

    mz_cut, abundance_cut = calc.cut_mz_domain_noise()

    assert np.array_equal(mz_cut, np.array([3.0, 2.0]))
    assert np.array_equal(abundance_cut, np.array([30.0, 20.0]))


def test_cut_mz_domain_noise_warns_and_returns_empty_for_invalid_roi():
    calc = _build_calc(
        mz_profile=[1.0, 2.0, 3.0, 4.0],
        abundance_profile=[10.0, 20.0, 30.0, 40.0],
        min_mz=10.0,
        max_mz=11.0,
    )

    with pytest.warns(UserWarning, match="Empty noise ROI"):
        mz_cut, abundance_cut = calc.cut_mz_domain_noise()

    assert mz_cut.size == 0
    assert abundance_cut.size == 0


def test_run_noise_threshold_calc_signal_noise_changes_with_roi():
    mz_profile = [1.0, 2.0, 3.0, 4.0, 5.0]
    abundance_profile = [1.0, 2.0, 100.0, 200.0, 300.0]

    calc_full_roi = _build_calc(
        mz_profile=mz_profile,
        abundance_profile=abundance_profile,
        min_mz=1.0,
        max_mz=5.0,
    )
    calc_narrow_roi = _build_calc(
        mz_profile=mz_profile,
        abundance_profile=abundance_profile,
        min_mz=1.0,
        max_mz=2.0,
    )

    avg_full, std_full = calc_full_roi.run_noise_threshold_calc()
    avg_narrow, std_narrow = calc_narrow_roi.run_noise_threshold_calc()

    # For signal_noise, run_noise_threshold_calc uses get_noise_average(abundance_cut),
    # so changing ROI should change median/std when ROI content changes.
    assert avg_full != avg_narrow
    assert std_full != std_narrow
