"""Absolute intensity floor before LC-MS persistent-homology peak picking."""

import pandas as pd

from corems.encapsulation.factory.parameters import LCMSParameters
from corems.encapsulation.factory.processingSetting import (
    LiquidChromatographSetting,
    settings_to_export_dict,
)
from corems.mass_spectra.calc.lc_calc import PHCalculations

# Max intensity is 1000. A relative floor of 0.1 keeps intensity > 100.
_INTENSITIES = (50.0, 200.0, 1000.0)


def _frame():
    data = pd.DataFrame(
        {
            "mz": [200.0, 201.0, 202.0],
            "scan": [1, 1, 1],
            "intensity": list(_INTENSITIES),
        }
    )
    scan_df = pd.DataFrame({"scan": [1], "scan_time": [1.0]})
    return data, scan_df


def _target_window():
    return {
        "target_mz_list": [201.0],
        "target_rt_list": [1.0],
        "mz_tolerance_ppm": 20_000,
        "rt_tolerance": 0.5,
    }


class _ProfileHost(PHCalculations):
    def __init__(self, abs_floor):
        data, scan_df = _frame()
        self.parameters = LCMSParameters()
        self.parameters.lc_ms.ph_inten_min_rel = 0.1
        self.parameters.lc_ms.ph_inten_min_abs = abs_floor
        self.parameters.lc_ms.verbose_processing = False
        self._ms_unprocessed = {1: data}
        self.scan_df = scan_df
        self.mass_features = {}
        self.kept = None

    def check_if_grid(self, data):
        return True

    def _record(self, data_thres):
        self.kept = sorted(float(value) for value in data_thres["intensity"])

    def _find_mass_features_ph_single(self, data_thres, *args, **kwargs):
        self._record(data_thres)

    def _find_mass_features_ph_partition(self, data_thres, *args, **kwargs):
        self._record(data_thres)


class _CentroidHost(PHCalculations):
    def __init__(self, abs_floor):
        data, scan_df = _frame()
        self.parameters = LCMSParameters()
        self.parameters.lc_ms.ph_inten_min_rel = 0.1
        self.parameters.lc_ms.ph_inten_min_abs = abs_floor
        self.parameters.lc_ms.verbose_processing = False
        self._ms_unprocessed = {1: data}
        self.scan_df = scan_df
        self.mass_features = {}
        self.kept = None

    def roll_up_dataframe(self, df, **kwargs):
        if len(df):
            self.kept = sorted(float(value) for value in df["intensity"])
        return df.iloc[0:0].copy()


def test_ph_inten_min_abs_defaults_to_zero_and_is_exported():
    settings = LiquidChromatographSetting()
    assert settings.ph_inten_min_abs == 0.0
    assert settings_to_export_dict(settings)["ph_inten_min_abs"] == 0.0


def test_zero_absolute_floor_keeps_the_relative_mask():
    profile = _ProfileHost(0.0)
    profile.find_mass_features_ph()
    centroid = _CentroidHost(0.0)
    centroid.find_mass_features_ph_centroid()
    assert profile.kept == [200.0, 1000.0]
    assert centroid.kept == [200.0, 1000.0]


def test_absolute_floor_drops_points_the_relative_floor_would_keep():
    profile = _ProfileHost(500.0)
    profile.find_mass_features_ph()
    centroid = _CentroidHost(500.0)
    centroid.find_mass_features_ph_centroid()
    assert profile.kept == [1000.0]
    assert centroid.kept == [1000.0]


def test_relative_floor_still_drops_points_above_the_absolute_floor():
    profile = _ProfileHost(10.0)
    profile.find_mass_features_ph()
    centroid = _CentroidHost(10.0)
    centroid.find_mass_features_ph_centroid()
    assert profile.kept == [200.0, 1000.0]
    assert centroid.kept == [200.0, 1000.0]


def test_targeted_search_ignores_the_absolute_floor():
    profile = _ProfileHost(500.0)
    profile.find_mass_features_ph(
        targeted_search=True, target_search_dict=_target_window()
    )
    centroid = _CentroidHost(500.0)
    centroid.find_mass_features_ph_centroid(
        targeted_search=True, target_search_dict=_target_window()
    )
    assert profile.kept == [50.0, 200.0, 1000.0]
    assert centroid.kept == [50.0, 200.0, 1000.0]
