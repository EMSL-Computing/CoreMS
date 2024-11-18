__author__ = "Yuri E. Corilo"
__date__ = "Feb 13, 2020"

from pathlib import Path
import numpy as np
from pandas import Series

from corems.mass_spectra.calc import SignalProcessing as sp


class GC_Calculations:
    def calibrate_ri(self, ref_dict, cal_file_path):
        if not self:
            self.process_chromatogram()

        for gcms_peak in self:
            gcms_peak.calc_ri(ref_dict)

        self.ri_pairs_ref = ref_dict
        if isinstance(cal_file_path, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            self.cal_file_path = Path(cal_file_path)
        else:
            self.cal_file_path = cal_file_path

    def smooth_tic(self, tic):
        implemented_smooth_method = self.chromatogram_settings.implemented_smooth_method

        pol_order = self.chromatogram_settings.savgol_pol_order

        window_len = self.chromatogram_settings.smooth_window

        window = self.chromatogram_settings.smooth_method

        return sp.smooth_signal(
            tic, window_len, window, pol_order, implemented_smooth_method
        )

    def centroid_detector(self, tic, rt):
        noise_std = self.chromatogram_settings.std_noise_threshold

        method = self.chromatogram_settings.noise_threshold_method

        # peak picking
        min_height = self.chromatogram_settings.peak_height_min_percent
        min_datapoints = self.chromatogram_settings.min_peak_datapoints

        # baseline detection
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent
        max_height = self.chromatogram_settings.peak_height_max_percent

        peak_indexes_generator = sp.peak_detector_generator(
            tic,
            noise_std,
            method,
            rt,
            max_height,
            min_height,
            max_prominence,
            min_datapoints,
        )

        return peak_indexes_generator

    def remove_outliers(self, data):
        from numpy import percentile

        q25, q75 = percentile(data, 25), percentile(data, 75)
        iqr = q75 - q25
        if self.parameters.verbose_processing:
            print("Percentiles: 25th=%.3f, 75th=%.3f, IQR=%.3f" % (q25, q75, iqr))
        # calculate the outlier cutoff
        cut_off = iqr * 1.5
        lower, upper = q25 - cut_off, q75 + cut_off
        # identify outliers
        outliers = [x for x in data if x < lower or x > upper]
        if self.parameters.verbose_processing:
            print("Identified outliers: %d" % len(outliers))
        # remove outliers
        nanfilled_outliers = Series(
            [x if lower <= x <= upper else np.nan for x in data]
        )

        return nanfilled_outliers
