import numpy as np
from bisect import bisect_left
from scipy.optimize import curve_fit


__author__ = "Yuri E. Corilo"
__date__ = "March 11, 2020"


class GCPeakCalculation(object):
    """
    Class for performing peak calculations in GC chromatography.

    Methods
    -------
    * `calc_area(self, tic: List[float], dx: float) -> None`: Calculate the area under the curve of the chromatogram.
    * `linear_ri(self, right_ri: float, left_ri: float, left_rt: float, right_rt: float) -> float`: Calculate the retention index using linear interpolation.
    * `calc_ri(self, rt_ri_pairs: List[Tuple[float, float]]) -> int`: Calculate the retention index based on the given retention time - retention index pairs.
    """

    def calc_area(self, tic: list[float], dx: float) -> None:
        """
        Calculate the area under the curve of the chromatogram.

        Parameters
        ----------
        tic : List[float]
            The total ion current (TIC) values.
        dx : float
            The spacing between data points.
        """
        yy = tic[self.start_scan : self.final_scan]
        self._area = np.trapz(yy, dx=dx)

    def linear_ri(
        self, right_ri: float, left_ri: float, left_rt: float, right_rt: float
    ) -> float:
        """
        Calculate the retention index using linear interpolation.

        Parameters
        ----------
        right_ri : float
            The retention index at the right reference point.
        left_ri : float
            The retention index at the left reference point.
        left_rt : float
            The retention time at the left reference point.
        right_rt : float
            The retention time at the right reference point.

        Returns
        -------
        float
            The calculated retention index.
        """
        return left_ri + (
            (right_ri - left_ri)
            * (self.retention_time - left_rt)
            / (right_rt - left_rt)
        )

    def calc_ri(self, rt_ri_pairs: list[tuple[float, float]]) -> None:
        """
        Calculate the retention index based on the given retention time - retention index pairs.

        Parameters
        ----------
        rt_ri_pairs : List[Tuple[float, float]]
            The list of retention time - retention index pairs.

        """
        current_rt = self.retention_time

        rts = [rt_ri[0] for rt_ri in rt_ri_pairs]
        index = bisect_left(rts, current_rt)

        if index >= len(rt_ri_pairs):
            index -= 1

        current_ref = rt_ri_pairs[index]

        if current_rt == current_ref[0]:
            self._ri = current_ref[1]

        else:
            if index == 0:
                index += 1

            left_rt = rt_ri_pairs[index - 1][0]
            left_ri = rt_ri_pairs[index - 1][1]

            right_rt = rt_ri_pairs[index][0]
            right_ri = rt_ri_pairs[index][1]

            self._ri = self.linear_ri(right_ri, left_ri, left_rt, right_rt)


class LCMSMassFeatureCalculation:
    """Class for performing peak calculations in LC-MS mass spectrometry.

    This class is intended to be used as a mixin class for the LCMSMassFeature class.
    """

    def calc_dispersity_index(self):
        """
        Calculate the dispersity index of the mass feature.

        This function calculates the dispersity index of the mass feature and
        stores the result in the `_dispersity_index` attribute. The dispersity index is calculated as the standard
        deviation of the retention times that account for 50% of the cummulative intensity, starting from the most
        intense point, as described in [1]. Note that this calculation is done on the full EIC data, not just the data within a mass feature's
        integration bounds.

        Returns
        -------
        None, stores the result in the `_dispersity_index` attribute of the class.

        Raises
        ------
        ValueError
            If the EIC data are not available.

        References
        ----------
        1) Boiteau, Rene M., et al. "Relating Molecular Properties to the Persistence of Marine Dissolved
        Organic Matter with Liquid Chromatography–Ultrahigh-Resolution Mass Spectrometry."
        Environmental Science & Technology 58.7 (2024): 3267-3277.
        """
        # Check if the EIC data is available
        if self.eic_list is None:
            raise ValueError(
                "EIC data are not available. Please add the EIC data first."
            )

        # Sort the EIC data and RT data by descending intensity
        eic_in = self._eic_data.eic.argsort()
        sorted_eic = self._eic_data.eic[eic_in[::-1]]
        sorted_rt = self._eic_data.time[eic_in[::-1]]

        # Calculate the dispersity index
        cum_sum = np.cumsum(sorted_eic) / np.sum(sorted_eic)
        rt_summ = sorted_rt[np.where(cum_sum < 0.5)]
        if len(rt_summ) > 1:
            d = np.std(rt_summ)
            self._dispersity_index = d
        elif len(rt_summ) == 1:
            self._dispersity_index = 0

    def calc_fraction_height_width(self, fraction: float):
        """
        Calculate the height width of the mass feature at a specfic fraction of the maximum intensity.

        This function returns a tuple with the minimum and maximum half-height width based on scan resolution.

        Parameters
        ----------
        fraction : float
            The fraction of the maximum intensity to calculate the height width.
            For example, 0.5 will calculate the half-height width.

        Returns
        -------
        Tuple[float, float, bool]
            The minimum and maximum half-height width based on scan resolution (in minutes), and a boolean indicating if the width was estimated.
        """

        # Pull out the EIC data
        eic = self._eic_data.eic_smoothed

        # Find the indices of the maximum intensity on either side
        max_index = np.where(self._eic_data.scans == self.apex_scan)[0][0]
        left_index = max_index
        right_index = max_index
        while eic[left_index] > eic[max_index] * fraction and left_index > 0:
            left_index -= 1
        while (
            eic[right_index] > eic[max_index] * fraction and right_index < len(eic) - 1
        ):
            right_index += 1

        # Get the retention times of the indexes just below the half height
        left_rt = self._eic_data.time[left_index]
        right_rt = self._eic_data.time[right_index]

        # If left_rt and right_rt are outside the bounds of the integration, set them to the bounds and set estimated to True
        estimated = False
        if left_rt < self.eic_rt_list[0]:
            left_rt = self.eic_rt_list[0]
            left_index = np.where(self._eic_data.scans == self._eic_data.apexes[0][0])[
                0
            ][0]
            estimated = True
        if right_rt > self.eic_rt_list[-1]:
            right_rt = self.eic_rt_list[-1]
            right_index = np.where(
                self._eic_data.scans == self._eic_data.apexes[0][-1]
            )[0][0]
            estimated = True
        half_height_width_max = right_rt - left_rt

        # Get the retention times of the indexes just above the half height
        left_rt = self._eic_data.time[left_index + 1]
        right_rt = self._eic_data.time[right_index - 1]
        half_height_width_min = right_rt - left_rt

        return half_height_width_min, half_height_width_max, estimated

    def calc_half_height_width(self, accept_estimated: bool = False):
        """
        Calculate the half-height width of the mass feature.

        This function calculates the half-height width of the mass feature and
        stores the result in the `_half_height_width` attribute

        Returns
        -------
        None, stores the result in the `_half_height_width` attribute of the class.
        """
        min_, max_, estimated = self.calc_fraction_height_width(0.5)
        if not estimated or accept_estimated:
            self._half_height_width = np.array([min_, max_])

    def calc_tailing_factor(self, accept_estimated: bool = False):
        """
        Calculate the peak asymmetry of the mass feature.

        This function calculates the peak asymmetry of the mass feature and
        stores the result in the `_tailing_factor` attribute.
        Calculations completed at 5% of the peak height in accordance with the USP tailing factor calculation.

        Returns
        -------
        None, stores the result in the `_tailing_factor` attribute of the class.

        References
        ----------
        1) JIS K0124:2011 General rules for high performance liquid chromatography
        2) JIS K0214:2013 Technical terms for analytical chemistry
        """
        # First calculate the width of the peak at 5% of the peak height
        width_min, width_max, estimated = self.calc_fraction_height_width(0.05)

        if not estimated or accept_estimated:
            # Next calculate the width of the peak at 95% of the peak height
            eic = self._eic_data.eic_smoothed
            max_index = np.where(self._eic_data.scans == self.apex_scan)[0][0]
            left_index = max_index
            while eic[left_index] > eic[max_index] * 0.05 and left_index > 0:
                left_index -= 1

            left_half_time_min = (
                self._eic_data.time[max_index] - self._eic_data.time[left_index]
            )
            left_half_time_max = (
                self._eic_data.time[max_index] - self._eic_data.time[left_index + 1]
            )

            tailing_factor = np.mean([width_min, width_max]) / (
                2 * np.mean([left_half_time_min, left_half_time_max])
            )

            self._tailing_factor = tailing_factor

    def calc_gaussian_similarity(self):
        """
        Calculate the Gaussian similarity score of the mass feature.
        
        This function fits a Gaussian curve to the EIC data and calculates 
        the R-squared correlation between the actual EIC and the fitted Gaussian.
        A higher score indicates better Gaussian shape (ideal chromatographic peak).
        
        Returns
        -------
        None, stores the result in the `_gaussian_similarity` attribute of the class.
        
        Raises
        ------
        ValueError
            If the EIC data are not available.
        """
        # Check if the EIC data is available
        if self.eic_list is None:
            raise ValueError(
                "EIC data are not available. Please add the EIC data first."
            )
        
        # Get EIC data within integration bounds
        time_data = np.array(self.eic_rt_list)
        intensity_data = np.array(self.eic_list)
        
        if len(time_data) < 4:  # Need minimum points for meaningful fit
            self._gaussian_similarity = 0.0
            return
        
        # Normalize intensities for fitting
        max_intensity = np.max(intensity_data)
        if max_intensity == 0:
            self._gaussian_similarity = 0.0
            return
        
        try:
            # Define Gaussian function
            def gaussian(x, amplitude, mean, stddev, baseline):
                return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2)) + baseline
            
            # Initial parameter estimates
            amplitude_init = max_intensity
            mean_init = time_data[np.argmax(intensity_data)]
            stddev_init = (time_data[-1] - time_data[0]) / 6  # Rough estimate
            baseline_init = np.min(intensity_data)
            
            # Fit Gaussian curve
            popt, _ = curve_fit(
                gaussian, 
                time_data, 
                intensity_data,
                p0=[amplitude_init, mean_init, stddev_init, baseline_init],
                maxfev=1000,
                bounds=(
                    [0, time_data[0], 0, 0],  # Lower bounds
                    [np.inf, time_data[-1], np.inf, max_intensity]  # Upper bounds
                )
            )
            
            # Calculate fitted values
            fitted_intensities = gaussian(time_data, *popt)
            
            # Calculate R-squared (coefficient of determination)
            ss_res = np.sum((intensity_data - fitted_intensities) ** 2)
            ss_tot = np.sum((intensity_data - np.mean(intensity_data)) ** 2)
            
            if ss_tot == 0:
                self._gaussian_similarity = 1.0 if ss_res == 0 else 0.0
            else:
                r_squared = 1 - (ss_res / ss_tot)
                # Ensure score is between 0 and 1
                self._gaussian_similarity = max(0.0, min(1.0, r_squared))
                
        except (RuntimeError, ValueError, TypeError):
            # Fitting failed, assign poor similarity
            self._gaussian_similarity = 0.0

    def calc_noise_score(self, noise_window_factor: float = 2.0):
        """
        Calculate the noise score of the mass feature.
        
        This function estimates the signal-to-noise ratio by comparing the peak
        intensity to the baseline noise level in surrounding regions.
        
        Parameters
        ----------
        noise_window_factor : float, default 2.0
            Factor to determine noise estimation window size relative to peak width.
            Larger values use wider windows for noise estimation.
        
        Returns
        -------
        None, stores the result in the `_noise_score` attribute of the class.
        
        Raises
        ------
        ValueError
            If the EIC data are not available.
        """
        # Check if the EIC data is available
        if self.eic_list is None:
            raise ValueError(
                "EIC data are not available. Please add the EIC data first."
            )
        
        # Get full EIC data (not just integration bounds)
        full_time = self._eic_data.time
        full_eic = self._eic_data.eic_smoothed
        
        # Get peak information
        apex_rt = self.retention_time
        peak_intensity = np.max(self.eic_list)
        
        # Estimate peak width (use half-height width if available, otherwise estimate)
        if hasattr(self, '_half_height_width') and self._half_height_width is not None:
            peak_width = np.mean(self._half_height_width)
        else:
            # Estimate width from integration bounds
            peak_width = self.eic_rt_list[-1] - self.eic_rt_list[0]

        # Define noise estimation windows
        noise_window_size = peak_width * noise_window_factor
        left_noise_start = apex_rt - peak_width - noise_window_size
        left_noise_end = apex_rt - peak_width
        right_noise_start = apex_rt + peak_width
        right_noise_end = apex_rt + peak_width + noise_window_size
        
        # Extract noise regions
        left_noise_mask = (full_time >= left_noise_start) & (full_time <= left_noise_end)
        right_noise_mask = (full_time >= right_noise_start) & (full_time <= right_noise_end)
        
        left_noise = full_eic[left_noise_mask]
        right_noise = full_eic[right_noise_mask]
        
        # Combine noise regions
        noise_data = np.concatenate([left_noise, right_noise])
        
        if len(noise_data) == 0:
            # No noise data available, use minimum intensity as baseline
            baseline_intensity = np.min(full_eic)
            noise_std = baseline_intensity * 0.1  # Rough estimate
        else:
            # Calculate baseline and noise statistics
            baseline_intensity = np.median(noise_data)  # Use median for robustness
            noise_std = np.std(noise_data)
            
            # Handle case where noise_std is very small
            if noise_std == 0:
                noise_std = baseline_intensity * 0.01 if baseline_intensity > 0 else 1
        
        # Calculate signal-to-noise ratio
        signal = peak_intensity - baseline_intensity
        if signal <= 0 or noise_std == 0:
            self._noise_score = 0.0
        else:
            snr = signal / noise_std
            # Convert to a 0-1 score using sigmoid-like transformation
            # Score approaches 1 for high S/N ratios, 0 for low ones
            self._noise_score = min(1.0, snr / (snr + 10.0))
