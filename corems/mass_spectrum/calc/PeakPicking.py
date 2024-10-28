"""
@author: Yuri E. Corilo
@date: Jun 27, 2019
"""

import warnings
from numpy import (
    hstack,
    inf,
    isnan,
    where,
    array,
    polyfit,
    nan,
    pad,
    zeros,
    searchsorted,
)
from corems.encapsulation.constant import Labels
from corems.mass_spectra.calc import SignalProcessing as sp


class PeakPicking:
    """Class for peak picking.

    Parameters
    ----------
    None

    Attributes
    ----------
    None

    Methods
    -------
    * prepare_peak_picking_data().
        Prepare the mz, abundance, and frequence data for peak picking.
    * cut_mz_domain_peak_picking().
        Cut the m/z domain for peak picking.
    * extrapolate_axes_for_pp(mz=None, abund=None, freq=None).
        Extrapolate the m/z axis and fill the abundance axis with 0s.
    * do_peak_picking().
        Perform peak picking.
    * find_minima(apex_index, abundance, len_abundance, right=True).
        Find the minima of a peak.
    * linear_fit_calc(intes, massa, index_term, index_sign).
        Algebraic solution to a linear fit.
    * calculate_resolving_power(intes, massa, current_index).
        Calculate the resolving power of a peak.
    * cal_minima(mass, abun).
        Calculate the minima of a peak.
    * calc_centroid(mass, abund, freq).
        Calculate the centroid of a peak.
    * get_threshold(intes).
        Get the intensity threshold for peak picking.
    * algebraic_quadratic(list_mass, list_y).
        Find the apex of a peak - algebraically.
    * find_apex_fit_quadratic(mass, abund, freq, current_index).
        Find the apex of a peak.
    * check_prominence(abun, current_index, len_abundance, peak_height_diff).
        Check the prominence of a peak.
    * use_the_max(mass, abund, current_index, len_abundance, peak_height_diff).
        Use the max peak height as the centroid.
    * calc_centroid_legacy(mass, abund, freq).
        Legacy centroid calculation. Deprecated - for deletion.

    """

    def prepare_peak_picking_data(self):
        """Prepare the data for peak picking.

        This function will prepare the m/z, abundance, and frequency data for peak picking according to the settings.

        Returns
        -------
        mz : ndarray
            The m/z axis.
        abundance : ndarray
            The abundance axis.
        freq : ndarray or None
            The frequency axis, if available.
        """
        # First apply cut_mz_domain_peak_picking
        mz, abundance, freq = self.cut_mz_domain_peak_picking()

        # Then extrapolate the axes for peak picking
        if self.settings.picking_point_extrapolate > 0:
            mz, abundance, freq = self.extrapolate_axes_for_pp(mz, abundance, freq)
        return mz, abundance, freq

    def cut_mz_domain_peak_picking(self):
        """
        Cut the m/z domain for peak picking.

        Simplified function

        Returns
        -------
        mz_domain_X_low_cutoff : ndarray
            The m/z values within the specified range.
        mz_domain_low_Y_cutoff : ndarray
            The abundance values within the specified range.
        freq_domain_low_Y_cutoff : ndarray or None
            The frequency values within the specified range, if available.

        """
        max_picking_mz = self.settings.max_picking_mz
        min_picking_mz = self.settings.min_picking_mz

        # min_start =  where(self.mz_exp_profile  > min_picking_mz)[0][0]
        # max_final =  where(self.mz_exp_profile < max_picking_mz)[-1][-1]
        min_start = searchsorted(a=self.mz_exp_profile, v=min_picking_mz)
        max_final = searchsorted(a=self.mz_exp_profile, v=max_picking_mz)

        if self.has_frequency:
            if self.freq_exp_profile.any():
                return (
                    self.mz_exp_profile[min_start:max_final],
                    self.abundance_profile[min_start:max_final],
                    self.freq_exp_profile[min_start:max_final],
                )

        else:
            return (
                self.mz_exp_profile[min_start:max_final],
                self.abundance_profile[min_start:max_final],
                None,
            )

    def legacy_cut_mz_domain_peak_picking(self):
        """
        Cut the m/z domain for peak picking.
        DEPRECATED
        Returns
        -------
        mz_domain_X_low_cutoff : ndarray
            The m/z values within the specified range.
        mz_domain_low_Y_cutoff : ndarray
            The abundance values within the specified range.
        freq_domain_low_Y_cutoff : ndarray or None
            The frequency values within the specified range, if available.

        """
        max_picking_mz = self.settings.max_picking_mz
        min_picking_mz = self.settings.min_picking_mz

        min_final = where(self.mz_exp_profile > min_picking_mz)[-1][-1]
        min_start = where(self.mz_exp_profile > min_picking_mz)[0][0]

        (
            mz_domain_X_low_cutoff,
            mz_domain_low_Y_cutoff,
        ) = (
            self.mz_exp_profile[min_start:min_final],
            self.abundance_profile[min_start:min_final],
        )

        max_final = where(self.mz_exp_profile < max_picking_mz)[-1][-1]
        max_start = where(self.mz_exp_profile < max_picking_mz)[0][0]

        if self.has_frequency:
            if self.freq_exp_profile.any():
                freq_domain_low_Y_cutoff = self.freq_exp_profile[min_start:min_final]

                return (
                    mz_domain_X_low_cutoff[max_start:max_final],
                    mz_domain_low_Y_cutoff[max_start:max_final],
                    freq_domain_low_Y_cutoff[max_start:max_final],
                )

        else:
            return (
                mz_domain_X_low_cutoff[max_start:max_final],
                mz_domain_low_Y_cutoff[max_start:max_final],
                None,
            )

    @staticmethod
    def extrapolate_axis(initial_array, pts):
        """
        This function will extrapolate an input array in both directions by N pts.

        Parameters
        ----------
        initial_array : ndarray
            The input array.
        pts : int
            The number of points to extrapolate.

        Returns
        -------
        ndarray
            The extrapolated array.

        Notes
        --------
        This is a static method.
        """
        initial_array_len = len(initial_array)
        right_delta = initial_array[-1] - initial_array[-2]
        left_delta = initial_array[1] - initial_array[0]

        # Create an array with extra space for extrapolation
        pad_array = zeros(initial_array_len + 2 * pts)

        # Copy original array into the middle of the padded array
        pad_array[pts : pts + initial_array_len] = initial_array

        # Extrapolate the right side
        for pt in range(pts):
            final_value = initial_array[-1]
            value_to_add = right_delta * (pt + 1)
            new_value = final_value + value_to_add
            pad_array[initial_array_len + pts + pt] = new_value

        # Extrapolate the left side
        for pt in range(pts):
            first_value = initial_array[0]
            value_to_subtract = left_delta * (pt + 1)
            new_value = first_value - value_to_subtract
            pad_array[pts - pt - 1] = new_value

        return pad_array

    def extrapolate_axes_for_pp(self, mz=None, abund=None, freq=None):
        """Extrapolate the m/z axis and fill the abundance axis with 0s.

        Parameters
        ----------
        mz : ndarray or None
            The m/z axis, if available. If None, the experimental m/z axis is used.
        abund : ndarray or None
            The abundance axis, if available. If None, the experimental abundance axis is used.
        freq : ndarray or None
            The frequency axis, if available. If None, the experimental frequency axis is used.

        Returns
        -------
        mz : ndarray
            The extrapolated m/z axis.
        abund : ndarray
            The abundance axis with 0s filled.
        freq : ndarray or None
            The extrapolated frequency axis, if available.

        Notes
        --------
        This function will extrapolate the mz axis by the number of datapoints specified in the settings,
        and fill the abundance axis with 0s.
        This should prevent peak picking issues at the spectrum edge.

        """
        # Check if the input arrays are provided
        if mz is None or abund is None:
            mz, abund = self.mz_exp_profile, self.abundance_profile
            if self.has_frequency:
                freq = self.freq_exp_profile
            else:
                freq = None
        pts = self.settings.picking_point_extrapolate
        if pts == 0:
            return mz, abund, freq

        mz = self.extrapolate_axis(mz, pts)
        abund = pad(abund, (pts, pts), mode="constant", constant_values=(0, 0))
        if freq is not None:
            freq = self.extrapolate_axis(freq, pts)
        return mz, abund, freq

    def do_peak_picking(self):
        """Perform peak picking."""
        mz, abundance, freq = self.prepare_peak_picking_data()

        if (
            self.label == Labels.bruker_frequency
            or self.label == Labels.midas_frequency
        ):
            self.calc_centroid(mz, abundance, freq)

        elif self.label == Labels.thermo_profile:
            self.calc_centroid(mz, abundance, freq)

        elif self.label == Labels.bruker_profile:
            self.calc_centroid(mz, abundance, freq)

        elif self.label == Labels.booster_profile:
            self.calc_centroid(mz, abundance, freq)

        elif self.label == Labels.simulated_profile:
            self.calc_centroid(mz, abundance, freq)

        else:
            raise Exception("Unknow mass spectrum type", self.label)

    def find_minima(self, apex_index, abundance, len_abundance, right=True):
        """Find the minima of a peak.

        Parameters
        ----------
        apex_index : int
            The index of the peak apex.
        abundance : ndarray
            The abundance values.
        len_abundance : int
            The length of the abundance array.
        right : bool, optional
            Flag indicating whether to search for minima to the right of the apex (default is True).

        Returns
        -------
        int
            The index of the minima.

        """
        j = apex_index

        if right:
            minima = abundance[j] > abundance[j + 1]
        else:
            minima = abundance[j] > abundance[j - 1]

        while minima:
            if j == 1 or j == len_abundance - 2:
                break

            if right:
                j += 1

                minima = abundance[j] >= abundance[j + 1]

            else:
                j -= 1
                minima = abundance[j] >= abundance[j - 1]

        if right:
            return j
        else:
            return j

    @staticmethod
    def linear_fit_calc(intes, massa, index_term, index_sign):
        """
        Algebraic solution to a linear fit - roughly 25-50x faster than numpy polyfit when passing only two vals and doing a 1st order fit

        Parameters
        ----------
        intes : ndarray
            The intensity values.
        massa : ndarray
            The mass values.
        index_term : int
            The index of the current term.
        index_sign : str
            The index sign

        Returns
        -------
        ndarray
            The coefficients of the linear fit.

        Notes
        --------
        This is a static method.
        """
        if index_sign == "+":
            x1, x2 = massa[index_term], massa[index_term + 1]
            y1, y2 = intes[index_term], intes[index_term + 1]
        elif index_sign == "-":
            x1, x2 = massa[index_term], massa[index_term - 1]
            y1, y2 = intes[index_term], intes[index_term - 1]
        else:
            warnings.warn("error in linear fit calc, unknown index sign")

        # Calculate the slope (m)
        slope = (y2 - y1) / (x2 - x1)

        # Calculate the intercept (b)
        intercept = y1 - slope * x1

        # The coefficients array would be [slope, intercept]
        coefficients = array([slope, intercept])
        return coefficients

    def calculate_resolving_power(self, intes, massa, current_index):
        """Calculate the resolving power of a peak.

        Parameters
        ----------
        intes : ndarray
            The intensity values.
        massa : ndarray
            The mass values.
        current_index : int
            The index of the current peak.

        Returns
        -------
        float
            The resolving power of the peak.

        Notes
        --------
        This is a conservative calculation of resolving power,
        the peak need to be resolved at least at the half-maximum magnitude,
        otherwise, the combined full width at half maximum is used to calculate resolving power.

        """

        peak_height = intes[current_index]
        target_peak_height = peak_height / 2

        peak_height_minus = peak_height
        peak_height_plus = peak_height

        # There are issues when a peak is at the high or low limit of a spectrum in finding its local minima and maxima
        # This solution will return nan for resolving power when a peak is possibly too close to an edge to avoid the issue

        if current_index < 5:
            warnings.warn("peak at low spectrum edge, returning no resolving power")
            return nan
        elif abs(current_index - len(intes)) < 5:
            warnings.warn("peak at high spectrum edge, returning no resolving power")
            return nan
        else:
            pass

        index_minus = current_index
        while peak_height_minus >= target_peak_height:
            index_minus = index_minus - 1
            if index_minus < 0:
                warnings.warn(
                    "Res. calc. warning - peak index minus adjacent to spectrum edge \n \
                        Zeroing the first 5 data points of abundance. Peaks at spectrum edge may be incorrectly reported \n \
                        Perhaps try to increase picking_point_extrapolate (e.g. to 3)"
                )
                # Pad the first 5 data points with zeros and restart the loop
                intes[:5] = 0
                peak_height_minus = target_peak_height
                index_minus = current_index
            else:
                peak_height_minus = intes[index_minus]

        if self.mspeaks_settings.legacy_centroid_polyfit:
            x = [massa[index_minus], massa[index_minus + 1]]
            y = [intes[index_minus], intes[index_minus + 1]]
            coefficients = polyfit(x, y, 1)
        else:
            coefficients = self.linear_fit_calc(
                intes, massa, index_minus, index_sign="+"
            )

        a = coefficients[0]
        b = coefficients[1]
        if self.mspeaks_settings.legacy_resolving_power:
            y_intercept = intes[index_minus] + (
                (intes[index_minus + 1] - intes[index_minus]) / 2
            )
        else:
            y_intercept = target_peak_height
        massa1 = (y_intercept - b) / a

        index_plus = current_index
        while peak_height_plus >= target_peak_height:
            index_plus = index_plus + 1

            try:
                peak_height_plus = intes[index_plus]
            except IndexError:
                warnings.warn(
                    "Res. calc. warning - peak index plus adjacent to spectrum edge \n \
                        Zeroing the last 5 data points of abundance. Peaks at spectrum edge may be incorrectly reported\
                        Perhaps try to increase picking_point_extrapolate (e.g. to 3)"
                )
                # Pad the first 5 data points with zeros and restart the loop
                intes[-5:] = 0
                peak_height_plus = target_peak_height
                index_plus = current_index

        if self.mspeaks_settings.legacy_centroid_polyfit:
            x = [massa[index_plus], massa[index_plus - 1]]
            y = [intes[index_plus], intes[index_plus - 1]]
            coefficients = polyfit(x, y, 1)
        else:
            coefficients = self.linear_fit_calc(
                intes, massa, index_plus, index_sign="-"
            )

        a = coefficients[0]
        b = coefficients[1]

        if self.mspeaks_settings.legacy_resolving_power:
            y_intercept = intes[index_plus - 1] + (
                (intes[index_plus] - intes[index_plus - 1]) / 2
            )
        else:
            y_intercept = target_peak_height

        massa2 = (y_intercept - b) / a

        if massa1 > massa2:
            resolvingpower = massa[current_index] / (massa1 - massa2)

        else:
            resolvingpower = massa[current_index] / (massa2 - massa1)

        return resolvingpower

    def cal_minima(self, mass, abun):
        """Calculate the minima of a peak.

        Parameters
        ----------
        mass : ndarray
            The mass values.
        abun : ndarray
            The abundance values.

        Returns
        -------
        ndarray or None
            The mass values at the minima, if found.

        """
        abun = -abun

        dy = abun[1:] - abun[:-1]

        # replaces nan for infinity
        indices_nan = where(isnan(abun))[0]

        if indices_nan.size:
            abun[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf

        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            return mass[indexes], abun[indexes]

    def calc_centroid(self, mass, abund, freq):
        """Calculate the centroid of a peak.

        Parameters
        ----------
        mass : ndarray
            The mass values.
        abund : ndarray
            The abundance values.
        freq : ndarray or None
            The frequency values, if available.

        Returns
        -------
        None

        """

        max_height = self.mspeaks_settings.peak_height_max_percent
        max_prominence = self.mspeaks_settings.peak_max_prominence_percent
        min_peak_datapoints = self.mspeaks_settings.min_peak_datapoints
        peak_derivative_threshold = self.mspeaks_settings.peak_derivative_threshold
        max_abun = max(abund)
        peak_height_diff = lambda hi, li: ((abund[hi] - abund[li]) / max_abun) * 100

        domain = mass
        signal = abund
        len_signal = len(signal)

        signal_threshold, factor = self.get_threshold(abund)
        max_signal = factor

        correct_baseline = False

        include_indexes = sp.peak_picking_first_derivative(
            domain,
            signal,
            max_height,
            max_prominence,
            max_signal,
            min_peak_datapoints,
            peak_derivative_threshold,
            signal_threshold=signal_threshold,
            correct_baseline=correct_baseline,
            abun_norm=1,
            plot_res=False,
        )

        for indexes_tuple in include_indexes:
            apex_index = indexes_tuple[1]

            peak_indexes = self.check_prominence(
                abund, apex_index, len_signal, peak_height_diff
            )

            if peak_indexes:
                mz_exp_centroid, freq_centr, intes_centr = self.find_apex_fit_quadratic(
                    mass, abund, freq, apex_index
                )

                if mz_exp_centroid:
                    peak_resolving_power = self.calculate_resolving_power(
                        abund, mass, apex_index
                    )
                    s2n = intes_centr / self.baseline_noise_std
                    self.add_mspeak(
                        self.polarity,
                        mz_exp_centroid,
                        abund[apex_index],
                        peak_resolving_power,
                        s2n,
                        indexes_tuple,
                        exp_freq=freq_centr,
                        ms_parent=self,
                    )
                # pyplot.plot(domain[start_index: final_index + 1], signal[start_index:final_index + 1], c='black')
                # pyplot.show()

    def get_threshold(self, intes):
        """Get the intensity threshold for peak picking.

        Parameters
        ----------
        intes : ndarray
            The intensity values.

        Returns
        -------
        float
            The intensity threshold.
        float
            The factor to multiply the intensity threshold by.
        """

        intes = array(intes).astype(float)

        noise_threshold_method = self.settings.noise_threshold_method

        if noise_threshold_method == "minima":
            if self.is_centroid:
                warnings.warn(
                    "Auto threshould is disabled for centroid data, returning 0"
                )
                factor = 1
                abundance_threshold = 1e-20
            # print(self.settings.noise_threshold_min_std)
            else:
                abundance_threshold = self.baseline_noise + (
                    self.settings.noise_threshold_min_std * self.baseline_noise_std
                )
                factor = 1

        elif noise_threshold_method == "signal_noise":
            abundance_threshold = self.settings.noise_threshold_min_s2n
            if self.is_centroid:
                factor = 1
            else:
                factor = self.baseline_noise_std

        elif noise_threshold_method == "relative_abundance":
            abundance_threshold = self.settings.noise_threshold_min_relative_abundance
            factor = intes.max() / 100

        elif noise_threshold_method == "absolute_abundance":
            abundance_threshold = self.settings.noise_threshold_absolute_abundance
            factor = 1

        elif noise_threshold_method == "log":
            if self.is_centroid:
                raise Exception("log noise Not tested for centroid data")
            abundance_threshold = self.settings.noise_threshold_log_nsigma
            factor = self.baseline_noise_std

        else:
            raise Exception(
                "%s method was not implemented, please refer to corems.mass_spectrum.calc.NoiseCalc Class"
                % noise_threshold_method
            )

        return abundance_threshold, factor

    @staticmethod
    def algebraic_quadratic(list_mass, list_y):
        """
        Find the apex of a peak - algebraically.
        Faster than using numpy polyfit by ~28x per fit.

        Parameters
        ----------
        list_mass : ndarray
            list of m/z values (3 points)
        list_y : ndarray
            list of abundance values (3 points)

        Returns
        -------
        a, b, c: float
            coefficients of the quadratic equation.

        Notes
        --------
        This is a static method.
        """
        x_1, x_2, x_3 = list_mass
        y_1, y_2, y_3 = list_y

        a = (
            y_1 / ((x_1 - x_2) * (x_1 - x_3))
            + y_2 / ((x_2 - x_1) * (x_2 - x_3))
            + y_3 / ((x_3 - x_1) * (x_3 - x_2))
        )

        b = (
            -y_1 * (x_2 + x_3) / ((x_1 - x_2) * (x_1 - x_3))
            - y_2 * (x_1 + x_3) / ((x_2 - x_1) * (x_2 - x_3))
            - y_3 * (x_1 + x_2) / ((x_3 - x_1) * (x_3 - x_2))
        )

        c = (
            y_1 * x_2 * x_3 / ((x_1 - x_2) * (x_1 - x_3))
            + y_2 * x_1 * x_3 / ((x_2 - x_1) * (x_2 - x_3))
            + y_3 * x_1 * x_2 / ((x_3 - x_1) * (x_3 - x_2))
        )
        return a, b, c

    def find_apex_fit_quadratic(self, mass, abund, freq, current_index):
        """
        Find the apex of a peak.

        Parameters
        ----------
        mass : ndarray
            The mass values.
        abund : ndarray
            The abundance values.
        freq : ndarray or None
            The frequency values, if available.
        current_index : int
            The index of the current peak.


        Returns
        -------
        float
            The m/z value of the peak apex.
        float
            The frequency value of the peak apex, if available.
        float
            The abundance value of the peak apex.

        """
        # calc prominence
        # peak_indexes = self.check_prominence(abund, current_index, len_abundance, peak_height_diff )

        # if not peak_indexes:

        #    return None, None, None, None

        # else:

        # fit parabola to three most abundant datapoints
        list_mass = [
            mass[current_index - 1],
            mass[current_index],
            mass[current_index + 1],
        ]
        list_y = [
            abund[current_index - 1],
            abund[current_index],
            abund[current_index + 1],
        ]

        if self.mspeaks_settings.legacy_centroid_polyfit:
            z = polyfit(list_mass, list_y, 2)
            a = z[0]
            b = z[1]
        else:
            a, b, c = self.algebraic_quadratic(list_mass, list_y)

        calculated = -b / (2 * a)

        if calculated < 1 or int(calculated) != int(list_mass[1]):
            mz_exp_centroid = list_mass[1]

        else:
            mz_exp_centroid = calculated

        if (
            self.label == Labels.bruker_frequency
            or self.label == Labels.midas_frequency
        ):
            # fit parabola to three most abundant frequency datapoints
            list_freq = [
                freq[current_index - 1],
                freq[current_index],
                freq[current_index + 1],
            ]
            if self.mspeaks_settings.legacy_centroid_polyfit:
                z = polyfit(list_mass, list_y, 2)
                a = z[0]
                b = z[1]
            else:
                a, b, c = self.algebraic_quadratic(list_mass, list_y)

            calculated_freq = -b / (2 * a)

            if calculated_freq < 1 or int(calculated_freq) != freq[current_index]:
                freq_centr = list_freq[1]

            else:
                freq_centr = calculated_freq

        else:
            freq_centr = None

        if self.mspeaks_settings.legacy_centroid_polyfit:
            abundance_centroid = abund[current_index]
        else:
            abundance_centroid = a * mz_exp_centroid**2 + b * mz_exp_centroid + c

        return mz_exp_centroid, freq_centr, abundance_centroid

    def check_prominence(
        self, abun, current_index, len_abundance, peak_height_diff
    ) -> tuple or False:
        """Check the prominence of a peak.

        Parameters
        ----------
        abun : ndarray
            The abundance values.
        current_index : int
            The index of the current peak.
        len_abundance : int
            The length of the abundance array.
        peak_height_diff : function
            The function to calculate the peak height difference.

        Returns
        -------
        tuple or False
            A tuple containing the indexes of the peak, if the prominence is above the threshold.
            Otherwise, False.

        """

        final_index = self.find_minima(current_index, abun, len_abundance, right=True)

        start_index = self.find_minima(current_index, abun, len_abundance, right=False)

        peak_indexes = (current_index - 1, current_index, current_index + 1)

        if (
            min(
                peak_height_diff(current_index, start_index),
                peak_height_diff(current_index, final_index),
            )
            > self.mspeaks_settings.peak_min_prominence_percent
        ):
            return peak_indexes

        else:
            return False

    def use_the_max(self, mass, abund, current_index, len_abundance, peak_height_diff):
        """Use the max peak height as the centroid

        Parameters
        ----------
        mass : ndarray
            The mass values.
        abund : ndarray
            The abundance values.
        current_index : int
            The index of the current peak.
        len_abundance : int
            The length of the abundance array.
        peak_height_diff : function
            The function to calculate the peak height difference.

        Returns
        -------
        float
            The m/z value of the peak apex.
        float
            The abundance value of the peak apex.
        tuple or None
            A tuple containing the indexes of the peak, if the prominence is above the threshold.
            Otherwise, None.
        """

        peak_indexes = self.check_prominence(
            abund, current_index, len_abundance, peak_height_diff
        )

        if not peak_indexes:
            return None, None, None

        else:
            return mass[current_index], abund[current_index], peak_indexes

    def calc_centroid_legacy(self, mass, abund, freq):
        """Legacy centroid calculation
        Deprecated - for deletion.

        """
        warnings.warn(
            "Legacy centroid calculation is deprecated. Please use the new centroid calculation method."
        )
        pass
        if False:
            len_abundance = len(abund)

            max_abundance = max(abund)

            peak_height_diff = (
                lambda hi, li: ((abund[hi] - abund[li]) / max_abundance) * 100
            )

            abundance_threshold, factor = self.get_threshold(abund)
            # print(abundance_threshold, factor)
            # find indices of all peaks
            dy = abund[1:] - abund[:-1]

            # replaces nan for infi nity
            indices_nan = where(isnan(abund))[0]

            if indices_nan.size:
                abund[indices_nan] = inf
                dy[where(isnan(dy))[0]] = inf

            indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

            # noise threshold
            if indexes.size and abundance_threshold is not None:
                indexes = indexes[abund[indexes] / factor >= abundance_threshold]
            # filter out 'peaks' within 3 points of the spectrum limits
            # remove entries within 3 points of upper limit
            indexes = [x for x in indexes if (len_abundance - x) > 3]
            # remove entries within 3 points of zero
            indexes = [x for x in indexes if x > 3]

            for current_index in indexes:
                if self.label == Labels.simulated_profile:
                    mz_exp_centroid, intes_centr, peak_indexes = self.use_the_max(
                        mass, abund, current_index, len_abundance, peak_height_diff
                    )
                    if mz_exp_centroid:
                        peak_resolving_power = self.calculate_resolving_power(
                            abund, mass, current_index
                        )
                        s2n = intes_centr / self.baseline_noise_std
                        freq_centr = None
                        self.add_mspeak(
                            self.polarity,
                            mz_exp_centroid,
                            abund[current_index],
                            peak_resolving_power,
                            s2n,
                            peak_indexes,
                            exp_freq=freq_centr,
                            ms_parent=self,
                        )

                else:
                    mz_exp_centroid, freq_centr, intes_centr, peak_indexes = (
                        self.find_apex_fit_quadratic(
                            mass,
                            abund,
                            freq,
                            current_index,
                            len_abundance,
                            peak_height_diff,
                        )
                    )
                    if mz_exp_centroid:
                        try:
                            peak_resolving_power = self.calculate_resolving_power(
                                abund, mass, current_index
                            )
                        except IndexError:
                            print("index error, skipping peak")
                            continue

                        s2n = intes_centr / self.baseline_noise_std
                        self.add_mspeak(
                            self.polarity,
                            mz_exp_centroid,
                            abund[current_index],
                            peak_resolving_power,
                            s2n,
                            peak_indexes,
                            exp_freq=freq_centr,
                            ms_parent=self,
                        )
