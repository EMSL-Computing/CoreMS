import numpy as np
from bisect import bisect_left

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
        intense point, as described in [1].

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
        eic_in = self.eic_list.argsort()
        sorted_eic = self.eic_list[eic_in[::-1]]
        sorted_rt = self.eic_rt_list[eic_in[::-1]]

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
