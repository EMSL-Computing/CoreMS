# from numpy import array, polyfit, poly1d, where, trapz
from numpy import trapz

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
        self._area = trapz(yy, dx=dx)

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
