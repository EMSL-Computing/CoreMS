__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"

from threading import Thread
from dataclasses import dataclass
from typing import List


@dataclass
class SearchResults:
    calculated_mz: float
    exp_mz: float
    error: float
    tolerance: float


class MZSearch(Thread):
    def __init__(
        self,
        exp_mzs: List[float],
        calculated_mzs: List[float],
        tolerance,
        method="ppm",
        average_target_mz=True,
    ):
        """
        Parameters
        ----------
        calculated_mzs: [float] calculated m/z
        exp_mzs: [float] experimental m/z
        method: string,
            ppm or ppb
        call run to trigger the m/z search algorithm
        or start if using it as thread
        """
        Thread.__init__(self)
        # placeholder for the results
        self._matched_mz = {}

        self._calculated_mzs = calculated_mzs

        self._matched_mz = {}

        self._averaged_target_mz = []

        self._exp_mzs = exp_mzs

        self._tolerance = tolerance
        self.method = method

        if average_target_mz:
            self.colapse_calculated()

    @property
    def results(self):
        """{calculated_mz: [SearchResults]}
        contains the results of the search
        """
        return self._matched_mz

    @property
    def averaged_target_mz(self):
        """[float]
        contains the average target m/z to be searched against
        """
        return self._averaged_target_mz

    @property
    def calculated_mzs(self):
        """[float]
        contains the mz target to be searched against
        """
        if self.averaged_target_mz:
            return sorted(self.averaged_target_mz)
        else:
            return sorted(list(self._calculated_mzs))

    @property
    def exp_mzs(self):
        """[float]
        contains the exp mz to be searched against
        """
        return self._exp_mzs

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, method):
        """
        method: string,
           ppm or ppb
        """
        if method not in ["ppm" or "ppb"]:
            raise ValueError("Method should be ppm or ppb")
        self._method = method

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, tolerance):
        """
        method: string,
           ppm or ppb
        """
        if tolerance < 0:
            raise ValueError("Tolerance needs to be a positive number")
        self._tolerance = tolerance

    def colapse_calculated(self):
        if len(self.calculated_mzs) > 1:
            all_mz = []
            subset = set()

            i = -1
            while True:
                i = i + 1

                if i == len(self.calculated_mzs) - 1:
                    all_mz.append({i})
                    # print(i, 'break1')
                    break

                if i >= len(self.calculated_mzs) - 1:
                    # print(i, 'break2')
                    break

                error = self.calc_mz_error(
                    self.calculated_mzs[i], self.calculated_mzs[i + 1]
                )

                # print(self.tolerance)

                check_error = self.check_ppm_error(self.tolerance, error)

                if not check_error:
                    start_list = {i}

                else:
                    start_list = set()

                while check_error:
                    start_list.add(i)
                    start_list.add(i + 1)

                    i = i + 1

                    if i == len(self.calculated_mzs) - 1:
                        start_list.add(i)
                        # print(i, 'break3')
                        break

                    error = self.calc_mz_error(
                        self.calculated_mzs[i], self.calculated_mzs[i + 1]
                    )
                    check_error = self.check_ppm_error(self.tolerance, error)

                if start_list:
                    all_mz.append(start_list)

            results = []
            for each in all_mz:
                # print(each)
                mzs = [self.calculated_mzs[i] for i in each]
                results.append(sum(mzs) / len(mzs))

            # print(results)
            self._averaged_target_mz = results

    def run(self):
        dict_nominal_exp_mz = self.get_nominal_exp(self.exp_mzs)

        for calculated_mz in self.calculated_mzs:
            nominal_selected_mz = int(calculated_mz)

            if nominal_selected_mz in dict_nominal_exp_mz.keys():
                self.search_mz(self.results, dict_nominal_exp_mz, calculated_mz, 0)

            elif nominal_selected_mz - 1 in dict_nominal_exp_mz.keys():
                self.search_mz(self.results, dict_nominal_exp_mz, calculated_mz, -1)

            elif nominal_selected_mz + 1 in dict_nominal_exp_mz.keys():
                self.search_mz(self.results, dict_nominal_exp_mz, calculated_mz, +1)

            else:
                continue

    @staticmethod
    def calc_mz_error(calculated_mz, exp_mz, method="ppm"):
        """
        Parameters
        ----------
        calculated_mz: float,
        exp_mz:float
        method: string,
            ppm or ppb
        """
        if method == "ppm":
            multi_factor = 1000000

        elif method == "ppb":
            multi_factor = 1000000

        else:
            raise Exception(
                "method needs to be ppm or ppb, \
                             you have entered %s"
                % method
            )

        return ((exp_mz - calculated_mz) / calculated_mz) * multi_factor

    @staticmethod
    def check_ppm_error(tolerance, error):
        return True if -tolerance <= error <= tolerance else False

    def get_nominal_exp(self, exp_mzs) -> dict:
        dict_nominal_exp_mz = {}

        for exp_mz in exp_mzs:
            nominal_mz = int(exp_mz)

            if nominal_mz not in dict_nominal_exp_mz.keys():
                dict_nominal_exp_mz[int(exp_mz)] = [exp_mz]
            else:
                dict_nominal_exp_mz[int(exp_mz)].append(exp_mz)

        return dict_nominal_exp_mz

    def search_mz(self, results, dict_nominal_exp_mz, calculated_mz, offset) -> None:
        nominal_calculated_mz = int(calculated_mz) + offset
        matched_n_precursors = dict_nominal_exp_mz.get(nominal_calculated_mz)

        for precursor_mz in matched_n_precursors:
            error = self.calc_mz_error(calculated_mz, precursor_mz, method=self.method)

            if self.check_ppm_error(self.tolerance, error):
                new_match = SearchResults(
                    calculated_mz, precursor_mz, error, self.tolerance
                )

                if calculated_mz not in results.keys():
                    results[calculated_mz] = [new_match]

                else:
                    results[calculated_mz].append(new_match)
