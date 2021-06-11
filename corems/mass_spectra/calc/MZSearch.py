__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"

from threading import Thread
from dataclasses import dataclass
from corems import Vector

@dataclass
class SearchResults:

    calculated_mz: float
    exp_mz: float
    error: float
    tolerance: float

class MZSearch(Thread):

    def __init__(self, exp_mzs: Vector, calculated_mzs: Vector, tolerance, method="ppm"):
        '''
        Parameters
        ----------
        calculated_mzs: [float] calculated m/z
        exp_mzs: [float] experimental m/z
        method: string,
            ppm or ppb
        call run to trigger the m/z search algorithm
        or start if using it as thread
        '''
        Thread.__init__(self)
        # placeholder for the results
        self._matched_mz = {}

        self._calculated_mzs = calculated_mzs
        self._exp_mzs = exp_mzs

        self.tolerance = tolerance
        self.method = method

    @property
    def results(self):
        ''' {calculated_mz: [SearchResults]}
            contains the results of the search
        '''
        return self._matched_mz

    @property
    def calculated_mzs(self):
        ''' [float]
            contains the mz target to be searched against
        '''
        return self._calculated_mzs

    @property
    def exp_mzs(self):
        ''' [float]
            contains the exp mz to be searched against
        '''
        return self._exp_mzs

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, method):
        '''
         method: string,
            ppm or ppb
        '''
        if method not in ['ppm' or 'ppb']:
            raise ValueError("Method should be ppm or ppb")
        self._method = method

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, tolerance):
        '''
         method: string,
            ppm or ppb
        '''
        if tolerance < 0:
            raise ValueError("Tolerance needs to be a positive number")
        self._tolerance = tolerance

    def run(self):

        dict_nominal_exp_mz = self.get_nominal_exp(self.exp_mzs)

        for calculated_mz in self.calculated_mzs:

            nominal_selected_mz = int(calculated_mz)

            if nominal_selected_mz in dict_nominal_exp_mz.keys():

                self.search_mz(dict_nominal_exp_mz, calculated_mz, 0)

            elif nominal_selected_mz - 1 in dict_nominal_exp_mz.keys():

                self.search_mz(dict_nominal_exp_mz, calculated_mz, -1)

            elif nominal_selected_mz + 1 in dict_nominal_exp_mz.keys():

                self.search_mz(dict_nominal_exp_mz, calculated_mz, +1)

            else:

                continue

    @staticmethod
    def calc_mz_error(calculated_mz, exp_mz, method='ppm'):
        '''
        Parameters
        ----------
        calculated_mz: float,
        exp_mz:float
        method: string,
            ppm or ppb
        '''
        if method == 'ppm':
            multi_factor = 1000000

        elif method == 'ppb':
            multi_factor = 1000000

        else:
            raise Exception("method needs to be ppm or ppb, \
                             you have entered %s" % method)

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

    def search_mz(self, dict_nominal_exp_mz, calculated_mz, offset) -> None:

        nominal_calculated_mz = int(calculated_mz) + offset
        matched_n_precursors = dict_nominal_exp_mz.get(nominal_calculated_mz)

        for precursor_mz in matched_n_precursors:
            error = self.calc_mz_error(calculated_mz, precursor_mz,
                                       method=self.method)

            if self.check_ppm_error(self.tolerance, error):

                new_match = SearchResults(calculated_mz,
                                          precursor_mz,
                                          error,
                                          self.tolerance)

                if calculated_mz not in self.results.keys():
                    self.results[calculated_mz] = [new_match]
                else:
                    self.results[calculated_mz].append(new_match)
