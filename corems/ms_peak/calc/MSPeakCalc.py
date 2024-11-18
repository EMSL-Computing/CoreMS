__author__ = "Yuri E. Corilo"
__date__ = "Jun 04, 2019"

import warnings

import pyswarm
from lmfit import models
from numpy import (
    ceil,
    exp,
    flip,
    floor,
    linspace,
    log,
    nan,
    pi,
    poly1d,
    polyfit,
    rint,
    sqrt,
    square,
    trapz,
)

from corems.encapsulation.constant import Atoms
from corems.encapsulation.factory.parameters import MSParameters


class MSPeakCalculation:
    """Class to perform calculations on MSPeak objects.

    This class provides methods to perform various calculations on MSPeak objects, such as calculating Kendrick Mass Defect (KMD) and Kendrick Mass (KM), calculating peak area, and fitting peak lineshape using different models.

    Parameters
    ----------
    None

    Attributes
    ----------
    _ms_parent : MSParent
        The parent MSParent object associated with the MSPeakCalculation object.
    mz_exp : float
        The experimental m/z value of the peak.
    peak_left_index : int
        The start scan index of the peak.
    peak_right_index : int
        The final scan index of the peak.
    resolving_power : float
        The resolving power of the peak.

    Methods
    -------
    * _calc_kmd(dict_base).
        Calculate the Kendrick Mass Defect (KMD) and Kendrick Mass (KM) for a given base formula.
    * calc_area().
        Calculate the peak area using numpy's trapezoidal fit.
    * fit_peak(mz_extend=6, delta_rp=0, model='Gaussian').
        Perform lineshape analysis on a peak using lmfit module.
    * voigt_pso(w, r, yoff, width, loc, a).
        Calculate the Voigt function for particle swarm optimization (PSO) fitting.
    * objective_pso(x, w, u).
        Calculate the objective function for PSO fitting.
    * minimize_pso(lower, upper, w, u).
        Minimize the objective function using the particle swarm optimization algorithm.
    * fit_peak_pso(mz_extend=6, upsample_multiplier=5).
        Perform lineshape analysis on a peak using particle swarm optimization (PSO) fitting.
    * voigt(oversample_multiplier=1, delta_rp=0, mz_overlay=1).
        [Legacy] Perform voigt lineshape analysis on a peak.
    * pseudovoigt(oversample_multiplier=1, delta_rp=0, mz_overlay=1, fraction=0.5).
        [Legacy] Perform pseudovoigt lineshape analysis on a peak.
    * lorentz(oversample_multiplier=1, delta_rp=0, mz_overlay=1).
        [Legacy] Perform lorentz lineshape analysis on a peak.
    * gaussian(oversample_multiplier=1, delta_rp=0, mz_overlay=1).
        [Legacy] Perform gaussian lineshape analysis on a peak.
    * get_mz_domain(oversample_multiplier, mz_overlay).
        [Legacy] Resample/interpolate datapoints for lineshape analysis.
    * number_possible_assignments().
        Return the number of possible molecular formula assignments for the peak.
    * molecular_formula_lowest_error().
        Return the molecular formula with the smallest absolute mz error.
    * molecular_formula_highest_prob_score().
        Return the molecular formula with the highest confidence score.
    * molecular_formula_earth_filter(lowest_error=True).
        Filter molecular formula using the 'Earth' filter.
    * molecular_formula_water_filter(lowest_error=True).
        Filter molecular formula using the 'Water' filter.
    * molecular_formula_air_filter(lowest_error=True).
        Filter molecular formula using the 'Air' filter.
    * cia_score_S_P_error().
        Compound Identification Algorithm SP Error - Assignment Filter.
    * cia_score_N_S_P_error().
        Compound Identification Algorithm NSP Error - Assignment Filter.

    """

    def _calc_kmd(self, dict_base):
        """Calculate the Kendrick Mass Defect (KMD) and Kendrick Mass (KM) for a given base formula

        Parameters
        ----------
        dict_base : dict
            dictionary with the base formula to be used in the calculation
            Default is CH2, e.g.
                dict_base = {"C": 1, "H": 2}
        """

        if self._ms_parent:
            # msPeak obj does have a ms object parent
            kendrick_rounding_method = (
                self._ms_parent.mspeaks_settings.kendrick_rounding_method
            )  # rounding method can be one of floor, ceil or round
            # msPeak obj does not have a ms object parent
        else:
            kendrick_rounding_method = MSParameters.ms_peak.kendrick_rounding_method

        mass = 0
        for atom in dict_base.keys():
            mass += Atoms.atomic_masses.get(atom) * dict_base.get(atom)

        kendrick_mass = (int(mass) / mass) * self.mz_exp

        if kendrick_rounding_method == "ceil":
            nominal_km = ceil(kendrick_mass)

        elif kendrick_rounding_method == "round":
            nominal_km = rint(kendrick_mass)

        elif kendrick_rounding_method == "floor":
            nominal_km = floor(kendrick_mass)

        else:
            raise Exception(
                "%s method was not implemented, please refer to corems.ms_peak.calc.MSPeakCalc Class"
                % kendrick_rounding_method
            )

        kmd = nominal_km - kendrick_mass

        # kmd = (nominal_km - km) * 1
        # kmd = round(kmd,0)

        return kmd, kendrick_mass, nominal_km

    def calc_area(self):
        """Calculate the peak area using numpy's trapezoidal fit

        uses provided mz_domain to accurately integrate areas independent of digital resolution

        Returns
        -------
        float
            peak area
        """
        if self.peak_right_index > self.peak_left_index:
            yy = self._ms_parent.abundance_profile[
                self.peak_left_index : self.peak_right_index
            ]
            xx = self._ms_parent.mz_exp_profile[
                self.peak_left_index : self.peak_right_index
            ]
            # check if the axis is high to low m/z or not. if its MSFromFreq its high mz first, if its from Profile, its low mz first
            if xx[0] > xx[-1]:
                xx = flip(xx)
                yy = flip(yy)
            return float(trapz(yy, xx))

        else:
            warnings.warn(
                "Peak Area Calculation for m/z {} has failed".format(self.mz_exp)
            )
            return nan

    def fit_peak(self, mz_extend=6, delta_rp=0, model="Gaussian"):
        """Lineshape analysis on a peak using lmfit module.

        Model and fit peak lineshape by defined function - using lmfit module
        Does not oversample/resample/interpolate data points
        Better to go back to time domain and perform more zero filling - if possible.

        Parameters
        ----------
        mz_extend : int
            extra points left and right of peak definition to include in fitting
        delta_rp : float
            delta resolving power to add to resolving power
        model : str
            Type of lineshape model to use.
            Models allowed: Gaussian, Lorentz, Voigt

        Returns
        -----
        mz_domain : ndarray
            x-axis domain for fit
        fit_peak : lmfit object
            fit results object from lmfit module

        Notes
        -----
        Returns the calculated mz domain, initial defined abundance profile, and the fit peak results object from lmfit module
        mz_extend here extends the x-axis domain so that we have sufficient points either side of the apex to fit.
        Takes about 10ms per peak
        """
        start_index = (
            self.peak_left_index - mz_extend if not self.peak_left_index == 0 else 0
        )
        final_index = (
            self.peak_right_index + mz_extend
            if not self.peak_right_index == len(self._ms_parent.mz_exp_profile)
            else self.peak_right_index
        )

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / (self.resolving_power + delta_rp)

            mz_domain = self._ms_parent.mz_exp_profile[start_index:final_index]
            abundance_domain = self._ms_parent.abundance_profile[
                start_index:final_index
            ]

            if model == "Gaussian":
                # stardard deviation
                sigma = self.fwhm / (2 * sqrt(2 * log(2)))
                amplitude = (sqrt(2 * pi) * sigma) * self.abundance
                model = models.GaussianModel()
                params = model.make_params(
                    center=self.mz_exp, amplitude=amplitude, sigma=sigma
                )

            elif model == "Lorentz":
                # stardard deviation
                sigma = self.fwhm / 2
                amplitude = sigma * pi * self.abundance
                model = models.LorentzianModel()
                params = model.make_params(
                    center=self.mz_exp, amplitude=amplitude, sigma=sigma
                )

            elif model == "Voigt":
                # stardard deviation
                sigma = self.fwhm / 3.6013
                amplitude = (sqrt(2 * pi) * sigma) * self.abundance
                model = models.VoigtModel()
                params = model.make_params(
                    center=self.mz_exp, amplitude=amplitude, sigma=sigma, gamma=sigma
                )
            else:
                raise LookupError("model lineshape not known or defined")

            # calc_abundance = model.eval(params=params, x=mz_domain) #Same as initial fit, returned in fit_peak object
            fit_peak = model.fit(abundance_domain, params=params, x=mz_domain)
            return mz_domain, fit_peak

        else:
            raise LookupError(
                "resolving power is not defined, try to use set_max_resolving_power()"
            )

    def voigt_pso(self, w, r, yoff, width, loc, a):
        """Voigt function for particle swarm optimisation (PSO) fitting

        From https://github.com/pnnl/nmrfit/blob/master/nmrfit/equations.py.
        Calculates a Voigt function over w based on the relevant properties of the distribution.

        Parameters
        ----------
        w : ndarray
            Array over which the Voigt function will be evaluated.
        r : float
            Ratio between the Guassian and Lorentzian functions.
        yoff : float
            Y-offset of the Voigt function.
        width : float
            The width of the Voigt function.
        loc : float
            Center of the Voigt function.
        a : float
            Area of the Voigt function.
        Returns
        -------
        V : ndarray
            Array defining the Voigt function over w.

        References
        ----------
        1. https://github.com/pnnl/nmrfit

        Notes
        -----
        Particle swarm optimisation (PSO) fitting function can be significantly more computationally expensive than lmfit, with more parameters to optimise.

        """
        # Lorentzian component
        L = (2 / (pi * width)) * 1 / (1 + ((w - loc) / (0.5 * width)) ** 2)

        # Gaussian component
        G = (
            (2 / width)
            * sqrt(log(2) / pi)
            * exp(-(((w - loc) / (width / (2 * sqrt(log(2))))) ** 2))
        )

        # Voigt body
        V = (yoff + a) * (r * L + (1 - r) * G)

        return V

    def objective_pso(self, x, w, u):
        """Objective function for particle swarm optimisation (PSO) fitting

        The objective function used to fit supplied data.  Evaluates sum of squared differences between the fit and the data.

        Parameters
        ----------
        x : list of floats
            Parameter vector.
        w : ndarray
            Array of frequency data.
        u : ndarray
            Array of data to be fit.

        Returns
        -------
        rmse : float
            Root mean square error between the data and fit.

        References
        ----------
        1. https://github.com/pnnl/nmrfit

        """
        # global parameters
        r, width, loc, a = x
        yoff = 0

        # calculate fit for V
        V_fit = self.voigt_pso(w, r, yoff, width, loc, a)

        # real component RMSE
        rmse = sqrt(square((u - V_fit)).mean(axis=None))

        # return the total RMSE
        return rmse

    def minimize_pso(self, lower, upper, w, u):
        """Minimization function for particle swarm optimisation (PSO) fitting

        Minimizes the objective function using the particle swarm optimization algorithm.
        Minimization function based on defined parameters


        Parameters
        ----------
        lower : list of floats
            Lower bounds for the parameters.
        upper : list of floats
            Upper bounds for the parameters.
        w : ndarray
            Array of frequency data.
        u : ndarray
            Array of data to be fit.

        Notes
        -----
        Particle swarm optimisation (PSO) fitting function can be significantly more computationally expensive than lmfit, with more parameters to optimise.
        Current parameters take ~2 seconds per peak.


        References
        ----------
        1. https://github.com/pnnl/nmrfit

        """
        # TODO - allow support to pass swarmsize, maxiter, omega, phip, phig parameters.
        # TODO - Refactor PSO fitting into its own class?

        xopt, fopt = pyswarm.pso(
            self.objective_pso,
            lower,
            upper,
            args=(w, u),
            swarmsize=1000,
            maxiter=5000,
            omega=-0.2134,
            phip=-0.3344,
            phig=2.3259,
        )
        return xopt, fopt

    def fit_peak_pso(self, mz_extend: int = 6, upsample_multiplier: int = 5):
        """Lineshape analysis on a peak using particle swarm optimisation (PSO) fitting

        Function to fit a Voigt peakshape using particle swarm optimisation (PSO).
        Should return better results than lmfit, but much more computationally expensive

        Parameters
        ----------
        mz_extend : int, optional
            extra points left and right of peak definition to include in fitting. Defaults to 6.
        upsample_multiplier : int, optional
            factor to increase x-axis points by for simulation of fitted lineshape function. Defaults to 5.

        Returns
        -------
        xopt : array
            variables describing the voigt function.
            G/L ratio, width (fwhm), apex (x-axis), area.
            y-axis offset is fixed at 0
        fopt : float
            objective score (rmse)
        psfit : array
            recalculated y values based on function and optimised fit
        psfit_hdp : tuple of arrays
            0 - linspace x-axis upsampled grid
            1 - recalculated y values based on function and upsampled x-axis grid
            Does not change results, but aids in visualisation of the 'true' voigt lineshape

        Notes
        -----
        Particle swarm optimisation (PSO) fitting function can be significantly more computationally expensive than lmfit, with more parameters to optimise.
        """
        # TODO - Add ability to pass pso args (i.e. swarm size, maxiter, omega, phig, etc)
        # TODO: fix xopt. Magnitude mode data through CoreMS/Bruker starts at 0 but is noise centered well above 0.
        # Thermo data is noise reduced by also noise subtracted, so starts at 0
        # Absorption mode/phased data will have positive and negative components and may not be baseline corrected

        start_index = (
            self.peak_left_index - mz_extend if not self.peak_left_index == 0 else 0
        )
        final_index = (
            self.peak_right_index + mz_extend
            if not self.peak_right_index == len(self._ms_parent.mz_exp_profile)
            else self.peak_right_index
        )

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / (self.resolving_power)

            mz_domain = self._ms_parent.mz_exp_profile[start_index:final_index]
            abundance_domain = self._ms_parent.abundance_profile[
                start_index:final_index
            ]
            lower = [0, self.fwhm * 0.8, (self.mz_exp - 0.0005), 0]
            upper = [
                1,
                self.fwhm * 1.2,
                (self.mz_exp + 0.0005),
                self.abundance / self.signal_to_noise,
            ]
            xopt, fopt = self.minimize_pso(lower, upper, mz_domain, abundance_domain)

            psfit = self.voigt_pso(mz_domain, xopt[0], 0, xopt[1], xopt[2], xopt[3])
            psfit_hdp_x = linspace(
                min(mz_domain), max(mz_domain), num=len(mz_domain) * upsample_multiplier
            )
            psfit_hdp = self.voigt_pso(
                psfit_hdp_x, xopt[0], 0, xopt[1], xopt[2], xopt[3]
            )
            return xopt, fopt, psfit, (psfit_hdp_x, psfit_hdp)
        else:
            raise LookupError(
                "resolving power is not defined, try to use set_max_resolving_power()"
            )

    def voigt(self, oversample_multiplier=1, delta_rp=0, mz_overlay=1):
        """[Legacy] Voigt lineshape analysis function
        Legacy function for voigt lineshape analysis

        Parameters
        ----------
        oversample_multiplier : int
            factor to increase x-axis points by for simulation of fitted lineshape function
        delta_rp : float
            delta resolving power to add to resolving power
        mz_overlay : int
            extra points left and right of peak definition to include in fitting

        Returns
        -------
        mz_domain : ndarray
            x-axis domain for fit
        calc_abundance : ndarray
            calculated abundance profile based on voigt function
        """

        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / (
                self.resolving_power + delta_rp
            )  # self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / 3.6013

            # half width baseline distance

            # mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)
            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)

            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )

            # TODO derive amplitude
            amplitude = (sqrt(2 * pi) * sigma) * self.abundance

            model = models.VoigtModel()

            params = model.make_params(
                center=self.mz_exp, amplitude=amplitude, sigma=sigma, gamma=sigma
            )

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance

        else:
            raise LookupError(
                "resolving power is not defined, try to use set_max_resolving_power()"
            )

    def pseudovoigt(
        self, oversample_multiplier=1, delta_rp=0, mz_overlay=1, fraction=0.5
    ):
        """[Legacy] pseudovoigt lineshape function

        Legacy function for pseudovoigt lineshape analysis.
        Note - Code may not be functional currently.

        Parameters
        ----------
        oversample_multiplier : int, optional
            factor to increase x-axis points by for simulation of fitted lineshape function. Defaults to 1.
        delta_rp : float, optional
            delta resolving power to add to resolving power. Defaults to 0.
        mz_overlay : int, optional
            extra points left and right of peak definition to include in fitting. Defaults to 1.
        fraction : float, optional
            fraction of gaussian component in pseudovoigt function. Defaults to 0.5.

        """
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / (
                self.resolving_power + delta_rp
            )  # self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / 2

            # half width baseline distance

            # mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)
            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)

            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            model = models.PseudoVoigtModel()

            # TODO derive amplitude
            gamma = sigma

            amplitude = (sqrt(2 * pi) * sigma) * self.abundance
            amplitude = (sqrt(pi / log(2)) * (pi * sigma * self.abundance)) / (
                (pi * (1 - gamma)) + (sqrt(pi * log(2)) * gamma)
            )

            params = model.make_params(center=self.mz_exp, sigma=sigma)

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance

        else:
            raise LookupError(
                "resolving power is not defined, try to use set_max_resolving_power()"
            )

    def lorentz(self, oversample_multiplier=1, delta_rp=0, mz_overlay=1):
        """[Legacy] Lorentz lineshape analysis function

        Legacy function for lorentz lineshape analysis

        Parameters
        ----------
        oversample_multiplier : int
            factor to increase x-axis points by for simulation of fitted lineshape function
        delta_rp : float
            delta resolving power to add to resolving power
        mz_overlay : int
            extra points left and right of peak definition to include in fitting

        Returns
        -------
        mz_domain : ndarray
            x-axis domain for fit
        calc_abundance : ndarray
            calculated abundance profile based on lorentz function

        """
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / (
                self.resolving_power + delta_rp
            )  # self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / 2

            # half width baseline distance
            hw_base_distance = 8 * sigma

            # mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)

            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)
            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            model = models.LorentzianModel()

            amplitude = sigma * pi * self.abundance

            params = model.make_params(
                center=self.mz_exp, amplitude=amplitude, sigma=sigma
            )

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance

        else:
            raise LookupError(
                "resolving power is not defined, try to use set_max_resolving_power()"
            )

    def gaussian(self, oversample_multiplier=1, delta_rp=0, mz_overlay=1):
        """[Legacy] Gaussian lineshape analysis function
        Legacy gaussian lineshape analysis function

        Parameters
        ----------
        oversample_multiplier : int
            factor to increase x-axis points by for simulation of fitted lineshape function
        delta_rp : float
            delta resolving power to add to resolving power
        mz_overlay : int
            extra points left and right of peak definition to include in fitting

        Returns
        -------
        mz_domain : ndarray
            x-axis domain for fit
        calc_abundance : ndarray
            calculated abundance profile based on gaussian function


        """

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / (
                self.resolving_power + delta_rp
            )  # self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / (2 * sqrt(2 * log(2)))

            # half width baseline distance
            # hw_base_distance = (3.2 * s)

            # match_loz_factor = 3

            # n_d = hw_base_distance * match_loz_factor

            # mz_domain = linspace(
            #    self.mz_exp - n_d, self.mz_exp + n_d, datapoint)

            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)

            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )

            # calc_abundance = norm.pdf(mz_domain, self.mz_exp, s)

            model = models.GaussianModel()

            amplitude = (sqrt(2 * pi) * sigma) * self.abundance

            params = model.make_params(
                center=self.mz_exp, amplitude=amplitude, sigma=sigma
            )

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance

        else:
            raise LookupError(
                "resolving power is not defined, try to use set_max_resolving_power()"
            )

    def get_mz_domain(self, oversample_multiplier, mz_overlay):
        """[Legacy] function to resample/interpolate datapoints for lineshape analysis

        This code is used for the legacy line fitting functions and not recommended.
        Legacy function to support expanding mz domain for legacy lineshape functions

        Parameters
        ----------
        oversample_multiplier : int
            factor to increase x-axis points by for simulation of fitted lineshape function
        mz_overlay : int
            extra points left and right of peak definition to include in fitting

        Returns
        -------
        mz_domain : ndarray
            x-axis domain for fit

        """
        start_index = (
            self.peak_left_index - mz_overlay if not self.peak_left_index == 0 else 0
        )
        final_index = (
            self.peak_right_index + mz_overlay
            if not self.peak_right_index == len(self._ms_parent.mz_exp_profile)
            else self.peak_right_index
        )

        if oversample_multiplier == 1:
            mz_domain = self._ms_parent.mz_exp_profile[start_index:final_index]

        else:
            # we assume a linear correlation for m/z and datapoits
            # which is only true if the m/z range in narrow (within 1 m/z unit)
            # this is not true for a wide m/z range

            indexes = range(start_index, final_index + 1)
            mz = self._ms_parent.mz_exp_profile[indexes]
            pol = poly1d(polyfit(indexes, mz, 1))
            oversampled_indexes = linspace(
                start_index,
                final_index,
                (final_index - start_index) * oversample_multiplier,
            )
            mz_domain = pol(oversampled_indexes)

        return mz_domain

    @property
    def number_possible_assignments(
        self,
    ):
        return len(self.molecular_formulas)

    def molecular_formula_lowest_error(self):
        """Return the molecular formula with the smallest absolute mz error"""

        return min(self.molecular_formulas, key=lambda m: abs(m.mz_error))

    def molecular_formula_highest_prob_score(self):
        """Return the molecular formula with the highest confidence score score"""

        return max(self.molecular_formulas, key=lambda m: abs(m.confidence_score))

    def molecular_formula_earth_filter(self, lowest_error=True):
        """Filter molecular formula using the 'Earth' filter

        This function applies the Formularity-esque 'Earth' filter to possible molecular formula assignments.
        Earth Filter:
            O > 0 AND N <= 3 AND P <= 2 AND 3P <= O

        If the lowest_error method is also used, it will return the single formula annotation with the smallest absolute error which also fits the Earth filter.
        Otherwise, it will return all Earth-filter compliant formulas.

        Parameters
        ----------
        lowest_error : bool, optional.
            Return only the lowest error formula which also fits the Earth filter.
            If False, return all Earth-filter compliant formulas. Default is True.

        Returns
        -------
        list
            List of molecular formula objects which fit the Earth filter

        References
        ----------
        1. Nikola Tolic et al., "Formularity: Software for Automated Formula Assignment of Natural and Other Organic Matter from Ultrahigh-Resolution Mass Spectra"
            Anal. Chem. 2017, 89, 23, 12659–12665
            doi: 10.1021/acs.analchem.7b03318
        """

        candidates = list(
            filter(
                lambda mf: mf.get("O") > 0
                and mf.get("N") <= 3
                and mf.get("P") <= 2
                and (3 * mf.get("P")) <= mf.get("O"),
                self.molecular_formulas,
            )
        )
        if len(candidates) > 0:
            if lowest_error:
                return min(candidates, key=lambda m: abs(m.mz_error))
            else:
                return candidates
        else:
            return candidates

    def molecular_formula_water_filter(self, lowest_error=True):
        """Filter molecular formula using the 'Water' filter

        This function applies the Formularity-esque 'Water' filter to possible molecular formula assignments.
        Water Filter:
            O > 0 AND N <= 3 AND S <= 2 AND P <= 2

        If the lowest_error method is also used, it will return the single formula annotation with the smallest absolute error which also fits the Water filter.
        Otherwise, it will return all Water-filter compliant formulas.

        Parameters
        ----------
        lowest_error : bool, optional
            Return only the lowest error formula which also fits the Water filter.
            If False, return all Water-filter compliant formulas. Defaults to 2

        Returns
        -------
        list
            List of molecular formula objects which fit the Water filter

        References
        ----------
        1. Nikola Tolic et al., "Formularity: Software for Automated Formula Assignment of Natural and Other Organic Matter from Ultrahigh-Resolution Mass Spectra"
            Anal. Chem. 2017, 89, 23, 12659–12665
            doi: 10.1021/acs.analchem.7b03318
        """

        candidates = list(
            filter(
                lambda mf: mf.get("O") > 0
                and mf.get("N") <= 3
                and mf.get("S") <= 2
                and mf.get("P") <= 2,
                self.molecular_formulas,
            )
        )
        if len(candidates) > 0:
            if lowest_error:
                return min(candidates, key=lambda m: abs(m.mz_error))
            else:
                return candidates
        else:
            return candidates

    def molecular_formula_air_filter(self, lowest_error=True):
        """Filter molecular formula using the 'Air' filter

        This function applies the Formularity-esque 'Air' filter to possible molecular formula assignments.
        Air Filter:
            O > 0 AND N <= 3 AND S <= 1 AND P = 0 AND 3(S+N) <= O

        If the lowest_error method is also used, it will return the single formula annotation with the smallest absolute error which also fits the Air filter.
        Otherwise, it will return all Air-filter compliant formulas.

        Parameters
        ----------
        lowest_error : bool, optional
            Return only the lowest error formula which also fits the Air filter.
            If False, return all Air-filter compliant formulas. Defaults to True.

        Returns
        -------
        list
            List of molecular formula objects which fit the Air filter

        References
        ----------
        1. Nikola Tolic et al., "Formularity: Software for Automated Formula Assignment of Natural and Other Organic Matter from Ultrahigh-Resolution Mass Spectra"
            Anal. Chem. 2017, 89, 23, 12659–12665
            doi: 10.1021/acs.analchem.7b03318
        """

        candidates = list(
            filter(
                lambda mf: mf.get("O") > 0
                and mf.get("N") <= 2
                and mf.get("S") <= 1
                and mf.get("P") == 0
                and 3 * (mf.get("S") + mf.get("N")) <= mf.get("O"),
                self.molecular_formulas,
            )
        )

        if len(candidates) > 0:
            if lowest_error:
                return min(candidates, key=lambda m: abs(m.mz_error))
            else:
                return candidates
        else:
            return candidates

    def cia_score_S_P_error(self):
        """Compound Identification Algorithm SP Error - Assignment Filter

        This function applies the Compound Identification Algorithm (CIA) SP Error filter to possible molecular formula assignments.

        It takes the molecular formula with the lowest S+P count, and returns the formula with the lowest absolute error from this subset.

        Returns
        -------
        MolecularFormula
            A single molecular formula which fits the rules of the CIA SP Error filter


        References
        ----------
        1. Elizabeth B. Kujawinski and Mark D. Behn, "Automated Analysis of Electrospray Ionization Fourier Transform Ion Cyclotron Resonance Mass Spectra of Natural Organic Matter"
            Anal. Chem. 2006, 78, 13, 4363–4373
            doi: 10.1021/ac0600306
        """
        # case EFormulaScore.HAcap:

        lowest_S_P_mf = min(
            self.molecular_formulas, key=lambda mf: mf.get("S") + mf.get("P")
        )
        lowest_S_P_count = lowest_S_P_mf.get("S") + lowest_S_P_mf.get("P")

        list_same_s_p = list(
            filter(
                lambda mf: mf.get("S") + mf.get("P") == lowest_S_P_count,
                self.molecular_formulas,
            )
        )

        # check if list is not empty
        if list_same_s_p:
            return min(list_same_s_p, key=lambda m: abs(m.mz_error))

        else:
            return lowest_S_P_mf

    def cia_score_N_S_P_error(self):
        """Compound Identification Algorithm NSP Error - Assignment Filter

        This function applies the Compound Identification Algorithm (CIA) NSP Error filter to possible molecular formula assignments.

        It takes the molecular formula with the lowest N+S+P count, and returns the formula with the lowest absolute error from this subset.

        Returns
        -------
        MolecularFormula
            A single molecular formula which fits the rules of the CIA NSP Error filter

        References
        ----------
        1. Elizabeth B. Kujawinski and Mark D. Behn, "Automated Analysis of Electrospray Ionization Fourier Transform Ion Cyclotron Resonance Mass Spectra of Natural Organic Matter"
            Anal. Chem. 2006, 78, 13, 4363–4373
            doi: 10.1021/ac0600306

        Raises
        -------
        Exception
            If no molecular formula are associated with mass spectrum peak.
        """
        # case EFormulaScore.HAcap:
        if self.molecular_formulas:
            lowest_N_S_P_mf = min(
                self.molecular_formulas,
                key=lambda mf: mf.get("N") + mf.get("S") + mf.get("P"),
            )
            lowest_N_S_P_count = (
                lowest_N_S_P_mf.get("N")
                + lowest_N_S_P_mf.get("S")
                + lowest_N_S_P_mf.get("P")
            )

            list_same_N_S_P = list(
                filter(
                    lambda mf: mf.get("N") + mf.get("S") + mf.get("P")
                    == lowest_N_S_P_count,
                    self.molecular_formulas,
                )
            )

            if list_same_N_S_P:
                SP_filtered_list = list(
                    filter(
                        lambda mf: (mf.get("S") <= 3) and (mf.get("P") <= 1),
                        list_same_N_S_P,
                    )
                )

                if SP_filtered_list:
                    return min(SP_filtered_list, key=lambda m: abs(m.mz_error))

                else:
                    return min(list_same_N_S_P, key=lambda m: abs(m.mz_error))

            else:
                return lowest_N_S_P_mf
        else:
            raise Exception(
                "No molecular formula associated with the mass spectrum peak at m/z: %.6f"
                % self.mz_exp
            )
