__author__ = "Yuri E. Corilo"
__date__ = "Jun 04, 2019"

import warnings


from scipy.stats import norm, cauchy
from numpy import linspace, sqrt, log, trapz, pi, log, poly1d, polyfit,flip, square,exp, nan, ceil, rint, floor
from corems.encapsulation.constant import Atoms
from corems.encapsulation.factory.parameters import MSParameters
from lmfit import models
import pyswarm

class MSPeakCalculation(object):

    '''
    classdocs
    '''

    def _calc_kdm(self, dict_base):
        '''dict_base = {"C": 1, "H": 2}
        '''
        kendrick_rounding_method = MSParameters.ms_peak.kendrick_rounding_method # rounding method can be one of floor, ceil or round

        mass = 0
        for atom in dict_base.keys():
            mass += Atoms.atomic_masses.get(atom) * dict_base.get(atom)

        kendrick_mass = (int(mass) / mass) * self.mz_exp

        if kendrick_rounding_method == 'ceil':

            nominal_km = ceil(kendrick_mass)

        elif kendrick_rounding_method == 'round': 

            nominal_km = rint(kendrick_mass)

        elif kendrick_rounding_method == 'floor':

            nominal_km = floor(kendrick_mass)

        else:
            raise  Exception("%s method was not implemented, please refer to corems.ms_peak.calc.MSPeakCalc Class" % kendrick_rounding_method)

        kmd = (nominal_km - kendrick_mass) 

        # kmd = (nominal_km - km) * 1
        #kdm = round(kmd,0)

        return kmd, kendrick_mass, nominal_km

    def calc_area(self):
        '''
        Calculate the peak area using numpy's trapezoidal fit
        uses provided mz_domain to accurately integrate areas independent of digital resolution
        '''
        if self.final_index > self.start_index:

            yy = self._ms_parent.abundance_profile[self.start_index:self.final_index]
            xx = self._ms_parent.mz_exp_profile[self.start_index:self.final_index]
            xx = flip(xx)
            return float(trapz(yy, xx))

        else:

            warnings.warn("Peak Area Calculation for m/z {} has failed".format(self.mz_exp))
            return nan

    def fit_peak(self,mz_extend=6, delta_rp = 0, model='Gaussian'):
        '''
        Model and fit peak lineshape by defined function - using lmfit module
        Do not oversample/resample/interpolate data points 
        Better to go back to time domain and perform more zero filling
        Models allowed: Gaussian, Lorentz, Voigt
        Returns the calculated mz domain, initial defined abundance profile, and the fit peak results object from lmfit module
        mz_extend here extends the x-axis domain so that we have sufficient points either side of the apex to fit.
        Takes about 10ms per peak
        '''
        start_index = self.start_index - mz_extend  if not self.start_index == 0 else 0
        final_index = self.final_index + mz_extend  if not self.final_index == len(self._ms_parent.mz_exp_profile) else self.final_index

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = (self.mz_exp / (self.resolving_power + delta_rp))

            mz_domain = self._ms_parent.mz_exp_profile[start_index:final_index]
            abundance_domain = self._ms_parent.abundance_profile[start_index:final_index]

            if model=='Gaussian':
                # stardard deviation
                sigma = self.fwhm / (2 * sqrt(2 * log(2)))
                amplitude = (sqrt(2*pi)*sigma) * self.abundance
                model = models.GaussianModel()
                params = model.make_params(center=self.mz_exp, amplitude=amplitude, sigma = sigma)

            elif model=='Lorentz':
                # stardard deviation
                sigma = self.fwhm / 2
                amplitude = sigma* pi * self.abundance
                model = models.LorentzianModel()
                params = model.make_params(center=self.mz_exp, amplitude=amplitude, sigma = sigma)

            elif model=='Voigt':
                # stardard deviation
                sigma = self.fwhm / 3.6013
                amplitude = (sqrt(2*pi)*sigma) * self.abundance
                model = models.VoigtModel()
                params = model.make_params(center=self.mz_exp, amplitude=amplitude, sigma = sigma, gamma = sigma)
            else:
                raise LookupError('model lineshape not known or defined')

            #calc_abundance = model.eval(params=params, x=mz_domain) #Same as initial fit, returned in fit_peak object
            fit_peak = model.fit(abundance_domain,params=params, x=mz_domain)
            return mz_domain, fit_peak

        else:
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')


    def voigt_pso(self,w, r, yoff, width, loc, a):
        """
        From https://github.com/pnnl/nmrfit/blob/master/nmrfit/equations.py
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
        """
        # Lorentzian component
        L = (2 / (pi * width)) * 1 / (1 + ((w - loc) / (0.5 * width))**2)

        # Gaussian component
        G = (2 / width) * sqrt(log(2) / pi) * exp(-((w - loc) / (width / (2 * sqrt(log(2)))))**2)

        # Voigt body
        V = (yoff + a) * (r * L + (1 - r) * G)

        return V


    def objective_pso(self,x, w, u):
        """
        The objective function used to fit supplied data.  Evaluates sum of squared differences
        between the fit and the data.
        Parameters
        ----------
        x : list of floats
            Parameter vector.
        w : ndarray
            Array of frequency data.
        Returns
        -------
        rmse : float
            Root mean square error between the data and fit.
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

    def minimize_pso(self,lower, upper, w, u):
        '''
        Minimization function based on defined parameters
        To Do - allow support to pass swarmsize, maxiter, omega, phip, phig parameters.
        Current parameters take ~2 seconds per peak.
        '''
        xopt, fopt = pyswarm.pso(self.objective_pso, lower, upper, args=(w, u),
                                    swarmsize=1000,
                                    maxiter=5000,
                                    omega=-0.2134,
                                    phip=-0.3344,
                                    phig=2.3259)
        return xopt, fopt

    def fit_peak_pso(self, mz_extend=6,upsample_multiplier=5):
        '''
        Function to fit a Voigt peakshape using particle swarm optimisation (PSO)
        Should return better results than lmfit, but much more computationally expensive
        # To Do - Add ability to pass pso args (i.e. swarm size, maxiter, omega, phig, etc)
        Parameters
        ----------
        mz_extend : int
            extra points left and right of peak definition to include in fitting
        upsample_multiplier : int
            factor to increase x-axis points by for simulation of fitted lineshape function
        Returns
        -------
        xopt : array
            variables describing the voigt function
            G/L ratio, width (fwhm), apex (x-axis), area
            y-axis offset is fixed at 0 
                # To do: fix this. Magnitude mode data through CoreMS/Bruker starts at 0 but is noise centered well above 0.
                # Thermo data is noise reduced by also noise subtracted, so starts at 0
                # Absorption mode/phased data will have positive and negative components and may not be baseline corrected
        fopt : float
            objective score (rmse)
        psfit : array
            recalculated y values based on function and optimised fit
        psfit_hdp : tuple of arrays
            0 - linspace x-axis upsampled grid
            1 - recalculated y values based on function and upsampled x-axis grid
            Does not change results, but aids in visualisation of the 'true' voigt lineshape
        '''
        start_index = self.start_index - mz_extend  if not self.start_index == 0 else 0
        final_index = self.final_index + mz_extend  if not self.final_index == len(self._ms_parent.mz_exp_profile) else self.final_index

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = (self.mz_exp / (self.resolving_power))

            mz_domain = self._ms_parent.mz_exp_profile[start_index:final_index]
            abundance_domain = self._ms_parent.abundance_profile[start_index:final_index]
            lower = [0, self.fwhm*0.8, (self.mz_exp-0.0005), 0]
            upper = [1, self.fwhm*1.2, (self.mz_exp+0.0005), self.abundance/self.signal_to_noise]
            xopt, fopt = self.minimize_pso(lower,upper,mz_domain,abundance_domain)
            
            psfit = self.voigt_pso(mz_domain,xopt[0],0,xopt[1],xopt[2],xopt[3])
            psfit_hdp_x = linspace(min(mz_domain),max(mz_domain),num=len(mz_domain)*upsample_multiplier)
            psfit_hdp = self.voigt_pso(psfit_hdp_x,xopt[0],0,xopt[1],xopt[2],xopt[3])
            return xopt, fopt, psfit, (psfit_hdp_x, psfit_hdp)
        else:
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')

             
    def voigt(self, oversample_multiplier=1, delta_rp = 0, mz_overlay=1):
        '''
        Legacy function for voigt lineshape analysis
        '''
        
        
        if self.resolving_power:

            # full width half maximum distance
            self.fwhm = (self.mz_exp / (self.resolving_power + delta_rp))#self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / 3.6013

            # half width baseline distance
            
            #mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)
            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)    
            
            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            
            #TODO derive amplitude
            amplitude = (sqrt(2*pi)*sigma) * self.abundance

            model = models.VoigtModel()

            params = model.make_params(center=self.mz_exp, amplitude=amplitude, sigma = sigma, gamma = sigma)

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance
        
        else:
            
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')

    def pseudovoigt(self, oversample_multiplier=1, delta_rp = 0, mz_overlay=1, fraction =0.5):
        '''
        Legacy pseudovoigt lineshape function
        '''
        if self.resolving_power:

            # full width half maximum distance
            self.fwhm = (self.mz_exp / (self.resolving_power + delta_rp))#self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / 2

            # half width baseline distance
            
            #mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)
            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)    
            
            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            model = models.PseudoVoigtModel()
            
            # TODO derive amplitude
            gamma = sigma
            
            amplitude = (sqrt(2*pi)*sigma) * self.abundance
            amplitude = (sqrt(pi/log(2)) * (pi*sigma*self.abundance)) /( (pi*(1-gamma)) + (sqrt(pi*log(2)) * gamma) )

            params = model.make_params(center=self.mz_exp, sigma = sigma)

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance
        
        else:
            
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')


    def lorentz(self, oversample_multiplier=1, delta_rp = 0, mz_overlay=1):
        '''
        Legacy lorentz lineshape analysis function
        '''
        if self.resolving_power:

            # full width half maximum distance
            self.fwhm = (self.mz_exp / (self.resolving_power + delta_rp))#self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / 2

            # half width baseline distance
            hw_base_distance = (8 * sigma)

            #mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)
            
            
            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)    
            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            model = models.LorentzianModel()
            
            amplitude = sigma* pi * self.abundance

            params = model.make_params(center=self.mz_exp, amplitude=amplitude, sigma = sigma)

            calc_abundance = model.eval(params=params, x=mz_domain)

            return mz_domain, calc_abundance
        
        else:
            
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')

    def gaussian(self, oversample_multiplier=1, delta_rp = 0, mz_overlay=1):
        '''
        Legacy gaussian lineshape analysis function
        '''

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = (self.mz_exp / (self.resolving_power + delta_rp))#self.resolving_power)

            # stardart deviation
            sigma = self.fwhm / (2 * sqrt(2 * log(2)))

            # half width baseline distance
            #hw_base_distance = (3.2 * s)

            #match_loz_factor = 3

            #n_d = hw_base_distance * match_loz_factor

            #mz_domain = linspace(
            #    self.mz_exp - n_d, self.mz_exp + n_d, datapoint)

            mz_domain = self.get_mz_domain(oversample_multiplier, mz_overlay)    
            
            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            
            #calc_abundance = norm.pdf(mz_domain, self.mz_exp, s)

            model = models.GaussianModel()
            
            amplitude = (sqrt(2*pi)*sigma) * self.abundance

            params = model.make_params(center=self.mz_exp, amplitude=amplitude, sigma = sigma)

            calc_abundance = model.eval(params=params, x=mz_domain)
            
            return mz_domain, calc_abundance 

        else:
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')

    def get_mz_domain(self, oversample_multiplier, mz_overlay):
        '''
        Legacy function to support expanding mz domain for legacy lineshape functions
        '''
        start_index = self.start_index - mz_overlay  if not self.start_index == 0 else 0
        final_index = self.final_index + mz_overlay  if not self.final_index == len(self._ms_parent.mz_exp_profile) else self.final_index

        if oversample_multiplier == 1:

            mz_domain = self._ms_parent.mz_exp_profile[start_index: final_index]
            
        else:
            # we assume a linear correlation for m/z and datapoits 
            # which is only true if the m/z range in narrow (within 1 m/z unit)
            # this is not true for a wide m/z range
                        
            indexes = range(start_index, final_index+1)
            mz = self._ms_parent.mz_exp_profile[indexes]
            pol = poly1d(polyfit(indexes, mz, 1))
            oversampled_indexes = linspace(start_index, final_index, (final_index-start_index) * oversample_multiplier)    
            mz_domain = pol(oversampled_indexes)

        return mz_domain
    
    @property
    def number_possible_assignments(self,):
        
        return len(self.molecular_formulas)

    def molecular_formula_lowest_error(self):
       
       return min(self.molecular_formulas, key=lambda m: abs(m.mz_error))

    def molecular_formula_highest_prob_score(self):
       
       return max(self.molecular_formulas, key=lambda m: abs(m.confidence_score))

    def molecular_formula_earth_filter(self, lowest_error=True):
        
        candidates = list(filter(lambda mf: mf.get("O") > 0 and mf.get("N") <=3 and mf.get("P") <= 2 and (3 * mf.get("P")) <= mf.get("O"), self.molecular_formulas))

        if lowest_error:
            return min(candidates, key=lambda m: abs(m.mz_error))
        else:
            return candidates

    def molecular_formula_water_filter(self, lowest_error=True):
       
        candidates = list(filter(lambda mf: mf.get("O") > 0 and mf.get("N") <=3 and mf.get("S") <=2 and  mf.get("P") <= 2, self.molecular_formulas))

        if lowest_error:
            return min(candidates, key=lambda m: abs(m.mz_error))
        else:
            return candidates
    
    def molecular_formula_air_filter(self, lowest_error=True):
       
        candidates = list(filter(lambda mf: mf.get("O") > 0 and mf.get("N") <=2 and mf.get("S") <=1 and  mf.get("P") == 0 and 3* (mf.get("S") + mf.get("N")) <= mf.get("O"), self.molecular_formulas))
        
        if lowest_error:
            return min(candidates, key=lambda m: abs(m.mz_error))
        else:
            return candidates

    def cia_score_S_P_error(self):
        #case EFormulaScore.HAcap:

        lowest_S_P_mf = min(self.molecular_formulas, key=lambda mf: mf.get('S') + mf.get('P'))
        lowest_S_P_count = lowest_S_P_mf.get("S") + lowest_S_P_mf.get("P")
        
        list_same_s_p = list(filter(lambda mf: mf.get('S') + mf.get('P') == lowest_S_P_count, self.molecular_formulas))

        #check if list is not empty
        if list_same_s_p:
        
            return min(list_same_s_p, key=lambda m: abs(m.mz_error))
        
        else:
        
            return lowest_S_P_mf
    
    def cia_score_N_S_P_error(self):
        #case EFormulaScore.HAcap:
        if self.molecular_formulas:

            lowest_N_S_P_mf = min(self.molecular_formulas, key=lambda mf: mf.get('N') + mf.get('S') + mf.get('P'))
            lowest_N_S_P_count = lowest_N_S_P_mf.get("N") + lowest_N_S_P_mf.get("S") + lowest_N_S_P_mf.get("P")

            list_same_N_S_P = list(filter(lambda mf: mf.get('N') + mf.get('S') + mf.get('P') == lowest_N_S_P_count, self.molecular_formulas))

            if list_same_N_S_P:

                SP_filtered_list =  list(filter(lambda mf: (mf.get("S") <= 3 ) and  (mf.get("P")  <= 1 ), list_same_N_S_P))
                
                if SP_filtered_list:
                    
                    return min(SP_filtered_list, key=lambda m: abs(m.mz_error)) 
                
                else:    
                    
                    return min(list_same_N_S_P, key=lambda m: abs(m.mz_error))            
            
            else:
                
                return lowest_N_S_P_mf 
        else:
            raise Exception("No molecular formula associated with the mass spectrum peak at m/z: %.6f" % self.mz_exp)
