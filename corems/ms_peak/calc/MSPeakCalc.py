
from scipy.stats import norm, cauchy
from numpy import linspace, sqrt, log, trapz, pi, log, poly1d, polyfit
from corems.encapsulation.constant import Atoms

__author__ = "Yuri E. Corilo"
__date__ = "Jun 04, 2019"


class MSPeakCalculation(object):

    '''
    classdocs
    '''

    def _calc_kdm(self, dict_base):
        '''dict_base = {"C": 1, "H": 2}
        '''
        
        mass = 0
        for atom in dict_base.keys():
            mass += Atoms.atomic_masses.get(atom) * dict_base.get(atom)

        kendrick_mass = (int(mass) / mass) * self.mz_exp

        nominal_km = int(kendrick_mass)

        kmd = (nominal_km - kendrick_mass) 

        # kmd = (nominal_km - km) * 1
        #kdm = round(kmd,0)

        return kmd, kendrick_mass, nominal_km

    def calc_area(self, dx=1):
        
        if self.final_index > self.start_index:
        
            yy = self.ms_parent.abundance_profile[self.start_index:self.final_index]
            
            return trapz(yy, dx = dx)
        
        else:
            
            return None
            
    def voigt(self, oversample_multiplier=1, delta_rp = 0, mz_overlay=1):
        
        from lmfit import models
        
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
        
        from lmfit import models
        
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

    def get_mz_domain(self, oversample_multiplier, mz_overlay):
        
        import matplotlib.pyplot as plt

        start_index = self.start_index - mz_overlay  if not self.start_index == 0 else 0
        final_index = self.final_index + mz_overlay  if not self.final_index == len(self.ms_parent.mz_exp_profile) else self.final_index

        if oversample_multiplier == 1:

            mz_domain = self.ms_parent.mz_exp_profile[start_index: final_index]
            
        else:
            # we assume a linear correlation for m/z and datapoits 
            # which is only true if the m/z range in narrow (within 1 m/z unit)
            # this is not true for a wide m/z range
                         
            indexes = range(start_index, final_index+1)
            mz = self.ms_parent.mz_exp_profile[indexes]
            pol = poly1d(polyfit(indexes, mz, 1))
            oversampled_indexes = linspace(start_index, final_index, (final_index-start_index) * oversample_multiplier)    
            mz_domain = pol(oversampled_indexes)

        return mz_domain

    def lorentz(self, oversample_multiplier=1, delta_rp = 0, mz_overlay=1):

        from lmfit import models

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
        from lmfit import models

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
