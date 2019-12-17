
from scipy.stats import norm, cauchy
from numpy import linspace, sqrt, log
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
            mass = mass + Atoms.atomic_masses.get(atom) * dict_base.get(atom)

        kendrick_mass = (int(mass) / mass) * self.mz_exp

        nominal_km = int(kendrick_mass)

        kmd = (nominal_km - kendrick_mass) * 100

        # kmd = (nominal_km - km) * 1
        kdm = round(kmd, 0)

        return kdm, kendrick_mass, nominal_km

    def lorentz_pdf(self, datapoint=10000):

        if self.resolving_power:

            # full width half maximum distance
            self.fwhm = (self.mz_exp / self.resolving_power)#self.resolving_power)

            # stardart deviation
            γ = self.fwhm / 2

            # half width baseline distance
            hw_base_distance = (8 * γ)

            #mz_domain = linspace(self.mz_exp - hw_base_distance,
            #                     self.mz_exp + hw_base_distance, datapoint)
            mz_domain = linspace(self.nominal_mz_exp - 0.1,
                                 self.nominal_mz_exp + 1.1, datapoint)
            
            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            calc_abundance = cauchy.pdf(mz_domain, self.mz_exp, γ)

            return mz_domain, (calc_abundance * self.abundance / max(calc_abundance))
        
        else:
            
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')

    def gaussian_pdf(self, datapoint=10000):

        # check if MSPeak contains the resolving power info
        if self.resolving_power:
            # full width half maximum distance
            self.fwhm = self.mz_exp / self.resolving_power

            # stardart deviation
            s = self.fwhm / (2 * sqrt(2 * log(2)))

            # half width baseline distance
            hw_base_distance = (3.2 * s)

            match_loz_factor = 3

            n_d = hw_base_distance * match_loz_factor

            mz_domain = linspace(
                self.mz_exp - n_d, self.mz_exp + n_d, datapoint)

            # gaussian_pdf = lambda x0, x, s: (1/ math.sqrt(2*math.pi*math.pow(s,2))) * math.exp(-1 * math.pow(x-x0,2) / 2*math.pow(s,2) )
            calc_abundance = norm.pdf(mz_domain, self.mz_exp, s)

            return mz_domain, (calc_abundance * self.abundance / max(calc_abundance))

        else:
            raise LookupError(
                'resolving power is not defined, try to use set_max_resolving_power()')
