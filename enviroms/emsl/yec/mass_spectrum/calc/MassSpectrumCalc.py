__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"

from numpy import power, multiply, sqrt
from enviroms.emsl.yec.mass_spectrum.calc.NoiseCalc import NoiseThreshouldCalc
from enviroms.emsl.yec.mass_spectrum.calc.PeakPicking import PeakPicking


class MassSpecCalc(PeakPicking, NoiseThreshouldCalc):
    '''
    classdocs
    '''
    
    def _f_to_mz(self):

        Aterm, Bterm, Cterm = self.Aterm, self.Bterm, self.Cterm
        # Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        if Cterm == 0:

            mz_domain = (Aterm / self.frequency_domain) + \
                (Bterm / power(self.frequency_domain, 2))

        # @will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:

            mz_domain = Aterm/(2*(Bterm + self.frequency_domain)) + sqrt(Aterm**2 + 4*Bterm *
                                                                         Cterm + 4*Cterm*self.frequency_domain)/(2*(Bterm + self.frequency_domain))

        return mz_domain
    
    def assign_molecular_formulas(self):
        '''call assigment algorithms here'''
        formula_dict = {'C': 5, 'H': 20, 'O': 4, "IonType": 'open_shell'}

        for mspeak in self.mspeaks:

            mspeak.molecular_formula = formula_dict

    def number_average_molecular_weight(self, profile=False):

        # mode is profile or centroid data
        if profile:
            a = multiply(self.exp_mz, self.abundance)
            b = self.abundance
            return a.sum()/b.sum()

        else:

            return sum(self.exp_mz_centroide*self.abundance_centroid)/sum(self.abundance_centroid)

    def weight_average_molecular_weight(self, profile=False):

        # implement from MassSpectralPeaks objs

        if profile:
            a = multiply(power(self.exp_mz, 2), self.abundance)
            b = self.exp_mz*self.abundance
            return a.sum() / b.sum()

        else:
            return sum(power(self.exp_mz_centroide, 2)*self.abundance_centroid)/sum(self.exp_mz_centroide*self.abundance_centroid)
