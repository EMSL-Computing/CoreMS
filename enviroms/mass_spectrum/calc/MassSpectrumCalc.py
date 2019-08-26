__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"

from numpy import power, multiply, sqrt, multiply
from enviroms.mass_spectrum.calc.NoiseCalc import NoiseThreshouldCalc
from enviroms.mass_spectrum.calc.PeakPicking import PeakPicking


class MassSpecCalc(PeakPicking, NoiseThreshouldCalc):
    '''
    classdocs
    '''

    def _f_to_mz(self):

        Aterm, Bterm, Cterm = self.Aterm, self.Bterm, self.Cterm
        # Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        
        if Cterm == 0:
            
            if Bterm == 0:
                #uncalibrated data
                mz_domain = Aterm / self.freq_exp 
                
            else:
                
                mz_domain = (Aterm / (self.freq_exp)) + (Bterm / power((self.freq_exp), 2))

        # @will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:

            mz_domain = (Aterm / self.freq_exp) + (Bterm / power(self.freq_exp, 2)) + Cterm

        return mz_domain

    def _f_to_mz_bruker(self):
        ''' burker equations for converting frequency (Hz) to m/z, 
        nOmega aquistion is not yet implemented here'''

        Aterm, Bterm, Cterm = self.Aterm, self.Bterm, self.Cterm
        # Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        print(Aterm, Bterm, Cterm)
        if Cterm == 0:
            
            if Bterm == 0:
                #uncalibrated data
                return Aterm / self.freq_exp 
            
            else:
                #calc2
                return Aterm / (self.freq_exp + Bterm)

        # @will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:
            diff = Aterm * Aterm
            
            #this sign(diff + 4) changes on older aquistion software
            diff = diff + 4 * Cterm * (self.freq_exp - Bterm)
            diff = sqrt(diff)
            diff = -Aterm+diff
            #calc3
            return (2*Cterm)/diff
            return diff/2* (self.freq_exp - Bterm)

    def assign_molecular_formulas(self):
        '''call assigment algorithms here'''
        formula_dict = {'C': 5, 'H': 20, 'O': 4, "IonType": 'open_shell'}

        for mspeak in self.mspeaks:

            mspeak.molecular_formula = formula_dict
   
    def number_average_molecular_weight(self, profile=False):

        # mode is profile or centroid data
        if profile:
            a = multiply(self.mz_exp, self.abundance)
            b = self.abundance
            return a.sum()/b.sum()

        else:

            return sum(self.mz_exp_centroide*self.abundance_centroid)/sum(self.abundance_centroid)
    
    def weight_average_molecular_weight(self, profile=False):

        # implement from MassSpectralPeaks objs

        if profile:
            a = multiply(power(self.mz_exp, 2), self.abundance)
            b = self.mz_exp*self.abundance
            return a.sum() / b.sum()

        else:
            return sum(power(self.mz_exp_centroide, 2)*self.abundance_centroid)/sum(self.mz_exp_centroide*self.abundance_centroid)
