'''
Created on Jun 14, 2019

@author: eber373
'''


'''comment I'm not yet 100% sure about this structure so you with caution. I might decouple this classes later if 
they need to be reused somewhere else then the the MassSpectrum'''

from numpy import power, sqrt

from enviroms.emsl.yec.calc.mass_spec.NoiseCalc import NoiseThreshouldCalc
from enviroms.emsl.yec.calc.mass_spec.PeakPicking import PeakPicking


class MassSpecCalc(NoiseThreshouldCalc, PeakPicking):
    
    '''
    classdocs
    '''
    def _f_to_mz(self):
        
        Aterm, Bterm, Cterm = self.Aterm, self.Bterm, self.Cterm 
        #Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        if Cterm == 0:
            
            mz_domain = (Aterm/ self.frequency_domain ) + (Bterm / power(self.frequency_domain, 2))
           
        #@will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:
            
            mz_domain =  Aterm/(2*(Bterm + self.frequency_domain)) + sqrt(Aterm**2 + 4*Bterm*Cterm + 4*Cterm*self.frequency_domain)/(2*(Bterm + self.frequency_domain))
        
        return mz_domain 
    
    def assign_molecular_formulas(self):
        
        '''call assigment algorithms here'''
        formula_dict = {'C':5,'H':20, 'O':4, "IonType": 'open_shell' }
        
        for mspeak in self.mspeaks:
            
            mspeak.molecular_formula = formula_dict
        
        
    
    def number_average_molecular_weight(self, profile=False):
        #mode is profile or centroid data 
        if profile: 
            return sum(self.exp_mz*self.abundance)/sum(self.abundance)
        else:
            return sum(self.exp_mz_centroide*self.abundance_centroid)/sum(self.abundance_centroid)
        
    def weight_average_molecular_weight(self, profile=False):
        #implement from MassSpectralPeaks objs
        if profile: 
            
            return  sum(power(self.exp_mz,2)*self.abundance)/sum(self.exp_mz*self.abundance)
        
        else:
            
            return sum(power(self.exp_mz_centroide,2)*self.abundance_centroid)/sum(self.exp_mz_centroide*self.abundance_centroid)
    
        

        