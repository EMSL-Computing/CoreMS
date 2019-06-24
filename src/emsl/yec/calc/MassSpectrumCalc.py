'''
Created on Jun 14, 2019

@author: eber373
'''
from numpy import power, sqrt, array


class MassSpecCalculations(object):
    
    '''
    classdocs
    '''
    
    def f_to_mz(self):
        
        Aterm, Bterm, Cterm = self.calibration_terms
        #Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        if Cterm == 0:
            
            self.mz_domain = (Aterm/ self.frequency_domain ) + (Bterm / power(self.frequency_domain, 2))
           
        #@will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:
            
            self.mz_domain =  Aterm/(2*(Bterm + self.frequency_domain)) + sqrt(Aterm**2 + 4*Bterm*Cterm + 4*Cterm*self.frequency_domain)/(2*(Bterm + self.frequency_domain))
    
    @property
    def number_average_molecular_weight(self):
        #implement from MassSpectralPeaks objs
        return 
    @property
    def weight_average_molecular_weight(self):
        #implement from MassSpectralPeaks objs
        return 
        