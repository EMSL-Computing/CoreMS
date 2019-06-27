'''
Created on Jun 14, 2019

@author: eber373
'''


from numpy import power, sqrt
from emsl.yec.calc.mass_spec.NoiseCalc import NoiseThreshouldCalc
'''comment I'm not yet 100% sure about this structure so you with caution. I might decouple this classes later if 
they need to be reused somewhere else then the the MassSpectrum'''
        
class MassSpecCalc(NoiseThreshouldCalc):
    
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
    
    def number_average_molecular_weight(self, profile_data=True):
        #mode is profile or centroid data 
        if profile_data: 
            return sum(self.exp_mz*self.magnitude)/sum(self.magnitude)
        else:
            return sum(self.exp_mz_centroide*self.magnitude_centroid)/sum(self.magnitude)
        
    def weight_average_molecular_weight(self, profile_data=True):
        #implement from MassSpectralPeaks objs
        if profile_data: 
            
            return  sum(power(self.exp_mz,2)*self.magnitude)/sum(self.exp_mz*self.magnitude)
        
        else:
            
            return power(sum(self.exp_mz_centroide*self.magnitude_centroid),2)/sum(self.exp_mz_centroide*self.magnitude_centroid)
    
        

        