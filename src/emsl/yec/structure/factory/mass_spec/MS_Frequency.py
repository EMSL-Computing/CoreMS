'''
Created on Jun 27, 2019

@author: corilo
'''
from numpy import flip

from emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase


class MassSpecfromFreq(MassSpecBase):
       
    def __init__(self,frequency_domain, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        super().__init__(None, flip(magnitude), d_params)
        
        self.frequency_domain = frequency_domain
        
        self._set_mz_domain()
        
        """implement here, code in MassSpecCalc"""
        #self.calc_noise_threshould
        #self.peak_picking
        #self.calc_resolving_power  
   
   
    def _set_mz_domain(self):
            
        self._exp_mz = flip(self._f_to_mz())
        
    def process_mass_spec(self):
        
        self.set_noise_treshould()
        self.find_peaks()    