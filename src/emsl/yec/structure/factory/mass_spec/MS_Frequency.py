from emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase
from numpy import flip

__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"

class MassSpecfromFreq(MassSpecBase):
       
    def __init__(self,frequency_domain, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        super().__init__(None, flip(magnitude), d_params)
        
        self.frequency_domain = frequency_domain
        
        self._set_mz_domain()
        
        ''' use this call to automatically process data as the object is created, Setting need to be changed before initiating the class to be in effect'''
        #self.process_mass_spec()
        
        """implement here, code in MassSpecCalc"""
        #self.calc_noise_threshould
        #self.peak_picking
        #self.calc_resolving_power  
   
    def _set_mz_domain(self):
            
        self._exp_mz = flip(self._f_to_mz())
        
    def cal_noise_treshould(self, auto=True):
        
        self._baselise_noise, self._baselise_noise_std  = self.run_noise_threshould_calc(auto)
        
    def process_mass_spec(self):
        
        self.cal_noise_treshould()
        
        self.find_peaks()    