'''
Created on Jun 27, 2019

@author: Yuri E. Corilo
'''
from emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase


class MassSpecProfile(MassSpecBase):
    
    '''
    classdocs
    '''
    
    def __init__(self, exp_mz, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        super().__init__(exp_mz, magnitude, d_params)
        
        self.calc_noise_threshould
        self.peak_picking
        self.calc_resolving_power  
        
        #for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value)           