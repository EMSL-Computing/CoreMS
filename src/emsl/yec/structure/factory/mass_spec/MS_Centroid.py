'''
Created on Jun 27, 2019

@author: eber373
'''
from emsl.yec.structure.factory.mass_spec.MassSpectrumClass import MassSpecBase

class MassSpecCentroid(MassSpecBase):
    
    '''
    classdocs
    '''
    def __init__(self, exp_mz, magnitude, dataframe, d_params, **kwargs): 
                 
        '''
        Constructor
        dataframe will contain l_base_noise, l_signal_to_noise, l_charge, l_resolving_power
        needs to simulate peak shape and pass as exp_mz and magnitude.
        '''
        super().__init__(exp_mz, magnitude, d_params)
        
        #self.__process__from__centroid(l_base_noise, l_signal_to_noise, l_charge, l_resolving_power)
    def __process__from__centroid(self, l_base_noise, l_signal_to_noise, l_charge, l_resolving_power):
        
        if l_base_noise.any():
            self.set_base_noise(l_base_noise)
        
        if l_signal_to_noise.any():
            self.set_noise(l_signal_to_noise)
        
        if l_charge.any():
            self.set_charge(l_charge)
        
        if l_resolving_power.any():
            self.set_resolving_power(l_resolving_power)  