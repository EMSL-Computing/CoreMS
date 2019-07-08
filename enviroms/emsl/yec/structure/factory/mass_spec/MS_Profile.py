from enviroms.emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase

__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"


class MassSpecProfile(MassSpecBase):
    
    '''
    classd
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