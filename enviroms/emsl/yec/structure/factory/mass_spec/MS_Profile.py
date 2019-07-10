from enviroms.emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase

__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"


class MassSpecProfile(MassSpecBase):
    
    '''
    classd
    '''
    
    def __init__(self, dataframe, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        self.label = 'Profile'
        exp_mz = dataframe['m/z'].values
        abundance = dataframe['Abundance'].values
        super().__init__(exp_mz, abundance, d_params)
        
        self.stn = dataframe["S/N"].values
        self.resolving_power = dataframe['Resolving Power'].values
        
    def process_mass_spec(self):
    
        self.find_peaks() 
        