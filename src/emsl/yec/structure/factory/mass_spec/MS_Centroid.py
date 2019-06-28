'''
Created on Jun 27, 2019

@author: eber373
'''
from emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase


class MassSpecCentroid(MassSpecBase):
    
    '''
    classdocs
    '''
    def __init__(self, dataframe, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        '''needs to simulate peak shape and pass as exp_mz and magnitude.'''
        exp_mz = dataframe['m/z'].values
        magnitude = dataframe['Abundance'].values
        super().__init__(exp_mz, magnitude, d_params)
        
        self._set_parameters_objects(d_params)
        self.__process__from__centroid(dataframe)
        
    def __simulate_profile__data__(self):
        
        #needs theoretical resolving power calculation
        #return exp_mz, magnitude
        peakshape = None#Gaussian
        
    def __process__from__centroid(self, dataframe):
        
        # this need to change after mass spec deconvolution
        ion_charge = self.polarity
        for index in range(dataframe['m/z'].size):
            
            exp_mz_centroid = dataframe['m/z'][index]
            intes_centr = dataframe['Abundance'][index]
            s_n = dataframe['S/N'][index]
            peak_resolving_power = dataframe['Resolving Power'][index]
            self.add_mspeak(ion_charge, exp_mz_centroid, intes_centr, peak_resolving_power, s_n, 0)    