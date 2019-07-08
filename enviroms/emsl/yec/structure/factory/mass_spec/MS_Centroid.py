from enviroms.emsl.yec.structure.factory.mass_spec.MS_Base import MassSpecBase

__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"

class MassSpecCentroid(MassSpecBase):
    
    '''
    classdocs
    '''
    def __init__(self, dataframe, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        '''needs to simulate peak shape and pass as exp_mz and magnitude.'''
        exp_mz_centroid = dataframe['m/z'].values
        magnitude_centroid = dataframe['Abundance'].values
        exp_mz, magnitude = self.__simulate_profile__data__(exp_mz_centroid, magnitude_centroid)
       
        #print( exp_mz)
        super().__init__(exp_mz, magnitude, d_params)
        
        self._set_parameters_objects(d_params)
        self.__process__from__centroid(dataframe)
        
    def __simulate_profile__data__(self, exp_mz_centroid, magnitude_centroid):
        
        '''needs theoretical resolving power calculation and define peak shape
        this is a quick fix to be able to plot as lines
        peakshape = #Gaussian'''
        
        x, y = [], []
        for i in range(len(exp_mz_centroid)):
            x.append(exp_mz_centroid[i]- 0.0000001)
            x.append(exp_mz_centroid[i])
            x.append(exp_mz_centroid[i]+ 0.0000001)
            y.append(0)
            y.append(magnitude_centroid[i])
            y.append(0)
        return x, y
         
    def __process__from__centroid(self, dataframe):
        
        # this need to change after mass spec deconvolution
        ion_charge = self.polarity
        for index in range(dataframe['m/z'].size):
            exp_mz_centroid = dataframe['m/z'][index]
            intes_centr = dataframe['Abundance'][index]
            s_n = dataframe['S/N'][index]
            peak_resolving_power = dataframe['Resolving Power'][index]
            self.add_mspeak(ion_charge, exp_mz_centroid, intes_centr, peak_resolving_power, s_n, index)  
            
              