
from pandas import DataFrame
import h5py
from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.input import numpyArray

class ReadHDF_Booster(MassListBaseClass):
    
    '''
    The MassSpectra object contains lots of MassSpectrum

    Parameters
    ----------
    arg : str
        The arg is used for ...
    *args
        The variable arguments are used for ...
    **kwargs
        The keyword arguments are used for ...

    Attributes
    ----------
    arg : str
        This is where we store arg,
    '''

    def __init__(self, file_location, polarity, delimiter="  ", isCentroid=False):
        
        super(ReadHDF_Booster, self).__init__(file_location, polarity, isCentroid=False)
        
    def get_mass_spectrum(self, auto_process=True, auto_noise=False):
        
        h5pydata = h5py.File(self.file_location, 'r')
        
        scans = list(h5pydata.keys())
        
        # only one mass spectrum
        if len(scans) == 1:
            
            booster_data = h5pydata[scans[0]]
            
        if self.isCentroid:
            
            raise NotImplementedError
                        
            #mz = booster_data[0]
            #abun = booster_data[1]
            #rp = ?
            #snt = ?
            #return numpyArray.ms_from_array_centroid(mz, abun, rp, snt, self.file_location.stem, polarity=self.polarity=,  auto_process=True, auto_noise=auto_noise)
            
        else:
            mz = booster_data[0]
            abun = booster_data[1]
            return numpyArray.ms_from_array_profile(mz, abun, self.file_location.stem, polarity=self.polarity, auto_process=True, auto_noise=auto_noise)

    


    
    
    