
from pandas import DataFrame
import h5py

from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile
from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.input import InputParameters

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
        
    def get_data_profile(self, mz, abundance, rt, auto_process, auto_noise):

        data_dict = {'m/z': mz, 'Abundance': abundance}
    
        df = DataFrame(data_dict)
        
        output_parameters = self.get_output_parameters(rt)
            
        return MassSpecProfile(df, output_parameters, auto_process=auto_process, auto_noise=auto_noise)
        
    def get_mass_spectrum(self, auto_process=True, auto_noise=False):
        
        h5pydata = h5py.File(self.file_location, 'r')
        
        scans = list(h5pydata.keys())
        
        # only one mass spectrum
        if len(scans) == 1:
            
            booster_data = h5pydata[scans[0]]
            
            if self.isCentroid:
                
                raise NotImplementedError
                #return numpyArray.ms_from_array_centroid(mz, abun, rp, snt, self.file_location.stem, polarity=self.polarity=,  auto_process=True, auto_noise=auto_noise)
                
            else:
                mz = booster_data[0]
                abun = booster_data[1]
                rt = float(scans[0])
                return self.get_data_profile(mz, abun, rt, auto_process, auto_noise)

    def get_output_parameters(self, rt):
        
        d_parms = InputParameters.d_parms(self.file_location)
        
        d_parms["polarity"] = self.polarity
        
        d_parms["filename_path"] = self.file_location
        
        d_parms["mobility_scan"] = 0
        
        d_parms["mobility_rt"] = 0
        
        d_parms["scan_number"] = 0
        
        d_parms["rt"] = rt

        d_parms['label'] = Labels.simulated_profile
        
        return d_parms


    
    
    