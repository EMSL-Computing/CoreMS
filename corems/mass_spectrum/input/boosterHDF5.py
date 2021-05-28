
from io import BytesIO

import h5py
from s3path import S3Path
from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile
from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters

class ReadHDF_BoosterMassSpectrum(MassListBaseClass):
    
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

    def __init__(self, file_location, isCentroid=False):
        
        self.polarity = self.get_polarity(file_location)

        super().__init__(file_location, isCentroid=False)
        
    def get_data_profile(self, mz, abundance, auto_process, auto_noise):

        data_dict = {Labels.mz: mz, Labels.abundance: abundance}
    
        output_parameters = self.get_output_parameters()
            
        return MassSpecProfile(data_dict, output_parameters, auto_process=auto_process, auto_noise=auto_noise)
    
    def get_attr_data(self, scan, attr_srt):

        return self.h5pydata[self.scans[scan]].attrs[attr_srt]

    def get_polarity(self, file_location):

        if isinstance(file_location, S3Path):
            data = BytesIO(file_location.open('rb').read())
        else:
            data = file_location
        
        self.h5pydata = h5py.File(data, 'r')

        self.scans = list(self.h5pydata.keys())
        
        polarity = self.get_attr_data(0,'r_h_polarity')
        
        if polarity == 'negative scan':
            
            return -1
        else:
            return +1    

    def get_mass_spectrum(self, auto_process=True, auto_noise=True):
        
        # only one mass spectrum
        if len(self.scans) == 1:
            
            booster_data = self.h5pydata[self.scans[0]]
            
            if self.isCentroid:
                
                raise NotImplementedError
                #return numpyArray.ms_from_array_centroid(mz, abun, rp, snt, self.file_location.stem, polarity=self.polarity,  auto_process=True, auto_noise=auto_noise)
                
            else:
                
                mz = booster_data[0]
                abun = booster_data[1]
                
                return self.get_data_profile(mz, abun, auto_process, auto_noise)

    def get_output_parameters(self):
        
        d_params = default_parameters(self.file_location)
        
        d_params["polarity"] = self.polarity
        
        d_params["filename_path"] = self.file_location
        
        d_params["mobility_scan"] = 0
        
        d_params["mobility_rt"] = 0
        
        d_params["scan_number"] = 0
        
        d_params["rt"] = self.get_attr_data(0, 'r_h_start_time')

        d_params['label'] = Labels.booster_profile

        d_params["Aterm"] = self.get_attr_data(0, 'r_cparams')[0]

        d_params["Bterm"] = self.get_attr_data(0, 'r_cparams')[1]

        return d_params


    
    
    