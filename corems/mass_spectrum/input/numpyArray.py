
__author__ = "Yuri E. Corilo"
__date__ = "Oct 23, 2019"

from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.encapsulation.settings.input import InputParameters
from corems.encapsulation.constant import Labels
from pandas import DataFrame

def ms_from_array_profile(mz, abundance,  dataname, polarity=-1,  auto_process=True, auto_noise=False):

    data_dict = {'m/z': mz, 'Abundance': abundance}
    
    df = DataFrame(data_dict)
    
    output_parameters = get_output_parameters(polarity, dataname)
        
    return MassSpecProfile(df, output_parameters, auto_process=auto_process, auto_noise=auto_noise)

def ms_from_array_centroid(mz, abundance, rp, s2n, dataname, polarity=-1, auto_process=True):

    data_dict = {'m/z': mz, 'Abundance': abundance, 'S/N' : s2n, 'Resolving Power' : rp}
    
    df = DataFrame(data_dict)

    output_parameters = get_output_parameters(polarity, dataname)
        
    return MassSpecCentroid(df, output_parameters, auto_process=auto_process)
    
def get_output_parameters(polarity, file_location):
        
        d_parms = InputParameters.d_parms(file_location)
        
        d_parms["polarity"] = polarity
        
        d_parms["filename_path"] = file_location
        
        d_parms["mobility_scan"] = 0
        
        d_parms["mobility_rt"] = 0
        
        d_parms["scan_number"] = 0
        
        d_parms["rt"] = 0

        d_parms['label'] = Labels.simulated_profile
        
        return d_parms
