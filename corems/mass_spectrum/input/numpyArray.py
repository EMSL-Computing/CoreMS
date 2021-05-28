
__author__ = "Yuri E. Corilo"
__date__ = "Oct 23, 2019"

from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.encapsulation.factory.parameters import default_parameters
from corems.encapsulation.constant import Labels

def ms_from_array_profile(mz, abundance,  dataname, polarity=-1,  auto_process=True, auto_noise=True, data_type=Labels.simulated_profile):

    data_dict = {Labels.mz: mz, Labels.abundance: abundance}
    
    output_parameters = get_output_parameters(polarity, dataname)
    
    output_parameters[Labels.label] = data_type

    ms = MassSpecProfile(data_dict, output_parameters, auto_process=auto_process, auto_noise=auto_noise)
    
    return ms

def ms_from_array_centroid(mz, abundance, rp, s2n, dataname, polarity=-1, auto_process=True):

    data_dict = {Labels.mz: mz, Labels.abundance: abundance, Labels.s2n : s2n, Labels.rp: rp}
    
    output_parameters = get_output_parameters(polarity, dataname)
        
    return MassSpecCentroid(data_dict, output_parameters)
    
def get_output_parameters(polarity, file_location):
        
        d_params = default_parameters(file_location)
        
        d_params['analyzer'] = 'Generic Simulated'

        d_params['instrument_label'] = 'Generic Simulated'

        d_params["polarity"] = polarity
        
        d_params["filename_path"] = file_location
        
        d_params["mobility_scan"] = 0
        
        d_params["mobility_rt"] = 0
        
        d_params["scan_number"] = 0
        
        d_params["rt"] = 0

        d_params[Labels.label] = Labels.simulated_profile
        
        return d_params
