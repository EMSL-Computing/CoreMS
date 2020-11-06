from corems.encapsulation.factory.processingSetting  import MolecularFormulaSearchSettings, TransientSetting, MassSpecPeakSetting, MassSpectrumSetting
from corems.encapsulation.factory.processingSetting  import CompoundSearchSettings, GasChromatographSetting
from corems.encapsulation.factory.processingSetting  import DataInputSetting

class MSParameters:

    molecular_search =  MolecularFormulaSearchSettings()
    transient =  TransientSetting()
    mass_spectrum = MassSpectrumSetting()
    ms_peak =  MassSpecPeakSetting()
    data_input =  DataInputSetting()
  
class GCMSParameters:

    molecular_search = CompoundSearchSettings()
    gc_ms = GasChromatographSetting()
    

def default_parameters(file_location): #pragma: no cover

        parameters = dict()

        parameters["Aterm"] = 0

        parameters["Bterm"] = 0

        parameters["Cterm"] = 0

        parameters["exc_high_freq"] = 0

        parameters["exc_low_freq"] = 0

        parameters["bandwidth"] = 0

        parameters['analyzer'] = 'Unknown'
        
        parameters['instrument_label'] = 'Unknown' 

        parameters['sample_name'] = 'Unknown'

        parameters["number_data_points"] = 0

        parameters["polarity"] = 'Unknown'

        parameters["filename_path"] = str(file_location)

        """scan_number and rt will be need to lc ms"""

        parameters["mobility_scan"] = 0

        parameters["mobility_rt"] = 0

        parameters["scan_number"] = 0

        parameters["rt"] = 0

        return parameters    