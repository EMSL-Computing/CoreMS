from corems.encapsulation.factory.processingSetting  import LiquidChromatographSetting, MolecularFormulaSearchSettings, TransientSetting, MassSpecPeakSetting, MassSpectrumSetting
from corems.encapsulation.factory.processingSetting  import CompoundSearchSettings, GasChromatographSetting
from corems.encapsulation.factory.processingSetting  import DataInputSetting

class MSParameters:
    """MSParameters class is used to store the parameters used for the processing of the mass spectrum
    
    Each attibute is a class that contains the parameters for the processing of the mass spectrum, see the corems.encapsulation.factory.processingSetting module for more details.

    Attributes
    -----------
    molecular_search: MolecularFormulaSearchSettings
        MolecularFormulaSearchSettings object
    transient: TransientSetting
        TransientSetting object
    mass_spectrum: MassSpectrumSetting
        MassSpectrumSetting object
    ms_peak: MassSpecPeakSetting
        MassSpecPeakSetting object
    data_input: DataInputSetting
        DataInputSetting object
    """

    molecular_search = MolecularFormulaSearchSettings()
    transient = TransientSetting()
    mass_spectrum = MassSpectrumSetting()
    ms_peak = MassSpecPeakSetting()
    data_input = DataInputSetting()

class GCMSParameters:
    """GCMSParameters class is used to store the parameters used for the processing of the gas chromatograph mass spectrum

    Each attibute is a class that contains the parameters for the processing of the data, see the corems.encapsulation.factory.processingSetting module for more details.

    Attributes
    -----------
    molecular_search: MolecularFormulaSearchSettings
        MolecularFormulaSearchSettings object
    gc_ms: GasChromatographSetting
        GasChromatographSetting object
    """

    molecular_search = CompoundSearchSettings()
    gc_ms = GasChromatographSetting()

class LCMSParameters:
    """LCMSParameters class is used to store the parameters used for the processing of the liquid chromatograph mass spectrum

    Each attibute is a class that contains the parameters for the processing of the data, see the corems.encapsulation.factory.processingSetting module for more details.

    Attributes
    -----------
    lc_ms: LiquidChromatographSetting
        LiquidChromatographSetting object
    mass_spectrum: MassSpectrumSetting
        MassSpectrumSetting object
    ms_peak: MassSpecPeakSetting
        MassSpecPeakSetting object
    ms1_molecular_search: MolecularFormulaSearchSettings
        MolecularFormulaSearchSettings object
    ms2_molecular_search: MolecularFormulaSearchSettings
        MolecularFormulaSearchSettings object
    """
    lc_ms = LiquidChromatographSetting()
    
    mass_spectrum = MassSpectrumSetting()

    ms_peak = MassSpecPeakSetting()

    ms1_molecular_search = MolecularFormulaSearchSettings()
    
    ms2_molecular_search = MolecularFormulaSearchSettings()

def default_parameters(file_location):  # pragma: no cover
    """Generate parameters dictionary with the default parameters for data processing
       To gather parameters from instrument data during the data parsing step, a parameters dictionary with the default parameters needs to be generated.
       This dictionary acts as a placeholder and is later used as an argument for all the class constructor methods during instantiation. 
       The data gathered from the instrument is added to the class properties.

    Parameters
    ----------
    file_location: str
        path to the file

    Returns
    -------
    parameters: dict
        dictionary with the default parameters for data processing    
    """

    parameters = dict()

    parameters["Aterm"] = 0

    parameters["Bterm"] = 0

    parameters["Cterm"] = 0

    parameters["exc_high_freq"] = 0

    parameters["exc_low_freq"] = 0

    parameters["mw_low"] = 0
        
    parameters["mw_high"] = 0

    parameters["qpd_enabled"] = 0

    parameters["bandwidth"] = 0

    parameters['analyzer'] = 'Unknown'

    parameters['acquisition_time'] = None

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
