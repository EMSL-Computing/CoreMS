import dataclasses

from corems.encapsulation.factory.processingSetting  import LiquidChromatographSetting, MolecularFormulaSearchSettings, TransientSetting, MassSpecPeakSetting, MassSpectrumSetting
from corems.encapsulation.factory.processingSetting  import CompoundSearchSettings, GasChromatographSetting
from corems.encapsulation.factory.processingSetting  import DataInputSetting

def reset_ms_parameters():
    """Reset the MSParameter class to the default values"""
    MSParameters.molecular_search = MolecularFormulaSearchSettings()
    MSParameters.transient = TransientSetting()
    MSParameters.mass_spectrum = MassSpectrumSetting()
    MSParameters.ms_peak = MassSpecPeakSetting()
    MSParameters.data_input = DataInputSetting()

def reset_gcms_parameters():
    """Reset the GCMSParameters class to the default values"""
    GCMSParameters.molecular_search = CompoundSearchSettings()
    GCMSParameters.gc_ms = GasChromatographSetting()

def reset_lcms_parameters():
    """Reset the LCMSParameters class to the default values"""
    LCMSParameters.lc_ms = LiquidChromatographSetting()
    LCMSParameters.mass_spectrum = MassSpectrumSetting()
    LCMSParameters.ms_peak = MassSpecPeakSetting()
    LCMSParameters.ms1_molecular_search = MolecularFormulaSearchSettings()
    LCMSParameters.ms2_molecular_search = MolecularFormulaSearchSettings()

class MSParameters:
    """MSParameters class is used to store the parameters used for the processing of the mass spectrum
    
    Each attibute is a class that contains the parameters for the processing of the mass spectrum, see the corems.encapsulation.factory.processingSetting module for more details.

    Parameters
    ----------
    use_defaults: bool, optional
        if True, the class will be instantiated with the default values, otherwise the current values will be used. Default is False.

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

    Notes
    -----
    One can use the use_defaults parameter to reset the parameters to the default values.
    Alternatively, to use the current values - modify the class's contents before instantiating the class.
    """

    molecular_search = MolecularFormulaSearchSettings()
    transient = TransientSetting()
    mass_spectrum = MassSpectrumSetting()
    ms_peak = MassSpecPeakSetting()
    data_input = DataInputSetting()

    def __init__(self, use_defaults = False) -> None:
        if not use_defaults:
            self.molecular_search = dataclasses.replace(MSParameters.molecular_search)
            self.transient = dataclasses.replace(MSParameters.transient)
            self.mass_spectrum = dataclasses.replace(MSParameters.mass_spectrum)
            self.ms_peak = dataclasses.replace(MSParameters.ms_peak)
            self.data_input = dataclasses.replace(MSParameters.data_input)
        else:
            self.molecular_search = MolecularFormulaSearchSettings()
            self.transient = TransientSetting()
            self.mass_spectrum = MassSpectrumSetting()
            self.ms_peak = MassSpecPeakSetting()
            self.data_input = DataInputSetting()

class GCMSParameters:
    """GCMSParameters class is used to store the parameters used for the processing of the gas chromatograph mass spectrum

    Each attibute is a class that contains the parameters for the processing of the data, see the corems.encapsulation.factory.processingSetting module for more details.

    Parameters
    ----------
    use_defaults: bool, optional
        if True, the class will be instantiated with the default values, otherwise the current values will be used. Default is False.

    Attributes
    -----------
    molecular_search: MolecularFormulaSearchSettings
        MolecularFormulaSearchSettings object
    gc_ms: GasChromatographSetting
        GasChromatographSetting object

    Notes
    -----
    One can use the use_defaults parameter to reset the parameters to the default values.
    Alternatively, to use the current values - modify the class's contents before instantiating the class.
    """

    molecular_search = CompoundSearchSettings()
    gc_ms = GasChromatographSetting()

    def __init__(self, use_defaults = False) -> None:
        if not use_defaults:
            self.molecular_search = dataclasses.replace(GCMSParameters.molecular_search)
            self.gc_ms = dataclasses.replace(GCMSParameters.gc_ms)
        else:
            self.molecular_search = CompoundSearchSettings()
            self.gc_ms = GasChromatographSetting()

class LCMSParameters:
    """LCMSParameters class is used to store the parameters used for the processing of the liquid chromatograph mass spectrum

    Each attibute is a class that contains the parameters for the processing of the data, see the corems.encapsulation.factory.processingSetting module for more details.

    Parameters
    ----------
    use_defaults: bool, optional
        if True, the class will be instantiated with the default values, otherwise the current values will be used. Default is False.

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

    Notes
    -----
    One can use the use_defaults parameter to reset the parameters to the default values.
    Alternatively, to use the current values - modify the class's contents before instantiating the class.
    """
    lc_ms = LiquidChromatographSetting()
    mass_spectrum = MassSpectrumSetting()
    ms_peak = MassSpecPeakSetting()
    ms1_molecular_search = MolecularFormulaSearchSettings()
    ms2_molecular_search = MolecularFormulaSearchSettings()

    def __init__(self, use_defaults = False) -> None:
        if not use_defaults:
            self.lc_ms = dataclasses.replace(LCMSParameters.lc_ms)
            self.mass_spectrum = dataclasses.replace(LCMSParameters.mass_spectrum)
            self.ms_peak = dataclasses.replace(LCMSParameters.ms_peak)
            self.ms1_molecular_search = dataclasses.replace(LCMSParameters.ms1_molecular_search)
            self.ms2_molecular_search = dataclasses.replace(LCMSParameters.ms2_molecular_search)
        else:
            self.lc_ms = LiquidChromatographSetting()
            self.mass_spectrum = MassSpectrumSetting()
            self.ms_peak = MassSpecPeakSetting()
            self.ms1_molecular_search = MolecularFormulaSearchSettings()
            self.ms2_molecular_search = MolecularFormulaSearchSettings()

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
