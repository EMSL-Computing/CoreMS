from corems.encapsulation.factory.processingSetting  import LiquidChromatographSetting, MolecularFormulaSearchSettings, TransientSetting, MassSpecPeakSetting, MassSpectrumSetting
from corems.encapsulation.factory.processingSetting  import CompoundSearchSettings, GasChromatographSetting
from corems.encapsulation.factory.processingSetting  import DataInputSetting

class MSParameters:

    molecular_search = MolecularFormulaSearchSettings()
    transient = TransientSetting()
    mass_spectrum = MassSpectrumSetting()
    ms_peak = MassSpecPeakSetting()
    data_input = DataInputSetting()

class GCMSParameters:

    molecular_search = CompoundSearchSettings()
    gc_ms = GasChromatographSetting()

class LCMSParameters:
    
    '''
    enforce_target_ms2: bool
            only perform EIC for target_mz if the m/z was selected as precursor for ms2
    scans: list or tuple
        list of select scan to average or a tuple containing the range to average
    peak_height_max_percent: float
        1-100 % used for baseline detection use 0.1 for second_derivative and 10 for other methods    
    peak_max_prominence_percent: float
        1-100 % used for baseline detection
    peak_height_min_percent: float
        0-100 % used for peak detection
    eic_signal_threshold: 
        0-100 % used for extracted ion chromatogram peak detection
    '''

    lc_ms = LiquidChromatographSetting()
    
    mass_spectrum = MassSpectrumSetting()

    ms_peak = MassSpecPeakSetting()

    ms1_molecular_search = MolecularFormulaSearchSettings()
    
    ms2_molecular_search = MolecularFormulaSearchSettings()

def default_parameters(file_location):  # pragma: no cover

    parameters = dict()

    parameters["Aterm"] = 0

    parameters["Bterm"] = 0

    parameters["Cterm"] = 0

    parameters["exc_high_freq"] = 0

    parameters["exc_low_freq"] = 0

    parameters["bandwidth"] = 0

    parameters['analyzer'] = 'Unknown'

    parameters['aquisition_time'] = None

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
