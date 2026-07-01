import dataclasses

from corems.encapsulation.factory.processingSetting import (
    LiquidChromatographSetting,
    MolecularFormulaSearchSettings,
    TransientSetting,
    MassSpecPeakSetting,
    MassSpectrumSetting,
)
from corems.encapsulation.factory.processingSetting import (
    CompoundSearchSettings,
    GasChromatographSetting,
)
from corems.encapsulation.factory.processingSetting import DataInputSetting

def hush_output():
    """Toggle all the verbose_processing flags to False on the MSParameters, GCMSParameters and LCMSParameters classes"""
    MSParameters.molecular_search.verbose_processing = False
    MSParameters.mass_spectrum.verbose_processing = False
    GCMSParameters.gc_ms.verbose_processing = False
    LCMSParameters.lc_ms.verbose_processing = False

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
    reset_ms_parameters()
    LCMSParameters.lc_ms = LiquidChromatographSetting()


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

    def __init__(self, use_defaults=False) -> None:
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

    def copy(self):
        """Create a copy of the MSParameters object"""
        new_ms_parameters = MSParameters()
        new_ms_parameters.molecular_search = dataclasses.replace(self.molecular_search)
        new_ms_parameters.transient = dataclasses.replace(self.transient)
        new_ms_parameters.mass_spectrum = dataclasses.replace(self.mass_spectrum)
        new_ms_parameters.ms_peak = dataclasses.replace(self.ms_peak)
        new_ms_parameters.data_input = dataclasses.replace(self.data_input)

        return new_ms_parameters

    def print(self):
        """Print the MSParameters object"""
        for k, v in self.__dict__.items():
            print(k, type(v).__name__)

            for k2, v2 in v.__dict__.items():
                print("    {}: {}".format(k2, v2))

    def __eq__(self, value: object) -> bool:
        # Check that the object is of the same type
        if not isinstance(value, MSParameters):
            return False
        equality_check = []
        equality_check.append(self.molecular_search == value.molecular_search)
        equality_check.append(self.transient == value.transient)
        equality_check.append(self.mass_spectrum == value.mass_spectrum)
        equality_check.append(self.ms_peak == value.ms_peak)
        equality_check.append(self.data_input == value.data_input)

        return all(equality_check)


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

    def __init__(self, use_defaults=False) -> None:
        if not use_defaults:
            self.molecular_search = dataclasses.replace(GCMSParameters.molecular_search)
            self.gc_ms = dataclasses.replace(GCMSParameters.gc_ms)
        else:
            self.molecular_search = CompoundSearchSettings()
            self.gc_ms = GasChromatographSetting()

    def copy(self):
        """Create a copy of the GCMSParameters object"""
        new_gcms_parameters = GCMSParameters()
        new_gcms_parameters.molecular_search = dataclasses.replace(
            self.molecular_search
        )
        new_gcms_parameters.gc_ms = dataclasses.replace(self.gc_ms)

        return new_gcms_parameters

    def __eq__(self, value: object) -> bool:
        # Check that the object is of the same type
        if not isinstance(value, GCMSParameters):
            return False
        equality_check = []
        equality_check.append(self.molecular_search == value.molecular_search)
        equality_check.append(self.gc_ms == value.gc_ms)

        return all(equality_check)

    def print(self):
        """Print the GCMSParameters object"""
        for k, v in self.__dict__.items():
            print(k, type(v).__name__)

            for k2, v2 in v.__dict__.items():
                print("    {}: {}".format(k2, v2))


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
    mass_spectrum: dict
        dictionary with the mass spectrum parameters for ms1 and ms2, each value is a MSParameters object

    Notes
    -----
    One can use the use_defaults parameter to reset the parameters to the default values.
    Alternatively, to use the current values - modify the class's contents before instantiating the class.
    """

    lc_ms = LiquidChromatographSetting()
    mass_spectrum = {"ms1": MSParameters(), "ms2": MSParameters()}

    def __init__(self, use_defaults=False) -> None:
        if not use_defaults:
            self.lc_ms = dataclasses.replace(LCMSParameters.lc_ms)
            self.mass_spectrum = {
                "ms1": MSParameters(use_defaults=False),
                "ms2": MSParameters(use_defaults=False),
            }
        else:
            self.lc_ms = LiquidChromatographSetting()
            self.mass_spectrum = {
                "ms1": MSParameters(use_defaults=True),
                "ms2": MSParameters(use_defaults=True),
            }

    def copy(self):
        """Create a copy of the LCMSParameters object"""
        new_lcms_parameters = LCMSParameters()
        new_lcms_parameters.lc_ms = dataclasses.replace(self.lc_ms)
        for key in self.mass_spectrum:
            new_lcms_parameters.mass_spectrum[key] = self.mass_spectrum[key].copy()

        return new_lcms_parameters

    def __eq__(self, value: object) -> bool:
        # Check that the object is of the same type
        if not isinstance(value, LCMSParameters):
            return False
        equality_check = []
        equality_check.append(self.lc_ms == value.lc_ms)

        # Check that the mass_spectrum dictionary has the same keys
        equality_check.append(self.mass_spectrum.keys() == value.mass_spectrum.keys())

        # Check that the values of the mass_spectrum dictionary are equal
        for key in self.mass_spectrum.keys():
            equality_check.append(
                self.mass_spectrum[key].mass_spectrum
                == value.mass_spectrum[key].mass_spectrum
            )
            equality_check.append(
                self.mass_spectrum[key].ms_peak == value.mass_spectrum[key].ms_peak
            )
            equality_check.append(
                self.mass_spectrum[key].molecular_search
                == value.mass_spectrum[key].molecular_search
            )
            equality_check.append(
                self.mass_spectrum[key].transient == value.mass_spectrum[key].transient
            )
            equality_check.append(
                self.mass_spectrum[key].data_input
                == value.mass_spectrum[key].data_input
            )

        return all(equality_check)

    def print(self):
        """Print the LCMSParameters object"""
        # Print the lcms paramters
        for k, v in self.__dict__.items():
            if k == "lc_ms":
                print(k, type(v).__name__)

        for k2, v2 in self.mass_spectrum.items():
            """Print the MSParameters object"""
            for k3, v3 in v2.__dict__.items():
                print("{} - {}: {}".format(k2, k3, type(v3).__name__))

                for k4, v4 in v3.__dict__.items():
                    print("    {}: {}".format(k4, v4))


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

    parameters["analyzer"] = "Unknown"

    parameters["acquisition_time"] = None

    parameters["instrument_label"] = "Unknown"

    parameters["sample_name"] = "Unknown"

    parameters["number_data_points"] = 0

    parameters["polarity"] = "Unknown"

    parameters["filename_path"] = str(file_location)

    """scan_number and rt will be need to lc ms"""

    parameters["mobility_scan"] = 0

    parameters["mobility_rt"] = 0

    parameters["scan_number"] = 0

    parameters["rt"] = 0

    return parameters
