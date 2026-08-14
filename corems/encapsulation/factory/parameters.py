import dataclasses

from corems.encapsulation.factory.processingSetting import (
    LiquidChromatographSetting,
    MolecularFormulaSearchSettings,
    TransientSetting,
    MassSpecPeakSetting,
    MassSpectrumSetting,
    LCMSCollectionSettings,
    SpectralSimilaritySearchSettings,
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
    MSParameters.spectral_similarity_search = SpectralSimilaritySearchSettings()


def reset_gcms_parameters():
    """Reset the GCMSParameters class to the default values"""
    GCMSParameters.molecular_search = CompoundSearchSettings()
    GCMSParameters.gc_ms = GasChromatographSetting()


def reset_lcms_parameters():
    """Reset the LCMSParameters class to the default values"""
    reset_ms_parameters()
    LCMSParameters.lc_ms = LiquidChromatographSetting()
    LCMSParameters.mass_spectrum = {
        "ms1": MSParameters(use_defaults=True),
        "ms2": MSParameters(use_defaults=True),
    }


def settings_from_lcms(lcms_obj, profile: str = "ms2") -> SpectralSimilaritySearchSettings:
    """Return spectral similarity search settings for an LCMS object profile key.

    Parameters
    ----------
    lcms_obj
        Object with ``parameters.mass_spectrum`` (e.g. LCMSBase).
    profile : str, optional
        Key in ``parameters.mass_spectrum`` (default ``"ms2"``).

    Returns
    -------
    SpectralSimilaritySearchSettings
        Live reference on that LCMS object (mutations affect only this sample).
    """
    return lcms_obj.parameters.mass_spectrum[profile].spectral_similarity_search


def settings_from_lcms_collection(lcms_collection, profile: str = "ms2") -> SpectralSimilaritySearchSettings:
    """Return spectral similarity settings from the **first** sample in a collection.

    This is a read convenience for library build when all samples already share
    the same settings.  Mutating the returned object updates **only the first
    sample**.  To change settings for every sample (and keep them equal for a
    later re-instantiation of the collection), configure a settings instance
    then call :func:`apply_spectral_similarity_search_to_collection`.

    Parameters
    ----------
    lcms_collection
        Sequence of LCMS objects (e.g. LCMSCollection).
    profile : str, optional
        Key in ``parameters.mass_spectrum`` (default ``"ms2"``).

    Returns
    -------
    SpectralSimilaritySearchSettings
    """
    first = lcms_collection[0]
    return settings_from_lcms(first, profile=profile)


def apply_spectral_similarity_search_to_collection(
    lcms_collection,
    settings: SpectralSimilaritySearchSettings,
    profile: str = "ms2",
) -> None:
    """Copy *settings* onto every LCMS sample's ``mass_spectrum[profile]`` bag.

    Use this after configuring spectral similarity / networking knobs so all
    members of a collection stay aligned (collection workflows assume shared
    processing parameters).

    Parameters
    ----------
    lcms_collection
        Sequence of LCMS objects (e.g. LCMSCollection).
    settings : SpectralSimilaritySearchSettings
        Settings to broadcast.  Each sample receives ``settings.copy()`` so
        later per-sample mutations do not cross-link.
    profile : str, optional
        Key in ``parameters.mass_spectrum`` (default ``"ms2"``).
    """
    for sample in lcms_collection:
        if profile not in sample.parameters.mass_spectrum:
            raise KeyError(
                f"profile={profile!r} not in sample parameters.mass_spectrum "
                f"(keys={list(sample.parameters.mass_spectrum)})"
            )
        sample.parameters.mass_spectrum[
            profile
        ].spectral_similarity_search = settings.copy()


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
    spectral_similarity_search: SpectralSimilaritySearchSettings
        Spectral similarity library search and molecular networking settings
        (not molecular formula search). Used for LCMS MS2 profiles; present
        but typically unused on MS1 mass spectrum parameters.

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
    spectral_similarity_search = SpectralSimilaritySearchSettings()

    def __init__(self, use_defaults=False) -> None:
        if not use_defaults:
            self.molecular_search = dataclasses.replace(MSParameters.molecular_search)
            self.transient = dataclasses.replace(MSParameters.transient)
            self.mass_spectrum = dataclasses.replace(MSParameters.mass_spectrum)
            self.ms_peak = dataclasses.replace(MSParameters.ms_peak)
            self.data_input = dataclasses.replace(MSParameters.data_input)
            self.spectral_similarity_search = MSParameters.spectral_similarity_search.copy()
        else:
            self.molecular_search = MolecularFormulaSearchSettings()
            self.transient = TransientSetting()
            self.mass_spectrum = MassSpectrumSetting()
            self.ms_peak = MassSpecPeakSetting()
            self.data_input = DataInputSetting()
            self.spectral_similarity_search = SpectralSimilaritySearchSettings()

    def copy(self):
        """Create a copy of the MSParameters object"""
        new_ms_parameters = MSParameters()
        new_ms_parameters.molecular_search = dataclasses.replace(self.molecular_search)
        new_ms_parameters.transient = dataclasses.replace(self.transient)
        new_ms_parameters.mass_spectrum = dataclasses.replace(self.mass_spectrum)
        new_ms_parameters.ms_peak = dataclasses.replace(self.ms_peak)
        new_ms_parameters.data_input = dataclasses.replace(self.data_input)
        new_ms_parameters.spectral_similarity_search = self.spectral_similarity_search.copy()

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
        equality_check.append(self.spectral_similarity_search == value.spectral_similarity_search)

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
        (includes ``spectral_similarity_search`` on each bag; use ``mass_spectrum["ms2"].spectral_similarity_search`` by default)

    Notes
    -----
    One can use the use_defaults parameter to reset the parameters to the default values.
    Alternatively, to use the current values - modify the class's contents before instantiating the class.

    Annotation fields on ``lc_ms`` (``ms2_min_fe_score``, etc.) are legacy aliases;
    prefer ``mass_spectrum["ms2"].spectral_similarity_search``. Use
    :meth:`sync_ms2_annotation_from_lc_ms` / :meth:`sync_ms2_annotation_to_lc_ms`.
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
        self.sync_ms2_annotation_from_lc_ms()

    def sync_ms2_annotation_from_lc_ms(self, profile: str = "ms2") -> None:
        """Copy legacy ``lc_ms`` annotation fields into ``mass_spectrum[profile].spectral_similarity_search``."""
        if profile not in self.mass_spectrum:
            return
        dest = self.mass_spectrum[profile].spectral_similarity_search
        dest.ms2_min_fe_score = self.lc_ms.ms2_min_fe_score
        dest.search_as_lipids = self.lc_ms.search_as_lipids
        dest.include_fragment_types = self.lc_ms.include_fragment_types

    def sync_ms2_annotation_to_lc_ms(self, profile: str = "ms2") -> None:
        """Copy ``mass_spectrum[profile].spectral_similarity_search`` annotation fields into legacy ``lc_ms``."""
        if profile not in self.mass_spectrum:
            return
        src = self.mass_spectrum[profile].spectral_similarity_search
        self.lc_ms.ms2_min_fe_score = src.ms2_min_fe_score
        self.lc_ms.search_as_lipids = src.search_as_lipids
        self.lc_ms.include_fragment_types = src.include_fragment_types

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
            equality_check.append(
                self.mass_spectrum[key].spectral_similarity_search
                == value.mass_spectrum[key].spectral_similarity_search
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


class LCMSCollectionParameters:
    """LCMSCollectionParameters class is used to store the parameters used for the processing of the LCMS collection

    Each attribute is a class that contains the parameters for the processing of the LCMS collection, 
    see the corems.encapsulation.factory.processingSetting module for more details.

    Parameters
    ----------
    use_defaults: bool, optional
        if True, the class will be instantiated with the default values, otherwise the current values will be used. 
        Default is False.

    Attributes
    -----------
    lcms_collection: LCMSCollectionSettings
        LCMSCollectionSettings object

    Notes
    -----
    One can use the use_defaults parameter to reset the parameters to the default values.
    Alternatively, to use the current values - modify the class's contents before instantiating the class.
    """
    
    lcms_collection = LCMSCollectionSettings()

    def __init__(self, use_defaults=False) -> None:
        if not use_defaults:
            self.lcms_collection = dataclasses.replace(LCMSCollectionParameters.lcms_collection)
        else:
            self.lcms_collection = LCMSCollectionSettings()

    def copy(self):
        """Create a copy of the LCMSCollectionParameters object"""
        new_lcms_collection_parameters = LCMSCollectionParameters()
        new_lcms_collection_parameters.lcms_collection = dataclasses.replace(self.lcms_collection)
        return new_lcms_collection_parameters

    def __eq__(self, value: object) -> bool:
        # Check that the object is of the same type
        if not isinstance(value, LCMSCollectionParameters):
            return False
        return self.lcms_collection == value.lcms_collection

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
