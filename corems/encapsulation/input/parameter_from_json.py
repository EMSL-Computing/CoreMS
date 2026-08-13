from pathlib import Path
import dataclasses
import json
import toml

from corems.encapsulation.factory.parameters import MSParameters, LCMSParameters
from corems.encapsulation.factory.processingSetting import (
    MolecularFormulaSearchSettings,
    TransientSetting,
    SpectralSimilaritySearchSettings,
)
from corems.encapsulation.factory.processingSetting import (
    MassSpectrumSetting,
    DataInputSetting,
)
from corems.encapsulation.factory.processingSetting import MassSpecPeakSetting, GasChromatographSetting, CompoundSearchSettings, LCMSCollectionSettings


def load_and_set_toml_parameters_ms(mass_spec_obj, parameters_path=False):
    """Load parameters from a toml file and set the parameters in the mass_spec_obj

    Parameters
    ----------
    mass_spec_obj : MassSpectrum
        corems MassSpectrum object

    parameters_path : str, optional
        path to the parameters file, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """

    if parameters_path:
        file_path = Path(parameters_path)

    else:
        filename = "SettingsCoreMS.toml"
        file_path = Path.cwd() / filename

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = toml.load(stream)
            _set_dict_data_ms(data_loaded, mass_spec_obj)
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def load_and_set_parameters_ms(mass_spec_obj, parameters_path=False):
    """Load parameters from a json file and set the parameters in the mass_spec_obj

    Parameters
    ----------
    mass_spec_obj : MassSpectrum
        corems MassSpectrum object
    parameters_path : str, optional
        path to the parameters file, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """

    if parameters_path:
        file_path = Path(parameters_path)

    else:
        filename = "SettingsCoreMS.json"
        file_path = Path.cwd() / filename

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = json.load(stream)
            _set_dict_data_ms(data_loaded, mass_spec_obj)
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def load_and_set_toml_parameters_gcms(gcms_obj, parameters_path=False):
    """Load parameters from a toml file and set the parameters in the GCMS object

    Parameters
    ----------
    gcms_obj : GCMSBase
        corems GCMSBase object
    parameters_path : str, optional
        path to the parameters file, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """

    if parameters_path:
        file_path = Path(parameters_path)

    else:
        filename = "SettingsCoreMS.toml"
        file_path = Path.cwd() / filename

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = toml.load(stream)
            _set_dict_data_gcms(data_loaded, gcms_obj)
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def load_and_set_parameters_gcms(gcms_obj, parameters_path=False):
    """Load parameters from a json file and set the parameters in the GCMS object

    Parameters
    ----------
    gcms_obj : GCMSBase
        corems GCMSBase object
    parameters_path : str, optional
        path to the parameters file, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """

    if parameters_path:
        file_path = Path(parameters_path)

    else:
        filename = "SettingsCoreMS.json"
        file_path = Path.cwd() / filename

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = json.load(stream)
            _set_dict_data_gcms(data_loaded, gcms_obj)
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def load_and_set_json_parameters_lcms(lcms_obj, parameters_path=False):
    """Load parameters from a json file and set the parameters in the LCMS object

    Parameters
    ----------
    lcms_obj : LCMSBase
        corems LCMSBase object
    parameters_path : str
        path to the parameters file saved as a .json, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """

    if parameters_path:
        file_path = Path(parameters_path)

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = json.load(stream)
            _set_dict_data_lcms(data_loaded, lcms_obj)
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def load_and_set_toml_parameters_lcms(lcms_obj, parameters_path=False):
    """Load parameters from a toml file and set the parameters in the LCMS object

    Parameters
    ----------
    lcms_obj : LCMSBase
        corems LCMSBase object
    parameters_path : str
        path to the parameters file saved as a .toml, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """

    if parameters_path:
        file_path = Path(parameters_path)

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = toml.load(stream)
            _set_dict_data_lcms(data_loaded, lcms_obj)
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def _set_dict_data_gcms(data_loaded, gcms_obj):
    """Set the parameters in the GCMS object from a dict

    This function is called by load_and_set_parameters_gcms and load_and_set_toml_parameters_gcms and should not be called directly.

    Parameters
    ----------
    data_loaded : dict
        dict with the parameters
    gcms_obj : GCMSBase
        corems GCMSBase object
    """

    classes = [
        GasChromatographSetting(),
        CompoundSearchSettings(),
    ]

    labels = ["GasChromatograph", "MolecularSearch"]

    label_class = zip(labels, classes)

    if data_loaded:
        for label, classe in label_class:
            class_data = data_loaded.get(label)
            # not always we will not all the settings
            # this allow a class data to be none and continue
            # to import the other classes
            if class_data:
                for item, value in class_data.items():
                    setattr(classe, item, value)

    gcms_obj.chromatogram_settings = classes[0]
    gcms_obj.molecular_search_settings = classes[1]


def _set_dict_data_lcms(data_loaded, lcms_obj):
    """Set the parameters on a LCMS object from a dict

    This function is called by load_and_set_parameters_lcms and load_and_set_toml_parameters_lcms and should not be called directly.

    Parameters
    ----------
    data_loaded : dict
        dict with the parameters
    lcms_obj : LCMSBase
        corems LCMSBase object
    """

    # Load the lcms parameters
    default_params = LCMSParameters(use_defaults=True)
    lcms_params = data_loaded.get("LiquidChromatograph")
    if lcms_params is not None:
        for item, value in lcms_params.items():
            # If the original value is a tuple but the new one is a list we need to convert the list to a tuple
            if isinstance(value, list) and isinstance(
                getattr(default_params.lc_ms, item), tuple
            ):
                setattr(lcms_obj.parameters.lc_ms, item, tuple(value))
            else:
                setattr(lcms_obj.parameters.lc_ms, item, value)

    def _apply_settings_dict(param_instance, class_data):
        if class_data is None:
            return param_instance
        writable = {f.name for f in dataclasses.fields(param_instance)}
        for item, value in class_data.items():
            if item not in writable:
                continue
            if item == "usedAtoms":
                for atom, atom_value in value.items():
                    value[atom] = tuple(atom_value)
            if item == "additional_similarities" and isinstance(value, list):
                setattr(param_instance, item, list(value))
            elif item == "similarity_thresholds" and isinstance(value, dict):
                setattr(param_instance, item, dict(value))
            elif isinstance(value, list) and isinstance(
                getattr(param_instance, item), tuple
            ):
                setattr(param_instance, item, tuple(value))
            else:
                setattr(param_instance, item, value)
        return param_instance

    def set_ms_params_by_key(ms_key):
        classes = [
            MassSpectrumSetting,
            MassSpecPeakSetting,
            MolecularFormulaSearchSettings,
            DataInputSetting,
            TransientSetting,
            SpectralSimilaritySearchSettings,
        ]

        labels = [
            "mass_spectrum",
            "ms_peak",
            "molecular_search",
            "data_input",
            "transient",
            "spectral_similarity_search",
        ]

        profile = data_loaded.get("mass_spectrum", {}).get(ms_key) or {}
        label_class = zip(labels, classes)

        for label, classe in label_class:
            class_data = profile.get(label)
            param_instance = classe()
            param_instance = _apply_settings_dict(param_instance, class_data)
            setattr(lcms_obj.parameters.mass_spectrum[ms_key], label, param_instance)

    # Load the mass spectrum parameters
    mass_spectrum_data = data_loaded.get("mass_spectrum") or {}
    for ms_key in mass_spectrum_data.keys():
        lcms_obj.parameters.mass_spectrum[ms_key] = MSParameters()
        set_ms_params_by_key(ms_key)

    # Legacy LiquidChromatograph annotation → default ms2 profile.
    # Only when the nested spectral_similarity_search section is absent; nested is
    # the source of truth and must not push back onto lc_ms (keeps export/import
    # round-trips equal for objects that only set nested fields).
    if lcms_params is not None and "ms2" in lcms_obj.parameters.mass_spectrum:
        profile_data = mass_spectrum_data.get("ms2") or {}
        if "spectral_similarity_search" not in profile_data:
            lcms_obj.parameters.sync_ms2_annotation_from_lc_ms(profile="ms2")


def _set_dict_data_ms(data_loaded, mass_spec_obj):
    """Set the parameters in the MassSpectrum object from a dict

    This function is called by load_and_set_parameters_ms and load_and_set_toml_parameters_ms and should not be called directly.

    Parameters
    ----------
    data_loaded : dict
        dict with the parameters
    mass_spec_obj : MassSpectrum
        corems MassSpectrum object
    """

    from copy import deepcopy

    classes = [
        MolecularFormulaSearchSettings(),
        TransientSetting(),
        MassSpectrumSetting(),
        MassSpecPeakSetting(),
    ]

    labels = ["MolecularFormulaSearch", "Transient", "MassSpectrum", "MassSpecPeak"]

    label_class = zip(labels, classes)

    if data_loaded:
        for label, classe in label_class:
            class_data = data_loaded.get(label)
            # not always we will have all the settings classes
            # this allow a class data to be none and continue
            # to import the other classes
            if class_data:
                for item, value in class_data.items():
                    setattr(classe, item, value)

    mass_spec_obj.molecular_search_settings = classes[0]
    mass_spec_obj.transient_settings = classes[1]
    mass_spec_obj.settings = classes[2]
    mass_spec_obj.mspeaks_settings = classes[3]


def load_and_set_toml_parameters_class(
    parameter_label, instance_parameters_class, parameters_path=False
):
    """Load parameters from a toml file and set the parameters in the instance_parameters_class

    Parameters
    ----------
    parameter_label : str
        label of the parameters in the toml file
    instance_parameters_class : object
        instance of the parameters class
    parameters_path : str, optional
        path to the parameters file, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found

    Returns
    -------
    object
        instance of the parameters class
    """

    if parameters_path:
        file_path = Path(parameters_path)

    else:
        file_path = Path.cwd() / "SettingsCoreMS.toml"

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = toml.load(stream)
            parameter_class = _set_dict_data(
                data_loaded, parameter_label, instance_parameters_class
            )

            return parameter_class
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def load_and_set_parameters_class(
    parameter_label, instance_parameters_class, parameters_path=False
):
    """Load parameters from a json file and set the parameters in the instance_parameters_class

    Parameters
    ----------
    parameter_label : str
        label of the parameters in the json file
    instance_parameters_class : object
        instance of the parameters class
    parameters_path : str, optional
        path to the parameters file, by default False

    Raises
    ------
    FileNotFoundError
        if the file is not found

    Returns
    -------
    object
        instance of the parameters class
    """

    if parameters_path:
        file_path = Path(parameters_path)

    else:
        file_path = Path.cwd() / "SettingsCoreMS.json"

    if file_path.exists():
        with open(
            file_path,
            "r",
            encoding="utf8",
        ) as stream:
            data_loaded = json.load(stream)
            parameter_class = _set_dict_data(
                data_loaded, parameter_label, instance_parameters_class
            )

            return parameter_class
    else:
        raise FileNotFoundError("Could not locate %s", file_path)


def _set_dict_data(data_loaded, parameter_label, instance_ParameterClass):
    """Set the parameters in an instance of a parameter class from a dict

    This function is called by load_and_set_parameters_class and load_and_set_toml_parameters_class and should not be called directly.

    Parameters
    ----------
    data_loaded : dict
        dict with the parameters
    parameter_label : str
        label of the parameters in the json file
    instance_ParameterClass : object
        instance of the parameters class

    Returns
    -------
    object
        instance of the parameters class
    """

    classes = [instance_ParameterClass]

    labels = [parameter_label]

    label_class = zip(labels, classes)

    if data_loaded:
        for label, classe in label_class:
            class_data = data_loaded.get(label)
            # not always we will have all the settings classes
            # this allow a class data to be none and continue
            # to import the other classes
            if class_data:
                for item, value in class_data.items():
                    setattr(classe, item, value)

    return classes[0]


def load_and_set_json_parameters_lcms_collection(lcms_collection, parameters_path):
    """Load parameters from a json file and set the parameters in the LCMS collection object

    Parameters
    ----------
    lcms_collection : LCMSCollection
        corems LCMSCollection object
    parameters_path : str or Path
        path to the parameters file saved as a .json

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """
    file_path = Path(parameters_path)

    if file_path.exists():
        with open(file_path, "r", encoding="utf8") as stream:
            data_loaded = json.load(stream)
            _set_dict_data_lcms_collection(data_loaded, lcms_collection)
    else:
        raise FileNotFoundError(f"Could not locate {file_path}")


def load_and_set_toml_parameters_lcms_collection(lcms_collection, parameters_path):
    """Load parameters from a toml file and set the parameters in the LCMS collection object

    Parameters
    ----------
    lcms_collection : LCMSCollection
        corems LCMSCollection object
    parameters_path : str or Path
        path to the parameters file saved as a .toml

    Raises
    ------
    FileNotFoundError
        if the file is not found
    """
    file_path = Path(parameters_path)

    if file_path.exists():
        with open(file_path, "r", encoding="utf8") as stream:
            data_loaded = toml.load(stream)
            _set_dict_data_lcms_collection(data_loaded, lcms_collection)
    else:
        raise FileNotFoundError(f"Could not locate {file_path}")


def _set_dict_data_lcms_collection(data_loaded, lcms_collection):
    """Set the parameters in the LCMS collection object from a dict

    This function is called by load_and_set_json_parameters_lcms_collection and 
    load_and_set_toml_parameters_lcms_collection and should not be called directly.

    Parameters
    ----------
    data_loaded : dict
        dict with the parameters
    lcms_collection : LCMSCollection
        corems LCMSCollection object
    """
    classes = [LCMSCollectionSettings()]
    labels = ["LCMSCollection"]

    label_class = zip(labels, classes)

    if data_loaded:
        for label, classe in label_class:
            class_data = data_loaded.get(label)
            # not always we will have all the settings
            # this allows a class data to be none and continue
            # to import the other classes
            if class_data:
                for attr, value in class_data.items():
                    if hasattr(classe, attr):
                        setattr(classe, attr, value)
        
        lcms_collection.parameters.lcms_collection = classes[0]
