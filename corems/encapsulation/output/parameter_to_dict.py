from corems.encapsulation.factory.parameters import (
    MSParameters,
    GCMSParameters,
    LCMSParameters,
)


def get_dict_all_default_data():
    """Return a dictionary with all default parameters for MS and GCMS"""
    ms_params = MSParameters(use_defaults=True)
    gcms_params = GCMSParameters(use_defaults=True)

    return {
        "MolecularFormulaSearch": ms_params.molecular_search.__dict__,
        "Transient": ms_params.transient.__dict__,
        "MassSpectrum": ms_params.mass_spectrum.__dict__,
        "MassSpecPeak": ms_params.ms_peak.__dict__,
        "DataInput": ms_params.data_input.__dict__,
        "MolecularSearch": gcms_params.molecular_search.__dict__,
        "GasChromatograph": gcms_params.gc_ms.__dict__,
    }


def get_dict_data_lcms(lcms_obj):
    """Return a dictionary with all parameters for LCMSBase object

    Parameters
    ----------
    lcms_obj: LCMSBase
        LCMSBase object

    Returns
    -------
    dict
        dictionary with all parameters for LCMSBase object
    """
    output_dict = {}
    output_dict["LiquidChromatograph"] = lcms_obj.parameters.lc_ms.__dict__
    output_dict["mass_spectrum"] = {}
    for key, value in lcms_obj.parameters.mass_spectrum.items():
        output_dict["mass_spectrum"][key] = {}
        for k, v in value.__dict__.items():
            output_dict["mass_spectrum"][key][k] = v.__dict__
    return output_dict


def get_dict_lcms_default_data():
    """Return a dictionary with all default parameters for LCMS"""
    default_params = LCMSParameters(use_defaults=True)

    output_dict = {}
    output_dict["LiquidChromatograph"] = default_params.lc_ms.__dict__
    output_dict["mass_spectrum"] = {}
    for key, value in default_params.mass_spectrum.items():
        output_dict["mass_spectrum"][key] = {}
        for k, v in value.__dict__.items():
            output_dict["mass_spectrum"][key][k] = v.__dict__
    return output_dict


def get_dict_data_ms(mass_spec):
    """Return a dictionary with all parameters for MassSpectrum object

    Parameters
    ----------
    mass_spec: MassSpectrum
        MassSpectrum object

    Returns
    -------
    dict
        dictionary with all parameters for MassSpectrum object
    """
    ms_params = mass_spec.parameters
    return {
        "MolecularFormulaSearch": ms_params.molecular_search.__dict__,
        "Transient": ms_params.transient.__dict__,
        "MassSpectrum": ms_params.mass_spectrum.__dict__,
        "MassSpecPeak": ms_params.ms_peak.__dict__,
        "DataInput": ms_params.data_input.__dict__,
    }


def get_dict_ms_default_data():
    """Return a dictionary with all default parameters for MS including data input"""
    ms_params = MSParameters(use_defaults=True)

    return {
        "MolecularFormulaSearch": ms_params.molecular_search.__dict__,
        "Transient": ms_params.transient.__dict__,
        "MassSpectrum": ms_params.mass_spectrum.__dict__,
        "MassSpecPeak": ms_params.ms_peak.__dict__,
        "DataInput": ms_params.data_input.__dict__,
    }


def get_dict_gcms_default_data():
    """Return a dictionary with all default parameters for GCMS"""
    default_gcms_params = GCMSParameters(use_defaults=True)

    return {
        "MolecularSearch": default_gcms_params.molecular_search.__dict__,
        "GasChromatograph": default_gcms_params.gc_ms.__dict__,
    }


def get_dict_data_gcms(gcms):
    """Return a dictionary with all parameters for GCMS"""

    return {
        "MolecularSearch": gcms.molecular_search_settings.__dict__,
        "GasChromatograph": gcms.chromatogram_settings.__dict__,
    }
