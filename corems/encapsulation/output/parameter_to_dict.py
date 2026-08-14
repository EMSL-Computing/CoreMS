from corems.encapsulation.factory.parameters import (
    MSParameters,
    GCMSParameters,
    LCMSParameters
    )
from corems.encapsulation.factory.processingSetting import settings_to_export_dict


def get_dict_all_default_data():
    """Return a dictionary with all default parameters for MS and GCMS"""
    ms_params = MSParameters(use_defaults=True)
    gcms_params = GCMSParameters(use_defaults=True)

    return {
        "MolecularFormulaSearch": settings_to_export_dict(ms_params.molecular_search),
        "Transient": settings_to_export_dict(ms_params.transient),
        "MassSpectrum": settings_to_export_dict(ms_params.mass_spectrum),
        "MassSpecPeak": settings_to_export_dict(ms_params.ms_peak),
        "DataInput": settings_to_export_dict(ms_params.data_input),
        "SpectralSimilaritySearch": settings_to_export_dict(ms_params.spectral_similarity_search),
        "MolecularSearch": settings_to_export_dict(gcms_params.molecular_search),
        "GasChromatograph": settings_to_export_dict(gcms_params.gc_ms),
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
        dictionary with all parameters for LCMSBase object.
        Legacy annotation fields on LiquidChromatograph are omitted;
        MS2 spectral search lives under each mass_spectrum profile's spectral_similarity_search.
    """
    output_dict = {}
    output_dict["LiquidChromatograph"] = settings_to_export_dict(
        lcms_obj.parameters.lc_ms
    )
    output_dict["mass_spectrum"] = {}
    for key, value in lcms_obj.parameters.mass_spectrum.items():
        output_dict["mass_spectrum"][key] = {}
        for k, v in value.__dict__.items():
            output_dict["mass_spectrum"][key][k] = settings_to_export_dict(v)
    return output_dict


def get_dict_lcms_default_data():
    """Return a dictionary with all default parameters for LCMS"""
    default_params = LCMSParameters(use_defaults=True)

    output_dict = {}
    output_dict["LiquidChromatograph"] = settings_to_export_dict(default_params.lc_ms)
    output_dict["mass_spectrum"] = {}
    for key, value in default_params.mass_spectrum.items():
        output_dict["mass_spectrum"][key] = {}
        for k, v in value.__dict__.items():
            output_dict["mass_spectrum"][key][k] = settings_to_export_dict(v)
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
        "MolecularFormulaSearch": settings_to_export_dict(ms_params.molecular_search),
        "Transient": settings_to_export_dict(ms_params.transient),
        "MassSpectrum": settings_to_export_dict(ms_params.mass_spectrum),
        "MassSpecPeak": settings_to_export_dict(ms_params.ms_peak),
        "DataInput": settings_to_export_dict(ms_params.data_input),
        "SpectralSimilaritySearch": settings_to_export_dict(ms_params.spectral_similarity_search),
    }


def get_dict_ms_default_data():
    """Return a dictionary with all default parameters for MS including data input"""
    ms_params = MSParameters(use_defaults=True)

    return {
        "MolecularFormulaSearch": settings_to_export_dict(ms_params.molecular_search),
        "Transient": settings_to_export_dict(ms_params.transient),
        "MassSpectrum": settings_to_export_dict(ms_params.mass_spectrum),
        "MassSpecPeak": settings_to_export_dict(ms_params.ms_peak),
        "DataInput": settings_to_export_dict(ms_params.data_input),
        "SpectralSimilaritySearch": settings_to_export_dict(ms_params.spectral_similarity_search),
    }


def get_dict_gcms_default_data():
    """Return a dictionary with all default parameters for GCMS"""
    default_gcms_params = GCMSParameters(use_defaults=True)

    return {
        "MolecularSearch": settings_to_export_dict(
            default_gcms_params.molecular_search
        ),
        "GasChromatograph": settings_to_export_dict(default_gcms_params.gc_ms),
    }


def get_dict_data_gcms(gcms):
    """Return a dictionary with all parameters for GCMS"""

    return {
        "MolecularSearch": settings_to_export_dict(gcms.molecular_search_settings),
        "GasChromatograph": settings_to_export_dict(gcms.chromatogram_settings),
    }


def get_dict_data_lcms_collection(lcms_collection):
    """Return a dictionary with all parameters for LCMSCollection object

    Parameters
    ----------
    lcms_collection: LCMSCollection
        LCMSCollection object

    Returns
    -------
    dict
        dictionary with all parameters for LCMSCollection object
    """
    output_dict = {}
    output_dict["LCMSCollection"] = settings_to_export_dict(
        lcms_collection.parameters.lcms_collection
    )
    return output_dict


def get_dict_lcms_collection_default_data():
    """Return a dictionary with all default parameters for LCMS Collection"""
    from corems.encapsulation.factory.processingSetting import LCMSCollectionSettings
    
    default_params = LCMSCollectionSettings()

    output_dict = {}
    output_dict["LCMSCollection"] = settings_to_export_dict(default_params)
    return output_dict
