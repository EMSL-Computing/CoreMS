from corems.encapsulation.factory.parameters import MSParameters, GCMSParameters, LCMSParameters

def get_dict_all_default_data():
    """ Return a dictionary with all default parameters for MS and GCMS
    
    """
    return { "MolecularFormulaSearch": MSParameters.molecular_search.__dict__,
             "Transient": MSParameters.transient.__dict__,
             "MassSpectrum": MSParameters.mass_spectrum.__dict__,
             "MassSpecPeak": MSParameters.ms_peak.__dict__,
             "DataInput": MSParameters.data_input.__dict__,
             "MolecularSearch": GCMSParameters.molecular_search.__dict__,
             "GasChromatograph": GCMSParameters.gc_ms.__dict__,
            }

def get_dict_data_lcms(lcms_obj):
    """ Return a dictionary with all parameters for LCMS
    
    """
    return { "LiquidChromatograph": lcms_obj.parameters.lc_ms.__dict__,
             "MassSpectrum": lcms_obj.parameters.mass_spectrum.__dict__,
             "MassSpecPeak": lcms_obj.parameters.ms_peak.__dict__,
             "MS1MolecularSearch": lcms_obj.parameters.ms1_molecular_search.__dict__, 
             "MS2MolecularSearch": lcms_obj.parameters.ms2_molecular_search.__dict__,
            }

def get_dict_lcms_default_data():
    """ Return a dictionary with all default parameters for LCMS
    
    """
    
    return { "LiquidChromatograph": LCMSParameters.lc_ms.__dict__,
             "MassSpectrum": LCMSParameters.mass_spectrum.__dict__,
             "MassSpecPeak": LCMSParameters.ms_peak.__dict__,
             "MS1MolecularSearch": LCMSParameters.ms1_molecular_search.__dict__, 
             "MS2MolecularSearch": LCMSParameters.ms2_molecular_search.__dict__,
            }

def get_dict_data_ms(mass_spec):
    """ Return a dictionary with all parameters for MS
    
    """

    if mass_spec._transient_settings:

        return { "MolecularFormulaSearch": mass_spec.molecular_search_settings.__dict__,
                "Transient": mass_spec.transient_settings.__dict__,
                "MassSpectrum": mass_spec.settings.__dict__,
                "MassSpecPeak": mass_spec.mspeaks_settings.__dict__
                }
    else:
        
        return { "MolecularFormulaSearch": mass_spec.molecular_search_settings.__dict__,
                "MassSpectrum": mass_spec.settings.__dict__,
                "MassSpecPeak": mass_spec.mspeaks_settings.__dict__
                }
                
def get_dict_ms_default_data():
    """ Return a dictionary with all default parameters for MS including data input
    
    """
    
    return { "MolecularFormulaSearch": MSParameters.molecular_search.__dict__,
             "Transient": MSParameters.transient.__dict__,
             "MassSpectrum": MSParameters.mass_spectrum.__dict__,
             "MassSpecPeak": MSParameters.ms_peak.__dict__,
             "DataInput": MSParameters.data_input.__dict__,
            }

def get_dict_gcms_default_data():
    """ Return a dictionary with all default parameters for GCMS
    
    """
    
    return { "MolecularSearch": GCMSParameters.molecular_search.__dict__,
             "GasChromatograph": GCMSParameters.gc_ms.__dict__,
            }

def get_dict_data_gcms(gcms):
    """ Return a dictionary with all parameters for GCMS
    
    """

    return { "MolecularSearch": gcms.molecular_search_settings.__dict__,
             "GasChromatograph":  gcms.chromatogram_settings.__dict__,
            }          