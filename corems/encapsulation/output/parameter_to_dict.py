from corems.encapsulation.factory.parameters import MSParameters, GCMSParameters

def get_dict_all_default_data():
    
    return { "MolecularFormulaSearch": MSParameters.molecular_search.__dict__,
             "Transient": MSParameters.transient.__dict__,
             "MassSpectrum": MSParameters.mass_spectrum.__dict__,
             "MassSpecPeak": MSParameters.ms_peak.__dict__,
             "DataInput": MSParameters.data_input.__dict__,
             "MolecularSearch": GCMSParameters.molecular_search.__dict__,
             "GasChromatograph": GCMSParameters.gc_ms.__dict__,
            }

def get_dict_data_ms(mass_spec):

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
    
    return { "MolecularFormulaSearch": MSParameters.molecular_search.__dict__,
             "Transient": MSParameters.transient.__dict__,
             "MassSpectrum": MSParameters.mass_spectrum.__dict__,
             "MassSpecPeak": MSParameters.ms_peak.__dict__,
             "DataInput": MSParameters.data_input.__dict__,
            }

def get_dict_gcms_default_data():
    
    return { "MolecularSearch": GCMSParameters.molecular_search.__dict__,
             "GasChromatograph": GCMSParameters.gc_ms.__dict__,
            }

def get_dict_data_gcms(gcms):

    return { "MolecularSearch": gcms.molecular_search_settings.__dict__,
             "GasChromatograph":  gcms.chromatogram_settings.__dict__,
            }          