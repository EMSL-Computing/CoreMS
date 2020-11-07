from pathlib import Path
import json

from corems.encapsulation.factory.processingSetting  import MolecularFormulaSearchSettings, TransientSetting
from corems.encapsulation.factory.processingSetting  import MassSpectrumSetting
from corems.encapsulation.factory.processingSetting  import MassSpecPeakSetting
from corems.encapsulation.factory.processingSetting  import GasChromatographSetting
from corems.encapsulation.factory.processingSetting import CompoundSearchSettings, DataInputSetting

def load_and_set_parameters_ms(mass_spec_obj, parameters_path=False):   
    
    if parameters_path:
        
        file_path = Path(parameters_path)

    else:
        
        filename='SettingsCoreMS.json'
        file_path = Path.cwd() / filename 

    if file_path.exists():  

            with open(file_path, 'r', encoding='utf8',) as stream:
                data_loaded = json.load(stream)
                _set_dict_data_ms(data_loaded, mass_spec_obj)
    else:
        
        raise FileNotFoundError("Could not locate %s", file_path)   

def load_and_set_parameters_gcms(gcms_obj, parameters_path=False):   
    
    if parameters_path:
        
        file_path = Path(parameters_path)

    else:
        
        filename='SettingsCoreMS.json'
        file_path = Path.cwd() / filename 

    if file_path.exists():  

            with open(file_path, 'r', encoding='utf8',) as stream:
                data_loaded = json.load(stream)
                _set_dict_data_gcms(data_loaded, gcms_obj)
    else:
        
        raise FileNotFoundError("Could not locate %s", file_path)   
    
def _set_dict_data_gcms(data_loaded, gcms_obj):
    
    classes = [GasChromatographSetting(),
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

def _set_dict_data_ms(data_loaded, mass_spec_obj):
    
    from copy import deepcopy

    classes = [MolecularFormulaSearchSettings(), 
               TransientSetting(),
               MassSpectrumSetting(),
               MassSpecPeakSetting()
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


def load_and_set_parameters_class(parameter_label, instance_parameters_class, parameters_path=False):   
    
    if parameters_path: file_path = Path(parameters_path)

    else: file_path = Path.cwd() / 'SettingsCoreMS.json' 
        
    if file_path.exists():
        
        with open(file_path, 'r', encoding='utf8',) as stream:
            
            data_loaded = json.load(stream)
            parameter_class = _set_dict_data(data_loaded, parameter_label, instance_parameters_class)
            
            return parameter_class
    else:
        
        raise FileNotFoundError("Could not locate %s", file_path)  
    
def _set_dict_data(data_loaded, parameter_label, instance_ParameterClass):
    
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