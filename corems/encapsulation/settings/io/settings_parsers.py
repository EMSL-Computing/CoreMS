import json
from pathlib import Path
from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings
from corems.encapsulation.settings.input.ProcessingSetting import TransientSetting
from corems.encapsulation.settings.input.ProcessingSetting import MassSpectrumSetting
from corems.encapsulation.settings.input.ProcessingSetting import MassSpecPeakSetting
from corems.encapsulation.settings.input.InputSetting import DataInputSetting

def get_dict_data():
    
    moleculaSearchSettings = {}
    for item, value in MoleculaSearchSettings.__dict__.items():
        if not item.startswith('__'):
            moleculaSearchSettings[item] =  value
    
    transientSetting = {}
    for item, value in TransientSetting.__dict__.items():
        if not item.startswith('__'):
            transientSetting[item] =  value
    
    massSpectrumSetting = {}
    for item, value in MassSpectrumSetting.__dict__.items():
        if not item.startswith('__'):
            massSpectrumSetting[item] =  value
    
    massSpecPeakSetting = {}
    for item, value in MassSpecPeakSetting.__dict__.items():
        if not item.startswith('__'):
            massSpecPeakSetting[item] =  value                        
    
    dataInputSetting = {}
    for item, value in DataInputSetting.__dict__.items():
        if not item.startswith('__'):
            dataInputSetting[item] =  value  

    return { "MoleculaSearch": moleculaSearchSettings,
             "Transient": transientSetting,
             "MassSpectrum": massSpectrumSetting,
             "MassSpecPeak": massSpecPeakSetting,
             "DataInput": dataInputSetting,
            }

def set_dict_data(data_loaded):
    
    labels = ["MoleculaSearch", "Transient", "MassSpectrum", "MassSpecPeak", "DataInput"]
    classes = [MoleculaSearchSettings, TransientSetting, MassSpectrumSetting, MassSpecPeakSetting]
    
    label_class = zip(labels, classes)
    
    if data_loaded:
    
        for label, classe in label_class:

            class_data = data_loaded.get(label)
            for item, value in class_data.items():
                setattr(classe, item, value)
        
    else:
        
        Warning("Could not load the settings, using the defaults values")    

def dump_search_settings_json( filename='SettingsCoreMS.json'):
    
    '''Write JSON file into current directory
    '''        
    
    data_dict = get_dict_data()

    file_path = Path.cwd() / filename 
    
    with open(file_path, 'w', encoding='utf8', ) as outfile:
            
        import re
        #pretty print 
        output = json.dumps(data_dict, sort_keys=True, indent=4, separators=(',', ': '))
        output = re.sub(r'",\s+', '", ', output)
        
        outfile.write(output)
        
def load_search_setting_json(settings_path=False):
    
    '''LOAD JSON file from current directory
        
        if setting path:  
            setting_path: PATH 
        else:
            setting_path: False
    '''        
    
    if settings_path:
        
        file_path = Path(settings_path)

    else:
        
        filename='SettingsCoreMS.json'
        file_path = Path.cwd() / filename 

    if Path.exists:  
        
        with open(file_path, 'r', encoding='utf8',) as stream:
            
            stream_lines = [n for n in stream.readlines() if not '//' in n.strip()]
            jdata = ''.join(stream_lines)
            data_loaded = json.loads(jdata)
            set_dict_data(data_loaded)
    else:
        
        raise FileNotFoundError("Could not locate %s", file_path)        
