import io, os
from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings, MoleculaLookupDictSettings
from ruamel import yaml
import json

def get_dict_data():
    
    data = {}
    for item, value in MoleculaSearchSettings.__dict__.items():
        if not item.startswith('__'):
            data[item] =  value
    return data

def set_dict_data(data_loaded):
    
    if data_loaded:
        for item, value in data_loaded.items():
           setattr(MoleculaSearchSettings, item, value)
    else:
        Warning("Could not load the settings, using the defaults values")    

def dump_search_settings_yaml(filename='SearchConfig'):
    '''Write YAML file into current directory
    '''  
    data_dict = get_dict_data()
    
    file_path = os.getcwd() + os.path.normcase('/' + filename+'.yml')
    
    with io.open(file_path, 'w', encoding='utf8') as outfile:
        yaml.dump(data_dict, outfile, default_flow_style=False, allow_unicode=True)
    
def load_search_setting_yaml(setting_path=False):
    '''LOAD YAML file from current directory
        
        if setting path:  
            setting_path: PATH 
        else:
            setting_path: False
    ''' 
    
    if setting_path:
        file_path = setting_path

    else:
        filename='searchSettings.yml'
        file_path = os.getcwd() + os.path.normcase('/' + filename)

    with open(file_path, 'r') as stream:
        data_loaded = yaml.safe_load(stream)
        set_dict_data(data_loaded)
    #MoleculaSearchSettings.__dict__ = data_loaded

def dump_search_settings_json( filename='SearchConfig'):
    '''Write JSON file into current directory
    '''        
    
    data_dict = get_dict_data()

    file_path = os.getcwd() + os.path.normcase('/' + filename+'.json')
    
    with io.open(file_path, 'w', encoding='utf8', ) as outfile:
            
        import re
        output = json.dumps(data_dict, sort_keys=True, indent=4, separators=(',', ': '))
        output = re.sub(r'",\s+', '", ', output)
        outfile.write(output)
        
def load_search_setting_json(setting_path=False):
    '''LOAD JSON file from current directory
        
        if setting path:  
            setting_path: PATH 
        else:
            setting_path: False
    '''        
    
    if setting_path:
        file_path = setting_path

    else:
        filename='searchSettings.json'
        file_path = os.getcwd() + os.path.normcase('/' + filename)
        
    with open(file_path, 'r', encoding='utf8',) as stream:
        
        stream_lines = [n for n in stream.readlines() if not n.startswith('    //')]
        jdata = ''.join(stream_lines)
       
        data_loaded = json.loads(jdata)
        set_dict_data(data_loaded)
