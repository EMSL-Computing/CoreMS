import sys
sys.path.append(".")

from corems.encapsulation.input import parameter_from_json
from corems.encapsulation.output import parameter_to_json, parameter_to_dict
from corems.encapsulation.factory.processingSetting  import MolecularLookupDictSettings

def test_toml():
      
    parameter_to_json.dump_all_settings_toml()
    parameter_to_json.dump_gcms_settings_toml()
    parameter_to_json.dump_ms_settings_toml()

def test_json():
      
    parameter_to_json.dump_all_settings_json()
    parameter_to_json.dump_gcms_settings_json()
    parameter_to_json.dump_ms_settings_json()
   
def test_data():
    
    parameter_to_dict.get_dict_ms_default_data()
    parameter_to_dict.get_dict_gcms_default_data()
    

def test_settings_search():

    test = MolecularLookupDictSettings()
    test.usedAtoms['C'] = (0,0)
    test.url_database = 'test'

if __name__ == "__main__":
    
    test_json()
    test_data()
    test_settings_search()