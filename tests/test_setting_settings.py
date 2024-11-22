import sys
import os

from corems.encapsulation.output import parameter_to_json, parameter_to_dict
from corems.encapsulation.factory.processingSetting  import MolecularLookupDictSettings

def test_toml():
      
    parameter_to_json.dump_all_settings_toml()
    assert os.path.exists('SettingsCoreMS.toml')
    os.remove('SettingsCoreMS.toml')

    parameter_to_json.dump_gcms_settings_toml()
    assert os.path.exists('SettingsCoreMS.toml')
    os.remove('SettingsCoreMS.toml')

    parameter_to_json.dump_ms_settings_toml()
    assert os.path.exists('SettingsCoreMS.toml')
    os.remove('SettingsCoreMS.toml')

def test_json():
      
    parameter_to_json.dump_all_settings_json()
    assert os.path.exists('SettingsCoreMS.json')
    os.remove('SettingsCoreMS.json')

    parameter_to_json.dump_gcms_settings_json()
    assert os.path.exists('SettingsCoreMS.json')
    os.remove('SettingsCoreMS.json')

    parameter_to_json.dump_ms_settings_json()
    assert os.path.exists('SettingsCoreMS.json')
    os.remove('SettingsCoreMS.json')
   
def test_data():
    
    param_dict = parameter_to_dict.get_dict_ms_default_data()
    assert len(param_dict) > 4
    param_dict = parameter_to_dict.get_dict_gcms_default_data()
    assert  len(param_dict) > 1
    

def test_settings_search():

    test = MolecularLookupDictSettings().__dict__
    assert len(test) > 10
    assert "usedAtoms" in test
    assert "url_database" in test

