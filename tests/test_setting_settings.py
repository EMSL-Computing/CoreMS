import sys
sys.path.append(".")

from corems.encapsulation.settings.io import settings_parsers
from corems.encapsulation.settings.processingSetting import MolecularLookupDictSettings

def test_json():
      
      settings_parsers.dump_search_settings_json()
      settings_parsers.load_search_setting_json()

def test_data():
    
    settings_parsers.get_dict_data()

def test_settings_search():

    test = MolecularLookupDictSettings()
    test.usedAtoms['C'] = (0,0)
    test.db_directory = 'test'

if __name__ == "__main__":
    
    test_json()
    test_data()
    test_settings_search()