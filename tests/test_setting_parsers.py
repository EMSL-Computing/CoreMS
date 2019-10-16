import sys
sys.path.append(".")

from corems.encapsulation.settings.io import settings_parsers

def test_dump():
      
      settings_parsers.dump_search_settings_json()
  
def test_load():
    
    settings_parsers.load_search_setting_json()

if __name__ == "__main__":
    test_load()
    #test_dump()