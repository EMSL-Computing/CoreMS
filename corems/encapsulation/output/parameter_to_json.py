from corems.encapsulation.output import parameter_to_dict 


def dump_all_settings_json(filename='SettingsCoreMS.json', file_path=None):
    
    from pathlib import Path
    import json
    '''Write JSON file into current directory
    '''        
    data_dict_all = parameter_to_dict.get_dict_all_default_data()
    
    if not file_path:
        file_path = Path.cwd() / filename 
    
    with open(file_path, 'w', encoding='utf8', ) as outfile:
            
        import re
        #pretty print 
        output = json.dumps(data_dict_all, sort_keys=False, indent=4, separators=(',', ': '))
        output = re.sub(r'",\s+', '", ', output)
        
        outfile.write(output)

def dump_ms_settings_json(filename='SettingsCoreMS.json', file_path=None):
    
    from pathlib import Path
    import json
    '''Write JSON file into current directory
    '''        
    data_dict = parameter_to_dict.get_dict_ms_default_data()

    if not file_path:
        
        file_path = Path.cwd() / filename 
    
    with open(file_path, 'w', encoding='utf8', ) as outfile:
            
        import re
        #pretty print 
        output = json.dumps(data_dict, sort_keys=False, indent=4, separators=(',', ': '))
        output = re.sub(r'",\s+', '", ', output)
        
        outfile.write(output)

def dump_gcms_settings_json(filename='SettingsCoreMS.json', file_path=None):
    '''Write JSON file into current directory
    '''        
    from pathlib import Path
    import json
    
    data_dict = parameter_to_dict.get_dict_gcms_default_data()

    if not file_path:
        
        file_path = Path.cwd() / filename 
    
    with open(file_path, 'w', encoding='utf8', ) as outfile:
            
        import re
        #pretty print 
        output = json.dumps(data_dict, sort_keys=False, indent=4, separators=(',', ': '))
        output = re.sub(r'",\s+', '", ', output)
        
        outfile.write(output)       