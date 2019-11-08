
from pandas import DataFrame
import h5py

from corems.encapsulation.settings.io.settings_parsers import set_dict_data
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroid
from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.input import InputSetting

class ReadCoreMSHDF_MassSpectrum(ReadCoremsMasslist):
    
    '''
    The MassSpectra object contains lots of MassSpectrum

    Parameters
    ----------
    arg : str
        The arg is used for ...
    *args
        The variable arguments are used for ...
    **kwargs
        The keyword arguments are used for ...

    Attributes
    ----------
    arg : str
        This is where we store arg,
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.h5pydata = h5py.File(self.file_location, 'r')

        self.scans = list(self.h5pydata.keys())

    #override baseclass  
    def load_settings(self, scan):

        loaded_settings = {}
        loaded_settings['MoleculaSearch'] = self.get_attr_data(scan, 'MoleculaSearchSetting')
        loaded_settings['MassSpecPeak'] = self.get_attr_data(scan, 'MassSpecPeakSetting')
        
        loaded_settings['MassSpectrum'] = self.get_high_level_attr_data('MassSpectrumSetting')
        loaded_settings['Transient'] = self.get_high_level_attr_data('TransientSetting')
        
        set_dict_data(loaded_settings)

    #override baseclass  
    def get_dataframe(self, scan):

        columnsLabels = self.get_attr_data(scan, 'ColumnsLabels')

        corems_table_data = self.h5pydata[self.scans[scan]]

        list_dict = []
        for row in corems_table_data:
            data_dict = {}
            for data_index, data in enumerate(row):
                label = columnsLabels[data_index]    
                data_dict[label] = data
            
            list_dict.append(data_dict)
        
        return DataFrame(list_dict)
    
    
    def get_high_level_attr_data(self, attr_group):
        
        import json

        return json.loads(self.h5pydata.attrs[attr_group])
   
    def get_attr_data(self, scan, attr_group, attr_srt=None):
        
        import json

        if attr_srt:
            
            return json.loads(self.h5pydata[self.scans[scan]].attrs[attr_group])[attr_srt]
        
        else:
             
             return json.loads(self.h5pydata[self.scans[scan]].attrs[attr_group])
   
    #override baseclass  
    def get_output_parameters(self, polarity, scan=0):
        
        d_parms = InputSetting.d_parms(self.file_location)
        
        d_parms["polarity"] = polarity
        
        d_parms["filename_path"] = self.file_location
        
        d_parms["mobility_scan"] = self.get_attr_data( scan, 'MassSpecAttrs', 'mobility_scan')
        
        d_parms["mobility_rt"] = self.get_attr_data( scan, 'MassSpecAttrs', 'mobility_rt')
        
        d_parms["scan_number"] = scan
        
        d_parms["rt"] = self.get_attr_data( scan, 'MassSpecAttrs', 'rt')

        d_parms['label'] = Labels.corems_centroid

        d_parms["Aterm"] = self.get_attr_data( scan, 'MassSpecAttrs','Aterm')

        d_parms["Bterm"] = self.get_attr_data( scan, 'MassSpecAttrs','Bterm')
            
        d_parms["Cterm"] = self.get_attr_data( scan, 'MassSpecAttrs','Cterm')

        return d_parms


    
    
    