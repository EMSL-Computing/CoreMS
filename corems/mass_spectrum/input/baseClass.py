from pathlib import Path
from corems.encapsulation.settings.input.InputSetting import DataInputSetting
class MassListBaseClass:
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

    def __init__(self, file_location, polarity, delimiter="  ", isCentroid=True):
        
        '''
        Constructor
        '''
        
        self.file_location = Path(file_location)

        if not self.file_location.exists():
            raise FileExistsError("File does not exist: " + file_location)

        #(newline="\n")
        
        self._expected_columns = {"m/z", "Abundance", "S/N", "Resolving Power"}
        
        self.polarity = polarity
        
        self.delimiter  = delimiter
        
        self.isCentroid  = isCentroid

    def get_output_parameters(self, polarity):
        
        output_parameters = dict()
        
        output_parameters["Aterm"] = None
        
        output_parameters["Bterm"] = None
        
        output_parameters["Cterm"] = None
        
        output_parameters["polarity"] = polarity
        
        output_parameters["filename_path"] = self.file_location
        
        '''scan_number and rt will be need to lc ms'''
         
        output_parameters["mobility_scan"] = 0
        
        output_parameters["mobility_rt"] = 0
        
        output_parameters["scan_number"] = 0
        
        output_parameters["rt"] = 0
        
        return output_parameters
    
    def clean_data_frame(self, dataframe):
       
        #print(dataframe.head())    
        
        for column_name in dataframe.columns:
            
            expected_column_name = DataInputSetting.header_translate.get(column_name)
            if expected_column_name not in self._expected_columns:
                
                #dataframe = dataframe.drop(column_name, axis=1)
                del dataframe[column_name]
        #return dataframe
        
    def check_columns(self, header_labels):
        
        #print(header_labels)
        #print( DataInputSetting.header_translate.keys())
        #inverted_name_dict = {value: key for key, value in DataInputSetting.header_translate.items()}
        
        #print(inverted_name_dict)
        
        found_label = set()

        for label in header_labels:
            
            if not label in self._expected_columns:
                
                user_column_name = DataInputSetting.header_translate.get(label)
                
                if user_column_name in self._expected_columns:
                    
                    found_label.add(user_column_name)
    
            else:   
                found_label.add(label)
        
        not_found = self._expected_columns - found_label

        if len(not_found) > 0:
            
            raise Exception("Please make sure to include the columns %s" % ', '.join(not_found))    

