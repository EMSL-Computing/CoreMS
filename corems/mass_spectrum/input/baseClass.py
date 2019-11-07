from pathlib import Path

from pandas import read_csv, read_pickle

from corems.encapsulation.settings.input.InputSetting import DataInputSetting

class MassListBaseClass:
    '''
    # The MassListBaseClass object reads mass list data types and returns the mass spectrum obj

    # Parameters
    ----------
    ## file_location : str
        file_location is full data path
    ## * delimiter: str
        Delimeter to read text based files ("," "\t", " ", "  ", etc)
    ## **data_type : str
        The keyword argument data_type is used to determine what type of file to read:
        pandas(pickle), txt(.csv, txt, asci, etc) or HDF5
    ## **isCentroid : bool
        The keyword argument isCentroid is used to determine the mass spectrum data structure:
        if set to False will assume profile mode and will attempt to peak pick
        
    # Attributes
    ----------
    ## _expected_columns : set
       The file has to have a least the values inside the set.
	
    For label translation see:
		corems.encapsulation.settings.input.InputSetting 
       	or add your labels to the SettingsCoreMS.json file and parse the settings

    '''
    def __init__(self, file_location, delimiter="  ", data_type='txt', isCentroid=True):
        
        self.file_location = Path(file_location)

        if not self.file_location.exists():
            raise FileExistsError("File does not exist: %s" %  file_location)

        #(newline="\n")
        
        self._expected_columns = {"m/z", "Abundance", "S/N", "Resolving Power"}
        
        self.delimiter  = delimiter
        
        self.isCentroid  = isCentroid

        self.data_type = data_type

    def get_dataframe(self):

        if self.data_type == 'txt':

            dataframe = read_csv(self.file_location, delimiter=self.delimiter, engine='python')
            
        elif self.data_type == 'dataframe':

            dataframe = read_pickle(self.file_location)   
        
        else:

            raise TypeError('Data type %s is not supported' % self.data_type)

        return  dataframe 

    def load_settings(self,):

        #this will load the setting from SettingCoreMS.json
        # coreMSHFD5 overrides this function to import the attrs stored in the h5 file
        
        pass


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

