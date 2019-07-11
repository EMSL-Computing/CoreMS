from pandas import read_csv

from enviroms.emsl.yec.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroid



__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

        
class Read_MassList(object):
    '''
    classdocs
    '''

    def __init__(self, file_location, polarity, delimiter="  ", centroid=True):
        
        '''
        Constructor
        '''
        #(newline="\n")
        
        #change this dict VALUES to match your labels, THE ORDER WON'T MATTER
        self.name_dict = {'m/z':'m/z', 'Res.':'Resolving Power', 'I':'Abundance' , "S/N":"S/N"}
        
        
        self.__expected_columns = ['m/z', 'Abundance', 'S/N', 'Resolving Power']
        
        self.file_location = file_location
        
        self.polarity = polarity
        
        self.delimiter  = delimiter
        
        self.centroid  = centroid
         
        
    
    def __new__(cls,file_location, polarity, delimiter ):
        
        cls.__init__(cls, file_location, polarity, delimiter=delimiter)
        return cls.read_file(cls)
    
    def read_file(self):
        
        #delimiter = "  " or " " or  "," or "\t" etc  
        
        dataframe = read_csv(self.file_location, delimiter=self.delimiter, engine='python')
        
        self.check_columns(self, dataframe.columns)
            
        self.clean_data_frame(self, dataframe)
        
        dataframe.rename(columns=self.name_dict, inplace=True)
 
        output_parameters = self.get_output_parameters(self, self.polarity)
            
        if self.centroid:
            
            return MassSpecCentroid(dataframe, output_parameters)

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
            
            expected_column_name = self.name_dict.get(column_name)
            if expected_column_name not in self.__expected_columns:
                
                #dataframe = dataframe.drop(column_name, axis=1)
                del dataframe[column_name]
        #return dataframe
        
    def check_columns(self, data_frame_columns_names):
        
        inverted_name_dict = {value: key for key, value in self.name_dict.items()}
        
        missing_columns = []
        
        for column_name in self.__expected_columns:
            
            user_column_name = inverted_name_dict.get(column_name)
            
            if user_column_name not in data_frame_columns_names:
                missing_columns.append(column_name) 
        
        if len(missing_columns) > 0:
            
            raise Exception("Please make sure to include the columns %s" % ', '.join(missing_columns))