__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from pandas import read_csv

from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import AssignedMassSpecCentroid, MassSpecProfile, MassSpecCentroid
from corems.encapsulation.constant import Labels

class ReadCoremsMasslist(MassListBaseClass):

    def __init__(self, file_location, delimiter=","):
        
        '''
        Constructor
        '''
        
        self.dataframe = read_csv(file_location, delimiter=delimiter, engine='python')

        if not 'Ion Charge' in self.dataframe:
            raise ValueError("%s it is not a valid CoreMS file" % str(file_location))
        
        polarity = self.dataframe['Ion Charge'].values[0]

        super().__init__(file_location, polarity, delimiter=delimiter,isCentroid=True)
        
        
    def get_mass_spectrum(self, auto_process=True):
        
        #delimiter = "  " or " " or  "," or "\t" etc  
        
        self.check_columns(self.dataframe.columns)
        
        self.dataframe.rename(columns=self.name_dict, inplace=True)
 
        output_parameters = self.get_output_parameters(self.polarity)

        return AssignedMassSpecCentroid(self.dataframe, output_parameters, auto_process=auto_process)


class ReadMassList(MassListBaseClass):
    '''
    The ReadMassList object contains lots of MassSpectrum objs

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

    def get_mass_spectrum(self, auto_process=True):
        
        #delimiter = "  " or " " or  "," or "\t" etc  
        
        dataframe = read_csv(self.file_location, delimiter=self.delimiter, engine='python')
        
        self.check_columns(dataframe.columns)
            
        self.clean_data_frame(dataframe)
        
        dataframe.rename(columns=self.name_dict, inplace=True)
 
        output_parameters = self.get_output_parameters(self.polarity)
            
        if self.isCentroid:
            
            return MassSpecCentroid(dataframe, output_parameters, auto_process=auto_process)

        else:
            
            output_parameters[Labels.label] = Labels.bruker_profile
    
            return MassSpecProfile(dataframe, output_parameters, auto_process=auto_process)
    
    