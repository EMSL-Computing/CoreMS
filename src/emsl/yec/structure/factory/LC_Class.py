'''
Created on Jun 12, 2019
'''

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


from emsl.yec.calc.LC_Calc import LC_Calculations

class LC(LC_Calculations):
    '''
    classdocs
    '''


    def __init__(self, s_file_name, s_file_location, params):
        '''
        Constructor
        '''
        self.file_name = s_file_name
        self.file_location = s_file_location
        
    def set_total_ion_chromatogram(self):
        
        pass    
    
    