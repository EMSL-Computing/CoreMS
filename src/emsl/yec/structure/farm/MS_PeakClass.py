'''
Created on Jun 12, 2019
'''

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

class MSPeak(object):
    
    '''
    classdocs
    '''

    def __init__(self, f_mz, i_charge, f_mgt, f_rp, f_s2n):
        
        #needed to create the object
        self.mz = f_mz
        self.mass = f_mz/i_charge
        self.magnitude = f_mgt
        self.resolving_power = f_rp 
        self.signal_to_noise = f_s2n 
        
        #updated after molecular formula ID
        self.is_isotopologue = False    
        self.molecular_formula = {}
        self.assigment_score = 0
        self.isotopologue_indexes = []
        self.heteroatom_class_label = None #string
        self.ion_type = None #close or open shell
        
    def set_label_isotopologue(self, is_isotopologue): 
        
        self.is_isotopologue = is_isotopologue   
       
    def set_molecular_formula(self, d_mf):
        
        self.molecular_formula = d_mf
    
    def get_is_assigned(self):
        
        return bool(self.molecular_formula)     
        
    def set_isotopologues_indexes(self, l_isoto_indexes):
        
        self.isotopologue_indexes = l_isoto_indexes
        
            