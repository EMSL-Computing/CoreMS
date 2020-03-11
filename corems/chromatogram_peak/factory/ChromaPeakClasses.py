

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


from corems.chromatogram_peak.calc.GCPeakCalc import GCPeakCalculation 
from dataclasses import dataclass

@dataclass
class LowResCompoundRef:

    def __init__(self, compounds_dict):
    
        self.name = compounds_dict.get("id")
        self.ri = compounds_dict.get("RI") 
        self.rt = compounds_dict.get("RT") 
        self.casno  = compounds_dict.get("CASNO") 
        self.comment  = compounds_dict.get("COMMENT") 
        self.peaks_count = compounds_dict.get("NUM PEAKS") 
        
        self.mz  = compounds_dict.get('mz') 
        self.abundance  = compounds_dict.get("abundance") 

        self.source_temp_c  = compounds_dict.get("SOURCE TEMP C") 
        self.ev  = compounds_dict.get("EV") 
        self.formula  = compounds_dict.get("FORMULA") 
        self.source = compounds_dict.get("SOURCE") 

        self.similarity_score = None    

class ChromaPeakBase():
    '''
    classdocs
    '''
    def __init__(self, mass_spectrum_obj, start_index, index, final_index ):
        
        self.start_index = start_index
        self.final_index = final_index
        self.index = int(index)

        self.mass_spectrum = mass_spectrum_obj
        'updated individual calculation'
        self.area = None
       
        self.compounds = []
        

    def __len__(self):
        
        return len(self.compounds)
        
    def __getitem__(self, position):
        
        return self.compounds[position]

    def add_compound(self, compounds_dict, similarity):
        #implemented in child class
        pass
    def remove_compound(self, compounds_obj):
        
        self.compounds.remove(compounds_obj)
        
    def clear_compounds(self):
        
        self.compounds = []
       
    @property
    def highest_score_compound(self):
        
        return max(self, key = lambda c: c.similarity_score)

class GCPeak(ChromaPeakBase, GCPeakCalculation):

    def __init__(self, mass_spectrum_obj, indexes):
    
        super().__init__(mass_spectrum_obj, *indexes)

    def add_compound(self, compounds_dict, similarity):
       
       compound_obj = LowResCompoundRef(compounds_dict)
       compound_obj.similarity_score = similarity
       
       self.compounds.append(compound_obj)
       