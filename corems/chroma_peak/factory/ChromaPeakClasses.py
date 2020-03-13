

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
from dataclasses import dataclass

from numpy import array
from corems.chroma_peak.calc.ChromaPeakCalc import GCPeakCalculation 


@dataclass
class LowResCompoundRef:

    def __init__(self, compounds_dict):
        
        self.name = compounds_dict.get("id")
        self.ri = compounds_dict.get("ri")
        self.rt = compounds_dict.get("rt")
        self.casno  = compounds_dict.get("casno")
        self.comment  = compounds_dict.get("comment")
        self.peaks_count = compounds_dict.get("peaks_count")
        
        self.mz  = compounds_dict.get('mz') 
        self.abundance  = compounds_dict.get("abundance") 

        self.source_temp_c  = compounds_dict.get("source_temp_c") 
        self.ev  = compounds_dict.get("ev") 
        self.formula  = compounds_dict.get("formula") 
        self.source = compounds_dict.get("source") 

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
        self._area = None
       
        self._compounds = []
   
    def __len__(self):
        
        return len(self._compounds)
        
    def __getitem__(self, position):
        
        return self._compounds[position]

    def add_compound(self, compounds_dict, similarity):
        #implemented in child class
        pass

    def remove_compound(self, compounds_obj):
        
        self._compounds.remove(compounds_obj)
        
    def clear_compounds(self):
        
        self._compounds = []

    @property   
    def rt(self):
        return self.mass_spectrum.rt

    @property   
    def tic(self):
        return self.mass_spectrum.tic    

    @property   
    def area(self):
        return self._area

class GCPeak(ChromaPeakBase, GCPeakCalculation):

    def __init__(self, mass_spectrum_obj, indexes):
    
        super().__init__(mass_spectrum_obj, *indexes)
        
        self._ri = None

    def add_compound(self, compounds_dict, similarity):
       
       compound_obj = LowResCompoundRef(compounds_dict)
       compound_obj.similarity_score = similarity
       
       self._compounds.append(compound_obj)

    @property
    def ri(self):

        return self._ri

    @property
    def highest_score_compound(self):
        
        return max(self, key = lambda c: c.similarity_score)

                   