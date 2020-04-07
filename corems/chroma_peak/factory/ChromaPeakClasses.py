

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
from numpy import array
from corems.chroma_peak.calc.ChromaPeakCalc import GCPeakCalculation 
from corems.molecular_id.factory.EI_SQL import LowResCompoundRef



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

    def add_compound(self, compounds_dict, spectral_similarity_score, ri_score=None, similarity_score=None):
       
        compound_obj = LowResCompoundRef(compounds_dict)

        compound_obj.spectral_similarity_score = spectral_similarity_score

        compound_obj.ri_score = ri_score

        compound_obj.similarity_score = similarity_score

        self._compounds.append(compound_obj)

    @property
    def ri(self): return self._ri

    @property
    def highest_ss_compound(self):
        if self:
            return max(self, key = lambda c: c.spectral_similarity_score)
        else:
            None
    @property
    def highest_score_compound(self):
        if self:
            return max(self, key = lambda c: c.similarity_score)
        else: 
            return None
                   