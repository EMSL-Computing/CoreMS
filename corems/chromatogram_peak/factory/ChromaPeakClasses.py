

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


from corems.chromatogram_peak.calc.GCPeakCalc import GCPeakCalculation 

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
        self.similarity_score = []

    def __len__(self):
        
        return len(self.compounds)
        
    def __setitem__(self, position, compounds_obj):
        
        self.compounds[position] = compounds_obj

    def __getitem__(self, position):
        
        return self.compounds[position]

    def add_compound(self, compounds_obj, similarity):
       
       self.compounds.append(compounds_obj)
       self.similarity_score.append(similarity)

    def remove_compound(self, compounds_obj):
        
        index = self.compounds.index(compounds_obj)
        del self.compounds[index]
        del self.similarity_score[index]

    def clear_compounds(self):
        
        self.compounds = []
        self.similarity_score = []

class GCPeak(ChromaPeakBase, GCPeakCalculation):

    def __init__(self, mass_spectrum_obj, indexes):
    
        super().__init__(mass_spectrum_obj, *indexes)