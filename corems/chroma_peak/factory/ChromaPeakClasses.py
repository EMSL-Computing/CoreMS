

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
from numpy import array, trapz    
from corems.chroma_peak.calc.ChromaPeakCalc import GCPeakCalculation
import corems.mass_spectra.factory.LC_Class as lcms
from corems.mass_spectra.factory.LC_Temp import EIC_Data
from corems.molecular_id.factory.EI_SQL import LowResCompoundRef


class ChromaPeakBase():
    '''
    classdocs
    '''
    def __init__(self, chromatogram_parent, mass_spectrum_obj, start_index, index, final_index):
        
        self.start_scan = start_index
        self.final_scan = final_index
        self.apex_scan = int(index)
        self.chromatogram_parent = chromatogram_parent
        self.mass_spectrum = mass_spectrum_obj
        'updated individual calculation'
        self._area = None
   
    @property   
    def retention_time(self):
        return self.mass_spectrum.retention_time

    @property   
    def tic(self):
        return self.mass_spectrum.tic    

    @property   
    def area(self):
        return self._area

    @property
    def rt_list(self):
        return [self.chromatogram_parent.retention_time[i] for i in range(self.start_index, self.final_index+1) ]

    @property   
    def tic_list(self):
        return [self.chromatogram_parent.tic[i] for i in range(self.start_index, self.final_index+1) ]

class DataDependentPeak(ChromaPeakBase):

    def __init__(self, chromatogram_parent, mass_spectrum_obj, peak_indexes, eic_data: EIC_Data, possible_molform = None):

        eic_data_apex_index = eic_data.scans.index(peak_indexes[1])
        
        retention_time = eic_data.time[eic_data_apex_index]

        mass_spectrum_obj.retention_time = retention_time

        super().__init__(chromatogram_parent, mass_spectrum_obj, *peak_indexes)    
        
        self._eic_data = eic_data

        self._possible_molecular_formulae = possible_molform if possible_molform else []

    def __len__(self):
        
        return len(self._dependent_mass_spectra)
        
    def __getitem__(self, position):
        
        return self._dependent_mass_spectra[position]

    def add_dependent_mass_spectrum(self, mass_spectrum):
        
        self._dependent_mass_spectra.append(mass_spectrum)

    def add_molecular_formula(self, molfform):
        
        self._possible_molecular_formulae.append(molfform)

    @property
    def eic_rt_list(self):
        return [self._eic_data.time[i] for i in range(self.start_index, self.final_index+1) ]

    @property   
    def eci_list(self):
        return [self._eic_data.eic[i] for i in range(self.start_index, self.final_index+1) ]


class GCPeak(ChromaPeakBase, GCPeakCalculation):

    def __init__(self, chromatogram_parent, mass_spectrum_obj, indexes):

        self._compounds = []

        self._ri = None

        super().__init__(chromatogram_parent, mass_spectrum_obj, *indexes)

    def __len__(self):
        
        return len(self._compounds)
        
    def __getitem__(self, position):
        
        return self._compounds[position]

    def remove_compound(self, compounds_obj):
        
        self._compounds.remove(compounds_obj)
        
    def clear_compounds(self):
        
        self._compounds = []

    def add_compound(self, compounds_dict, spectral_similarity_scores, ri_score=None, similarity_score=None):

        compound_obj = LowResCompoundRef(compounds_dict)

        # add all spectral similarities methods as a dict
        compound_obj.spectral_similarity_scores = spectral_similarity_scores
        # TODO need to add spectral similarity score label in the options in the parameters encapsulation class
        compound_obj.spectral_similarity_score = spectral_similarity_scores.get("cosine_correlation")

        compound_obj.ri_score = ri_score

        compound_obj.similarity_score = similarity_score

        self._compounds.append(compound_obj)

        if similarity_score:
            self._compounds.sort(key=lambda c: c.similarity_score, reverse=True)
        else:
            self._compounds.sort(key=lambda c: c.spectral_similarity_score, reverse=True)

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

    @property
    def compound_names(self):
        if self:
            return [c.name for c in self]
        else: 
            return []

class GCPeakDeconvolved(GCPeak):
    
    def __init__(self, chromatogram_parent, mass_spectra, apex_index, rt_list, tic_list ):
        
        self._ri = None   

        self._rt_list = list(rt_list)
        
        self._tic_list = list(tic_list)

        self.mass_spectra = list(mass_spectra)

        super().__init__(chromatogram_parent, self.mass_spectra[apex_index], (0, apex_index, len(self.mass_spectra)-1))

    @property
    def rt_list(self):
        return self._rt_list
    
    @property   
    def tic_list(self):   
        return self._tic_list