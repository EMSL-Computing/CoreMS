

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
from numpy import array, trapz    
from corems.chroma_peak.calc.ChromaPeakCalc import GCPeakCalculation
#import corems.mass_spectra.factory.LC_Class as lcms
from corems.mass_spectra.factory.LC_Temp import EIC_Data
from corems.molecular_id.factory.EI_SQL import LowResCompoundRef

class ChromaPeakBase():
    """ Base class for chromatographic peak (ChromaPeak) objects.

    Parameters
    -------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    start_index : int
        The start index of the peak.
    index : int
        The index of the peak.
    final_index : int
        The final index of the peak.

    Attributes
    --------
    start_scan : int
        The start scan of the peak.
    final_scan : int 
        The final scan of the peak.
    apex_scan : int
        The apex scan of the peak.
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum : MassSpectrum
        The mass spectrum object.
    _area : float 
        The area of the peak.

    Properties
    --------
    * retention_time : float. 
        The retention time of the peak.
    * tic : float. 
        The total ion current of the peak.
    * area : float. 
        The area of the peak.
    * rt_list : list. 
        The list of retention times within the peak.
    * tic_list : list.
        The list of total ion currents within the peak.

    Methods
    --------   
    * None
    """
    
    def __init__(self, chromatogram_parent, mass_spectrum_obj, start_index, index, final_index):
        self.start_scan = start_index
        self.final_scan = final_index
        self.apex_scan = int(index)
        self.chromatogram_parent = chromatogram_parent
        self.mass_spectrum = mass_spectrum_obj
        self._area = None

    @property   
    def retention_time(self):
        """Retention Time"""
        return self.mass_spectrum.retention_time

    @property   
    def tic(self):
        """Total Ion Current"""
        return self.mass_spectrum.tic    

    @property   
    def area(self):
        """Peak Area"""
        return self._area

    @property
    def rt_list(self):
        """Retention Time List"""
        return [self.chromatogram_parent.retention_time[i] for i in range(self.start_scan, self.final_scan+1) ]

    @property   
    def tic_list(self):
        """Total Ion Current List"""
        return [self.chromatogram_parent.tic[i] for i in range(self.start_scan, self.final_scan+1) ]

class DataDependentPeak(ChromaPeakBase):
    """ Data dependent peak (DataDependentPeak) object.

    Parameters
    -------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    peak_indexes : tuple
        The start, apex and final indexes of the peak.
    eic_data : EIC_Data
        The EIC data object.
    
    Attributes
    --------    
    _eic_data : EIC_Data
        The EIC data object.    
    _possible_molecular_formulae : list 
        The list of possible molecular formulas.
    
    Properties
    --------
    * eic_rt_list : list. 
        EIC retention time list.
    * eic_list : list.
        EIC list.
    * targeted_molecular_formulas : list.
        The list of possible molecular formulas.
    
    Methods
    --------
    * add_dependent_mass_spectrum(mass_spectrum).  
        Add a dependent mass spectrum to the peak.  
    * add_molecular_formula(molfform).  
        Add a molecular formula to the peak.  
    
    
    """
    def __init__(self, chromatogram_parent, mass_spectrum_obj, peak_indexes, 
                 eic_data: EIC_Data, possible_molform = None):
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
        """Add a dependent mass spectrum to the peak.
        
        Parameters
        ----------
        mass_spectrum : MassSpectrum
            The mass spectrum object.
        """

        self._dependent_mass_spectra.append(mass_spectrum)

    def add_molecular_formula(self, molfform):
        """Add a molecular formula to the peak.
        
        Parameters
        ----------
        molfform : str
            The molecular formula.
        """
        self._possible_molecular_formulae.extend(molfform)

    @property
    def eic_rt_list(self):
        """EIC Retention Time List"""
        return [self._eic_data.time[i] for i in range(self.start_scan, self.final_scan+1) ]

    @property   
    def eic_list(self):
        """EIC List"""
        return [self._eic_data.eic[i] for i in range(self.start_scan, self.final_scan+1) ]

    @property
    def targeted_molecular_formulas(self):
        """Targeted Molecular Formulas""" #Is this correct?
        return self._possible_molecular_formulae

class GCPeak(ChromaPeakBase, GCPeakCalculation):
    """ Class representing a peak in a gas chromatography (GC) chromatogram.

    Parameters
    ----------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object associated with the peak.
    indexes : tuple
        The indexes of the peak in the chromatogram.

    Attributes
    ----------
    _compounds : list
        List of compounds associated with the peak.
    _ri : float or None
        Retention index of the peak.

    Methods
    -------
    * __len__(). Returns the number of compounds associated with the peak.  
    * __getitem__(position).  Returns the compound at the specified position.  
    * remove_compound(compounds_obj). Removes the specified compound from the peak.  
    * clear_compounds(). Removes all compounds from the peak.  
    * add_compound(compounds_dict, spectral_similarity_scores, ri_score=None, similarity_score=None). Adds a compound to the peak with the specified attributes.  
    * ri().  Returns the retention index of the peak.  
    * highest_ss_compound(). Returns the compound with the highest spectral similarity score.  
    * highest_score_compound(). Returns the compound with the highest similarity score.  
    * compound_names(). Returns a list of names of compounds associated with the peak.  
    """
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
        """ Adds a compound to the peak with the specified attributes.

        Parameters
        ----------
        compounds_dict : dict
            Dictionary containing the compound information.
        spectral_similarity_scores : dict
            Dictionary containing the spectral similarity scores.
        ri_score : float or None, optional
            The retention index score of the compound. Default is None.
        similarity_score : float or None, optional
            The similarity score of the compound. Default is None.
        """
        compound_obj = LowResCompoundRef(compounds_dict)
        compound_obj.spectral_similarity_scores = spectral_similarity_scores
        compound_obj.spectral_similarity_score = spectral_similarity_scores.get("cosine_correlation")
        #TODO check is the above line correct?
        compound_obj.ri_score = ri_score
        compound_obj.similarity_score = similarity_score
        self._compounds.append(compound_obj)
        if similarity_score:
            self._compounds.sort(key=lambda c: c.similarity_score, reverse=True)
        else:
            self._compounds.sort(key=lambda c: c.spectral_similarity_score, reverse=True)

    @property
    def ri(self):
        """Returns the retention index of the peak.

        Returns
        -------
        float or None
            The retention index of the peak.
        """
        return self._ri

    @property
    def highest_ss_compound(self):
        """ Returns the compound with the highest spectral similarity score.

        Returns
        -------
        LowResCompoundRef or None
            The compound with the highest spectral similarity score.
        """
        if self:
            return max(self, key=lambda c: c.spectral_similarity_score)
        else:
            return None
    
    @property
    def highest_score_compound(self):
        """ Returns the compound with the highest similarity score.

        Returns
        -------
        LowResCompoundRef or None
            The compound with the highest similarity score.
        """
        if self:
            return max(self, key=lambda c: c.similarity_score)
        else: 
            return None

    @property
    def compound_names(self):
        """ Returns a list of names of compounds associated with the peak.

        Returns
        -------
        list
            List of names of compounds associated with the peak.
        """
        if self:
            return [c.name for c in self]
        else: 
            return []

class GCPeakDeconvolved(GCPeak):
    """ Represents a deconvolved peak in a chromatogram.
    
    Parameters
    ----------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectra : list
        List of mass spectra associated with the peak.
    apex_index : int
        Index of the apex mass spectrum in the `mass_spectra` list.
    rt_list : list
        List of retention times.
    tic_list : list
        List of total ion currents.
    """
    
    def __init__(self, chromatogram_parent, mass_spectra, apex_index, rt_list, tic_list ):
        self._ri = None   
        self._rt_list = list(rt_list)
        self._tic_list = list(tic_list)
        self.mass_spectra = list(mass_spectra)
        super().__init__(chromatogram_parent, self.mass_spectra[apex_index], (0, apex_index, len(self.mass_spectra)-1))

    @property
    def rt_list(self):
        """ Get the list of retention times.
        
        Returns
        -------
        list
            The list of retention times.
        """
        return self._rt_list
    
    @property   
    def tic_list(self):   
        """ Get the list of total ion currents.
        
        Returns
        -------
        list
            The list of total ion currents.
        """
        return self._tic_list