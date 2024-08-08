__author__ = 'Yuri E. Corilo'
__date__ = 'Jul 02, 2019'

import dataclasses
import os
from typing import List, Dict

from corems.encapsulation.constant import Atoms, Labels

@dataclasses.dataclass
class TransientSetting:
    """Transient processing settings class

    Attributes
    ----------
    implemented_apodization_function : tuple
        Available apodization functions
    apodization_method : str
        Apodization function to use. Hanning is a good default for Fourier transform magnitude mode. For absorption mode processing, Half-Sine or Half-Kaiser may be more appropriate.
    number_of_truncations : int
        How many times to truncate the transient prior to Fourier transform
    number_of_zero_fills : int
        How many times to zero fille the transient prior to Fourier transform.
    next_power_of_two : bool
        If True, zero fill to the next power of two after the new length of len(transient)+(number_of_zero_fills*len(transient)). 
    kaiser_beta : float
        Beta parameter for Kaiser or Half-Kaiser apodisation function. 0 is rectangular,  5 is similar to Hamming,
        6 is similar to hanning, and 8.6 is similar to Blackman (from numpy docs)

    """
    implemented_apodization_function: tuple = ('Hamming', 'Hanning', 'Blackman','Full-Sine','Half-Sine','Kaiser','Half-Kaiser')
    apodization_method: str = 'Hanning'
    number_of_truncations: int = 0
    number_of_zero_fills: int = 1
    next_power_of_two: bool = False
    kaiser_beta: float = 8.6

    def __post_init__(self):

        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class DataInputSetting:
    """Data input settings class
    
    Attributes
    ----------
    header_translate : dict
        Dictionary with the header labels to be translated to the corems labels. For example, {'m/z':'m/z', 'Resolving Power':'Resolving Power', 'Abundance':'Abundance' , 'S/N':'S/N'}
    """
    # add to this dict the VALUES to match your labels, THE ORDER WON"T MATTER
    # "column_translate" : {"m/z":"m/z", "Resolving Power":"Resolving Power", "Abundance":"Abundance" , "S/N":"S/N"}
    header_translate: dict = dataclasses.field(default_factory=dict)

    def __post_init__(self):

        self.header_translate = {'m/z': Labels.mz, 
                                 'mOz': Labels.mz,
                                 'Mass': Labels.mz,
                                 'Resolving Power': Labels.rp,
                                 'Res.': Labels.rp,
                                 'resolution': Labels.rp,
                                 'Intensity': Labels.abundance,
                                 'I': Labels.abundance,
                                 'Abundance': Labels.abundance,
                                 'abs_abu': Labels.abundance,
                                 'Signal/Noise': Labels.s2n,
                                 'S/N': Labels.s2n,
                                 'sn': Labels.s2n}

    def add_mz_label(self, label):
        """Add a label to the header_translate dictionary to be translated to the corems label for mz."""
        self.header_translate[label] = Labels.mz

    def add_peak_height_label(self, label):
        """Add a label to the header_translate dictionary to be translated to the corems label for peak height."""

        self.header_translate[label] = Labels.abundance

    def add_sn_label(self, label):
        """Add a label to the header_translate dictionary to be translated to the corems label for signal to noise."""
        self.header_translate[label] = Labels.s2n

    def add_resolving_power_label(self, label):
        """Add a label to the header_translate dictionary to be translated to the corems label for resolving power."""
        self.header_translate[label] = Labels.rp

@dataclasses.dataclass
class LiquidChromatographSetting:
    """Liquid chromatograph processing settings class

    Attributes
    ----------
    scans : list or tuple, optional
        List of select scan to average or a tuple containing the range to average. Default is (0, 1).
    eic_tolerance_ppm : float, optional
        Mass tolerance in ppm for extracted ion chromatogram peak detection. Default is 5.
    correct_eic_baseline : bool, optional
        If True, correct the baseline of the extracted ion chromatogram. Default is True.
    enforce_target_ms2 : bool, optional
        Only extract EICs for target_mz if the m/z was selected as precursor for ms2. Default is True.
    average_target_mz : bool, optional
        If True, average the target_mz EICs. Default is True.
    smooth_window : int, optional
        Window size for smoothing the ion chromatogram. Default is 5.
    smooth_method : str, optional
        Smoothing method to use. Default is 'savgol'. Other options are 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar'.
    implemented_smooth_method : tuple, optional
        Smoothing methods that can be implemented. Values are ('savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar').
    savgol_pol_order : int, optional
        Polynomial order for Savitzky-Golay smoothing. Default is 2.
    peak_height_max_percent : float, optional
        1-100 % used for baseline detection use 0.1 for second_derivative and 10 for other methods. Default is 10.
    peak_max_prominence_percent : float, optional
        1-100 % used for baseline detection. Default is 1.
    peak_derivative_threshold : float, optional
        Threshold for defining derivative crossing. Default is 0.0005.
    min_peak_datapoints : float, optional
        minimum data point to define a chromatografic peak. Default is 5.
    noise_threshold_method : str, optional
        Method for detecting noise threshold. Default is 'manual_relative_abundance'.
    noise_threshold_methods_implemented : tuple, optional
        Methods for detected noise threshold that can be implemented. Default is ('auto_relative_abundance', 'manual_relative_abundance', 'second_derivative').
    std_noise_threshold : int, optional
        the amount of standard deviations used to calculate noise thresould, average + (std_noise_threshold * noise). Default is 3.
    peak_height_min_percent : float, optional
        0-100 % used for peak detection. Default is 0.1.
    eic_signal_threshold : float, optional
        0-100 % used for extracted ion chromatogram peak detection. Default is 0.01.
    """
   
    scans: list | tuple = (-1,-1)
    
    eic_tolerance_ppm: float = 5
    
    correct_eic_baseline = True

    enforce_target_ms2: bool = True

    average_target_mz: bool = True

    smooth_window: int = 5

    smooth_method: str = 'savgol'

    implemented_smooth_method: tuple = ('savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar')

    savgol_pol_order: int = 2

    peak_height_max_percent: float = 10  

    peak_max_prominence_percent: float = 1

    peak_derivative_threshold:float = 0.0005

    min_peak_datapoints: float = 5

    noise_threshold_method: str = 'manual_relative_abundance'

    noise_threshold_methods_implemented: tuple = ('auto_relative_abundance', 'manual_relative_abundance', 'second_derivative')

    std_noise_threshold: int = 3

    peak_height_min_percent: float = 0.1  

    eic_signal_threshold: float = 0.01  

    def __post_init__(self):
        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class MassSpectrumSetting:
    """Mass spectrum processing settings class

    Attributes
    ----------
    noise_threshold_method : str, optional
        Method for detecting noise threshold. Default is 'log'.
    noise_threshold_methods_implemented : tuple, optional
        Methods for detected noise threshold that can be implemented. Default is ('minima', 'signal_noise', 'relative_abundance', 'absolute_abundance', 'log').
    noise_threshold_min_std : int, optional
        Minumum value for noise thresholding when using 'minima' noise threshold method. Default is 6.
    noise_threshold_min_s2n : float, optional
        Minimum value for noise thresholding when using 'signal_noise' noise threshold method. Default is 4.
    noise_threshold_min_relative_abundance : float, optional
        Minimum value for noise thresholding when using 'relative_abundance' noise threshold method. Default is 6.
    noise_threshold_absolute_abundance : float, optional
        Minimum value for noise thresholding when using 'absolute_abundance' noise threshold method. Default is 1_000_000.
    noise_threshold_log_nsigma : int, optional
        Number of standard deviations to use when using 'log' noise threshold method. Default is 6.
    noise_threshold_log_nsigma_corr_factor : float, optional
        Correction factor for log noise threshold method. Default is 0.463.
    noise_threshold_log_nsigma_bins : int, optional
        Number of bins to use for histogram when using 'log' noise threshold method. Default is 500.
    noise_min_mz : float, optional
        Minimum m/z to use for noise thresholding. Default is 50.0.
    noise_max_mz : float, optional
        Maximum m/z to use for noise thresholding. Default is 1200.0.
    min_picking_mz : float, optional
        Minimum m/z to use for peak picking. Default is 50.0.
    max_picking_mz : float, optional
        Maximum m/z to use for peak picking. Default is 1200.0.
    picking_point_extrapolate : int, optional
        How many data points (in each direction) to extrapolate the mz axis and 0 pad the abundance axis. Default is 3.
    calib_minimize_method : str, optional
        Minimization method to use for calibration. Default is 'Powell'.
    calib_pol_order : int, optional
        Polynomial order to use for calibration. Default is 2.
    max_calib_ppm_error : float, optional
        Maximum ppm error to use for calibration. Default is 1.0.
    min_calib_ppm_error : float, optional
        Minimum ppm error to use for calibration. Default is -1.0.
    calib_sn_threshold : float, optional
        Signal to noise threshold to use for calibration. Default is 2.0.
    calibration_ref_match_method: string, optional
        Method for matching reference masses with measured masses for recalibration. Default is 'legacy'. 
    calibration_ref_match_tolerance: float, optional
        If using the new method for calibration reference mass matching, this tolerance is the initial matching tolerance. Default is 0.003
    do_calibration : bool, optional
        If True, perform calibration. Default is True.    
    """
    noise_threshold_method: str = 'log'

    noise_threshold_methods_implemented: tuple = ('minima', 'signal_noise', 'relative_abundance', 'absolute_abundance', 'log')

    noise_threshold_min_std: int = 6 # when using 'minima' method

    noise_threshold_min_s2n: float = 4 # when using 'signal_noise' method

    noise_threshold_min_relative_abundance: float = 6  # from 0-100, when using 'relative_abundance' method

    noise_threshold_absolute_abundance: float = 1_000_000 # when using 'absolute_abundance' method

    noise_threshold_log_nsigma: int = 6 # when using 'log' method
    noise_threshold_log_nsigma_corr_factor: float = 0.463 #mFT is 0.463, aFT is 1.0
    noise_threshold_log_nsigma_bins: int = 500 # bins for the histogram for the noise

    noise_min_mz: float = 50.0
    noise_max_mz: float = 1200.0

    min_picking_mz: float = 50.0
    max_picking_mz: float = 1200.0

    # How many data points (in each direction) to extrapolate the mz axis and 0 pad the abundance axis
    # This will fix peak picking at spectrum limit issues
    #  0 to keep normal behaviour, typical value 3 to fix
    picking_point_extrapolate: int = 0 

    calib_minimize_method: str = 'Powell'
    calib_pol_order: int = 2
    max_calib_ppm_error: float = 1.0
    min_calib_ppm_error: float = -1.0
    calib_sn_threshold: float = 2.0
    calibration_ref_match_method: str = 'legacy'
    calibration_ref_match_method_implemented: tuple = ('legacy', 'merged')
    calibration_ref_match_tolerance: float = 0.003
    calibration_ref_match_std_raw_error_limit: float = 1.5
    #calib_ref_mzs: list = [0]

    do_calibration: bool = True

    def __post_init__(self):
        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class MassSpecPeakSetting:
    """Mass spectrum peak processing settings class

    Attributes
    ----------
    kendrick_base : Dict, optional
        Dictionary specifying the elements and their counts in the Kendrick base.
        Defaults to {'C': 1, 'H': 2}.
    kendrick_rounding_method : str, optional
        Method for calculating the nominal Kendrick mass. Valid values are 'floor', 'ceil', or 'round'.
        Defaults to 'floor'.
    implemented_kendrick_rounding_methods : tuple
        Tuple of valid rounding methods for calculating the nominal Kendrick mass.
        Defaults to ('floor', 'ceil', 'round').
    peak_derivative_threshold : float, optional
        Threshold for defining derivative crossing. Should be a value between 0 and 1.
        Defaults to 0.0.
    peak_min_prominence_percent : float, optional
        Minimum prominence percentage used for peak detection. Should be a value between 1 and 100.
        Defaults to 0.1.
    min_peak_datapoints : float, optional
        Minimum number of data points used for peak detection. Should be a value between 0 and infinity.
        Defaults to 5.
    peak_max_prominence_percent : float, optional
        Maximum prominence percentage used for baseline detection. Should be a value between 1 and 100.
        Defaults to 0.1.
    peak_height_max_percent : float, optional
        Maximum height percentage used for baseline detection. Should be a value between 1 and 100.
        Defaults to 10.
    legacy_resolving_power : bool, optional
        Flag indicating whether to use the legacy (CoreMS v1) resolving power calculation.
        Defaults to True.
    """

    kendrick_base: Dict = dataclasses.field(default_factory=dict)
    
    kendrick_rounding_method: str = 'floor' # 'floor', 'ceil' or 'round' are valid methods for calculating nominal kendrick mass
    
    implemented_kendrick_rounding_methods : tuple = ('floor','ceil','round')

    peak_derivative_threshold: float =  0.0 # define derivative crossing threshould 0-1

    peak_min_prominence_percent: float = 0.1  # 1-100 % used for peak detection

    min_peak_datapoints: float = 5 # 0-inf used for peak detection

    peak_max_prominence_percent: float = 0.1  # 1-100 % used for baseline detection

    peak_height_max_percent: float = 10  # 1-100 % used for baseline detection

    legacy_resolving_power: bool = True # Use the legacy (CoreMS v1) resolving power calculation (True)
    
    def __post_init__(self):

        # default to CH2
        if not self.kendrick_base:
            self.kendrick_base = {'C': 1, 'H': 2}
        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class GasChromatographSetting:
    """Gas chromatograph processing settings class

    Attributes
    ----------
    use_deconvolution : bool, optional
        If True, use deconvolution. Default is False.
    implemented_smooth_method : tuple, optional
        Smoothing methods that can be implemented. Default is ('savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar').
    smooth_window : int, optional
        Window size for smoothing the ion chromatogram. Default is 5.
    smooth_method : str, optional
        Smoothing method to use. Default is 'savgol'. Other options are 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar'.
    savgol_pol_order : int, optional
        Polynomial order for Savitzky-Golay smoothing. Default is 2.
    peak_derivative_threshold : float, optional
        Threshold for defining derivative crossing. Should be a value between 0 and 1.
        Defaults to 0.0005.
    peak_height_max_percent : float, optional
        Maximum height percentage used for baseline detection. Should be a value between 1 and 100.
        Defaults to 10.
    peak_max_prominence_percent : float, optional
        Maximum prominence percentage used for baseline detection. Should be a value between 1 and 100.
        Defaults to 1.
    min_peak_datapoints : float, optional
        Minimum number of data points used for peak detection. Should be a value between 0 and infinity.
        Defaults to 5.
    max_peak_width : float, optional
        Maximum peak width used for peak detection. Should be a value between 0 and infinity.
        Defaults to 0.1.
    noise_threshold_method : str, optional
        Method for detecting noise threshold. Default is 'manual_relative_abundance'.
    noise_threshold_methods_implemented : tuple, optional
        Methods for detected noise threshold that can be implemented. Default is ('auto_relative_abundance', 'manual_relative_abundance', 'second_derivative').
    std_noise_threshold : int, optional
        Default is 3.
    peak_height_min_percent : float, optional
        0-100 % used for peak detection. Default is 0.1.
    peak_min_prominence_percent : float, optional
        0-100 % used for peak detection. Default is 0.1.
    eic_signal_threshold : float, optional
        0-100 % used for extracted ion chromatogram peak detection. Default is 0.01.
    max_rt_distance : float, optional
        Maximum distance allowance for hierarchical cluster, in minutes. Default is 0.025.
    """
    use_deconvolution: bool = False

    implemented_smooth_method: tuple = ('savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar')

    smooth_window: int = 5

    smooth_method: str = 'savgol'

    savgol_pol_order: int = 2

    peak_derivative_threshold:float = 0.0005

    peak_height_max_percent: float = 10  # 1-100 % used for baseline detection use 0.1 for second_derivative and 10 for other methods

    peak_max_prominence_percent: float = 1  # 1-100 % used for baseline detection

    min_peak_datapoints: float = 5

    max_peak_width: float = 0.1

    noise_threshold_method: str = 'manual_relative_abundance'

    noise_threshold_methods_implemented: tuple = ('auto_relative_abundance', 'manual_relative_abundance', 
                                                  'second_derivative')

    std_noise_threshold: int = 3

    peak_height_min_percent: float = 0.1  # 0-100 % used for peak detection

    peak_min_prominence_percent: float = 0.1  # 0-100 % used for peak detection

    eic_signal_threshold: float = 0.01  # 0-100 % used for extracted ion chromatogram peak detection

    max_rt_distance: float = 0.025  # minutes, max distance allowance hierarchical clutter

    def __post_init__(self):

        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class CompoundSearchSettings:
    """Settings for compound search

    Attributes
    ----------
    url_database : str, optional
        URL for the database. Default is 'sqlite:///db/pnnl_lowres_gcms_compounds.sqlite'.
    ri_search_range : float, optional
        Retention index search range. Default is 35.
    rt_search_range : float, optional
        Retention time search range, in minutes. Default is 1.0.
    correlation_threshold : float, optional
        Threshold for correlation for spectral similarity. Default is 0.5.
    score_threshold : float, optional
        Threshold for compsite score. Default is 0.0.
    ri_spacing : float, optional
        Retention index spacing. Default is 200.
    ri_std : float, optional
        Retention index standard deviation. Default is 3.
    ri_calibration_compound_names : list, optional
        List of compound names to use for retention index calibration. Default is ['Methyl Caprylate', 'Methyl Caprate', 'Methyl Pelargonate', 'Methyl Laurate', 'Methyl Myristate', 'Methyl Palmitate', 'Methyl Stearate', 'Methyl Eicosanoate', 'Methyl Docosanoate', 'Methyl Linocerate', 'Methyl Hexacosanoate', 'Methyl Octacosanoate', 'Methyl Triacontanoate'].
    
    """
    url_database: str = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/lowres" # 'postgresql://postgres:labthomson0102@172.22.113.27:5432/GCMS' # 'sqlite:///db/pnnl_lowres_gcms_compounds.sqlite'

    ri_search_range: float = 35

    rt_search_range: float = 1.0  # used for retention index calibration

    correlation_threshold: float = 0.5  # used for calibration, spectral similarity 

    score_threshold: float = 0.0

    ri_spacing: float = 200

    ri_std: float = 3  # in standard deviation

    ri_calibration_compound_names: List = dataclasses.field(default_factory=list)

    # calculates and export all spectral similarity methods
    exploratory_mode: bool = False

    score_methods: tuple = ('highest_sim_score', 'highest_ss')

    output_score_method: str = 'All'

    

    def __post_init__(self):
        # enforce datatype
        self.url_database = os.getenv('SPECTRAL_GCMS_DATABASE_URL', 'sqlite:///db/pnnl_lowres_gcms_compounds.sqlite')

        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

        self.ri_calibration_compound_names = ['Methyl Caprylate', 
                                              'Methyl Caprate', 
                                              'Methyl Pelargonate', 
                                              'Methyl Laurate', 
                                              'Methyl Myristate', 
                                              'Methyl Palmitate', 
                                              'Methyl Stearate', 
                                              'Methyl Eicosanoate', 
                                              'Methyl Docosanoate', 
                                              'Methyl Linocerate', 
                                              'Methyl Hexacosanoate', 
                                              'Methyl Octacosanoate', 
                                              'Methyl Triacontanoate']

class MolecularLookupDictSettings:
    """Settings for molecular searching

    These are used to generate the database entries, do not change.

    Attributes
    ----------
    usedAtoms : dict, optional
        Dictionary of atoms and ranges. Default is {'C': (1, 90), 'H': (4, 200), 'O': (0, 12), 'N': (0, 0), 'S': (0, 0), 'P': (0, 0), 'Cl': (0, 0)}.
    min_mz : float, optional
        Minimum m/z to use for searching. Default is 50.0.
    max_mz : float, optional
        Maximum m/z to use for searching. Default is 1200.0.
    min_dbe : float, optional
        Minimum double bond equivalent to use for searching. Default is 0.
    max_dbe : float, optional
        Maximum double bond equivalent to use for searching. Default is 50.
    use_pah_line_rule : bool, optional
        If True, use the PAH line rule. Default is False.
    isRadical : bool, optional
        If True, search for radical ions. Default is True.
    isProtonated : bool, optional
        If True, search for protonated ions. Default is True.
    min_op_filter : float, optional
        Minimum oxygen to phosphorous ratio. Default is 2.
    min_hc_filter : float, optional
        Minimum hydrogen to carbon ratio. Default is 0.3.
    min_oc_filter : float, optional
        Minimum oxygen to carbon ratio. Default is 1.2.
    max_oc_filter : float, optional
        Maximum oxygen to carbon ratio. Default is 0.0.
    url_database : str, optional
        URL for the database. Default is None.
    db_jobs : int, optional
        Number of jobs to use for database queries. Default is 1.
    used_atom_valences : dict, optional
        Dictionary of atoms and valences. Default is {'C': 4, '13C': 4, 'H': 1, 'O': 2, '18O': 2, 'N': 3, 'S': 2, '34S': 2, 'P': 3, 'Cl': 1, '37Cl': 1, 'Br': 1, 'Na': 1, 'F': 1, 'K': 0}.
        
    """
    ### DO NOT CHANGE IT! These are used to generate the database entries 

    ### DO change when creating a new application database 

    ### FOR search settings runtime and database query check use the MolecularFormulaSearchSettings class below

    ### C, H, N, O, S and P atoms are ALWAYS needed at usedAtoms
    ### if you don't want to include one of those atoms set the max and min at 0
    ### you can include any atom listed at Atoms class inside encapsulation.settings.constants module
    ### make sure to include the selected covalence at the used_atoms_valences when adding new atoms 
    ### NOTE : Adducts atoms have zero covalence
    ### NOTE : Not using static variable because this class is distributed using multiprocessing
    def __init__(self):

        self.usedAtoms = {'C': (1, 90),
                          'H': (4, 200),
                          'O': (0, 12),
                          'N': (0, 0),
                          'S': (0, 0),
                          'P': (0, 0),
                          'Cl': (0, 0),
                          }

        self.min_mz = 50

        self.max_mz = 1200

        self.min_dbe = 0

        self.max_dbe = 50

        # overwrites the dbe limits above to DBE = (C + heteroatoms) * 0.9
        self.use_pah_line_rule = False

        self.isRadical = True

        self.isProtonated = True

        self.min_op_filter = 2

        self.min_hc_filter = 0.3

        self.min_oc_filter = 1.2

        self.max_oc_filter = 0.0

        self.url_database = None

        self.db_jobs = 1

        self.used_atom_valences = {'C': 4,
                                   '13C': 4,
                                   'H': 1,
                                   'O': 2,
                                   '18O': 2,
                                   'N': 3,
                                   'S': 2,
                                   '34S': 2,
                                   'P': 3,
                                   'Cl': 1,
                                   '37Cl': 1,
                                   'Br': 1,
                                   'Na': 1,
                                   'F': 1,
                                   'K': 0,
                                   }

@dataclasses.dataclass
class MolecularFormulaSearchSettings:
    """Settings for molecular searching

    Attributes
    ----------
    use_isotopologue_filter : bool, optional
        If True, use isotopologue filter. Default is False.
    isotopologue_filter_threshold : float, optional
        Threshold for isotopologue filter. Default is 33.
    isotopologue_filter_atoms : tuple, optional
        Tuple of atoms to use for isotopologue filter. Default is ('Cl', 'Br').
    use_runtime_kendrick_filter : bool, optional
        If True, use runtime Kendrick filter. Default is False.
    use_min_peaks_filter : bool, optional
        If True, use minimum peaks filter. Default is True.
    min_peaks_per_class : int, optional
        Minimum number of peaks per class. Default is 15.
    url_database : str, optional
        URL for the database. Default is 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'.
    db_jobs : int, optional
        Number of jobs to use for database queries. Default is 3.
    db_chunk_size : int, optional
        Chunk size to use for database queries. Default is 300.
    ion_charge : int, optional
        Ion charge. Default is -1.
    min_hc_filter : float, optional
        Minimum hydrogen to carbon ratio. Default is 0.3.
    max_hc_filter : float, optional
        Maximum hydrogen to carbon ratio. Default is 3.
    min_oc_filter : float, optional
        Minimum oxygen to carbon ratio. Default is 0.0.
    max_oc_filter : float, optional
        Maximum oxygen to carbon ratio. Default is 1.2.
    min_op_filter : float, optional
        Minimum oxygen to phosphorous ratio. Default is 2.
    use_pah_line_rule : bool, optional
        If True, use the PAH line rule. Default is False.
    min_dbe : float, optional
        Minimum double bond equivalent to use for searching. Default is 0.
    max_dbe : float, optional
        Maximum double bond equivalent to use for searching. Default is 40.
    mz_error_score_weight : float, optional
        Weight for m/z error score to contribute to composite score. Default is 0.6.
    isotopologue_score_weight : float, optional
        Weight for isotopologue score to contribute to composite score. Default is 0.4.
    adduct_atoms_neg : tuple, optional
        Tuple of atoms to use in negative polarity. Default is ('Cl', 'Br').
    adduct_atoms_pos : tuple, optional
        Tuple of atoms to use in positive polarity. Default is ('Na', 'K').
    score_methods : tuple, optional
        Tuple of score method that can be implemented. 
        Default is ('S_P_lowest_error', 'N_S_P_lowest_error', 'lowest_error', 'prob_score', 'air_filter_error', 'water_filter_error', 'earth_filter_error').
    score_method : str, optional
        Score method to use. Default is 'prob_score'. Options are 'S_P_lowest_error', 'N_S_P_lowest_error', 'lowest_error', 'prob_score', 'air_filter_error', 'water_filter_error', 'earth_filter_error'.   
    output_min_score : float, optional
        Minimum score for output. Default is 0.1.
    output_score_method : str, optional
        Score method to use for output. Default is 'All Candidates'.
    output_min_score : float, optional
        Minimum score for output. Default is 0.1.
    isRadical : bool, optional
        If True, search for radical ions. Default is False.
    isProtonated : bool, optional
        If True, search for protonated ions. Default is True.
    isAdduct : bool, optional
        If True, search for adduct ions. Default is False.
    usedAtoms : dict, optional
        Dictionary of atoms and ranges. Default is {'C': (1, 90), 'H': (4, 200), 'O': (0, 12), 'N': (0, 0), 'S': (0, 0), 'P': (0, 0), 'Cl': (0, 0)}.
    ionization_type : str, optional
        Ionization type. Default is 'ESI'.
    min_ppm_error : float, optional
        Minimum ppm error. Default is -10.0.
    max_ppm_error : float, optional
        Maximum ppm error. Default is 10.0.
    min_abun_error : float, optional
        Minimum abundance error for isotolopologue search. Default is -100.0.
    max_abun_error : float, optional
        Maximum abundance error for isotolopologue search. Default is 100.0.
    mz_error_range : float, optional
        m/z error range. Default is 1.5.
    error_method : str, optional
        Error method. Default is 'None'. Options are 'distance', 'lowest', 'symmetrical','average' 'None'.
    mz_error_average : float, optional
        m/z error average. Default is 0.0.
    used_atom_valences : dict, optional
        Dictionary of atoms and valences. Default is {'C': 4, '13C': 4, 'H': 1, 'O': 2, '18O': 2, 'N': 3, 'S': 2, '34S': 2, 'P': 3, 'Cl': 1, '37Cl': 1, 'Br': 1, 'Na': 1, 'F': 1, 'K': 0}.
    """
    use_isotopologue_filter: bool = False

    isotopologue_filter_threshold: float = 33

    isotopologue_filter_atoms: tuple = ('Cl', 'Br')

    use_runtime_kendrick_filter: bool = False

    use_min_peaks_filter: bool = True

    min_peaks_per_class: int = 15

    url_database: str = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'

    db_jobs: int = 3

    db_chunk_size: int = 300

    #query setting========
    ion_charge: int = -1

    min_hc_filter: float = 0.3

    max_hc_filter: float = 3

    min_oc_filter: float = 0.0

    max_oc_filter: float = 1.2

    min_op_filter: float = 2

    use_pah_line_rule: bool = False

    min_dbe: float = 0

    max_dbe: float = 40

    mz_error_score_weight: float = 0.6

    isotopologue_score_weight: float = 0.4

    # look for close shell ions [M + Adduct]+ only considers metal set in the list adduct_atoms  
    adduct_atoms_neg: tuple = ('Cl', 'Br')

    adduct_atoms_pos: tuple = ('Na', 'K')

    score_methods: tuple = ('S_P_lowest_error', 'N_S_P_lowest_error', 'lowest_error', 'prob_score',
                            'air_filter_error', 'water_filter_error', 'earth_filter_error')

    score_method: str = 'prob_score'

    output_min_score: float = 0.1

    output_score_method: str = 'All Candidates'

    output_min_score: float = 0.1

    # depending on the polarity mode it looks for [M].+ , [M].-
    # query and automatically compile add entry if it doesn't exist

    isRadical: bool = False

    # depending on the polarity mode it looks for [M + H]+ , [M - H]+
    # query and automatically compile and push options if it doesn't exist
    isProtonated: bool = True

    isAdduct: bool = False

    usedAtoms: dict = dataclasses.field(default_factory=dict)

    # search setting ========

    ionization_type: str = 'ESI'

    # empirically set / needs optimization
    min_ppm_error: float = -10.0  # ppm

    # empirically set / needs optimization
    max_ppm_error: float = 10.0  # ppm

    # empirically set / needs optimization set for isotopologue search
    min_abun_error: float = -100.0  # percentage

    # empirically set / needs optimization set for isotopologue search
    max_abun_error: float = 100.0  # percentage

    # empirically set / needs optimization
    mz_error_range: float = 1.5

    # 'distance', 'lowest', 'symmetrical','average' 'None'
    error_method: str = 'None'

    mz_error_average: float = 0.0

    # used_atom_valences: {'C': 4, 'H':1, etc} = dataclasses.field(default_factory=dict)
    used_atom_valences: dict = dataclasses.field(default_factory=dict)

    def __post_init__(self):

        self.url_database = os.getenv('COREMS_DATABASE_URL', 'sqlite:///db/molformula.db')
        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

        # enforce C and H if either do not exists
        if 'C' not in self.usedAtoms.keys():
            self.usedAtoms['C'] = (1, 100)
        if 'H' not in self.usedAtoms.keys():
            self.usedAtoms['H'] = (1, 200)

        # add cummon values
        current_used_atoms = self.used_atom_valences.keys()
        
        for atom in Atoms.atoms_covalence.keys():
            
            if atom not in current_used_atoms:
                
                covalence = Atoms.atoms_covalence.get(atom)
                
                if isinstance(covalence , int):
                    self.used_atom_valences[atom] = covalence
                
                else:
                    #will get the first number of all possible covalances, which should be the most commum 
                    self.used_atom_valences[atom] = covalence[0]

