__author__ = 'Yuri E. Corilo'
__date__ = 'Jul 02, 2019'

from dataclasses import dataclass,field
from typing import List, Dict
from corems.encapsulation.constant import Labels

@dataclass
class TransientSetting:
    
    implemented_apodization_function: tuple = ('Hamming', 'Hanning', 'Blackman')
    apodization_method: str = 'Hanning'
    number_of_truncations: int = 0
    number_of_zero_fills: int = 1

@dataclass            
class MassSpectrumSetting:
    
    threshold_method: str = 'auto'
    
    implemented_noise_threshold_methods: tuple = ('auto', 'signal_noise', 'relative_abundance')
    
    noise_threshold_std: int = 6
    
    s2n_threshold: float = 4
    
    relative_abundance_threshold:float = 6 # from 1-100
    
    min_noise_mz: float = 100.0
    max_noise_mz:float = 1200.0
    
    min_picking_mz:float = 100.0
    max_picking_mz:float = 1200.0

@dataclass    
class MassSpecPeakSetting:
    
    kendrick_base: Dict = field(default_factory=dict)
    #kendrick_base : Dict =  {'C': 1, 'H':2}
    
    peak_min_prominence_percent :float = 1 #1-100 % used for peak detection

    peak_max_prominence_percent :float = 0.1 #1-100 % used for baseline detection

    def __post_init__(self):
        
        self.kendrick_base = {'C': 1, 'H':2}
    
@dataclass 
class GasChromatographSetting:
    
    implemented_smooth_method: tuple = ('savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar')
    
    smooth_window: int = 5
    
    smooth_method: str = 'savgol'
    
    savgol_pol_order: int = 2
    
    peak_height_max_percent: float = 10 #1-100 % used for baseline detection use 0.1 for second_derivative and 10 for other methods

    peak_max_prominence_percent:float = 1 #1-100 % used for baseline detection

    min_peak_datapoints:float = 3
   
    max_peak_width:float = 0.1

    noise_threshold_method:str = 'manual_relative_abundance'
    
    implemented_noise_threshold_methods: tuple = ('auto_relative_abundance', 'manual_relative_abundance', 'second_derivative')
    
    std_noise_threshold: int = 3

    peak_height_min_abun:float = 0.2 #1-100 % used for peak detection

    peak_min_prominence_percent:float = 1 #1-100 % used for peak detection


ri_calibration_compound_names = (" [C8] Methyl Caprylate [7.812]",
" [C10] Methyl Caprate [10.647]",
" [C9] Methyl Pelargonate [9.248]",
" [C12] Methyl Laurate [13.250]",
" [C14] Methyl Myristate [15.597]",
" [C16] Methyl Palmitate [17.723]",
" [C18] Methyl Stearate [19.663]",
" [C20] Methyl Eicosanoate [21.441]",
" [C22] Methyl Docosanoate [23.082]",
" [C24] Methyl Linocerate [24.603]",
" [C26] Methyl Hexacosanoate [26.023]",
" [C28] Methyl Octacosanoate [27.349]",
" [C30] Methyl Triacontanoate [28.723]")

@dataclass 
class CompoundSearchSettings:

    ri_search_range:float = 20

    rt_search_range:float = 0.5
    
    correlation_threshold:float = 0.5 # used for calibration, spectral similarity 
    
    score_threshold:float = 0.0

    ri_spacing:float = 200

    ri_std:float = 3 # in standard deviation

@dataclass 
class MolecularLookupDictSettings:
    
    '''
    ### DO NOT CHANGE IT! These are used to generate the database entries 

    ### DO change when creating a new application database 
    
    ### FOR search settings runtime and database query check use the MolecularFormulaSearchSettings class below

    ### C, H, N, O, S and P atoms are ALWAYS needed at usedAtoms
    ### if you don't want to include one of those atoms set the max and min at 0
    ### you can include any atom listed at Atoms class inside encapsulation.settings.constants module
    ### make sure to include the selected covalence at the used_atoms_valences when adding new atoms 
    ### NOTE : Adducts atoms have zero covalence
    ### NOTE : Not using static variable because this class is distributed using multiprocessing
    '''
    def __init__(self):
        
        
        self.usedAtoms = {'C': (1, 100),
                    'H': (4, 200),
                    'O': (0, 22),
                    'N': (0, 0),
                    'S': (0, 0),
                    'P': (0, 0),
                    'Cl': (0, 0),
                    }
        
        #min_mz changes automatically with mass spectrum
        self.min_mz = 100

        #max_mz changes automatically with mass spectrum
        self.max_mz = 1200

        self.min_dbe = 0

        self.max_dbe = 50

        #overwrites the dbe limits above to DBE = (C + heteroatoms) * 0.9
        self.use_pah_line_rule = False

        self.isRadical = True

        self.isProtonated = True

        self.op_filter = 2

        self.hc_filter = 0.3

        self.oc_filter = 1.2
    
        self.db_directory = None

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
        
@dataclass 
class MolecularFormulaSearchSettings:
    #corilo@molform
    #url_database = 'postgresql://postgres:labthomson0102@172.22.113.27:5432/' 
    url_database: str = None#'sqlite://'

    use_isotopologue_filter: bool = False

    isotopologue_filter_threshold:float = 33 # percentile predicted and found 

    isotopologue_filter_atoms:tuple = ('Cl', 'Br')

    use_runtime_kendrick_filter:bool = False

    use_min_peaks_filter:bool = True

    min_peaks_per_class:int = 15

    db_directory:str = False

    '''query setting'''
    ion_charge:int = -1

    hc_filter:float = 0.3

    oc_filter:float = 1.2

    op_filter:float = 2
    
    use_pah_line_rule:bool = False

    min_dbe:float = 0

    max_dbe:float = 50

    # look for close shell ions [M + Adduct]+ only considers metal set in the list adduct_atoms  
    adduct_atoms_neg:tuple = ('Cl', 'Br')
    
    adduct_atoms_pos:tuple = ('Na', 'K')

    score_methods:tuple = ('S_P_lowest_error', 'N_S_P_lowest_error', 'lowest_error', 'prob_score',
                     'air_filter_error', 'water_filter_error', 'earth_filter_error' )
    
    score_method:str = 'prob_score'

    # depending on the polarity mode it looks for [M].+ , [M].-
    # query and automatically compile add entry if it doesn't exist
    
    isRadical:bool = True
    
    # depending on the polarity mode it looks for [M + H]+ , [M - H]+
     # query and automatically compile and push options if it doesn't exist
    isProtonated:bool = True
    
    #kendrick_base: {'C': 1, 'H':2} = field(default_factory=dict)
    usedAtoms: dict = field(default_factory=dict)
    
    ''' search setting '''
    
    isAdduct:bool = True

    ionization_type:str = 'ESI'

    # empirically set / needs optimization
    min_ppm_error:float   = -5 #ppm

    # empirically set / needs optimization    
    max_ppm_error:float = 5 #ppm

    # empirically set / needs optimization set for isotopologue search
    min_abun_error:float = -100 # percentage 
    
    # empirically set / needs optimization set for isotopologue search
    max_abun_error:float = 100 # percentage 

    # empirically set / needs optimization
    mz_error_range:float = 1.5

    # 'distance', 'lowest', 'symmetrical','average' 'None'
    error_method:str = 'None'

    mz_error_average:float = 0

    #kendrick_base: {'C': 1, 'H':2} = field(default_factory=dict)
    used_atom_valences: dict = field(default_factory=dict)

    def __post_init__(self):
        
        self.usedAtoms = {   'C': (1, 100),
                    'H': (4, 200),
                    'O': (1, 26),
                    'N': (0, 0),
                    'S': (0, 1),
                    'P': (0, 0),
                    'Cl': (0, 0),
                }

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
                            'K': 1,
                            }

@dataclass
class DataInputSetting:
    
    #add to this dict the VALUES to match your labels, THE ORDER WON"T MATTER
    #"column_translate" : {"m/z":"m/z", "Resolving Power":"Resolving Power", "Abundance":"Abundance" , "S/N":"S/N"}
    header_translate: dict = field(default_factory=dict)

    def __post_init__(self):
        
        self.header_translate = {'m/z': Labels.mz, 
                        "Resolving Power":"Resolving Power",
                        "Res.":Labels.rp, 
                        'I':Labels.abundance,
                        "Abundance":"Abundance",
                        "Signal/Noise":"S/N",
                        "S/N":"S/N"}

                            
if __name__ == "__main__":
    a = DataInputSetting()
    print(a.__dict__)
    
    