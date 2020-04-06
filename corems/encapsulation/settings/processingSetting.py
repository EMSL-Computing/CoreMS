__author__ = 'Yuri E. Corilo'
__date__ = 'Jul 02, 2019'

class TransientSetting:
    
    implemented_apodization_function = ('Hamming', 'Hanning', 'Blackman')
    apodization_method = 'Hanning'
    number_of_truncations = 0
    number_of_zero_fills = 1
            
class MassSpectrumSetting:
    
    threshold_method = 'auto'
    
    implemented_noise_threshold_methods = ('auto', 'signal_noise', 'relative_abundance')
    
    noise_threshold_std = 6
    
    s2n_threshold = 4
    
    relative_abundance_threshold = 6 # from 1-100
    
    min_noise_mz = 100.0
    max_noise_mz = 1200.0
    
    min_picking_mz = 100.0
    max_picking_mz = 1200.0
    
class MassSpecPeakSetting:
    
    kendrick_base =  {'C': 1, 'H':2}
    
    peak_min_prominence_percent = 1 #1-100 % used for peak detection

    peak_max_prominence_percent = 0.1 #1-100 % used for baseline detection

class GasChromatographSetting:
    
    implemented_smooth_method = ['savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar']
    
    smooth_window = 5
    
    smooth_method = 'savgol'
    
    savgol_pol_order = 2
    
    peak_height_max_percent = 10 #1-100 % used for baseline detection

    peak_max_prominence_percent = 10 #1-100 % used for baseline detection

    min_peak_datapoints = 5
   
    max_peak_width = 0.1

    noise_threshold_method = 'auto'
    
    implemented_noise_threshold_methods = ('auto', 'relative_abundance')
    
    std_noise_threshold = 3

    peak_height_min_abun = 5 #1-100 % used for peak detection

    peak_min_prominence_percent = 1 #1-100 % used for peak detection
    

class CompoundSearchSettings:

    ri_search_range = 20

    rt_search_range = 0.5
    
    correlation_threshold = 0.95 # used for calibration, spectral similarity 
    
    score_threshold = 0.0

    ri_spacing = 200

    ri_window = 3 # in standard deviation

class MolecularLookup0DictSettings:
    
    '''
    ### DO NOT CHANGE IT! These are used to generate the database entries 

    ### DO change when creating a new application database 
    
    ### FOR search settings runtime and database query check use the MolecularSearchSettings class below

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
                    'O': (0, 4),
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

        #ion_charge changes automatically with mass spectrum
        self.ion_charge = -1

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
        

class MolecularSearchSettings:
    
    url_database = 'postgresql://postgres:labthomson0102@172.22.113.27:5432/' 
    #url_database = None#'sqlite://'

    use_isotopologue_filter = False

    isotopologue_filter_threshold = 33 # percentile predicted and found 

    isotopologue_filter_atoms = ['Cl', 'Br']

    use_runtime_kendrick_filter = False

    use_min_peaks_filter = True

    min_peaks_per_class = 15

    db_directory = False

    '''query setting'''
    ion_charge  = -1

    hc_filter = 0.3

    oc_filter = 1.2

    op_filter = 2
    
    use_pah_line_rule = False

    min_dbe = 0

    max_dbe = 50

    # look for close shell ions [M + Adduct]+ only considers metal set in the list adduct_atoms  
    adduct_atoms_neg = ['Cl', 'Br', 'F']
    
    adduct_atoms_pos = ['Na', 'K']

    score_methods = ['S_P_lowest_error', 'N_S_P_lowest_error', 'lowest_error', 'prob_score',
                     'air_filter_error', 'water_filter_error', 'earth_filter_error' ]
    
    score_method = 'N_S_P_lowest_error'

    # depending on the polarity mode it looks for [M].+ , [M].-
    # query and automatically compile add entry if it doesn't exist
    
    isRadical = True
    
    # depending on the polarity mode it looks for [M + H]+ , [M - H]+
     # query and automatically compile and push options if it doesn't exist
    isProtonated = True

    usedAtoms = {   'C': (1, 100),
                    'H': (4, 200),
                    'O': (1, 22),
                    'N': (0, 0),
                    'S': (0, 0),
                    'P': (0, 0),
                    'Cl': (0, 0),
                }
    
    ''' search setting '''
    
    isAdduct = True

    ionization_type = 'ESI'

    # empirically set / needs optimization
    min_ppm_error  = -5 #ppm

    # empirically set / needs optimization    
    max_ppm_error = 5 #ppm

    # empirically set / needs optimization set for isotopologue search
    min_abun_error = -100 # percentage 
    
    # empirically set / needs optimization set for isotopologue search
    max_abun_error = 100 # percentage 

    # empirically set / needs optimization
    mz_error_range = 1.5

    # 'distance', 'lowest', 'symmetrical','average' 'None'
    error_method = 'None'

    mz_error_average = 0

    used_atom_valences = {'C': 4,
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
