__author__ = 'Yuri E. Corilo'
__date__ = 'Jul 02, 2019'

import dataclasses
import os
from typing import List, Dict

from corems.encapsulation.constant import Atoms, Labels

@dataclasses.dataclass
class TransientSetting:

    implemented_apodization_function: tuple = ('Hamming', 'Hanning', 'Blackman')
    apodization_method: str = 'Hanning'
    number_of_truncations: int = 0
    number_of_zero_fills: int = 1

    def __post_init__(self):

        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class DataInputSetting:

    # add to this dict the VALUES to match your labels, THE ORDER WON'T MATTER
    # 'column_translate' : {'m/z':'m/z', 'Resolving Power':'Resolving Power', 'Abundance':'Abundance' , 'S/N':'S/N'}
    header_translate: dict = dataclasses.field(default_factory=dict)

    def __post_init__(self):

        self.header_translate = {'m/z': Labels.mz, 
                                 'mOz': Labels.mz,
                                 'Mass': Labels.mz,
                                 'Resolving Power': 'Resolving Power',
                                 'Res.': Labels.rp,
                                 'Intensity': Labels.abundance,
                                 'I': Labels.abundance,
                                 'Abundance': Labels.abundance,
                                 'Signal/Noise': Labels.s2n,
                                 'S/N': Labels.s2n,
                                 'abs_abu': Labels.abundance,
                                 'sn': Labels.s2n,
                                 'resolution': Labels.rp}

    def add_mz_label(self, label):

        self.header_translate[label] = Labels.mz

    def add_peak_height_label(self, label):

        self.header_translate[label] = Labels.abundance

    def add_sn_label(self, label):

        self.header_translate[label] = Labels.s2n

    def add_resolving_power_label(self, label):

        self.header_translate[label] = Labels.rp

@dataclasses.dataclass
class LiqChromatographSetting:
    start_scan: int = 1
    final_scan: int = 7

    def __post_init__(self):
        # enforce datatype
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, field.type):

                value = field.type(value)
                setattr(self, field.name, value)

@dataclasses.dataclass
class MassSpectrumSetting:

    threshold_method: str = 'relative_abundance'

    implemented_noise_threshold_methods: tuple = ('auto', 'signal_noise', 'relative_abundance')

    noise_threshold_std: int = 6

    s2n_threshold: float = 4

    relative_abundance_threshold: float = 6  # from 0-100

    min_noise_mz: float = 100.0
    max_noise_mz: float = 1200.0

    min_picking_mz: float = 100.0
    max_picking_mz: float = 1200.0

    calib_minimize_method: str = 'Powell'
    calib_pol_order: int = 2
    max_calib_ppm_error: float = 1
    min_calib_ppm_error: float = -1
    calib_sn_threshold: float = 10

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

    kendrick_base: Dict = dataclasses.field(default_factory=dict)

    # kendrick_base : Dict =  {'C': 1, 'H':2}

    peak_min_prominence_percent: float = 1  # 1-100 % used for peak detection

    peak_max_prominence_percent: float = 0.1  # 1-100 % used for baseline detection

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

    use_deconvolution: bool = False

    implemented_smooth_method: tuple = ('savgol', 'hanning', 'blackman', 'bartlett', 'flat', 'boxcar')

    smooth_window: int = 5

    smooth_method: str = 'savgol'

    savgol_pol_order: int = 2

    peak_height_max_percent: float = 10  # 1-100 % used for baseline detection use 0.1 for second_derivative and 10 for other methods

    peak_max_prominence_percent: float = 1  # 1-100 % used for baseline detection

    min_peak_datapoints: float = 5

    max_peak_width: float = 0.1

    noise_threshold_method: str = 'manual_relative_abundance'

    implemented_noise_threshold_methods: tuple = ('auto_relative_abundance', 'manual_relative_abundance', 'second_derivative')

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

    url_database: str = 'sqlite:///db/pnnl_lowres_gcms_compounds.sqlite'

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

        self.ri_calibration_compound_names = (' [C8] Methyl Caprylate [7.812]',
                                ' [C10] Methyl Caprate [10.647]',
                                ' [C9] Methyl Pelargonate [9.248]',
                                ' [C12] Methyl Laurate [13.250]',
                                ' [C14] Methyl Myristate [15.597]',
                                ' [C16] Methyl Palmitate [17.723]',
                                ' [C18] Methyl Stearate [19.663]',
                                ' [C20] Methyl Eicosanoate [21.441]',
                                ' [C22] Methyl Docosanoate [23.082]',
                                ' [C24] Methyl Linocerate [24.603]',
                                ' [C26] Methyl Hexacosanoate [26.023]',
                                ' [C28] Methyl Octacosanoate [27.349]',
                                ' [C30] Methyl Triacontanoate [28.72]')

@dataclasses.dataclass
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

        self.usedAtoms = {'C': (1, 90),
                          'H': (4, 200),
                          'O': (0, 12),
                          'N': (0, 0),
                          'S': (0, 0),
                          'P': (0, 0),
                          'Cl': (0, 0),
                          }

        self.min_mz = 100

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

    use_isotopologue_filter: bool = False

    isotopologue_filter_threshold: float = 33

    isotopologue_filter_atoms: tuple = ('Cl', 'Br')

    use_runtime_kendrick_filter: bool = False

    use_min_peaks_filter: bool = True

    min_peaks_per_class: int = 15

    url_database: str = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'

    db_jobs: int = 3

    '''query setting'''
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

    # depending on the polarity mode it looks for [M].+ , [M].-
    # query and automatically compile add entry if it doesn't exist

    isRadical: bool = False

    # depending on the polarity mode it looks for [M + H]+ , [M - H]+
    # query and automatically compile and push options if it doesn't exist
    isProtonated: bool = True

    isAdduct: bool = False

    usedAtoms: dict = dataclasses.field(default_factory=dict)

    ''' search setting '''

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
                self.used_atom_valences[atom] = Atoms.atoms_covalence.get(atom)

        '''
        self.used_atom_valences.update({'13C': 4,
                                        '18O': 2,
                                        '34S': 2,
                                        'Cl': 1,
                                        '37Cl': 1,
                                        'Br': 1,
                                        'Na': 1,
                                        'F': 1,
                                        'K': 1,
                                        })
        '''

if __name__ == '__main__':
    a = DataInputSetting()
    print(a.__dict__)
