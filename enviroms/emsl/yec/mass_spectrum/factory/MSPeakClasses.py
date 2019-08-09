
from enviroms.emsl.yec.encapsulation.settings.input.ProcessingSetting import \
    MassSpecPeakSetting
from enviroms.emsl.yec.mass_spectrum.calc.MSPeakCalc import \
    MassSpecPeakCalculation
from enviroms.emsl.yec.molecular_id.factory.MolecularFormulaFactory import \
    MolecularFormula

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class MSPeak(MassSpecPeakCalculation):
    '''
    classdocs
    '''
    def __init__(self, ion_charge, mz_exp, abundance, resolving_power, signal_to_noise, massspec_index, index, exp_freq=None):

        # needed to create the object
        self.ion_charge = ion_charge
        self.mz_exp = mz_exp
        self.mass = float(mz_exp) / float(ion_charge)
        self.freq_exp = exp_freq
        self.abundance = abundance
        self.resolving_power = resolving_power
        self.signal_to_noise = signal_to_noise
        self.mass_spec_index = massspec_index
        self.index = index
        
        kendrick_dict_base = MassSpecPeakSetting.kendrick_base
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)

        self.baseline_noise = None

        'updated after calibration'
        self.mz_recal = None

        'updated after molecular formula ID'

        self.molecular_formulas = []
        self._confidence_score = None

        self.isotopologue_indexes = []
        self.found_isotopologues = {}

    def __len__(self):
        
        return len(self.molecular_formulas)
        
    def __setitem__(self, position, molecular_formula_obj):
        
        self.molecular_formulas[position] = molecular_formula_obj

    def __getitem__(self, position):
        
        return self.molecular_formulas[position]

    def change_kendrick_base(self, kendrick_dict_base):
        '''kendrick_dict_base = {"C": 1, "H": 2}'''
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)

    def add_molecular_formula(self, molecular_formula_obj):

        self.molecular_formulas.append(molecular_formula_obj)
    
    def clear_molecular_formulas(self):
        
        self.molecular_formulas= []

    @property
    def nominal_mz_exp(self): return int(self.mz_exp)

    @property
    def kmd(self): return self._kdm

    @property
    def kendrick_mass(self): return self._kendrick_mass

    @property
    def knm(self): return self._nominal_km
    
    @property
    def is_assigned(self):

        return bool(self.molecular_formulas)
    
    @property
    def number_possible_assigments(self,):
        
        return len(self.molecular_formulas)
    
    @property
    def molecular_formula_lowest_error(self):
        
         return min(self.molecular_formulas, key=lambda m: abs(m._calc_assigment_mass_error(self.mz_exp)))



class ICRMassPeak(MSPeak):

    def __init__(self, *args):

        super().__init__(*args)

    def threoretical_resolving_power(self):
        return 0

class TOFMassPeak(MSPeak):

    def __init__(self, *args):

        super().__init__(*args)

    def threoretical_resolving_power(self):
        return 0

class OrbiMassPeak(MSPeak):

    def __init__(self, *args):

        super().__init__(*args)

    def threoretical_resolving_power(self):
        return 0       
