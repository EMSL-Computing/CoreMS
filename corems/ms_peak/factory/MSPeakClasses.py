
from corems.encapsulation.settings.input.ProcessingSetting import MassSpecPeakSetting
from corems.ms_peak.calc.MSPeakCalc import MSPeakCalculation


__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class _MSPeak(MSPeakCalculation):
    '''
    classdocs
    '''
    def __init__(self, ion_charge, mz_exp, abundance, resolving_power, signal_to_noise, massspec_index, index, exp_freq=None):

        # needed to create the object
        self.ion_charge = int(ion_charge)
        self._mz_exp = float(mz_exp)
        self.mass = float(mz_exp) / float(ion_charge)
        self.abundance = float(abundance)
        self.resolving_power = float(resolving_power)
        self.signal_to_noise = float(signal_to_noise)
        self.mass_spec_index = int(massspec_index)
        self.index = int(index)
        'updated after calibration'
        self.mz_cal = None
        'updated individual calculation'
        self.baseline_noise = None
        
        if exp_freq:
            self.freq_exp = float(exp_freq)

        kendrick_dict_base = MassSpecPeakSetting.kendrick_base
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)
 
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
    def mz_exp(self):
        if self.mz_cal:
            return self.mz_cal
        else:
            return self._mz_exp
    
    @mz_exp.setter
    def mz_exp(self, mz_exp):
        self._mz_exp = mz_exp

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

class ICRMassPeak(_MSPeak):

    def __init__(self, *args, exp_freq=None):

        super().__init__(*args,exp_freq=exp_freq)

    def resolving_power_calc(self, B, T):
        
        '''
        low pressure limits, 
        T: float 
            transient time
        B: float
            Magnetic Filed Strength (Tesla)    
        
        reference
        Marshall et al. (Mass Spectrom Rev. 1998 Jan-Feb;17(1):1-35.)
        DOI: 10.1002/(SICI)1098-2787(1998)17:1<1::AID-MAS1>3.0.CO;2-K
        '''
        return (1.274e7 * self.ion_charge * B * T)/ (self.mz_exp*self.ion_charge)

    def set_threoretical_resolving_power(self, B, T):

        self.resolving_power = self.resolving_power_calc(B, T) 
        
class TOFMassPeak(_MSPeak):

    def __init__(self, *args, exp_freq=None):

        super().__init__(*args,exp_freq=exp_freq)

    def set_threoretical_resolving_power(self):
        return 0

class OrbiMassPeak(_MSPeak):

    def __init__(self, *args, exp_freq=None):

        super().__init__(*args,exp_freq=exp_freq)

    def set_threoretical_resolving_power(self):
        return 0       
