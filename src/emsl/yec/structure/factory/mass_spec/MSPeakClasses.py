'''
Created on Jun 12, 2019
'''


from emsl.yec.calc.mass_spec.MSPeakCalc import MassSpecPeakCalculation
from emsl.yec.encapsulation.settings.ProcessingSetting import MassSpecPeakSetting
from emsl.yec.structure.factory.mass_spec.MolecularFormulaClasses import MolecularFormula


__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class MSPeak(MassSpecPeakCalculation):
    '''
    classdocs
    '''
    def __init__(self, ion_charge, exp_mz, abundance, resolving_power, signal_to_noise, massspec_index, exp_freq=None):
        
        #needed to create the object
        self.ion_charge = ion_charge
        self.exp_mz = exp_mz
        self.mass = float(exp_mz)/float(ion_charge)
        self.exp_freq = exp_freq
        self.abundance = abundance
        self.resolving_power = resolving_power 
        self.signal_to_noise = signal_to_noise 
        self.mass_spec_index = massspec_index
        
        kendrick_dict_base = MassSpecPeakSetting.kendrick_base
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(kendrick_dict_base)
        
        self.baseline_noise = None
        
        'updated after calibration'
        self.recal_mz = None
        
        'updated after molecular formula ID'
        
        self._molecular_formula= None
        self._confidence_score = None
        
        self.isotopologue_indexes = []
        self.found_isotopologues = {}
        
    def change_kendrick_base(self, kendrick_dict_base):
        '''kendrick_dict_base = {"C": 1, "H": 2}'''
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(kendrick_dict_base)
        
    @property    
    def kdm(self): return self._kdm
        
    @property    
    def kendrick_mass(self): return self._kendrick_mass
        
    @property    
    def nominal_km(self): return self._nominal_km
    
    @property
    def molecular_formula(self): return self._molecular_formula
    
    @property
    def confidence_score(self): return self._confidence_score
    
    @property    
    def is_assigned(self):
        
        return bool(self.molecular_formula)    
    
    @molecular_formula.setter
    def molecular_formula(self, formula_dict):
        
        self._molecular_formula = MolecularFormula(formula_dict, self.ion_charge, self.exp_mz)
        
        #self._calc_theoretical_mz(formula_dict)
        
        #self._calc_dbe(formula_dict)
        
        #self._set_ion_type()
    @confidence_score.setter
    def confidence_score(self): return self._calc_confidence_score() 
        
    def add_found_isotopologue(self, object_found):
        
        self.isotopologue_indexes.append(object_found.mass_spec_index)
        
        print(object_found)
        
        if self.found_isotopologues.has_key(object_found.qnt_c13):
            
            self.found_isotopologues[object_found.qnt_c13].append([object_found])
            
        else:    
            
            self.founded_isotopologos[object_found.qnt_c13] = [object_found] 
        
            