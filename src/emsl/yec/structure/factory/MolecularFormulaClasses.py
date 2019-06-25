'''
Created on Jun 24, 2019

@author: Yuri E. corilo 
'''
from emsl.yec.calc.MolecularFormulaCalc import MolecularFormulaCalc


class MolecularFormula(MolecularFormulaCalc):
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, ion_charge, exp_mz):
        
        self._d_molecular_formula = _d_molecular_formula
        self._dbe = self._calc_dbe()
        self._ion_type = self._set_ion_type()
        self._ion_chage = ion_charge
        self._theoretical_mz = self._calc_theoretical_mz()
        self._assigment_mass_error = self._set_assigment_mass_error(exp_mz)
        self.is_isotopologue = False    
        
        
    @property
    def dbe(self): return self._dbe
    
    @property
    def ion_charge(self): return self._ion_chage
         
    @property
    def theoretical_mz(self): return self._theoretical_mz
    
    @property    
    def assigment_mass_error(self): return self._assigment_mass_error
        
    @property
    def ion_type(self): return self._ion_type    
    
    @property
    def d_molecular_formula(self): return self._d_molecular_formula         
    
    @ion_type.setter
    def ion_type(self): return self._set_ion_type()
    
    @theoretical_mz.setter     
    def theoretical_mz(self): return self._calc_theoretical_mz()
    
    @dbe.setter
    def dbe(self): return self._calc_dbe()
    
    @assigment_mass_error.setter     
    def assigment_mass_error(self, exp_mz): return self._set_assigment_mass_error(exp_mz)