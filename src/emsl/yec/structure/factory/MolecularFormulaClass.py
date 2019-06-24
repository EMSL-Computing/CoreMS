'''
Created on Jun 24, 2019

@author: Yuri E. corilo 
'''
from emsl.yec.calc import MolecularFormulaCalc


class MolecularFormula(MolecularFormulaCalc):
    '''
    classdocs
    '''

    def __init__(self, _d_molecular_formula):
        
        self._d_molecular_formula = _d_molecular_formula
        self._theoretical_mz = self._calc_theoretical_mz(_d_molecular_formula)
        self._assigment_mass_error = self._set_assigment_mass_error()
        self._dbe = self._calc_dbe()
        self._ion_type = self._set_ion_type()
        
        self.is_isotopologue = False    
    
    @property
    def dbe(self): return self._dbe
         
    @property
    def theoretical_mz(self): return self._theoretical_mz
    
    @property    
    def assigment_mass_error(self): return self._assigment_mass_error
        
    @property
    def ion_type(self): return self._ion_type    
    
    @property
    def d_molecular_formula(self): return self._d_molecular_formula         
