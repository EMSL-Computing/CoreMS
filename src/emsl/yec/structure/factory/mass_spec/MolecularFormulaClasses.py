'''
Created on Jun 24, 2019

@author: Yuri E. corilo 
'''
from emsl.yec.calc.mass_spec.MolecularFormulaCalc import MolecularFormulaCalc


class MolecularFormula(MolecularFormulaCalc):
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, ion_charge, exp_mz):
        
        self._d_molecular_formula = _d_molecular_formula
        self._ion_chage = ion_charge
        self._assigment_mass_error = self._set_assigment_mass_error(exp_mz)
        self.is_isotopologue = False    
        
    @property
    def O_C(self): return self._d_molecular_formula.get("O")/self._d_molecular_formula.get("C")
    
    @property
    def H_C(self): return self._d_molecular_formula.get("H")/self._d_molecular_formula.get("C")
    
    @property
    def dbe(self): return self._calc_dbe()
    
    @property
    def ion_charge(self): return self._ion_chage
         
    @property
    def theoretical_mz(self): return self._calc_theoretical_mz()
    
    @property    
    def assigment_mass_error(self): return self._assigment_mass_error
        
    @property
    def ion_type(self): return self._set_ion_type()
    
    @property
    def atoms(self): return self._d_molecular_formula.keys()
             
    @property
    def atoms_qnt(self,atom): 
        if atom in self._d_molecular_formula:
            return self._d_molecular_formula.get(atom)
        else:
            raise Exception('Could not find %s in this Molecular Formula object'%str(atom))
    
    def _set_ion_type(self):
        
        if self.d_molecular_formula:
            
            return self._d_molecular_formula.get("IonType")
        
        else:
            
            raise Exception("Please set Molecular_fromula_first")  
   