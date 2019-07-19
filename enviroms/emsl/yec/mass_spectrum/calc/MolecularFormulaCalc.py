'''
Created on Jun 24, 2019

@author: eber373
'''
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms

class MolecularFormulaCalc:
    
    def _calc_theoretical_mz(self):
        
        if self.ion_charge:
            
            mass = 0
            
            for each_atom in self._d_molecular_formula.keys() :
                
                if each_atom != "IonType" and each_atom != 'HC':
                    
                    mass = mass + Atoms.atomic_masses[each_atom]  *  self._d_molecular_formula.get(each_atom)
                    
            return mass + ((-self.ion_charge) * Atoms.electron_mass)
        
        else:
            
            raise Exception("Please set ion charge first")
         
    def _calc_assigment_mass_error(self, exp_mz, method='ppm'):
        '''methodo should be ppm or ppb'''
        
        Hum_Milhao = 1000000
        Hum_Bilhao = 1000000000        
        
        if method == 'ppm':
            mult_factor = Hum_Milhao
        
        elif method == 'ppb':
            mult_factor = Hum_Bilhao
        
        else:
            raise Exception("method needs to be ppm or ppb, you have entered %s" % method)
              
        if self.exp_mz:
            #self.parent need to have a MassSpecPeak associated with the MolecularFormula class
            return ((self.theoretical_mz - self.exp_mz)/self.theoretical_mz)*mult_factor
        
        else:
            
            raise Exception("Please set theoretical_mz first")    
        
    
    def _calc_dbe(self):
            
            individual_dbe = 0
            
            for atom in self._d_molecular_formula.keys():
                if atom != "IonType":
                    n_atom = int(self._d_molecular_formula.get(atom))
                    valencia = Atoms.atoms_valence.get(atom)
                    
                    if valencia is not None:
                        if valencia > 0 :
                            #print atom, valencia, n_atom, individual_dbe
                            individual_dbe = individual_dbe + (n_atom * (valencia - 2))
                    else:
                        continue
            
            return 1 + (0.5 * individual_dbe)
    
    
    