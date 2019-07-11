'''
Created on Jun 24, 2019

@author: eber373
'''
from enviroms.emsl.yec.encapsulation.constant.Constants import Constants


class MolecularFormulaCalc(object):
    
    def _calc_theoretical_mz(self):
        
        if self.ion_charge:
            
            mass = 0
            
            for each_atom in self._d_molecular_formula.keys() :
                
                if each_atom != "IonType" and each_atom != 'HC':
                    
                    mass = mass + Constants.atomic_masses[each_atom]  *  self._d_molecular_formula.get(each_atom)
                    
            return mass + ((-self.ion_charge) * Constants.electron_mass)
                    
                    
        else:
            
            raise Exception("Please set ion charge first")
         
    def _set_assigment_mass_error(self, exp_mz, method='ppm'):
        '''methodo should be ppm or ppb'''
        
        Hum_Milhao = 1000000
        Hum_Bilhao = 1000000000        
        
        if method == 'ppm':
            mult_factor = Hum_Milhao
        
        elif method == 'ppb':
            mult_factor = Hum_Bilhao
        
        else:
            raise Exception("method needs to be ppm or ppb, you have entered %s" % method)
              
        if self.theoretical_mz:
            #self.parent need to be the MassSpecPeak associated with the MF class
            return ((self.theoretical_mz - exp_mz)/self.theoretical_mz)*mult_factor
        
        else:
            
            raise Exception("Please set theoretical_mz first")    
        
    
    def _calc_dbe(self):
            
            individual_dbe = 0
            
            for atom in self._d_molecular_formula.keys():
                if atom != "IonType":
                    n_atom = int(self._d_molecular_formula.get(atom))
                    valencia = Constants.atoms_valence.get(atom)
                    if valencia > 0 and valencia is not None:
                        #print atom, valencia, n_atom, individual_dbe
                        individual_dbe = individual_dbe + (n_atom * (valencia - 2))
                    else:
                        continue
            
            return 1 + (0.5 * individual_dbe)
           
    def to_string(self):
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list()
            
            formulastring = "" 
            
            for each in range(0, len(formulalist),2):
                
                formulastring = formulastring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
                
            return formulastring[0:-1]   
       
        else:
            
            raise Exception("Molecular formula identification not performed yet")    
    
    def to_list(self):
        
        if self._d_molecular_formula:
            atoms_in_ordem = ["C", "H", "N", "O", "S", "P"]
    
            formula_list = []
            
            for atomo in atoms_in_ordem:
    
                numero_atom = self._d_molecular_formula.get(atomo)
    
                if numero_atom:
                    
                    formula_list.append(atomo)
                    formula_list.append(numero_atom)
                #else:
                #    formula_list_zero_filled.append(atomo)
                #    formula_list_zero_filled.append(0)
    
            atomos_in_dict = self._d_molecular_formula.keys()
            for atomo in atomos_in_dict:
    
                if atomo not in atoms_in_ordem and atomo != "IonType" and atomo != "HC":
                    
                    formula_list.append(atomo)
                    formula_list.append(self._d_molecular_formula.get(atomo))
    
            return formula_list
        else:
            raise Exception("Molecular formula identification not performed yet")
        
    @property
    def heteroatomic_class_label(self):
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list()
            classstring = "" 
            
            for each in range(0, len(formulalist),2):
                
                if formulalist[each] != 'C' and formulalist[each] != 'H':
                     
                    classstring = classstring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
            
            if classstring == "":
                
                classstring = "HC "
                
            if self._d_molecular_formula.get("IonType") == 'RADICAL':    
                
                return classstring[0:-1] + " " + "-R"
            
            else:
                
                return classstring[0:-1] + " "
            
            'dict, tuple or string'
        else:
            raise Exception("Molecular formula identification not performed yet")    