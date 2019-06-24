'''
Created on Jun 12, 2019
'''

from emsl.yec.calc.MSPeakCalc import MassSpecPeakCalculation

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class MSPeak(MassSpecPeakCalculation):
    '''
    classdocs
    '''
    def __init__(self, ion_charge, exp_mz, magnitude, resolving_power, signal_to_noise, massspec_index, exp_freq=None):
        
        #needed to create the object
        self.ion_charge = ion_charge
        self.exp_mz = exp_mz
        self.mass = exp_mz/ion_charge
        self.exp_freq = exp_freq
        self.magnitude = magnitude
        self.resolving_power = resolving_power 
        self.signal_to_noise = signal_to_noise 
        self.mass_spec_index = massspec_index
        
        self._kdm = None
        self._kendrick_mass = None
        self._nominal_km = None
       
        self.set_kmd_and_nominal()
        
        self.baseline_noise = None
        self.recal_mz = None
        self._ion_type = None
        
        self._theoretical_mz = None
        self._d_molecular_formula = None
        self._dbe = None
        self._assigment_mass_error = None 
        self.confidence_score = None
        
        #updated after molecular formula ID
        self.is_isotopologue = False    
        self.isotopologue_indexes = []
        self.found_isotopologues = {}
    
    @property    
    def kdm(self): return self._kdm
        
    @property    
    def kendrick_mass(self): return self._kendrick_mass
        
    @property    
    def nominal_km(self): return self._nominal_km
        
    @property    
    def mass_error(self): return self._mass_error
        
    @property
    def ion_type(self): return self._ion_type
        
    @property
    def d_molecular_formula(self): return self._d_molecular_formula
        
    @property
    def theoretical_mz(self): return self._theoretical_mz
        
    @property
    def dbe(self): return self._dbe
    
    @property    
    def is_assigned(self):
        
        return bool(self.d_molecular_formula)    
    
    @property.setter
    def set_ion_type(self):
        
        if self.d_molecular_formula:
            
            self._ion_type = self.d_molecular_formula.get("IonType")
        
        else:
            
            raise Exception("Please set ion type first")            
    
    
    @property.setter
    def set_molecular_formula(self, formula_dict):
        
        self._d_molecular_formula = formula_dict
        
        self.__calc_theoretical_mz(formula_dict)
        
        self.__calc_dbe(formula_dict)
        
        self.set_ion_type()
    
    @property
    def molecular_formula_list(self):
        
        atoms_in_ordem = ["C", "H", "N", "O", "S"]

        formula_list = []
        
        for atomo in atoms_in_ordem:

            numero_atom = self.d_molecular_formula.get(atomo)

            if numero_atom:
                
                formula_list.append(atomo)
                formula_list.append(numero_atom)
            #else:
            #    formula_list_zero_filled.append(atomo)
            #    formula_list_zero_filled.append(0)

        atomos_in_dict = self.d_molecular_formula.keys()
        for atomo in atomos_in_dict:

            if atomo not in atoms_in_ordem and atomo != "IonType" and atomo != "HC":
                
                formula_list.append(atomo)
                formula_list.append(self.d_molecular_formula.get(atomo))

        return formula_list
     
    @property 
    def molecular_formula_string(self):
        
        formulalist = self.get_molecular_formula_list()
        formulastring = "" 
        
        for each in range(0, len(formulalist),2):
            
            formulastring = formulastring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
            
        return formulastring[0:-1]   
    
    @property
    def heteroatomic_class_label(self):
        
        if self.d_molecular_formula:
            
            formulalist = self.get_molecular_formula_list()
            classstring = "" 
            
            for each in range(0, len(formulalist),2):
                
                if formulalist[each] != 'C' and formulalist[each] != 'H':
                     
                    classstring = classstring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
            
            if classstring == "":
                
                classstring = "HC "
                
            if self.d_molecular_formula.get("IonType") == 'RADICAL':    
                
                return classstring[0:-1] + " " + "-R"
            
            else:
                
                return classstring[0:-1] + " "
            
            'dict, tuple or string'
        else:
            return None    
        
    def add_found_isotopologue(self, object_found):
        
        self.isotopologue_indexes.append(object_found.mass_spec_index)
        
        print(object_found)
        
        if self.found_isotopologues.has_key(object_found.qnt_c13):
            
            self.found_isotopologues[object_found.qnt_c13].append([object_found])
            
        else:    
            
            self.founded_isotopologos[object_found.qnt_c13] = [object_found] 
        
            