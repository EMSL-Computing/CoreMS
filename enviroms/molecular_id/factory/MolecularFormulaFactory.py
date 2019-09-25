from copy import deepcopy
from collections import OrderedDict
from enviroms.mass_spectrum.calc.MolecularFormulaCalc import MolecularFormulaCalc
from enviroms.encapsulation.settings.input.ProcessingSetting import MassSpecPeakSetting
from enviroms.encapsulation.constant import Atoms, Labels



__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

class MolecularFormula(MolecularFormulaCalc):
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, ion_charge, exp_mz=None):
        
        #clear dictionary of atoms with 0 value
        
        self._d_molecular_formula = {key:val for key, val in _d_molecular_formula.items() if val != 0}

        self._ion_charge = ion_charge
        self._assigment_mass_error = None
        self._confidence_score = None
        self.is_isotopologue = False
        self.mspeak_indexes_isotopologues = []
        
        kendrick_dict_base = MassSpecPeakSetting.kendrick_base
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)  

        if exp_mz:
            self._assigment_mass_error = self._calc_assigment_mass_error(exp_mz)
            #self._confidence_score = self._calc_confidence_score()     
    
    def __len__(self):
        
        #crash if keys are not ordered
            return len(self._d_molecular_formula.keys())
        
    def __getitem__(self, atom):
        
            #atom = list(self._d_molecular_formula.keys())[position]
            return self._d_molecular_formula[atom]

    @property
    def O_C(self): 
            
            if 'O' in self._d_molecular_formula.keys():
                return self._d_molecular_formula.get("O")/self._d_molecular_formula.get("C")
            else:
                return None    
    
    @property
    def H_C(self): return self._d_molecular_formula.get("H")/self._d_molecular_formula.get("C")
    
    @property
    def dbe(self): return self._calc_dbe()
    
    @property
    def mz_nominal_theo(self): return int(self._calc_mz_theor())

    @property    
    def mz_error(self): return self._assigment_mass_error
    
    @property
    def mz_theor(self): return self._calc_mz_theor()

    @property
    def ion_type(self): 
        ion_type = self._d_molecular_formula.get(Labels.ion_type)
        if ion_type == Labels.protonated_de_ion:
            if self.ion_charge > 0: 
                return Labels.protonaded
            else: 
                return Labels.de_protonated    
        else:
            return ion_type

    @ion_type.setter
    def ion_type(self, ion_type):
        if  ion_type in [Labels.protonated_de_ion, Labels.adduct_ion, Labels.radical_ion]:
            self._d_molecular_formula[Labels.ion_type] = ion_type
        else:
            raise TypeError("Ion type can only be: 'DE_OR_PROTONATED', 'RADICAL' or  'ADDUCT', not %s"%ion_type)   

    @property
    def ion_charge(self): return self._ion_charge
    
    @property
    def atoms(self): return [key for key in self._d_molecular_formula.keys() if key != Labels.ion_type]
    
    @property
    def confidence_score(self): return self._calc_confidence_score() 
        
    @property
    def kmd(self): return self._kdm

    @property
    def kendrick_mass(self): return self._kendrick_mass

    @property
    def knm(self): return self._nominal_km

    def change_kendrick_base(self, kendrick_dict_base):
        '''kendrick_dict_base = {"C": 1, "H": 2}'''
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)
                
    def isotopologues(self, min_abudance, current_abundance): 
        
        for mf in self._cal_isotopologues(self._d_molecular_formula, min_abudance, current_abundance ):

            yield MolecularFormulaIsotopologue(*mf, self.ion_charge)
    
    def atoms_qnt(self,atom): 
        if atom in self._d_molecular_formula:
            return self._d_molecular_formula.get(atom)
        else:
            raise Warning('Could not find %s in this Molecular Formula object'%str(atom))
    
    def atoms_symbol(self, atom): 
        '''return the atom symbol without the mass number'''
        return ''.join([i for i in atom if not i.isdigit()])

    @property       
    def to_string(self):
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list
            
            formulastring = "" 
            
            for each in range(0, len(formulalist),2):
                
                formulastring = formulastring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
                
            return formulastring[0:-1]   
       
        else:
            
            raise Exception("Molecular formula identification not performed yet")    
    
    def to_string_formated(self):
        
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

        if self._d_molecular_formula:
            formula_srt = ''
            for atom in Atoms.atoms_order:
                if atom in self.to_dict.keys():
                    formula_srt += atom.translate(SUP) + str(int(self.to_dict.get(atom))).translate(SUB)
            return formula_srt
        
        else:
            raise Exception("Molecular formula identification not performed yet")    

    @property
    def to_dict(self):
        return self._d_molecular_formula
    
    @property   
    def to_list(self):
        '''TODO ensure self._d_molecular_formula is a orderedDict'''
        
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
    
                if atomo not in atoms_in_ordem and atomo != Labels.ion_type:
                    
                    formula_list.append(atomo)
                    formula_list.append(self._d_molecular_formula.get(atomo))
    
            return formula_list
        else:
            raise Exception("Molecular formula identification not performed yet")
        
    @property
    def class_label(self):
        
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list
            classstring = "" 
            
            for each in range(0, len(formulalist),2):
                
                if formulalist[each] != 'C' and formulalist[each] != 'H':
                     
                    classstring = classstring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
            
            if classstring == "": classstring = "HC "
                
            classstring = classstring.strip()
            if self._d_molecular_formula.get(Labels.ion_type) == 'RADICAL':    
                
                return classstring + " " + "-R"
            
            else: return classstring
            
            'dict, tuple or string'
        else:
            raise Exception("Molecular formula identification not performed yet")        
    
    @property
    def class_dict(self):
        
        if self._d_molecular_formula:
            
            class_dict = {}
            for atom, qnt in self._d_molecular_formula.items():
    
                if atom != Labels.ion_type and atom !='C' and atom !='H':
                    
                    class_dict[atom] = qnt
                    
            return class_dict
        
        raise Exception("Molecular formula identification not performed yet")           
    

class MolecularFormulaIsotopologue(MolecularFormula):
        
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, prob_ratio, ion_charge, exp_mz=None):
        
        super().__init__(_d_molecular_formula,  ion_charge)
        #prob_ratio is relative to the monoisotopic peak p_isotopologue/p_mono_isotopic
        self.prob_ratio = prob_ratio
        
        self.is_isotopologue = True
        self.mspeak_index_mono_isotopic = None

        if exp_mz:
            
            self._assigment_mass_error = self._calc_assigment_mass_error(exp_mz)
            self._confidence_score = False