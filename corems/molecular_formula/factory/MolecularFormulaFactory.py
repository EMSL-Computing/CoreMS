from copy import deepcopy
from collections import OrderedDict
from corems.molecular_formula.calc.MolecularFormulaCalc import MolecularFormulaCalc
from corems.encapsulation.constant import Atoms, Labels
from numpy import nan
import re

__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

class MolecularFormula(MolecularFormulaCalc):
    '''
    classdocs
    '''
    def __init__(self, molecular_formula, ion_charge, ion_type=None, adduct_atom=None, mspeak_parent=None):
        
        #clear dictionary of atoms with 0 value
        
        if   type(molecular_formula) is dict:
                self._from_dict(molecular_formula, ion_type, adduct_atom)   
        
        elif type(molecular_formula) is list:
                self._from_list(molecular_formula, ion_type, adduct_atom)   
        
        elif type(molecular_formula) is str:
                self._from_str(molecular_formula, ion_type, adduct_atom)   

        
        self._ion_charge = ion_charge

        self._confidence_score = None        
        self._isotopologue_similarity = None
        self._mz_error_score = None
        self._mass_error_average_score = None

        self.is_isotopologue = False

        # parent mass spectrum peak obj instance
        self._mspeak_parent = mspeak_parent

        self.expected_isotopologues = []
        self.mspeak_mf_isotopologues_indexes = []
        
        if self._mspeak_parent:
            kendrick_dict_base = self._mspeak_parent._ms_parent.kendrick_base
        else:
            kendrick_dict_base = {'C':1, 'H':2}
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)  
        
    def __repr__(self):

        return "MolecularFormula({0},{1},ion type = {2}".format(self._d_molecular_formula, self.ion_charge, self.ion_type)
    
    def __str__(self):

        return "MolecularFormula {0}, ion_charge:{1}, ion type:{2}, m/z:{3} ".format(self.string, self.ion_charge, self.ion_type, self.mz_calc)
    
    def __len__(self):
        
        #crash if keys are not ordered
            return len(self._d_molecular_formula.keys())
        
    def __getitem__(self, atom):
        
            #atom = list(self._d_molecular_formula.keys())[position]
            if atom in self._d_molecular_formula.keys():
                return self._d_molecular_formula[atom]
            else:
                return 0
    
    def get(self, atom):
        
            #atom = list(self._d_molecular_formula.keys())[position]
            if atom in self._d_molecular_formula.keys():
                return self._d_molecular_formula[atom]
            else:
                return 0
                
    def _from_dict(self, molecular_formula, ion_type, adduct_atom):
        
        self._d_molecular_formula = {key:val for key, val in molecular_formula.items() if val != 0}
        if ion_type:
            self._d_molecular_formula[Labels.ion_type] = ion_type
        if adduct_atom:
            if adduct_atom in self._d_molecular_formula:
                self._d_molecular_formula[adduct_atom] += 1 
            else: self._d_molecular_formula[adduct_atom] = 1 

        
    @property
    def isotopologue_count_percentile(self, ):
        if not len(self.expected_isotopologues) == 0:
            return (len(self.mspeak_mf_isotopologues_indexes)/len(self.expected_isotopologues))*100
        else: 
            return 100

    def _from_list(self, molecular_formula_list, ion_type, adduct_atom):
        # list has to be in the format 
        #['C', 10, 'H', 21, '13C', 1, 'Cl', 1, etc]  
        self._d_molecular_formula = {}
        for each in range(0, len(molecular_formula_list),2):
            
            atoms_label =  molecular_formula_list[each]
            atoms_count = int(molecular_formula_list[each+1])
            if atoms_count > 0:
                self._d_molecular_formula[atoms_label] = int(atoms_count)
        
        self._d_molecular_formula[Labels.ion_type] = ion_type
        if adduct_atom:
            if adduct_atom in self._d_molecular_formula:
                self._d_molecular_formula[adduct_atom] += 1 
            else: self._d_molecular_formula[adduct_atom] = 1 

    def _from_str(self, molecular_formula_str,  ion_type, adduct_atom):
        # string has to be in the format 
        #'C10 H21 13C1 Cl1 37Cl1 etc'
        molecular_formula_list = molecular_formula_str.split(' ')
        final_formula = []
        for i in molecular_formula_list:
            atoms_count = self.split(Atoms.atoms_order, i)
            final_formula.extend(atoms_count)
        print(final_formula)
        self._from_list(final_formula, ion_type, adduct_atom)

    def split(self, delimiters, string, maxsplit=0): #pragma: no cover
        
        ''' does not work when formula has atoms with same caracaters:
            i.e - C10H21NNa
        '''
        regexPattern = '|'.join(map(re.escape, delimiters)) #pragma: no cover
        isotopes = re.findall(regexPattern, string) #pragma: no cover
        counts = re.split(regexPattern, string, maxsplit)  #pragma: no cover
       
        return [isotopes[0], int(counts[1])]

    @property
    def O_C(self): 
            
            if 'O' in self._d_molecular_formula.keys():
                return self._d_molecular_formula.get("O")/self._d_molecular_formula.get("C")
            else:
                return 0    
    
    @property
    def H_C(self): return self._d_molecular_formula.get("H")/self._d_molecular_formula.get("C")
    
    @property
    def dbe(self): return self._calc_dbe()
    
    @property
    def mz_nominal_calc(self): return int(self._calc_mz())

    @property    
    def mz_error(self): return self._calc_assignment_mass_error()

    @property
    def mz_calc(self): return self._calc_mz()

    @property
    def ion_type(self): 
        
        ion_type = self._d_molecular_formula.get(Labels.ion_type)
        if ion_type == Labels.protonated_de_ion:
            if self.ion_charge > 0: 
                return Labels.protonated
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
    def confidence_score(self): 
        
        if not self._confidence_score:
            
            self._confidence_score = self._calc_confidence_score()
        
        return self._confidence_score

    @property
    def isotopologue_similarity(self): 
        
        if not self._isotopologue_similarity:
           
           self._isotopologue_similarity = self._calc_isotopologue_confidence()  
       
        return self._isotopologue_similarity  
    
    @property
    def average_mz_error_score(self): 
        
        ''' includes the isotopologues'''
        
        if not self._mass_error_average_score:
           
           self._mass_error_average_score = self._calc_average_mz_score()  
        
        return self._mass_error_average_score

    @property
    def mz_error_score(self): 
        if not self._mz_error_score:
           
           self._mz_error_score = self._calc_mz_confidence()  
        
        return self._mz_error_score
    
    @property
    def kmd(self): return self._kdm

    @property
    def kendrick_mass(self): return self._kendrick_mass

    @property
    def knm(self): return self._nominal_km

    def change_kendrick_base(self, kendrick_dict_base):
        '''kendrick_dict_base = {"C": 1, "H": 2}'''
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(kendrick_dict_base)
                
    def isotopologues(self, min_abundance, current_mono_abundance, dynamic_range): 
        
        # this calculation ignores the hydrogen
        for mf in self._cal_isotopologues(self._d_molecular_formula, min_abundance, current_mono_abundance, dynamic_range ):
             
            yield MolecularFormulaIsotopologue(*mf, current_mono_abundance, self.ion_charge)
    
    def atoms_qnt(self,atom): 
        if atom in self._d_molecular_formula:
            return self._d_molecular_formula.get(atom)
        else:
            raise Warning('Could not find %s in this Molecular Formula object'%str(atom))
    
    def atoms_symbol(self, atom): 
        '''return the atom symbol without the mass number'''
        return ''.join([i for i in atom if not i.isdigit()])

    @property       
    def string(self):
        
        if self._d_molecular_formula:
            formula_srt = ''
            for atom in Atoms.atoms_order:
                if atom in self.to_dict().keys():
                    formula_srt += atom + str(int(self.to_dict().get(atom))) + ' '
            return formula_srt.strip()
        
        else:
            raise Exception("Molecular formula identification not performed yet")    
    
    @property
    def string_formated(self):
        
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        
        if self._d_molecular_formula:
            formula_srt = ''
            for atom in Atoms.atoms_order:
                if atom in self.to_dict().keys():
                    formula_srt += atom.translate(SUP) + str(int(self.to_dict().get(atom))).translate(SUB)
            return formula_srt
        
        else:
            raise Exception("Molecular formula identification not performed yet")    

    def to_dict(self):
        return self._d_molecular_formula

    def to_list(self):
        '''TODO ensure self._d_molecular_formula is a orderedDict'''
        
        if self._d_molecular_formula:
            formula_list = []    
            
            for atom, atom_number in self._d_molecular_formula.items():
    
                if atom != Labels.ion_type:
                    
                    formula_list.append(atom)
                    formula_list.append(atom_number)
    
            return formula_list
        else:
            raise Exception("Molecular formula identification not performed yet")
    
    @property
    def class_label(self):
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list()
            classstring = '' 
            
            for each in range(0, len(formulalist),2):
                
                if formulalist[each] != 'C' and formulalist[each] != 'H' and formulalist[each] != 'HC':
                     
                    classstring = classstring + str(formulalist[each]) + str(formulalist[each+1]) + ' '    
            
            if classstring == '': classstring = 'HC'
                
            classstring = classstring.strip()
            
            if self._d_molecular_formula.get(Labels.ion_type) == Labels.radical_ion:    
                
                return classstring + ' -R'
            
            #elif self._d_molecular_formula.get(Labels.ion_type) == Labels.adduct_ion:    
                
            #    return classstring + ' -A'

            else: return classstring
            
            #'dict, tuple or string'
        
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
    def __init__(self, _d_molecular_formula, prob_ratio, mono_abundance, ion_charge, mspeak_parent=None):
        
        super().__init__(_d_molecular_formula,  ion_charge)
        #prob_ratio is relative to the monoisotopic peak p_isotopologue/p_mono_isotopic
        
        self.prob_ratio = prob_ratio
        
        self.abundance_calc = mono_abundance * prob_ratio

        self.is_isotopologue = True
        
        self.mspeak_index_mono_isotopic = None

        self.mono_isotopic_formula_index = None
        # parent mass spectrum peak obj instance
        self._mspeak_parent = mspeak_parent

    
    @property
    def area_error(self):
        return self._calc_area_error()

    @property
    def abundance_error(self):
        return self._calc_abundance_error()

        