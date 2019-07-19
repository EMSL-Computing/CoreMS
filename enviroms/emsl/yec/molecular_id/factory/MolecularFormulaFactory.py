from scipy.stats import binom
from copy import deepcopy
from enviroms.emsl.yec.mass_spectrum.calc.MolecularFormulaCalc import MolecularFormulaCalc
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

class MolecularFormula(MolecularFormulaCalc):
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, ion_charge, exp_mz=None):
        
        self._d_molecular_formula = _d_molecular_formula
        self._ion_chage = ion_charge
        
        if exp_mz:
            
            self._assigment_mass_error = self._calc_assigment_mass_error(exp_mz)
            self._confidence_score = False
        
    @property
    def O_C(self): return self._d_molecular_formula.get("O")/self._d_molecular_formula.get("C")
    
    @property
    def H_C(self): return self._d_molecular_formula.get("H")/self._d_molecular_formula.get("C")
    
    @property
    def dbe(self): return self._calc_dbe()
    
    @property
    def ion_charge(self): return self._ion_chage
         
    @property
    def mz_theor(self): return self._calc_theoretical_mz()
    
    @property
    def mz_nominal_theo(self): return int(self._calc_theoretical_mz())

    @property    
    def assigment_mass_error(self): return self._assigment_mass_error
        
    @property
    def ion_type(self): return self._d_molecular_formula.get("IonType")
    
    @property
    def atoms(self): return [key for key in self._d_molecular_formula.keys() if key != 'IonType']
    
    @property
    def confidence_score(self): return self._confidence_score
             
    @confidence_score.setter
    def confidence_score(self): return self._calc_confidence_score() 
        
    @property
    def isotopologues(self): return [MolecularFormulaIsotopologue(mf) for mf in self._cal_isotopologues(self._d_molecular_formula)]
    
    def atoms_qnt(self,atom): 
        if atom in self._d_molecular_formula:
            return self._d_molecular_formula.get(atom)
        else:
            raise Warning('Could not find %s in this Molecular Formula object'%str(atom))
    
    def atoms_symbol(self, atom): 
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
    @property   
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
    def class_label(self):
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list
            classstring = "" 
            
            for each in range(0, len(formulalist),2):
                
                if formulalist[each] != 'C' and formulalist[each] != 'H':
                     
                    classstring = classstring + str(formulalist[each]) + str(formulalist[each+1]) + " "    
            
            if classstring == "": classstring = "HC "
                
            if self._d_molecular_formula.get("IonType") == 'RADICAL':    
                
                return classstring[0:-1] + " " + "-R"
            
            else: return classstring[0:-1] + " "
            
            'dict, tuple or string'
        else:
            raise Exception("Molecular formula identification not performed yet")        
    
    def _cal_isotopologues(self, formula_dict):
        
        """
        primary function to look for isotopologues based on a monoisotopic molecular formula
        INPUT {'C':10, 'H', 20, 'O', 2, etc} Atomic labels need to follow Atoms class 
        
        This function needs to be expanded to include the calculation of resolving power
        and plot the results.
        The above description is a project by itself, but one could work on a web app pet project.
        
        *   use this function at runtime during the molecular identification algorithm
            only when a positive ID is observed to the monoisotopic ion
        
        *   use this function to simulate mass spectrum 
            (needs resolving power calculation to be fully operational)        
            last update on 07-19-2019, Yuri E. Corilo 
        
        *   This version is still beta, the ration heavy_iso/mono_iso still needs work when dealing with two or more heavy isotopes
            please double-check your results;
            expect it to break when adding non-conventional atoms (not yet tested)
            it needs speed optimization 
        """
        
        iteration_limit = 10
        iso_threo_ratio = 1000
        list_isotopologues = []
        #formula_dict = self._d_molecular_formula
        
        # this list set a limit for the isotopolgue iterations
        list_of_possible_carbon_13 = range(1,iteration_limit,1)
        #print(formula_dict.keys())
        for atom in formula_dict.keys():
            
            #remove hidrogen and iontype label
            if  atom != 'IonType' and atom != "H":
                number_atoms = formula_dict.get(atom)
                #remove atomic mass label
                atom_symbol = ''.join([i for i in atom if not i.isdigit()])
                #get all isotopes labels;
                # most abundance is the mono isotope with the atomic symbol(i.e., no atomic mass label);
                #dict key returning 'C': [13C, 14C, etc]
                isotopes_name_labels = Atoms.isotopes.get(atom_symbol)
                #check for the existence of isotopes at the Atoms class #this way is not ideal;it needs to be changed
                
                if isotopes_name_labels:
                    if len(isotopes_name_labels) > 1:
                        #isotopes_name_labels[0] = FullName
                        #isotopes_name_labels[0] = Tuple with Full Label i.e, 34S, 13C etc
                        #print(isotopes_name_labels)
                        
                        for isotopo_label in isotopes_name_labels[1]:
                            #print("isotopo_label", isotopo_label)
                            index = 0
                            
                            while iso_threo_ratio > 0.01:
                                
                                iso_proba = Atoms.isotopic_abundance.get(isotopo_label)
                                qnt_iso_atom = list_of_possible_carbon_13[index] 
                                index = index + 1
                                
                                p12 = binom.pmf(0, number_atoms, iso_proba)
                                p13 = binom.pmf(qnt_iso_atom, number_atoms, iso_proba)
                                
                                iso_threo_ratio = p13 / p12
                                
                                #print("iso_threo_ratio", iso_threo_ratio)
                                if iso_threo_ratio > 0.01:
                                    
                                    #creates the new molecular formula dict
                                    new_formula_dict = deepcopy(formula_dict)
                                    #exlude atom so it wont run again on the recursive call
                                    del new_formula_dict[atom] 
                                    
                                    #recursive call
                                    recursive_isotopologues_list = self._cal_isotopologues(new_formula_dict)
                                    
                                    #add both mono and isotopologue atoms 
                                    new_formula_dict[isotopo_label] = qnt_iso_atom
                                    new_formula_dict[atom] = formula_dict.get(atom) - qnt_iso_atom
                                    if new_formula_dict[atom] == 0:
                                        del new_formula_dict[atom]
                                    
                                    #add both mono and isotopologue atoms for recursive call
                                    for item_dict in recursive_isotopologues_list:
                                        item_dict[isotopo_label] = qnt_iso_atom
                                        item_dict[atom] = formula_dict.get(atom) - qnt_iso_atom
                                        if item_dict[atom] == 0:
                                            del item_dict[atom]

                                    list_isotopologues.append(new_formula_dict)
                                    #sum up results
                                    recursive_isotopologues_list
                                    #list_isotopologues = list_isotopologues + recursive_isotopologues_list
                                    list_isotopologues.extend(x for x in recursive_isotopologues_list if x not in list_isotopologues)
                
                    iso_threo_ratio = 1000
        
        return  list_isotopologues
                
class MolecularFormulaIsotopologue(MolecularFormula):
        
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, iso_threo_ratio=1, exp_mz=None):
        
        self._d_molecular_formula = _d_molecular_formula
        self.isotopologue_ratio = iso_threo_ratio
        
        if exp_mz:
            
            self._assigment_mass_error = self._calc_assigment_mass_error(exp_mz)
            self._confidence_score = False