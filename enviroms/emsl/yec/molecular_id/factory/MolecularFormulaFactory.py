from copy import deepcopy
from enviroms.emsl.yec.mass_spectrum.calc.MolecularFormulaCalc import MolecularFormulaCalc
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
from IsoSpecPy import IsoSpecPy
from numpy import exp

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
    def ion_charge(self): return self._ion_charge
         
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
    def isotopologues(self): 
        
        return [MolecularFormulaIsotopologue(*mf, self.ion_charge) for mf in self._cal_isotopologues(self._d_molecular_formula)]
    
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
        INPUT {'C':10, 'H', 20, 'O', 2, etc} Atomic labels need to follow Atoms class atoms labels
        
        This function needs to be expanded to include the calculation of resolving power
        and plot the results.
        
        *   use this function at runtime during the molecular identification algorithm
            only when a positive ID is observed to the monoisotopic ion
        
        *   use this function to simulate mass spectrum 
            (needs resolving power calculation to be fully operational)        
            last update on 07-19-2019, Yuri E. Corilo 
        
        *   expect it to break when adding non-conventional atoms (not yet tested)
        
        *   it needs speed optimization update: (Using IsoSpeccPy Library (faster)) 
            https://github.com/MatteoLacki/IsoSpec
        """
        #this value needs to be calculated from the mass spec dinamic range
        min_relative_abundance = 0.01
        
        cut_off_to_IsoSpeccPy = 1 - min_relative_abundance
        isotopologue_and_pro_ratio_tuples = []
        atoms_labels = (atom for atom in formula_dict.keys() if atom != 'IonType' and atom != "H")
        atoms_count = []
        masses_list_tuples = []
        props_list_tuples = []
        all_atoms_list = []
        for atom_label in atoms_labels:
            
            if not len(Atoms.isotopes.get(atom_label))>1:
                'This atom_label has no heavy isotope'
                atoms_count.append(formula_dict.get(atom_label))
                mass = Atoms.atomic_masses.get(atom_label)
                prop = Atoms.isotopic_abundance.get(atom_label)
                masses_list_tuples.append([mass])
                props_list_tuples.append([prop])
                all_atoms_list.append(atom_label)
                
            else:
                
                isotopoes_label_list = Atoms.isotopes.get(atom_label)[1]
            
                if len(isotopoes_label_list) > 1:
                    isotopoes_labels = [i for i in isotopoes_label_list]
                else:
                    isotopoes_labels = [isotopoes_label_list[0]]
                
                #all_atoms_list.extend(isotopoes_labels) 
                isotopoes_labels = [atom_label] + isotopoes_labels
                
                all_atoms_list.extend(isotopoes_labels)
                #print(all_labels)
                masses = [Atoms.atomic_masses.get(atom_label) for atom_label in isotopoes_labels]
                props = [Atoms.isotopic_abundance.get(atom_label) for atom_label in isotopoes_labels]
                
                atoms_count.append(formula_dict.get(atom_label))
                masses_list_tuples.append(masses)
                props_list_tuples.append(props)
        
        iso = IsoSpecPy.IsoSpec(atoms_count,masses_list_tuples,props_list_tuples, cut_off_to_IsoSpeccPy )
        confs = iso.getConfs()
        
        masses = confs[0]
        probs = exp(confs[1])
        molecular_formulas = confs[2]
        
        #print(all_atoms_list)    
        #print("mass",masses)
        #print("probs",exp(probs))
        #print("molecular_formula",molecular_formulas)
        
        for isotopologue_index in range(1,len(iso),1):
            #skip_mono_isotopic 
            
            formula_list = molecular_formulas[isotopologue_index]
            new_formula_dict = dict(zip(all_atoms_list, formula_list))
            new_formula_dict["IonType"] = formula_dict.get("IonType")
            
            prop_mono_iso = probs[0]
            prop_ratio = probs[isotopologue_index]/prop_mono_iso
            
            isotopologue_and_pro_ratio_tuples.append((new_formula_dict,prop_ratio))
        
        return isotopologue_and_pro_ratio_tuples

class MolecularFormulaIsotopologue(MolecularFormula):
        
    '''
    classdocs
    '''
    def __init__(self, _d_molecular_formula, prop_ratio, ion_charge, exp_mz=None):
        
        super().__init__(_d_molecular_formula,  ion_charge)
        #prop_ratio is relative to the monoisotopic peak p_isotopologue/p_mono_isotopic
        self.prop_ratio = prop_ratio
        
        if exp_mz:
            
            self._assigment_mass_error = self._calc_assigment_mass_error(exp_mz)
            self._confidence_score = False