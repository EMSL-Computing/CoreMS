__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
#from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MolecularSpaceTableSetting
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSetting
from IsoSpecPy import IsoSpecPy
from numpy import exp

class MolecularFormulaCalc:
    
    def _calc_confidence_score(self):
        raise NotImplementedError
        return 0
    
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
              
        if exp_mz:
            #self.parent need to have a MassSpecPeak associated with the MolecularFormula class
            return ((self.theoretical_mz - self.exp_mz)/self.theoretical_mz)*mult_factor
        
        else:
            
            raise Exception("Please set theoretical_mz first")    
    
    def _calc_dbe(self):
            
            individual_dbe = 0
            
            for atom in self._d_molecular_formula.keys():
                
                if atom != "IonType":
                    
                    atom = ''.join([i for i in atom if not i.isdigit()]) 
                    
                    n_atom = int(self._d_molecular_formula.get(atom))
                    
                    valencia = MoleculaLookupTableSetting.used_atom_valences.get(atom)
                    #valencia = Atoms.atoms_valence.get(atom)
                    
                    if valencia and valencia > 0:
                        #print atom, valencia, n_atom, individual_dbe
                        individual_dbe = individual_dbe + (n_atom * (valencia - 2))
                    else:
                        continue
            
            return 1 + (0.5 * individual_dbe)
    
    @staticmethod
    def _cal_isotopologues(formula_dict):
        
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
            
        *   it needs speed optimization update: (Using IsoSpeccPy, a C Library (fast and acurate)) 
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
            #print( atom_label)
            
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
                    'This atom_label has two or more heavy isotope'
                    isotopoes_labels = [i for i in isotopoes_label_list]
                else:
                    'This atom_label only has one heavy isotope'
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
        
        for isotopologue_index in range(1,len(iso),1):
            #skip_mono_isotopic 
            
            formula_list = molecular_formulas[isotopologue_index]
            new_formula_dict = dict(zip(all_atoms_list, formula_list))
            new_formula_dict["IonType"] = formula_dict.get("IonType")
            
            prop_mono_iso = probs[0]
            prop_ratio = probs[isotopologue_index]/prop_mono_iso
            
            isotopologue_and_pro_ratio_tuples.append((new_formula_dict,prop_ratio))
        
        return isotopologue_and_pro_ratio_tuples        
    
    
    