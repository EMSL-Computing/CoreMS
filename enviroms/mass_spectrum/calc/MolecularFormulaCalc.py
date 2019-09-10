__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

from enviroms.encapsulation.Constants import Atoms
from enviroms.encapsulation.Constants import Labels
from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSettings
from IsoSpecPy import IsoSpecPy
from numpy import exp

class MolecularFormulaCalc:
    
    def _calc_resolving_power_low_pressure(self, B, T):
        #T << collisional damping time 
        #T = transient time (seconds)
        #B = Magnetic Strenght (Testa)
        return (1.274 * 10000000 * B * T) *(1/self.mz_theor)    

    def _calc_resolving_power_high_pressure(self, B, t):
        #T << collisional damping time 
        # t= collisional dumping contanst
        #T = transient time (seconds)
        #B = Magnetic Strenght (Testa)
        return (2.758 * 10000000 * B * T) *(1/self.mz_theor)    

    def _calc_confidence_score(self):
        raise NotImplementedError
        
    
    def _calc_mz_theor(self):
        
        if self.ion_charge:
            
            mass = 0
            
            for each_atom in self._d_molecular_formula.keys() :
                
                if each_atom != Labels.ion_type and each_atom != 'HC':
                    
                    try:
                        mass = mass + Atoms.atomic_masses[each_atom]  *  self._d_molecular_formula.get(each_atom)
                    except: print(Labels.ion_type, each_atom) 
            return mass + ((-self.ion_charge) * Atoms.electron_mass)
        
        else:
            
            raise Exception("Please set ion charge first")
         
    def _calc_assigment_mass_error(self, mz_exp, method='ppm'):
        '''methodo should be ppm or ppb'''
        
        Hum_Milhao = 1000000
        Hum_Bilhao = 1000000000        
        
        if method == 'ppm':
            mult_factor = Hum_Milhao
        
        elif method == 'ppb':
            mult_factor = Hum_Bilhao
        
        else:
            raise Exception("method needs to be ppm or ppb, you have entered %s" % method)
              
        if mz_exp:
            
            self._assigment_mass_error = ((self.mz_theor - mz_exp)/self.mz_theor)*mult_factor
            #self.parent need to have a MassSpecPeak associated with the MolecularFormula class
            return self._assigment_mass_error
        
        else:
            
            raise Exception("Please set mz_theor first")    
    
    def _calc_abundance_error(self, mono_abundance, iso_abundance, method='perc'):
        '''methodo should be ppm or ppb'''
        
        mult_factor = 100
        if self.prop_ratio:
            
            theor_abundance = mono_abundance* self.prop_ratio
            #self.parent need to have a MassSpecPeak associated with the MolecularFormula class
            return ((theor_abundance - iso_abundance )/theor_abundance)*mult_factor
        
        else:
            
            raise Exception("Please calc_isotopologues")    

    @property
    def dbe_ai(self):
            
        carbons =  self._d_molecular_formula.get('C')
        hydrogens = self._d_molecular_formula.get('H')
        oxygens = self._d_molecular_formula.get('O')
        return 1 + (((2*carbons) - hydrogens - (2*oxygens))*0.5)

    def _calc_dbe(self):
            
            individual_dbe = 0
            
            for atom in self._d_molecular_formula.keys():
                
                if atom != Labels.ion_type:
                    
                    n_atom = int(self._d_molecular_formula.get(atom))
                    
                    clean_atom = ''.join([i for i in atom if not i.isdigit()]) 
                    
                    valencia = MoleculaLookupTableSettings.used_atom_valences.get(clean_atom)
                    #valencia = Atoms.atoms_valence.get(atom)
                    
                    if valencia and valencia > 0:
                        #print atom, valencia, n_atom, individual_dbe
                        individual_dbe = individual_dbe + (n_atom * (valencia - 2))
                    else:
                        continue
            
            return 1 + (0.5 * individual_dbe)

    def _calc_kdm(self, dict_base):
        '''dict_base = {"C": 1, "H": 2}
        '''
        mass = 0
        for atom in dict_base.keys():
            mass = mass + Atoms.atomic_masses.get(atom) * dict_base.get(atom)
        
        kendrick_mass = (int(mass)/mass)*self.mz_theor
        
        nominal_km =int(kendrick_mass)
       
        kmd = (nominal_km - kendrick_mass) * 100
        
        #kmd = (nominal_km - km) * 1
        kdm  = round(kmd,0)
        
        return kdm, kendrick_mass, nominal_km

    
    @staticmethod
    def _cal_isotopologues(formula_dict, min_abudance, current_abundance):
        
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
        
        *   it might break when adding non-conventional atoms (not yet tested)
            
        *   it needs speed optimization; update: (Using IsoSpeccPy, a C Library (fast and acurate)) 
            https://github.com/MatteoLacki/IsoSpec
        """
        #this value needs to be calculated from the mass spec dinamic range
        #updated it to reflect min possible mass peak abundance
        #min_abudance/current_abundance will only work if the isotopologue abundance is less than the monoisotopic abundance
        #needs to be checked
        min_relative_abundance = min_abudance/current_abundance
        
        cut_off_to_IsoSpeccPy = 1 - min_relative_abundance
        isotopologue_and_pro_ratio_tuples = []
        atoms_labels = (atom for atom in formula_dict.keys() if atom != Labels.ion_type) #and atom != "H")
        
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
            new_formula_dict[Labels.ion_type] = formula_dict.get(Labels.ion_type)
            #new_formula_dict['H'] = formula_dict.get('H')

            prop_mono_iso = probs[0]
            prop_ratio = probs[isotopologue_index]/prop_mono_iso
            
            isotopologue_and_pro_ratio_tuples.append((new_formula_dict,prop_ratio))
        
        return isotopologue_and_pro_ratio_tuples        
    
    
    