__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

from IsoSpecPy import IsoSpecPy
from numpy import exp, dot, power
from numpy.linalg import norm
from pandas import DataFrame
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr, spearmanr, kendalltau
from corems.encapsulation.constant import Atoms
from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.processingSetting import MolecularSearchSettings


class MolecularFormulaCalc:
    
    def _calc_resolving_power_low_pressure(self, B, T):
        '''
        ## Parameters
        ----------
        #### T << collisional damping time 
        #### T = transient time (seconds)
        #### B = Magnetic Strength (Testa)
        '''
        return (1.274 * 10000000 * B * T) *(1/self.mz_calc)    

    def _calc_resolving_power_high_pressure(self, B, t):
        '''
        ## Parameters
        ----------
        #### T << collisional damping time 
        #### t= collisional dumping constant
        #### T = transient time (seconds)
        #### B = Magnetic Strength (Testa)
        '''
        return (2.758 * 10000000 * B * T) *(1/self.mz_calc)    

    def _calc_mz(self):
        
        if self.ion_charge:
            
            mass = 0
            
            for each_atom in self._d_molecular_formula.keys() :
                
                if each_atom != Labels.ion_type and each_atom != 'HC':
                    
                    try:
                        mass = mass + Atoms.atomic_masses[each_atom]  *  self._d_molecular_formula.get(each_atom)
                    except: print(Labels.ion_type, each_atom) 
            
            return mass + ((-1*self.ion_charge) * Atoms.electron_mass)
        
        else:
            
            raise Exception("Please set ion charge first")
         
    def _calc_assignment_mass_error(self, method='ppm'):
        
        '''
        ## Parameters
        ----------
        #### mz_exp: float, 
        ####       Experimental m/z 
        #### method: string, 
        ####       ppm or ppb
        '''
         
        if method == 'ppm':
            multi_factor = 1000000
        
        elif method == 'ppb':
            multi_factor = 1000000
        
        else:
            raise Exception("method needs to be ppm or ppb, you have entered %s" % method)
              
        if self._mspeak_parent.mz_exp:
            
            self._assignment_mass_error = ((self.mz_calc - self._mspeak_parent.mz_exp)/self.mz_calc)*multi_factor

            return ((self.mz_calc - self._mspeak_parent.mz_exp)/self.mz_calc)*multi_factor
        
        else:
            
            raise Exception("Please set mz_calc first")    
    
    def _calc_fine_isotopic_similarity(self):
        pass
    
    def _calc_mz_confidence(self):

        return  exp(-((self._mspeak_parent.mz_exp - self.mz_calc)**2 ) / (2 * self._mspeak_parent.predicted_std**2))
    def _calc_confidence_score(self):
        '''
        ### Assumes random mass error, i.e, spectrum has to be calibrated and with zero mean
        #### TODO: Add spectral similarity 

        ## Parameters
        ----------
        #### mz_exp:
        ####    Experimental m/z 
        #### predicted_std:
        ####    Standart deviation calculated from Resolving power optimization or constant set by User 
        
        '''
        abundance_weight, b = 1.2, 0.5

        if not self._mspeak_parent.predicted_std: return -1
           
        else:

            if self.is_isotopologue:
                
                return self._calc_mz_confidence()
        
            else:
                
                dict_mz_abund_ref = {}
                
                for mf in self.expected_isotopologues:
                    dict_mz_abund_ref[mf.mz_calc] = mf.abundance_calc#power(mf.abundance_calc, abundance_weight) *  power(mf.mz_calc, b)  
                
                dict_mz_abund_exp = {}
                
                accumulated_mz_score = []
                for mf in self.expected_isotopologues:
                    
                    if mf._mspeak_parent:
                        
                        dict_mz_abund_exp[mf.mz_calc] = mf._mspeak_parent.abundance#power(mf._mspeak_parent.abundance, abundance_weight) * power(mf.mz_calc, b)   
                        accumulated_mz_score.append(mf._calc_mz_confidence())
                    
                    else:
                        # fill missing mz with abundance 0
                        dict_mz_abund_exp[mf.mz_calc] = 0.0
                        accumulated_mz_score.append(0.0)

                df = DataFrame([dict_mz_abund_exp, dict_mz_abund_ref])
                #print(df.head())
                #calculate cosine correlation, 
                x = df.T[0].values
                y = df.T[1].values

                correlation = kendalltau(x, y)[0]
                #correlation = (1 - cosine(x, y))
                
                #correlation = dot(x, y)/(norm(x)*norm(y))
                accumulated_mz_score.append(self._calc_mz_confidence())
                
                average_mz_score = sum(accumulated_mz_score)/len(accumulated_mz_score)
                
                score = ((correlation) * (average_mz_score**2))**(1/3)
                
                print("correlation",correlation)
                print("average_mz_score",average_mz_score)
                print("mz_score",self._calc_mz_confidence())
                print("score",score)
                return score  


    def _calc_abundance_error(self, method='percentile'):
        '''method should be ppm, ppb or percentile'''
        
        mult_factor = 100

        iso_abundance = self._mspeak_parent.abundance
        mono_abundance =self._mspeak_parent._ms_parent[self.mspeak_index_mono_isotopic].abundance

        if self.prob_ratio:
            
            theor_abundance = mono_abundance* self.prob_ratio
            #self.parent need to have a MassSpecPeak associated with the MolecularFormula class
            return ((theor_abundance - iso_abundance )/theor_abundance)*mult_factor
        
        else:
            
            raise Exception("Please calc_isotopologues")    

    def _calc_area_error(self, method='percentile'):
        '''method should be ppm, ppb or percentile'''
        
        mult_factor = 100
        
        iso_area = self._mspeak_parent.area
        mono_area =self._mspeak_parent._ms_parent[self.mspeak_index_mono_isotopic].area

        if self.prob_ratio:
            
            if mono_area and iso_area: 

                #exp_ratio = iso_area/mono_area  
                print(mono_area, iso_area)            
                area_calc = mono_area* self.prob_ratio

                #self.parent need to have a MassSpecPeak associated with the MolecularFormula class
                return ((area_calc - iso_area )/area_calc)*mult_factor
                #return ((self.prob_ratio - exp_ratio )/self.prob_ratio)*mult_factor
            
            else:
                
                #centroid mass spectrum
                return 0
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
                    
                    valencia = MolecularSearchSettings.used_atom_valences.get(clean_atom)
                    #valencia = Atoms.atoms_covalence.get(atom)
                    
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
        
        kendrick_mass = (int(mass)/mass)*self.mz_calc
        
        nominal_km =int(kendrick_mass)
       
        kmd = (nominal_km - kendrick_mass) * 100
        
        #kmd = (nominal_km - km) * 1
        kdm  = round(kmd,0)
        
        return kdm, kendrick_mass, nominal_km

    
    @staticmethod
    def _cal_isotopologues(formula_dict, min_abundance, current_abundance):
        
        '''
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
            
        *   it needs speed optimization; update: (Using IsoSpeccPy, a C Library (fast and accurate)) 
            https://github.com/MatteoLacki/IsoSpec
        '''
        # updated it to reflect min possible mass peak abundance
        
        min_relative_abundance = min_abundance/current_abundance
        cut_off_to_IsoSpeccPy = 1 - min_relative_abundance
        
        #print(min_relative_abundance, min_abundance, current_abundance, cut_off_to_IsoSpeccPy)
        
        atoms_labels = (atom for atom in formula_dict.keys() if atom != Labels.ion_type and atom != 'H')
       
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
                
                isotopes_label_list = Atoms.isotopes.get(atom_label)[1]
            
                if len(isotopes_label_list) > 1:
                    'This atom_label has two or more heavy isotope'
                    isotopoes_labels = [i for i in isotopes_label_list]
                else:
                    'This atom_label only has one heavy isotope'
                    isotopoes_labels = [isotopes_label_list[0]]
                
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
        
        conf = iso.getConfs()
        
        masses = conf[0]
        probs = exp(conf[1])
        molecular_formulas = conf[2]
        
        new_formulas = []
        
        for isotopologue_index in range(0,len(iso),1):
            #skip_mono_isotopic 
            
            formula_list = molecular_formulas[isotopologue_index]
            new_formula_dict = dict(zip(all_atoms_list, formula_list))
            new_formula_dict[Labels.ion_type] = formula_dict.get(Labels.ion_type)
            if formula_dict.get('H'):
                new_formula_dict['H'] = formula_dict.get('H')

            new_formulas.append({x:y for x,y in new_formula_dict.items() if y!=0})

        if (new_formulas):
            # find where monoisotopic is
            index_mono = new_formulas.index(formula_dict)   
            # calculate ratio iso/mono
            probs = list(probs/probs[index_mono])
            
            # delete the monoisotopic
            del probs[index_mono]
            del new_formulas[index_mono]
            
        return zip(new_formulas, probs )
    
    
    