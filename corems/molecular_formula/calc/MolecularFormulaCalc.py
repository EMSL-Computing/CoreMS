__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

from IsoSpecPy import IsoSpecPy
from numpy import isnan, power, exp, nextafter
from pandas import DataFrame
from scipy.stats import pearsonr, spearmanr, kendalltau

from corems.encapsulation.constant import Atoms
from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import MSParameters
from corems.molecular_id.calc.SpectralSimilarity import SpectralSimilarity

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

    def _calc_resolving_power_high_pressure(self, B, T):
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
        
        def protonated_mass(mass, ion_charge):
            return (mass + (ion_charge * Atoms.atomic_masses.get("H")) + (ion_charge * -1 * Atoms.electron_mass))/abs(ion_charge)
        
        def radical_mass(mass, ion_charge):
            return (mass + (ion_charge * -1 * Atoms.electron_mass))/ abs(ion_charge)

        if self.ion_charge:
            
            mass = 0
            
            ion_type = self._d_molecular_formula.get(Labels.ion_type)

            for each_atom in self._d_molecular_formula.keys() :
                
                if each_atom != Labels.ion_type and each_atom != 'HC':
                    
                    try:
                        mass = mass + Atoms.atomic_masses[each_atom]  *  self._d_molecular_formula.get(each_atom)
                    except: print(Labels.ion_type, each_atom) 
            
            if ion_type == Labels.protonated_de_ion:
                return protonated_mass(mass, self.ion_charge)
            elif ion_type == Labels.radical_ion or ion_type == Labels.adduct_ion:   
                return radical_mass(mass, self.ion_charge)
            else:
                #formula is probably ion form used for bruker ref list
                return radical_mass(mass, self.ion_charge)
                
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
            
            self._assignment_mass_error = ((self._mspeak_parent.mz_exp - self.mz_calc )/self.mz_calc)*multi_factor

            return ((self._mspeak_parent.mz_exp - self.mz_calc)/self.mz_calc)*multi_factor
        
        else:
            
            raise Exception("No ms peak associated with the molecular formula instance %s", self)    
    
    
    def _calc_mz_confidence(self, mean=0):
        
        # predicted std not set, using 0.3
        if not self._mspeak_parent.predicted_std: self._mspeak_parent.predicted_std = 0.3
        
        #print( self._mspeak_parent.predicted_std)
        
        return  exp( -1 * (power((self.mz_error -  mean),2)  / (2 * power(self._mspeak_parent.predicted_std,2)) ))
    
    def _calc_isotopologue_confidence(self):

        if self.is_isotopologue:
            # confidence of isotopologue is pure mz error 
            # TODO add more features here 
            
            mformula_index = self.mono_isotopic_formula_index
            mspeak_index = self.mspeak_index_mono_isotopic

            mspeak = self._mspeak_parent._ms_parent[mspeak_index]
            
            expected_isotopologues = mspeak[mformula_index].expected_isotopologues
            
            mono_mz = mspeak[mformula_index].mz_calc
            mono_abundance = mspeak.abundance

        else:

            mono_mz = self.mz_calc
            mono_abundance = self._mspeak_parent.abundance

            expected_isotopologues = self.expected_isotopologues
            # has isotopologues based on current dinamic range
            
        if expected_isotopologues:
            
            dict_mz_abund_ref = {'mz':[mono_mz], 'abundance':[mono_abundance]}
            
            # get reference data
            for mf in expected_isotopologues:
                dict_mz_abund_ref['abundance'].append(mf.abundance_calc)
                dict_mz_abund_ref['mz'].append(mf.mz_calc)

            dict_mz_abund_exp = {mono_mz:mono_abundance}
            
            # get experimental data
            for mf in expected_isotopologues:
                
                # molecular formula has been assigned to a peak
                if mf._mspeak_parent:
                    #stores mspeak abundance
                    dict_mz_abund_exp[mf.mz_calc] = mf._mspeak_parent.abundance
                    
                else:
                    # fill missing mz with abundance 0 and mz error score of 0
                    dict_mz_abund_exp[mf.mz_calc] = nextafter(0, 1)
            
            distance = SpectralSimilarity(dict_mz_abund_exp, dict_mz_abund_ref).manhattan_distance()
            correlation = 1 - self.normalize_distance(distance, [0, 2])
            #correlation = dwt_correlation(dict_mz_abund_exp, dict_mz_abund_ref)
            #correlation = cosine_correlation(dict_mz_abund_exp, dict_mz_abund_ref)
            
            if correlation == 1:
                print(dict_mz_abund_exp,dict_mz_abund_ref)
            if isnan(correlation):
                #print(dict_mz_abund_exp,dict_mz_abund_ref)
                correlation = 0.00001
        
        else:
            
            # no isotopologue expected giving a correlation score of 0.0 but it needs optimization
            correlation = 0.0

        return correlation

    def normalize_distance(self, dist, dist_range):
    
        result = (dist - dist_range[0]) / (dist_range[1] - dist_range[0])

        if result < 0:
            result = 0.
        elif result > 1:
            result = 1.

        return result

    def _calc_average_mz_score(self):
        
        if self.is_isotopologue:
            # confidence of isotopologue is pure mz error 
            # TODO add more features here 
            
            mformula_index = self.mono_isotopic_formula_index
            mspeak_index = self.mspeak_index_mono_isotopic

            mspeak = self._mspeak_parent._ms_parent[mspeak_index]
            
            expected_isotopologues = mspeak[mformula_index].expected_isotopologues

        else:
            
            expected_isotopologues = self.expected_isotopologues
            # has isotopologues based on current dinamic range
        
        accumulated_mz_score = [self.mz_error_score]
        
        if expected_isotopologues:
            
            for mf in expected_isotopologues:
                # molecular formula has been assigned to a peak
                if mf._mspeak_parent:
                    #stores mspeak abundance
                    accumulated_mz_score.append(mf.mz_error_score)
                else:
                    # fill missing mz with abundance 0 and mz error score of 0
                    accumulated_mz_score.append(0.0)

        average_mz_score = sum(accumulated_mz_score)/len(accumulated_mz_score)
        
        if isnan(average_mz_score):
                average_mz_score = 0.0

        return average_mz_score       

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
        
        isotopologue_correlation = self.isotopologue_similarity
        average_mz_score = self.average_mz_error_score
        # add monoisotopic peak mz error score
        
        # calculate score with higher weight for mass error
        #score = power(((isotopologue_correlation) * (power(average_mz_score,3))),1/4)
        a = self._mspeak_parent._ms_parent.molecular_search_settings.mz_error_score_weight
        b = self._mspeak_parent._ms_parent.molecular_search_settings.isotopologue_score_weight
        
        score = (isotopologue_correlation*b) + (average_mz_score*a)
        
        #if round(average_mz_score,2) == 0.00:
        #    print(a,b, average_mz_score, isotopologue_correlation, score, isotopologue_correlation*b)
        

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
                    
                    if self._mspeak_parent:
                        valencia = self._mspeak_parent._ms_parent.molecular_search_settings.used_atom_valences.get(clean_atom)
                    else:
                        valencia = MSParameters.molecular_search.used_atom_valences.get(clean_atom)
                    #valencia = Atoms.atoms_covalence.get(atom)
                    
                    if type(valencia) is tuple:
                        valencia = valencia[0]
                    if valencia > 0:
                        #print atom, valencia, n_atom, individual_dbe
                        individual_dbe = individual_dbe + (n_atom * (valencia - 2))
                    else:
                        continue
            
            dbe = 1 + (0.5 * individual_dbe)
            
            
            
            if self.ion_type == Labels.adduct_ion:
                dbe = dbe + 0.5
            
            return dbe

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

    def _cal_isotopologues(self, formula_dict, min_abundance, current_abundance, ms_dynamic_range):
        
        '''
        primary function to look for isotopologues based on a monoisotopic molecular formula
        INPUT {'C':10, 'H', 20, 'O', 2, etc} Atomic labels need to follow Atoms class atoms labels
        
        This function needs to be expanded to include the calculation of resolving power
        and plot the results.
        
        *   use this function at runtime during the molecular identification algorithm
            only when a positive ID is observed to the monoisotopic ion
        
        *   use this function to simulate mass spectrum 
            (needs resolving power calculation to be fully operational)        
            last update on 05-26-2020, Yuri E. Corilo 
        
        *   it might break when adding non-conventional atoms (not yet tested)
            
        *   it needs speed optimization; update: (Using IsoSpeccPy, a C Library (fast and accurate)) 
            https://github.com/MatteoLacki/IsoSpec
        '''
        # updated it to reflect min possible mass peak abundance
        cut_off_to_IsoSpeccPy = 1-(1/ms_dynamic_range)
        
        #print("cut_off_to_IsoSpeccPy", cut_off_to_IsoSpeccPy, current_abundance, min_abundance, ms_dynamic_range)
        #print(cut_off_to_IsoSpeccPy)
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
                    isotopos_labels = [i for i in isotopes_label_list]
                else:
                    'This atom_label only has one heavy isotope'
                    isotopos_labels = [isotopes_label_list[0]]
                
                #all_atoms_list.extend(isotopos_labels) 
                isotopos_labels = [atom_label] + isotopos_labels
                
                all_atoms_list.extend(isotopos_labels)
                
                masses = [Atoms.atomic_masses.get(atom_label) for atom_label in isotopos_labels]
                props = [Atoms.isotopic_abundance.get(atom_label) for atom_label in isotopos_labels]
                
                atoms_count.append(formula_dict.get(atom_label))
                masses_list_tuples.append(masses)
                props_list_tuples.append(props)
        
        iso = IsoSpecPy.IsoSpec(atoms_count,masses_list_tuples,props_list_tuples, cut_off_to_IsoSpeccPy)
        
        conf = iso.getConfs()
        masses = conf[0]
        probs = exp(conf[1])
        molecular_formulas = conf[2]
        #print('conf', conf)
        #print('probs', conf[1])
        
        new_formulas = []
        
        for isotopologue_index in range(len(iso)):
            #skip_mono_isotopic 
            
            formula_list = molecular_formulas[isotopologue_index]
            new_formula_dict = dict(zip(all_atoms_list, formula_list))
            new_formula_dict[Labels.ion_type] = formula_dict.get(Labels.ion_type)
            if formula_dict.get('H'):
                new_formula_dict['H'] = formula_dict.get('H')

            new_formulas.append({x:y for x,y in new_formula_dict.items() if y!=0})
        
        # formula_dict in new_formulas check if monoisotopic is being returned
        if new_formulas:# and formula_dict in new_formulas:
            
            #print(conf)    
            #print(new_formulas)    
            #print(atoms_count)
            #print(all_atoms_list)
            #print(masses_list_tuples)
            #print(props_list_tuples)
            # find where monoisotopic is
            index_mono = new_formulas.index(formula_dict)   
            # calculate ratio iso/mono
            probs = list(probs/probs[index_mono])
            
            # delete the monoisotopic
            del probs[index_mono]
            del new_formulas[index_mono]
            
            #print('probs_exp', probs)
            for formulas, prob in zip(new_formulas, probs):
                
                theor_abundance = current_abundance* prob
                if theor_abundance > min_abundance:
                    #print(prob, theor_abundance, current_abundance)
                    yield (formulas, prob)
            #return zip(new_formulas, probs )
    
        #else:
        #    return []    
    