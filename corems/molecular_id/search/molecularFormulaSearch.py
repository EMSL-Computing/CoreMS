__author__ = "Yuri E. Corilo"
__date__ = "Jul 29, 2019"

import os, time
from os.path import join
from copy import deepcopy

from corems.encapsulation.constant import Labels
from corems.molecular_id.factory.MolecularLookupTable import  MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL as molform_db
#from corems.molecular_id.factory.molecularMongo import MolForm_Mongo as molform_db

last_error = 0
last_dif = 0
closest_error = 0
error_average = 0
nbValues = 0

class SearchMolecularFormulas:
     
    '''
    runworker()
    '''
    def __init__(self, first_hit=False, find_isotopologues=True):
        
        self.first_hit = first_hit
        self.find_isotopologues = find_isotopologues

    def __enter__(self):

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        
        return False

    def run_search(self, possible_formulas_dict, mass_spectrum_obj, min_abundance, is_adduct=False):
            
        all_assigned_indexes = list()

        for ms_peak in mass_spectrum_obj.sort_by_abundance():

            #already assigned a molecular formula
            if self.first_hit: 
                
                if ms_peak.is_assigned: continue
        
            nominal_mz  = ms_peak.nominal_mz_exp

            #get mono isotopic peaks that was added a molecular formula obj
            #TODO update error variables

            possible_formulas_dict_nm = possible_formulas_dict.get(nominal_mz)
            
            if possible_formulas_dict_nm:

                if is_adduct:
                    # change molecular formula from radical to adduct since the nitrogen number is the same 
                    for m_formula in possible_formulas_dict: m_formula.ion_type = Labels.adduct_ion

                ms_peak_indexes = SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).find_formulas(possible_formulas_dict_nm, min_abundance, mass_spectrum_obj, ms_peak)    

                all_assigned_indexes.extend(ms_peak_indexes)
        
        #filter per min peaks per mono isotopic class
        self.check_min_peaks(all_assigned_indexes, mass_spectrum_obj)
    
    def check_min_peaks(self, ms_peak_indexes, mass_spectrum_obj):
            
            if  mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter:

                if not len(ms_peak_indexes) >= mass_spectrum_obj.molecular_search_settings.min_peaks_per_class:
                    
                    for index in ms_peak_indexes: 
                        
                        mass_spectrum_obj[index].clear_molecular_formulas()

    def check_adduct_class(self, classe_dict, mass_spectrum_obj):
            
            return any([key in classe_dict.keys() for key in mass_spectrum_obj.molecular_search_settings.adduct_atoms_neg + mass_spectrum_obj.molecular_search_settings.adduct_atoms_pos])

    def run(self, classes, nominal_mz, min_abundance, 
            mass_spectrum_obj, ms_peak, dict_res):
       
        for classe_str, class_dict in classes:
            
            #is_adduct = self.check_adduct_class(classe_dict)        
            #print("classe", classe_str)
           
            #we might need to increase the search space to -+1 m_z 
            if mass_spectrum_obj.molecular_search_settings.isRadical or mass_spectrum_obj.molecular_search_settings.isAdduct:
                
                ion_type = Labels.radical_ion
                
                classes_formulas = dict_res.get(ion_type).get(classe_str)
                
                if classes_formulas: 
                    
                    possible_formulas = classes_formulas.get(nominal_mz)
                
                    if possible_formulas:
                        
                        is_adduct = self.check_adduct_class(possible_formulas[0].class_dict, mass_spectrum_obj)
                        
                        if not is_adduct and mass_spectrum_obj.molecular_search_settings.isRadical:
                            
                            if possible_formulas:
               
                                ms_peak_indexes = SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak)    
                                
                                self.check_min_peaks(ms_peak_indexes, mass_spectrum_obj)

                        elif is_adduct and mass_spectrum_obj.molecular_search_settings.isAdduct:
                           
                            #replace ion_type in the molecular_formula object
                            for m_formula in possible_formulas: m_formula.ion_type = Labels.adduct_ion

                            if possible_formulas:
               
                                ms_peak_indexes = SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak)    
                                
                                self.check_min_peaks(ms_peak_indexes, mass_spectrum_obj)

            if mass_spectrum_obj.molecular_search_settings.isProtonated:# and not is_adduct:
            
                ion_type = Labels.protonated_de_ion
                
                classes_formulas = dict_res.get(ion_type).get(classe_str)
                
                if classes_formulas: 
                    
                    possible_formulas = classes_formulas.get(nominal_mz)
                
                    if possible_formulas:
                        
                        is_adduct = self.check_adduct_class(possible_formulas[0].class_dict, mass_spectrum_obj)
                        
                        if not is_adduct:
                            
                            ms_peak_indexes = SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak)

                            self.check_min_peaks(ms_peak_indexes, mass_spectrum_obj)

    def run_worker_ms_peaks(self, ms_peaks, mass_spectrum_obj):

        
        #save initial settings min peaks per class filter 
        initial_min_peak_bool = deepcopy(mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter)

        #deactivate the usage of min peaks per class filter
        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False

        SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).reset_error(mass_spectrum_obj)

        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker(mass_spectrum_obj.molecular_search_settings)

        classes_str = [class_tuple[0] for class_tuple in classes]

        nominal_mzs = [ms_peak.nominal_mz_exp for ms_peak in  ms_peaks]

        dict_res = self.get_dict_molecular_database(classes_str, nominal_mzs, mass_spectrum_obj.molecular_search_settings)

        for ms_peak in  ms_peaks:

            if self.first_hit: 
                if ms_peak.is_assigned: continue
            
            nominal_mz  = ms_peak.nominal_mz_exp
           
            self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, dict_res)

        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = initial_min_peak_bool                
                        
    
    def run_worker_ms_peak(self, ms_peak, mass_spectrum_obj):
        
        #save initial settings min peaks per class filter 
        initial_min_peak_bool = deepcopy(mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter)

        #deactivate the usage of min peaks per class filter
        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False

        SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).reset_error(mass_spectrum_obj)
        
        min_abundance = 0.01

        classes = MolecularCombinations().runworker(mass_spectrum_obj.molecular_search_settings)
        
        nominal_mz = ms_peak.nominal_mz_exp

        classes_str = [class_tuple[0] for class_tuple in classes]

        dict_res = self.get_dict_molecular_database(classes_str, [nominal_mz],  mass_spectrum_obj.molecular_search_settings)
        
        self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, dict_res)

        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = initial_min_peak_bool
            
    def run_worker_mass_spectrum(self, mass_spectrum_obj):

        #number_of_process = multiprocessing.cpu_count()

        '''loading this on a shared memory would be better than having to serialize it for every process
            waiting for python 3.8 release'''
       
        SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).reset_error(mass_spectrum_obj)

        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker(mass_spectrum_obj.molecular_search_settings)
        
        nominal_mzs = mass_spectrum_obj.nominal_mz

        classes_str = [class_tuple[0] for class_tuple in classes]

        #query database
        dict_res = self.get_dict_molecular_database(classes_str, nominal_mzs, mass_spectrum_obj.molecular_search_settings)
        
        for classe_tuple in classes:
                
            classe_str  = classe_tuple[0]
            classe_dict = classe_tuple[1]

            is_adduct = self.check_adduct_class(classe_dict, mass_spectrum_obj)    
            
            if mass_spectrum_obj.molecular_search_settings.isProtonated and not is_adduct:
        
                    ion_type = Labels.protonated_de_ion

                    possible_formulas = dict_res.get(ion_type).get(classe_str)
                    
                    if possible_formulas:

                        self.run_search(possible_formulas, mass_spectrum_obj, min_abundance)    

            if mass_spectrum_obj.molecular_search_settings.isRadical and not is_adduct:
                
                    ion_type = Labels.radical_ion
                    
                    possible_formulas = dict_res.get(ion_type).get(classe_str)
                   
                    if possible_formulas:

                        self.run_search(possible_formulas, mass_spectrum_obj, min_abundance)    

            # looks for adduct, used_atom_valences should be 0 
            # this code does not support H exchance by halogen atoms
            if mass_spectrum_obj.molecular_search_settings.isAdduct and is_adduct:
                
                ion_type = Labels.radical_ion
                
                possible_formulas = dict_res.get(ion_type).get(classe_str)
                
                if possible_formulas:
                    
                    #replace ion_type in the molecular_formula object
                    self.run_search(possible_formulas, mass_spectrum_obj, min_abundance, is_adduct=is_adduct)          

    def search_mol_formulas(self, mass_spectrum_obj, possible_formulas_list, find_isotopologues=True):

        SearchMolecularFormulaWorker(find_isotopologues=find_isotopologues).reset_error(mass_spectrum_obj)

        initial_min_peak_bool = mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter
        
        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
        
        possible_formulas_dict_nm =  {}
        
        for mf in possible_formulas_list:
            
            nm = mf.mz_nominal_theo
            
            if nm in possible_formulas_dict_nm.keys():
                possible_formulas_dict_nm[nm].append(mf)
            
            else:    
                possible_formulas_dict_nm[nm] = [mf]

        min_abundance = mass_spectrum_obj.min_abundance

        self.run_search(possible_formulas_dict_nm, mass_spectrum_obj, min_abundance)          

        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = initial_min_peak_bool
        mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = initial_min_peak_bool

        mspeaks = [mspeak for mspeak in mass_spectrum_obj if mspeak.is_assigned]
        
        return mspeaks

    def get_dict_molecular_database(self, classes_str, nominal_mzs, molecular_search_settings):
            
        dict_res = {}
        
        #print (classes_str)
        if molecular_search_settings.isProtonated:
            
            ion_type = Labels.protonated_de_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs, molecular_search_settings)

        if molecular_search_settings.isRadical:

            ion_type = Labels.radical_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs,  molecular_search_settings)
        
        return dict_res
                
            
class SearchMolecularFormulaWorker:
    
    #TODO add reset error function
    # needs this warper to pass the class to multiprocessing
    
    def __init__(self, find_isotopologues=True):
        self.find_isotopologues = find_isotopologues
    
    def __call__(self, args):

        return self.find_formulas(*args)  # ,args[1]

    def reset_error(self, mass_spectrum_obj):
        global last_error, last_dif, closest_error, error_average, nbValues  
        last_error, last_dif, closest_error, nbValues  = 0.0, 0.0, 0.0, 0.0
        
        error_average = 0

    def set_last_error(self, error, mass_spectrum_obj ):
        
        #set the changes to the global variables, not internal ones
        global last_error, last_dif, closest_error, error_average, nbValues  
        
        if mass_spectrum_obj.molecular_search_settings.error_method == 'distance':
            
            dif = error - last_error
            if dif < last_dif:
                last_dif = dif
                closest_error = error
                mass_spectrum_obj.molecular_search_settings.min_mz_error = closest_error - mass_spectrum_obj.molecular_search_settings.mz_error_range
                mass_spectrum_obj.molecular_search_settings.max_mz_error = closest_error + mass_spectrum_obj.molecular_search_settings.mz_error_range

        elif mass_spectrum_obj.molecular_search_settings.error_method == 'lowest':
            
            if error < last_error:
                mass_spectrum_obj.molecular_search_settings.min_mz_error = error - mass_spectrum_obj.molecular_search_settings.mz_error_range
                mass_spectrum_obj.molecular_search_settings.max_mz_error = error + mass_spectrum_obj.molecular_search_settings.mz_error_range
                last_error = error
                
        
        elif mass_spectrum_obj.molecular_search_settings.error_method == 'symmetrical':
               
               mass_spectrum_obj.molecular_search_settings.min_mz_error = mass_spectrum_obj.molecular_search_settings.mz_error_average - mass_spectrum_obj.molecular_search_settings.mz_error_range
               mass_spectrum_obj.molecular_search_settings.max_mz_error = mass_spectrum_obj.molecular_search_settings.mz_error_average + mass_spectrum_obj.molecular_search_settings.mz_error_range
        
        elif mass_spectrum_obj.molecular_search_settings.error_method == 'average':

                nbValues += 1
                error_average = error_average + ((error - error_average) / nbValues)
                mass_spectrum_obj.molecular_search_settings.min_mz_error =  error_average - mass_spectrum_obj.molecular_search_settings.mz_error_range
                mass_spectrum_obj.molecular_search_settings.max_mz_error =  error_average + mass_spectrum_obj.molecular_search_settings.mz_error_range    
                
                
        else:
            #using set mass_spectrum_obj.molecular_search_settings.min_mz_error and max_mz_error range
            pass

        '''returns the error based on the selected method at mass_spectrum_obj.molecular_search_settings.method
        '''    
        
        
    def find_formulas(self, possible_formulas, min_abundance, 
                      mass_spectrum_obj, ms_peak):
        '''
        # uses the closest error the next search (this is not ideal, it needs to use confidence
        # metric to choose the right candidate then propagate the error using the error from the best candidate
        # it needs to add s/n to the equation
        # it need optimization to define the mz_error_range within a m/z unit since it is directly 
        # proportional with the mass, and inversely proportional to the rp. 
        # It's not linear, i.e., sigma ∝ mass 
        # the idea it to correlate sigma to resolving power, signal to noise and sample complexity per mz unit
        # method='distance'
        '''

        
        mspeak_assigned_index = list()

        min_mz_error = mass_spectrum_obj.molecular_search_settings.min_mz_error
        max_mz_error = mass_spectrum_obj.molecular_search_settings.max_mz_error
        
        min_abun_error = mass_spectrum_obj.molecular_search_settings.min_abun_error
        max_abun_error = mass_spectrum_obj.molecular_search_settings.max_abun_error
        
        #f = open("abundance_error.txt", "a+")    
        ms_peak_mz_exp, ms_peak_abundance = ms_peak.mz_exp, ms_peak.abundance
        #min_error = min([pmf._calc_assignment_mass_error(ms_peak_mz_exp) for pmf in possible_formulas])
        
        for possible_formula in possible_formulas:
            
            if possible_formula:
                
                error = possible_formula._calc_assignment_mass_error(ms_peak_mz_exp)
               
                if  min_mz_error <= error <= max_mz_error:
                    
                    #update the error
                    
                    self.set_last_error(error, mass_spectrum_obj)    
                   
                    #add molecular formula match to ms_peak
                    
                    possible_formula.set_assignment_mass_error(error)

                    ms_peak.add_molecular_formula(possible_formula)
                    
                    mspeak_assigned_index.append(ms_peak.index)
                    
                    if self.find_isotopologues:
                        #calculates and look for isotopologues
                        isotopologues = possible_formula.isotopologues(min_abundance, ms_peak_abundance)
                        
                        for isotopologue_formula in isotopologues:
                           
                            possible_formula.expected_isotopologues.append(isotopologue_formula)
                            #move this outside to improve preformace
                            #we need to increase the search space to -+1 m_z 
                            first_index, last_index = mass_spectrum_obj.get_nominal_mz_first_last_indexes(isotopologue_formula.mz_nominal_theo)
                            
                            for ms_peak_iso in mass_spectrum_obj[first_index:last_index]:
                                
                                error = isotopologue_formula._calc_assignment_mass_error(ms_peak_iso.mz_exp)    
                                
                                if  min_mz_error <= error <= max_mz_error:
                                    
                                    #need to define error distribution for abundance measurements
                                    
                                    abundance_error = isotopologue_formula._calc_abundance_error(ms_peak_abundance,ms_peak_iso.abundance )            
                                    # margin of error was set empirically/ needs statistical calculation
                                    #  of margin of error for the measurement of the abundances
                                    if min_abun_error <= abundance_error <= max_abun_error:
                                        
                                        #update the error   
                                        
                                        self.set_last_error(error, mass_spectrum_obj)    
                                        
                                        isotopologue_formula.set_assignment_mass_error(error)

                                        #add molecular formula match to ms_peak
                                        ms_peak_iso.add_molecular_formula(isotopologue_formula)
                                        
                                        #add mspeaks mono isotopic index to the isotopologue MolecularFormula obj
                                        isotopologue_formula.mspeak_index_mono_isotopic = ms_peak.index
                                        
                                        #add mspeaks isotopologue index to the mono isotopic MolecularFormula obj
                                        possible_formula.mspeak_indexes_isotopologues.append(ms_peak_iso.index)

        return mspeak_assigned_index

