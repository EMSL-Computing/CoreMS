__author__ = "Yuri E. Corilo"
__date__ = "Jul 29, 2019"

import os, time
from os.path import join
from copy import deepcopy

from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupDictSettings, MoleculaSearchSettings
from corems.molecular_id.factory.MolecularLookupTable import  MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL as molform_db
#from corems.molecular_id.factory.molecularMongo import MolForm_Mongo as molform_db


global last_error 
global last_dif 
global closest_error 
global error_average
global nbValues 

class SearchMolecularFormulas:
     
    '''
    runworker()
    '''
    def __init__(self, first_hit=False):
        
        self.first_hit = first_hit

    def __enter__(self):

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        
        return False

    def run_search(self, possible_formulas_dict, mass_spectrum_obj, min_abundance, is_adduct=False):
            
        all_assigned_indexes = list()

        for ms_peak in mass_spectrum_obj.sort_by_abundance():

            #already assinged a molecular formula
            if self.first_hit: 
                if ms_peak.is_assigned: continue
        
            nominal_mz  = ms_peak.nominal_mz_exp

            #get mono isotopic peaks that was added a molecular formula obj
            #TODO update error variables

            possible_formulas_nominal = possible_formulas_dict.get(nominal_mz)
            
            if possible_formulas_nominal:

                if is_adduct:
                    
                    for m_formula in possible_formulas_dict: m_formula.ion_type = Labels.adduct_ion

                ms_peak_indexes = SearchMolecularFormulaWorker().find_formulas(possible_formulas_nominal, min_abundance, mass_spectrum_obj, ms_peak)    

                all_assigned_indexes.extend(ms_peak_indexes)
        
        #filter per min peaks per mono isotopic class
        self.check_min_peaks(all_assigned_indexes, mass_spectrum_obj)
    
    def check_min_peaks(self, ms_peak_indexes, mass_spectrum_obj):
            
            if  MoleculaSearchSettings.use_min_peaks_filter:

                if not len(ms_peak_indexes) >= MoleculaSearchSettings.min_peaks_per_class:
                    
                    for index in ms_peak_indexes: 
                        
                        mass_spectrum_obj[index].clear_molecular_formulas()

    def check_adduct_class(self, classe_dict):
            
            return any([key in classe_dict.keys() for key in MoleculaSearchSettings.adduct_atoms_neg + MoleculaSearchSettings.adduct_atoms_pos])

    def run(self, classes, nominal_mz, min_abundance, 
            mass_spectrum_obj, ms_peak, dict_res):
       
        for classe_str, class_dict in classes:
            
            #is_adduct = self.check_adduct_class(classe_dict)        
            #print("classe", classe_str)
           
            #we might need to increase the search space to -+1 m_z 
            if MoleculaSearchSettings.isRadical or MoleculaSearchSettings.isAdduct:
                
                ion_type = Labels.radical_ion
                
                classes_formulas = dict_res.get(ion_type).get(classe_str)
                
                if classes_formulas: 
                    
                    possible_formulas = classes_formulas.get(nominal_mz)
                
                    if possible_formulas:
                        
                        is_adduct = self.check_adduct_class(possible_formulas[0].class_dict)
                        
                        if not is_adduct and MoleculaSearchSettings.isRadical:
                            
                            if possible_formulas:
               
                                ms_peak_indexes = SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak)    
                                
                                self.check_min_peaks(ms_peak_indexes, mass_spectrum_obj)

                        elif is_adduct and MoleculaSearchSettings.isAdduct:
                           
                            #replace ion_type in the molecular_formula object
                            for m_formula in possible_formulas: m_formula.ion_type = Labels.adduct_ion

                            if possible_formulas:
               
                                ms_peak_indexes = SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak)    
                                
                                self.check_min_peaks(ms_peak_indexes, mass_spectrum_obj)

            if MoleculaSearchSettings.isProtonated:# and not is_adduct:
            
                ion_type = Labels.protonated_de_ion
                
                classes_formulas = dict_res.get(ion_type).get(classe_str)
                
                if classes_formulas: 
                    
                    possible_formulas = classes_formulas.get(nominal_mz)
                
                    if possible_formulas:
                        
                        is_adduct = self.check_adduct_class(possible_formulas[0].class_dict)
                        
                        if not is_adduct:
                            
                            ms_peak_indexes = SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak)

                            self.check_min_peaks(ms_peak_indexes, mass_spectrum_obj)

    def run_worker_ms_peaks(self, ms_peaks, mass_spectrum_obj):

        
        #save initial settings min peaks per class filter 
        initial_min_peak_bool = deepcopy(MoleculaSearchSettings.use_min_peaks_filter)

        #deactivate the usage of min peaks per class filter
        MoleculaSearchSettings.use_min_peaks_filter = False

        SearchMolecularFormulaWorker().reset_error()

        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker()

        classes_str = [class_tuple[0] for class_tuple in classes]

        nominal_mzs = [ms_peak.nominal_mz_exp for ms_peak in  ms_peaks]

        dict_res = self.get_dict_molecular_database(classes_str, nominal_mzs)

        for ms_peak in  ms_peaks:

            if self.first_hit: 
                if ms_peak.is_assigned: continue
            
            nominal_mz  = ms_peak.nominal_mz_exp
           
            self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, dict_res)

        MoleculaSearchSettings.use_min_peaks_filter = initial_min_peak_bool                
                        
    def run_worker_ms_peak(self, ms_peak, mass_spectrum_obj):
        
        #save initial settings min peaks per class filter 
        initial_min_peak_bool = deepcopy(MoleculaSearchSettings.use_min_peaks_filter)

        #deactivate the usage of min peaks per class filter
        MoleculaSearchSettings.use_min_peaks_filter = False

        SearchMolecularFormulaWorker().reset_error()
        
        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker()
        
        nominal_mz = ms_peak.nominal_mz_exp

        classes_str = [class_tuple[0] for class_tuple in classes]

        dict_res = self.get_dict_molecular_database(classes_str, [nominal_mz])
        
        self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, dict_res)

        MoleculaSearchSettings.use_min_peaks_filter = initial_min_peak_bool
            
    def run_worker_mass_spectrum(self, mass_spectrum_obj):

        #number_of_process = multiprocessing.cpu_count()

        '''loading this on a shared memory would be better than having to serialize it for every process
            waiting for python 3.8 release'''
       
        SearchMolecularFormulaWorker().reset_error()

        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker()
        
        nominal_mzs = mass_spectrum_obj.nominal_mz

        classes_str = [class_tuple[0] for class_tuple in classes]

        #query database
        dict_res = self.get_dict_molecular_database(classes_str, nominal_mzs)
        
        for classe_tuple in classes:
                
            classe_str  = classe_tuple[0]
            classe_dict = classe_tuple[1]

            is_adduct = self.check_adduct_class(classe_dict)    
            
            if MoleculaSearchSettings.isProtonated and not is_adduct:
        
                    ion_type = Labels.protonated_de_ion

                    possible_formulas = dict_res.get(ion_type).get(classe_str)
                    
                    if possible_formulas:

                        self.run_search(possible_formulas, mass_spectrum_obj, min_abundance)    

            if MoleculaSearchSettings.isRadical and not is_adduct:
                
                    ion_type = Labels.radical_ion
                    
                    possible_formulas = dict_res.get(ion_type).get(classe_str)
                   
                    if possible_formulas:

                        self.run_search(possible_formulas, mass_spectrum_obj, min_abundance)    

            # looks for adduct, used_atom_valences should be 0 
            # this code does not support H exchance by halogen atoms
            if MoleculaSearchSettings.isAdduct and is_adduct:
                
                ion_type = Labels.radical_ion
                
                possible_formulas = dict_res.get(ion_type).get(classe_str)
                
                if possible_formulas:
                    
                    #replace ion_type in the molecular_formula object
                    self.run_search(possible_formulas, mass_spectrum_obj, min_abundance, is_adduct=is_adduct)          

    def search_mol_formulas(self, mass_spectrum_obj, possible_formulas_list, find_isotopologues=True):

        
        SearchMolecularFormulaWorker(find_isotopologues=find_isotopologues).reset_error()

        initial_min_peak_bool = deepcopy(MoleculaSearchSettings.use_min_peaks_filter)
        
        MoleculaSearchSettings.use_min_peaks_filter = False

        possible_formulas_dict_nm =  {}
        for mf in possible_formulas_list:
            nm = mf.mz_nominal_theo
            if nm in possible_formulas_dict_nm.keys():
                possible_formulas_dict_nm[nm].append(mf)
            else:    
                possible_formulas_dict_nm[nm] = [mf]

        min_abundance = mass_spectrum_obj.min_abundance

        self.run_search(possible_formulas_dict_nm, mass_spectrum_obj, min_abundance)          

        MoleculaSearchSettings.use_min_peaks_filter = initial_min_peak_bool

        mspeaks = [mspeak for mspeak in mass_spectrum_obj if mspeak.is_assigned]
        
        return mspeaks

    def get_dict_molecular_database(self, classes_str, nominal_mzs):
            
        dict_res = {}
        
        #print (classes_str)
        if MoleculaSearchSettings.isProtonated:
            
            ion_type = Labels.protonated_de_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs)

        if MoleculaSearchSettings.isRadical:

            ion_type = Labels.radical_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs)
        
        return dict_res
                
            
class SearchMolecularFormulaWorker:
    
    #TODO add reset erro function
    # needs this wraper to pass the class to multiprocessing
    
    def __init__(self,find_isotopologues=True):
        self.find_isotopologues = find_isotopologues
    
    def __call__(self, args):

        return self.find_formulas(*args)  # ,args[1]

    def reset_error(self):
        global last_error, last_dif, closest_error, error_average, nbValues  
        last_error, last_dif, closest_error, nbValues  = 0.0, 0.0, 0.0, 0.0
        error_average = MoleculaSearchSettings.mz_error_average

    def set_last_error(self, error ):
        
        #set the changes to the global variables, not internal ones
        global last_error, last_dif, closest_error, error_average, nbValues  
        
        if MoleculaSearchSettings.error_method == 'distance':
            
            dif = error - last_error
            if dif < last_dif:
                last_dif = dif
                closest_error = error
                MoleculaSearchSettings.min_mz_error = closest_error - MoleculaSearchSettings.mz_error_range
                MoleculaSearchSettings.max_mz_error = closest_error + MoleculaSearchSettings.mz_error_range

        elif MoleculaSearchSettings.error_method == 'lowest':
            
            if error < last_error:
                MoleculaSearchSettings.min_mz_error = error - MoleculaSearchSettings.mz_error_range
                MoleculaSearchSettings.max_mz_error = error + MoleculaSearchSettings.mz_error_range
                last_error = error
                
        
        elif MoleculaSearchSettings.error_method == 'symmetrical':
               
               MoleculaSearchSettings.min_mz_error = MoleculaSearchSettings.mz_error_average - MoleculaSearchSettings.mz_error_range
               MoleculaSearchSettings.max_mz_error = MoleculaSearchSettings.mz_error_average + MoleculaSearchSettings.mz_error_range
        
        elif MoleculaSearchSettings.error_method == 'average':

                nbValues += 1
                error_average = error_average + ((error - error_average) / nbValues)
                MoleculaSearchSettings.min_mz_error =  error_average - MoleculaSearchSettings.mz_error_range
                MoleculaSearchSettings.max_mz_error =  error_average + MoleculaSearchSettings.mz_error_range    
                
                
        else:
            #using set MoleculaSearchSettings.min_mz_error and max_mz_error range
            pass

        '''returns the error based on the selected method at MoleculaSearchSettings.method
        '''    
        
        
    def find_formulas(self, possible_formulas, min_abundance, 
                      mass_spectrum_obj, ms_peak):
        '''
        # uses the closest error the next search (this is not ideal, it needs to use confidence
        # metric to choose the right candidate then propagate the error using the error from the best candidate
        # it needs to add s/n to the equation
        # it need optimization to define the mz_error_range within a m/z unit since it is directly 
        # proportional with the mass, and inversially proportinal to the rp. 
        # It's not linear, i.e., sigma ∝ mass 
        # the idea it to correlate sigma to resolving power, signal to noise and sample complexity per mz unit
        # method='distance'
        '''

        
        mspeak_assigned_index = list()

        min_mz_error = MoleculaSearchSettings.min_mz_error
        max_mz_error = MoleculaSearchSettings.max_mz_error
        
        min_abun_error = MoleculaSearchSettings.min_abun_error
        max_abun_error = MoleculaSearchSettings.max_abun_error
        
        #f = open("abundance_error.txt", "a+")    
        ms_peak_mz_exp, ms_peak_abundance = ms_peak.mz_exp, ms_peak.abundance
        #min_error = min([pmf._calc_assigment_mass_error(ms_peak_mz_exp) for pmf in possible_formulas])
        
        for possible_formula in possible_formulas:
            
            if possible_formula:
                
                error = possible_formula._calc_assigment_mass_error(ms_peak_mz_exp)
                
                if  min_mz_error <= error <= max_mz_error:
                    
                    #update the error
                    
                    self.set_last_error(error)    
                   
                    #add molecular formula match to ms_peak
                    ms_peak.add_molecular_formula(possible_formula)
                    
                    mspeak_assigned_index.append(ms_peak.index)
                    
                    if self.find_isotopologues:
                        #calculates and look for isotopologues
                        isotopologues = possible_formula.isotopologues(min_abundance, ms_peak_abundance)
                        
                        for isotopologue_formula in isotopologues:
                            
                            #move this outside to impove preformace
                            #we need to increase the search space to -+1 m_z 
                            first_index, last_index = mass_spectrum_obj.get_nominal_mz_frist_last_indexes(isotopologue_formula.mz_nominal_theo)
                            
                            for ms_peak_iso in mass_spectrum_obj[first_index:last_index]:
                                
                                error = isotopologue_formula._calc_assigment_mass_error(ms_peak_iso.mz_exp)    
                                
                                #need to define error distribution for abundance measurements
                                
                                abundance_error = isotopologue_formula._calc_abundance_error(ms_peak_abundance,ms_peak_iso.abundance )            
                                # margin of error was set empirically/ needs statistical calculation
                                #  of margin of error for the measurement of the abundances
                                if min_abun_error <= abundance_error <= max_abun_error:
                                    
                                    #update the error   
                                    
                                    self.set_last_error(error)    
                                    
                                    #add molecular formula match to ms_peak
                                    ms_peak_iso.add_molecular_formula(isotopologue_formula)
                                    
                                    #add mspeaks mono isotopic index to the isotopologue MolecularFormula obj
                                    isotopologue_formula.mspeak_index_mono_isotopic = ms_peak.index
                                    
                                    #add mspeaks isotopologue index to the mono isotopic MolecularFormula obj
                                    possible_formula.mspeak_indexes_isotopologues.append(ms_peak_iso.index)

        return mspeak_assigned_index

