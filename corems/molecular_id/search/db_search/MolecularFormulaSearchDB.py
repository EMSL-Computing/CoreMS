__author__ = "Yuri E. Corilo"
__date__ = "Jul 29, 2019"


import os, time
from os.path import join

from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupDictSettings, MolecularSearchSettings
from corems.mass_spectrum.input.massList import Read_MassList
from corems.molecular_id.factory.db_search.MolecularLookupTableDB import  MolecularCombinations
from corems.molecular_id.factory.db_search.molecularSQL import MolForm_SQL as molform_db
#from corems.molecular_id.factory.db_search.molecularMongo import MolForm_Mongo as molform_db

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

    def run(self, classes, nominal_mz, min_abundance, 
            mass_spectrum_obj, ms_peak, last_error, last_dif, 
            closest_error, error_average, nbValues, dict_res):
        
        def check_adduct_class(classe_dict):
            
            return any([key in classe_dict.keys() for key in MolecularSearchSettings.adduct_atoms_neg])
       
        for classe_str, class_dict in classes:
            
            #is_adduct = check_adduct_class(classe_dict)        
            #print("classe", classe_str)
            possible_formulas = list()    
            
            #we might need to increase the search space to -+1 m_z 
            if MolecularSearchSettings.isRadical or MolecularSearchSettings.isAdduct:
            
                ion_type = Labels.radical_ion
                
                formulas = dict_res.get(ion_type).get(classe_str).get(nominal_mz)
                #print('query took %i seconds' % (time.time() - time1))
                
                if formulas:
                    
                    is_adduct = check_adduct_class(formulas[0].class_dict)
                    
                    if is_adduct and MolecularSearchSettings.isAdduct:
                        
                        #replace ion_type in the molecular_formula object
                        for m_formula in formulas: m_formula.ion_type = Labels.adduct_ion

                        possible_formulas.extend(formulas)
                    
                    elif not is_adduct and MolecularSearchSettings.isRadical:
                        
                        possible_formulas.extend(formulas)

            if MolecularSearchSettings.isProtonated:# and not is_adduct:
            
                ion_type = Labels.protonated_de_ion
                
                formulas = dict_res.get(ion_type).get(classe_str).get(nominal_mz)
                
                if formulas:
                    
                    is_adduct = check_adduct_class(formulas[0].class_dict)
                    
                    if not is_adduct:
                        
                        possible_formulas.extend(formulas)
            #print(formulas)
            if possible_formulas:
                
               
                SearchMolecularFormulaWorker().find_formulas(possible_formulas, min_abundance, mass_spectrum_obj, ms_peak, last_error, last_dif, closest_error, error_average, nbValues)
    
    def run_worker_ms_peaks(self, ms_peaks, mass_spectrum_obj):

        last_dif = 0
        
        last_error = 0
        
        closest_error = 0

        error_average = MolecularSearchSettings.mz_error_average
        
        nbValues = 0

        MolecularSearchSettings.min_mz =  min(ms_peaks, key=lambda m: m.mz_exp).mz_exp
    
        MolecularSearchSettings.max_mz = max(ms_peaks, key=lambda m: m.mz_exp).mz_exp
        
        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker()

        classes_str = [class_tuple[0] for class_tuple in classes]

        nominal_mzs = [ms_peak.nominal_mz_exp for ms_peak in  ms_peaks]

        dict_res = self.get_dict_molecular_database(classes_str, nominal_mzs)

        for ms_peak in  ms_peaks:

            if self.first_hit:
                #print('hell yeah')
                if ms_peak.is_assigned: continue
            
            nominal_mz  = ms_peak.nominal_mz_exp
            '''
            waiting for python 3.8 release to set mass_spectrum_obj and dict_molecular_lookup_table on share memory
            pool = multiprocessing.Pool(number_of_process)
            args = [ (dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz), min_abundance, mass_spectrum_obj, ms_peak_mz_exp, ms_peak_abundance)  for classe in classes ]
            pool.map(SearchMolecularFormulaWorker(), args)

            pool.close()
            pool.join()
            '''
            self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, last_error, last_dif, 
                        closest_error, error_average, nbValues, dict_res)
                        
    def run_worker_ms_peak(self, ms_peak, mass_spectrum_obj):
        
        '''
        waiting for python 3.8 release to set mass_spectrum_obj and dict_molecular_lookup_table on share memory (redis?)
        pool = multiprocessing.Pool(number_of_process)
        args = [ (dict_molecular_lookup_table.get(classe).get(ion_type).get(nominal_mz), min_abundance, mass_spectrum_obj, ms_peak_mz_exp, ms_peak_abundance)  for classe in classes ]
        pool.map(SearchMolecularFormulaWorker(), args)

        pool.close()
        pool.join()
        '''

        last_dif = 0
        
        last_error = 0
        
        closest_error = 0

        error_average = MolecularSearchSettings.mz_error_average
        
        nbValues = 0

        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker()
        
        nominal_mz  = [ms_peak.nominal_mz_exp]  
        
        classes_str = [class_tuple[0] for class_tuple in classes]

        dict_res = self.get_dict_molecular_database(classes_str, nominal_mz)
        
        self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, last_error, last_dif, 
                        closest_error, error_average, nbValues, dict_res)
            
    def run_worker_mass_spectrum(self, mass_spectrum_obj):

        #number_of_process = multiprocessing.cpu_count()

        '''loading this on a shared memory would be better than having to serialize it for every process
            waiting for python 3.8 release'''
        last_dif = 0
        
        last_error = 0
        
        closest_error = 0

        error_average = MolecularSearchSettings.mz_error_average
        
        nbValues = 0

        min_abundance = mass_spectrum_obj.min_abundance

        classes = MolecularCombinations().runworker()
        
        nominal_mzs = mass_spectrum_obj.nominal_mz

        classes_str = [class_tuple[0] for class_tuple in classes]

        dict_res = self.get_dict_molecular_database(classes_str, nominal_mzs)
        #print(dict_res)
        for ms_peak in mass_spectrum_obj.sort_by_abundance():
            
            if self.first_hit:
                    #print('hell yeah')
                    if ms_peak.is_assigned: continue
                
            nominal_mz  = ms_peak.nominal_mz_exp
                        
            self.run(classes, nominal_mz, min_abundance, 
                        mass_spectrum_obj, ms_peak, last_error, last_dif, 
                        closest_error, error_average, nbValues, dict_res)

    def get_dict_molecular_database(self, classes_str, nominal_mzs):
            
        dict_res = {}
        
        #print (classes_str)
        if MolecularSearchSettings.isProtonated:
            
            ion_type = Labels.protonated_de_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs)

        if MolecularSearchSettings.isRadical:

            ion_type = Labels.radical_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs)
        
        return dict_res
                
            
class SearchMolecularFormulaWorker:
    
    # needs this wraper to pass the class to multiprocessing
    def __call__(self, args):

        return self.find_formulas(*args)  # ,args[1]

    def set_last_error(self, error, last_error, last_dif, closest_error, error_average, nbValues ):
        
        
        if MolecularSearchSettings.error_method == 'distance':
            
            dif = error - last_error
            if dif < last_dif:
                last_dif = dif
                closest_error = error
                MolecularSearchSettings.min_mz_error = closest_error - MolecularSearchSettings.mz_error_range
                MolecularSearchSettings.max_mz_error = closest_error + MolecularSearchSettings.mz_error_range

        elif MolecularSearchSettings.error_method == 'lowest':
            
            if error < last_error:
                MolecularSearchSettings.min_mz_error = error - MolecularSearchSettings.mz_error_range
                MolecularSearchSettings.max_mz_error = error + MolecularSearchSettings.mz_error_range
                last_error = error
                
        
        elif MolecularSearchSettings.error_method == 'symmetrical':
               
               MolecularSearchSettings.min_mz_error = MolecularSearchSettings.mz_error_average - MolecularSearchSettings.mz_error_range
               MolecularSearchSettings.max_mz_error = MolecularSearchSettings.mz_error_average + MolecularSearchSettings.mz_error_range
        
        elif MolecularSearchSettings.error_method == 'average':

                nbValues += 1
                error_average = error_average + ((error - error_average) / nbValues)
                MolecularSearchSettings.min_mz_error =  error_average - MolecularSearchSettings.mz_error_range
                MolecularSearchSettings.max_mz_error =  error_average + MolecularSearchSettings.mz_error_range    
                
                
        else:
            #using set MolecularSearchSettings.min_mz_error and max_mz_error range
            pass

        '''returns the error based on the selected method at MolecularSearchSettings.method
        '''    
        return last_error, last_dif, closest_error, error_average, nbValues        
        
    def find_formulas(self, possible_formulas, min_abundance, 
                      mass_spectrum_obj, ms_peak, 
                      last_error, last_dif, 
                      closest_error, error_average, 
                      nbValues):
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

        min_mz_error = MolecularSearchSettings.min_mz_error
        max_mz_error = MolecularSearchSettings.max_mz_error
        min_abun_error = MolecularSearchSettings.min_abun_error
        max_abun_error = MolecularSearchSettings.max_abun_error
        max_dbe = MolecularSearchSettings.max_dbe
        min_dbe = MolecularSearchSettings.min_dbe
        
        #f = open("abundance_error.txt", "a+")    
        ms_peak_mz_exp, ms_peak_abundance = ms_peak.mz_exp, ms_peak.abundance
        #min_error = min([pmf._calc_assigment_mass_error(ms_peak_mz_exp) for pmf in possible_formulas])
        
        for possible_formula in possible_formulas:
            
            if possible_formula:
                
                error = possible_formula._calc_assigment_mass_error(ms_peak_mz_exp)
                
                if  min_mz_error <= error <= max_mz_error and min_dbe <= possible_formula.dbe <= max_dbe:
                    
                    #update the error
                    last_error, last_dif, closest_error, error_average, nbValues  = self.set_last_error(error, last_error, last_dif, closest_error, error_average, nbValues)    
                    
                    #add molecular formula match to ms_peak
                    ms_peak.add_molecular_formula(possible_formula)
                    
                   
                    #calculates and look for isotopologues
                    isotopologues = possible_formula.isotopologues(min_abundance, ms_peak_abundance)
                    
                    for isotopologue_formula in isotopologues:
                        
                        #move this outside to impove preformace
                        #we need to increase the search space to -+1 m_z 
                        first_index, last_index = mass_spectrum_obj.get_nominal_mz_frist_last_indexes(isotopologue_formula.mz_nominal_theo)
                        
                        for ms_peak_iso in mass_spectrum_obj[first_index:last_index]:
                            
                            error = isotopologue_formula._calc_assigment_mass_error(ms_peak_iso.mz_exp)    
                            
                            #need to define error distribution for abundance measurements
                            if  min_mz_error <= error <= max_mz_error and min_dbe <= isotopologue_formula.dbe <= max_dbe:
                                    
                                    abundance_error = isotopologue_formula._calc_abundance_error(ms_peak_abundance,ms_peak_iso.abundance )            
                                    # margin of error was set empirically/ needs statistical calculation
                                    #  of margin of error for the measurement of the abundances
                                    if min_abun_error <= abundance_error <= max_abun_error:
                                        
                                        #update the error   
                                        last_error, last_dif, closest_error, error_average, nbValues  = self.set_last_error(error, last_error, last_dif, closest_error, error_average, nbValues)    
                                        
                                        #add molecular formula match to ms_peak
                                        ms_peak_iso.add_molecular_formula(isotopologue_formula)
                                        
                                        #add mspeaks mono isotopic index to the isotopologue MolecularFormula obj
                                        isotopologue_formula.mspeak_index_mono_isotopic = ms_peak.index
                                        
                                        #add mspeaks isotopologue index to the mono isotopic MolecularFormula obj
                                        possible_formula.mspeak_indexes_isotopologues.append(ms_peak_iso.index)
                                        

                                    