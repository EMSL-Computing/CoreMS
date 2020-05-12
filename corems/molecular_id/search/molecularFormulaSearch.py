__author__ = "Yuri E. Corilo"
__date__ = "Jul 29, 2019"

import os, time
from os.path import join
from copy import deepcopy
import json 

import tqdm
from sqlalchemy.types import Binary

from corems.encapsulation.constant import Atoms, Labels
from corems.molecular_id.factory.CreateMolFormDB import NewMolecularCombinations
from corems.molecular_id.factory.molformSQL import MolForm_SQL, MolecularFormulaLink
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.molecular_id.calc.MolecularFilter import MolecularFormulaSearchFilters
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems import chunks, timeit
from sqlalchemy.sql.sqltypes import Integer
import multiprocessing

last_error = 0
last_dif = 0
closest_error = 0
error_average = 0
nbValues = 0

class SearchMolecularFormulas:
     
    '''
    runworker()
    '''
    def __init__(self, mass_spectrum_obj, sql_db=None, first_hit=False, find_isotopologues=True):

        self.first_hit = first_hit
        
        self.find_isotopologues = find_isotopologues

        self.mass_spectrum_obj = mass_spectrum_obj
        
        if not sql_db:

            self.sql_db = MolForm_SQL(url=mass_spectrum_obj.molecular_search_settings.url_database)
        
        else:
            
            self.sql_db = sql_db

    def __enter__(self):

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        
        self.sql_db.close()
        
        return False

    def run_search(self, query, min_abundance, ion_type):

        def get_formulas(nominal_overlay=0.1):
            
            nominal_mz  = ms_peak.nominal_mz_exp

            defect_mass = ms_peak.mz_exp - nominal_mz
            nominal_masses = [nominal_mz]

            if (defect_mass) >= 1-nominal_overlay:
                nominal_masses.append(nominal_mz + 1)
            elif (defect_mass) <= nominal_overlay:
                nominal_masses.append(nominal_mz - 1)
            
            list_formulas_candidates = []
            
            for nominal_mass in nominal_masses:
                if nominal_mass in query.keys():   
                    list_formulas_candidates.extend(query.get(nominal_mass))

            return list_formulas_candidates
            
        all_assigned_indexes = list()

        molecular_search_settings = self.mass_spectrum_obj.molecular_search_settings

        for ms_peak in self.mass_spectrum_obj.sort_by_abundance():

            #already assigned a molecular formula
            if self.first_hit: 
                
                if ms_peak.is_assigned: continue
        
            ms_peak_indexes = SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).find_formulas(get_formulas(), min_abundance, self.mass_spectrum_obj, ms_peak, ion_type)    

            all_assigned_indexes.extend(ms_peak_indexes)
        
        #all_assigned_indexes = MolecularFormulaSearchFilters().filter_isotopologue(all_assigned_indexes, self.mass_spectrum_obj)

        #all_assigned_indexes = MolecularFormulaSearchFilters().filter_kendrick(all_assigned_indexes, self.mass_spectrum_obj)

        #MolecularFormulaSearchFilters().check_min_peaks(all_assigned_indexes, self.mass_spectrum_obj)
        #filter per min peaks per mono isotopic class
    
    def run(self, classe, possible_formulas_query, min_abundance, mass_spectrum_obj, ion_type):
       
        all_assigned_indexes = list()
        molecular_search_settings = mass_spectrum_obj.molecular_search_settings    
        for ms_peak in mass_spectrum_obj:
           
            #already assigned a molecular formula
            if self.first_hit: 
                
                if ms_peak.is_assigned: continue
        
            max_mz = ms_peak.mz_exp + 0.1
            min_mz = ms_peak.mz_exp - 0.1

            self.session.query.filter(
                MolecularFormulaLink.classe == classe, 
                MolecularFormulaLink.protonated_mass(-1) <= max_mz,
                MolecularFormulaLink.protonated_mass(-1) >= min_mz,
                MolecularFormulaLink.DBE >= molecular_search_settings.min_dbe, 
                MolecularFormulaLink.DBE <= molecular_search_settings.max_dbe, 
                MolecularFormulaLink.C >= molecular_search_settings.usedAtoms.get("C")[0],
                MolecularFormulaLink.C <= molecular_search_settings.usedAtoms.get("C")[1], 
                MolecularFormulaLink.H >= molecular_search_settings.usedAtoms.get("H")[0],
                MolecularFormulaLink.H <= molecular_search_settings.usedAtoms.get("H")[1], 
                )
            #if formulas:
                
            #    ms_peak_indexes = SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).find_formulas(formulas, min_abundance, self.mass_spectrum_obj, ms_peak, ion_type)    

                #all_assigned_indexes.extend(ms_peak_indexes)
        
        #all_assigned_indexes = MolecularFormulaSearchFilters().filter_isotopologue(all_assigned_indexes, self.mass_spectrum_obj)

        #all_assigned_indexes = MolecularFormulaSearchFilters().filter_kendrick(all_assigned_indexes, self.mass_spectrum_obj)

        #MolecularFormulaSearchFilters().check_min_peaks(all_assigned_indexes, self.mass_spectrum_obj)

                            
    def run_worker_ms_peaks(self, ms_peaks):

        SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).reset_error(self.mass_spectrum_obj)

        min_abundance = self.mass_spectrum_obj.min_abundance

        classes = NewMolecularCombinations(self.sql_db).runworker(self.mass_spectrum_obj.molecular_search_settings)
        
        classes_str = [class_tuple[0] for class_tuple in classes]

        #query database
        pbar = tqdm.tqdm(classes)
        
        dict_res = self.sql_db.get_dict_entries(classes, self.mass_spectrum_obj.molecular_search_settings)
         
        for classe_tuple in pbar:

            #add filter here and get indexes of class assigned 
            classe_str  = classe_tuple[0]
            classe_dict = classe_tuple[1]
            #query for the classes molecular formulas
            if self.mass_spectrum_obj.molecular_search_settings.isProtonated:
                    
                    pbar.set_description_str(desc="Started molecular formula search for class %s, (de)protonated " % classe_str, refresh=True)

                    ion_type = Labels.protonated_de_ion

                    #classe_query = self.sql_db.session.query(MolecularFormulaLink).filter(MolecularFormulaLink.classe == classe_str)
        
                    self.run(classe_query, min_abundance, ms_peaks, ion_type)    

            if self.mass_spectrum_obj.molecular_search_settings.isRadical:
                    
                    pbar.set_description_str(desc="Started molecular formula search for class %s, radical " % classe_str, refresh=True)
                    
                    ion_type = Labels.radical_ion
                    
                    classe_query = self.sql_db.get_by_classe(classe_str, self.mass_spectrum_obj.molecular_search_settings)
                    
                    self.run(classe_query, min_abundance, ms_peaks, ion_type)    

            # looks for adduct, used_atom_valences should be 0 
            # this code does not support H exchance by halogen atoms
            if self.mass_spectrum_obj.molecular_search_settings.isAdduct:
                
                pbar.set_description_str(desc="Started molecular formula search for class %s, adduct " % classe_str, refresh=True)
                
                ion_type = Labels.radical_ion
                
                classe_query = self.sql_db.get_by_classe(classe_str, self.mass_spectrum_obj.molecular_search_settings)
                
                self.run(classe_query, min_abundance, ms_peaks, ion_type)
        
        self.sql_db.close()                
    
    @staticmethod
    def database_to_dict(classe_str_list, nominal_mzs, mf_search_settings):
        
        sql_db = MolForm_SQL(url=mf_search_settings.url_database)
        
        dict_res = {}

        if mf_search_settings.isProtonated:
            dict_res[Labels.protonated_de_ion] = sql_db.get_dict_by_classes(classe_str_list, Labels.protonated_de_ion, nominal_mzs, mf_search_settings)    
            
        elif mf_search_settings.isRadical:
            dict_res[Labels.isRadical] = sql_db.get_dict_by_classes(classe_str_list, Labels.isRadical, nominal_mzs, mf_search_settings)    

        elif mf_search_settings.isAdduct:
            
            adduct_list = mf_search_settings.adduct_atoms_neg if mf_search_settings.molecular_search_settings.polarity < 0 else mf_search_settings.adduct_atoms_pos
            dict_res[Labels.isAdduct] = sql_db.get_dict_by_classes(classe_str_list, Labels.isAdduct, nominal_mzs, mf_search_settings, adduct_list)    

        return dict_res

    @staticmethod
    def get_database_into_dict(args):
        
        classe_str, nominal_mzs, molecular_search_settings = args

        sql_db = MolForm_SQL(url=molecular_search_settings.url_database)
        
        dict_res = {}

        if molecular_search_settings.isProtonated:
            dict_res[classe_str] = sql_db.get_dict_by_class(classe_str, Labels.protonated_de_ion, nominal_mzs, molecular_search_settings)    
            
        elif molecular_search_settings.isRadical:
            dict_res[classe_str] = sql_db.get_dict_by_class(classe_str, Labels.isRadical, nominal_mzs, molecular_search_settings)    

        return dict_res, classe_str

    @timeit                    
    def run_worker_mass_spectrum(self):

        #number_of_process = multiprocessing.cpu_count()

        '''loading this on a shared memory would be better than having to serialize it for every process
            waiting for python 3.8 release'''
       
        SearchMolecularFormulaWorker(find_isotopologues=self.find_isotopologues).reset_error(self.mass_spectrum_obj)

        min_abundance = self.mass_spectrum_obj.min_abundance

        classes = NewMolecularCombinations(self.sql_db).runworker(self.mass_spectrum_obj.molecular_search_settings)
        
        nominal_mzs = self.mass_spectrum_obj.nominal_mz

        #query database
        
        ion_type = Labels.protonated_de_ion

        #one GB for 500 classes   
        for classe_chunk in chunks(classes, 300): 
            
            classes_str_list = [class_tuple[0] for class_tuple in classe_chunk]
            
            worker_args = [(classe_str, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings) for classe_str in classes_str_list] 
            
            #dict_res = {}
            #p = multiprocessing.Pool(4)
            #for dict_class, classe_str in p.imap_unordered(self.get_database_into_dict, worker_args):
            #    dict_res[classe_str] = dict_class[classe_str]
                
            #p.close()
            #p.join()
                
            dict_res = self.database_to_dict(classes_str_list, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings)
            
            pbar = tqdm.tqdm(classe_chunk)
            
            for classe_tuple in pbar:
                
                #add filter here and get indexes of class assigned 
                classe_str  = classe_tuple[0]
                classe_dict = classe_tuple[1]
                #query for the classes molecular formulas
                
                if self.mass_spectrum_obj.molecular_search_settings.isProtonated:
                        
                        ion_type = Labels.protonated_de_ion

                        pbar.set_description_str(desc="Started molecular formula search for class %s, (de)protonated " % classe_str, refresh=True)
                        
                        self.run_search(dict_res.get(ion_type).get(classe_str), min_abundance, ion_type)

                if self.mass_spectrum_obj.molecular_search_settings.isRadical:
                        
                        pbar.set_description_str(desc="Started molecular formula search for class %s, radical " % classe_str, refresh=True)
                        
                        ion_type = Labels.radical_ion
                        
                        self.run_search(dict_res.get(ion_type).get(classe_str), min_abundance, ion_type) 

                # looks for adduct, used_atom_valences should be 0 
                # this code does not support H exchance by halogen atoms
                if self.mass_spectrum_obj.molecular_search_settings.isAdduct:
                    
                    pbar.set_description_str(desc="Started molecular formula search for class %s, adduct " % classe_str, refresh=True)
                    
                    ion_type = Labels.radical_ion
                    
                    #classe_query = self.sql_db.get_by_classe(classe_str, self.mass_spectrum_obj.molecular_search_settings)
                    #dict_res = self.sql_db.get_dict_entries([classe_str], ion_type, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings)

                    #replace ion_type in the molecular_formula object
                    self.run_search(classe_str, min_abundance, ion_type)  

        self.sql_db.close()     

    def add_adducts(self, possible_formulas):
        
        ion_type = Labels.adduct_ion

        if self.mass_spectrum_obj.polarity < 0:
            adduct_atoms = self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_neg
            

        else:
            adduct_atoms = self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_pos
          
    def search_mol_formulas(self,  possible_formulas_list, find_isotopologues=True):

        SearchMolecularFormulaWorker(find_isotopologues=find_isotopologues).reset_error(self.mass_spectrum_obj)

        initial_min_peak_bool = self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter
        initial_runtime_kendrick_filter = self.mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter
        
        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
        self.mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter = False

        possible_formulas_dict_nm =  {}
        
        for mf in possible_formulas_list:
            
            nm = mf.nominal_mz
            
            if nm in possible_formulas_dict_nm.keys():
                
                possible_formulas_dict_nm[nm].append(mf)
            
            else:    
                possible_formulas_dict_nm[nm] = [mf]

        min_abundance = self.mass_spectrum_obj.min_abundance

        self.run_search(possible_formulas_dict_nm, min_abundance)          

        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = initial_min_peak_bool
        self.mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter = initial_runtime_kendrick_filter

        mspeaks = [mspeak for mspeak in self.mass_spectrum_obj if mspeak.is_assigned]
        
        self.sql_db.close()
        
        return mspeaks
    
            
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
                mass_spectrum_obj.molecular_search_settings.min_ppm_error  = closest_error - mass_spectrum_obj.molecular_search_settings.mz_error_range
                mass_spectrum_obj.molecular_search_settings.max_ppm_error = closest_error + mass_spectrum_obj.molecular_search_settings.mz_error_range

        elif mass_spectrum_obj.molecular_search_settings.error_method == 'lowest':
            
            if error < last_error:
                mass_spectrum_obj.molecular_search_settings.min_ppm_error  = error - mass_spectrum_obj.molecular_search_settings.mz_error_range
                mass_spectrum_obj.molecular_search_settings.max_ppm_error = error + mass_spectrum_obj.molecular_search_settings.mz_error_range
                last_error = error
                
        
        elif mass_spectrum_obj.molecular_search_settings.error_method == 'symmetrical':
               
               mass_spectrum_obj.molecular_search_settings.min_ppm_error  = mass_spectrum_obj.molecular_search_settings.mz_error_average - mass_spectrum_obj.molecular_search_settings.mz_error_range
               mass_spectrum_obj.molecular_search_settings.max_ppm_error = mass_spectrum_obj.molecular_search_settings.mz_error_average + mass_spectrum_obj.molecular_search_settings.mz_error_range
        
        elif mass_spectrum_obj.molecular_search_settings.error_method == 'average':

                nbValues += 1
                error_average = error_average + ((error - error_average) / nbValues)
                mass_spectrum_obj.molecular_search_settings.min_ppm_error  =  error_average - mass_spectrum_obj.molecular_search_settings.mz_error_range
                mass_spectrum_obj.molecular_search_settings.max_ppm_error =  error_average + mass_spectrum_obj.molecular_search_settings.mz_error_range    
                
                
        else:
            #using set mass_spectrum_obj.molecular_search_settings.min_ppm_error  and max_ppm_error range
            pass

        '''returns the error based on the selected method at mass_spectrum_obj.molecular_search_settings.method
        '''    
        
    @staticmethod
    def calc_error(mz_exp, mz_calc, method='ppm'):
        
        '''method should be ppm or ppb'''
         
        if method == 'ppm':
            multi_factor = 1000000
        
        elif method == 'ppb':
            multi_factor = 1000000
        
        elif method == 'perc':
            multi_factor = 100

        else:
            raise Exception("method needs to be ppm or ppb, you have entered %s" % method)

        if mz_exp:
            
            return ((mz_calc - mz_exp)/mz_calc)*multi_factor        
        
        else:
            
            raise Exception("Please set mz_calc first")    

    def find_formulas(self, formulas, min_abundance, 
                      mass_spectrum_obj, ms_peak, ion_type):
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

        min_ppm_error  = mass_spectrum_obj.molecular_search_settings.min_ppm_error 
        max_ppm_error = mass_spectrum_obj.molecular_search_settings.max_ppm_error
        
        min_abun_error = mass_spectrum_obj.molecular_search_settings.min_abun_error
        max_abun_error = mass_spectrum_obj.molecular_search_settings.max_abun_error
        
        #f = open("abundance_error.txt", "a+")    
        ms_peak_mz_exp, ms_peak_abundance = ms_peak.mz_exp, ms_peak.abundance
        #min_error = min([pmf.mz_error for pmf in possible_formulas])
        
        for possible_formula in formulas:
            
            if possible_formula:
                
                error = self.calc_error(ms_peak_mz_exp, possible_formula.protonated_mass(-1))
                
                #error = possible_formula.mz_error
               
                if  min_ppm_error  <= error <= max_ppm_error:
                    
                    #update the error
                    
                    self.set_last_error(error, mass_spectrum_obj)    
                   
                    #add molecular formula match to ms_peak
                    
                    # get molecular formula dict from sql obj
                    #formula_dict = pickle.loads(possible_formula.mol_formula)
                    formula_dict = possible_formula.formula_dict
                    # create the molecular formula obj to be stored
                    molecular_formula = MolecularFormula(formula_dict, ms_peak.ion_charge, ion_type=ion_type)

                    # add the molecular formula obj to the mspeak obj
                    # add the mspeak obj and it's index for tracking next assignment step
                    
                    if self.find_isotopologues:
                        
                        # calculates isotopologues
                        isotopologues = molecular_formula.isotopologues(min_abundance, ms_peak_abundance, mass_spectrum_obj.dynamic_range)
                        
                        # search for isotopologues
                        for isotopologue_formula in isotopologues:
                            
                            isotopologue_formula = isotopologue_formula

                            molecular_formula.expected_isotopologues.append(isotopologue_formula)
                            #move this outside to improve preformace
                            #we need to increase the search space to -+1 m_z 
                            first_index, last_index = mass_spectrum_obj.get_nominal_mz_first_last_indexes(isotopologue_formula.mz_nominal_calc)
                            
                            for ms_peak_iso in mass_spectrum_obj[first_index:last_index]:
                                
                                error = self.calc_error(ms_peak_iso.mz_exp, isotopologue_formula.mz_calc)
                                
                                if  min_ppm_error  <= error <= max_ppm_error:
                                    
                                    #need to define error distribution for abundance measurements
                                    
                                    if mass_spectrum_obj.is_centroid:
                                        
                                        abundance_error = self.calc_error(isotopologue_formula.abundance_calc, ms_peak_iso.abundance,method='perc')            
                                        
                                    else:
                                        
                                        abundance_error = self.calc_error(ms_peak.area, ms_peak_iso.area, method='perc')
                                        
                                    #area_error = self.calc_error(ms_peak.area, ms_peak_iso.area, method='perc')

                                    # margin of error was set empirically/ needs statistical calculation
                                    #  of margin of error for the measurement of the abundances
                                    if min_abun_error <= abundance_error <= max_abun_error:
                                        
                                        #update the error   
                                        
                                        self.set_last_error(error, mass_spectrum_obj)    
                                        
                                        #isotopologue_formula.mz_error = error

                                        #isotopologue_formula.area_error = area_error

                                        #isotopologue_formula.abundance_error = abundance_error

                                        isotopologue_formula.mspeak_index_mono_isotopic = ms_peak.index
                                        
                                        #add mspeaks isotopologue index to the mono isotopic MolecularFormula obj and the respective formula position  
                                        
                                        #add molecular formula match to ms_peak
                                        x = ms_peak_iso.add_molecular_formula(isotopologue_formula)
                                        
                                        molecular_formula.mspeak_mf_isotopologues_indexes.append((ms_peak_iso.index, x))
                                        #add mspeaks mono isotopic index to the isotopologue MolecularFormula obj
                                        

                    y = ms_peak.add_molecular_formula(molecular_formula)            

                    mspeak_assigned_index.append((ms_peak.index, y))


        return mspeak_assigned_index

