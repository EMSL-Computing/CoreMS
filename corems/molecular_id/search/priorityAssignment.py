
import os,  sys
sys.path.append('.')
from copy import deepcopy
from threading import Thread
from itertools import product

from corems.encapsulation.constant import Labels, Atoms
from corems.molecular_id.calc.MolecularFilter import MolecularFormulaSearchFilters
from corems.molecular_id.factory.MolecularLookupTable import MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL as molform_db
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulaWorker
from corems.molecular_id.factory.molecularSQL import MolForm_SQL 
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter


class OxygenPriorityAssignment(Thread):

    def __init__(self, mass_spectrum_obj):
        '''TODO:- add support for other atoms and adducts: Done
                - add dbe range on search runtime : Done
                - add docs
                - improve performace : Done 
        '''
        Thread.__init__(self)
        self.mass_spectrum_obj = mass_spectrum_obj
        #  initiated at create_molecular_database()
        #self.dict_molecular_lookup_table = None
        
        self.sql_db = MolForm_SQL(mass_spectrum_obj.polarity)

    def run(self):
        
        # get Oxygen classes dict and the associate mspeak class 
        # list_of_classes_min_max_dbe = self.class_and_dbes_in_order()
        # create database separated to give the user the chance to use mass spec filters
             
        assign_classes_order_str_dict_tuple_list = self.create_data_base()
        
        if assign_classes_order_str_dict_tuple_list:
            
            self.run_worker_mass_spectrum(assign_classes_order_str_dict_tuple_list)
        
        else:
            
            raise RuntimeError('call create_data_base() first')
    
    def create_data_base(self):
        
        def create_molecular_database():
            
            # checks and creates the database entries for the specified heteroatomic classes 
            
            min_o = min(self.mass_spectrum_obj, key=lambda msp: msp[0]['O'])[0]['O'] - 2
            
            if min_o <= 0:
                min_o = 1

            max_o = max(self.mass_spectrum_obj, key=lambda msp: msp[0]['O'])[0]['O'] + 2

            #min_dbe = min(self.mass_spectrum_obj, key=lambda msp: msp[0].dbe)[0].dbe

            #max_dbe = max(self.mass_spectrum_obj, key=lambda msp: msp[0].dbe)[0].dbe

            #self.lookupTableSettings.use_pah_line_rule = False
            
            #self.lookupTableSettings.min_dbe = min_dbe/2#min_dbe - 7 if  (min_dbe - 7) > 0 else 0
            
            #self.lookupTableSettings.max_dbe = max_dbe * 2 #max_dbe + 7
            
            self.mass_spectrum_obj.reset_indexes()

            #initial_ox = deepcopy(self.mass_spectrum_obj.molecular_search_settings.usedAtoms)

            self.mass_spectrum_obj.molecular_search_settings.usedAtoms['O'] = (min_o, max_o)

            classes = MolecularCombinations(self.sql_db).runworker(self.mass_spectrum_obj.molecular_search_settings)
            
            classes_str = [class_tuple[0] for class_tuple in classes]

            nominal_mzs = self.mass_spectrum_obj.nominal_mz

            self.dict_molecular_lookup_table = self.get_dict_molecular_database(classes_str, nominal_mzs)
        
        # get the most abundant peak and them every 14Da, only allow Ox and its derivatives
       
        find_formula_thread = FindOxygenPeaks(self.mass_spectrum_obj, self.sql_db)
        find_formula_thread.run()
        
        #mass spec obj indexes are set to interate over only the peaks with a molecular formula candidate
        find_formula_thread.set_mass_spec_indexes_by_found_peaks()
        
        #get the Ox class and the DBE for the lowest error molecular formula candidate
        dict_ox_class_and_ms_peak = self.ox_classes_and_peaks_in_order_()
                      
        # sort the classes by abundance
        
        assign_classes_order_str_dict_tuple_list = self.get_classes_in_order(dict_ox_class_and_ms_peak)
        
        create_molecular_database()
                
        return assign_classes_order_str_dict_tuple_list
        
    def run_worker_mass_spectrum(self, assign_classes_order_tuples):
        
        def check_adduct_class(classe_dict):

            return any([key in classe_dict.keys() for key in self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_neg])
        
        def set_min_max_dbe_by_oxygen(classe_dict):
            # calculates min and max DBE based on the Oxygen number
            # ref :https://pubs.acs.org/doi/full/10.1021/ac200464q
            # if class does not has O it use the pha rule
            # ref : Vlad Lobodin manuscript to be include here
            '''
            atoms_exchanges = ['N']
            if 'O' in classe_dict.keys():
                
                Oxygen_number = classe_dict.get("O")
                for atom in atoms_exchanges:
                    if atom in classe_dict.keys():
                        Oxygen_number += classe_dict.get(atom)

                self.mass_spectrum_obj.molecular_search_settings.min_dbe = (Oxygen_number/3) - 0.5 
                self.mass_spectrum_obj.molecular_search_settings.max_dbe = Oxygen_number*3 + 0.5 + 2
            
            else:
            '''    
            self.mass_spectrum_obj.molecular_search_settings.use_pah_line_rule = True

        def run_search(possible_formulas_dict, mass_spectrum_obj, min_abundance, is_adduct=False):
            
            all_assigned_indexes = list()
            
            for ms_peak in mass_spectrum_obj.sort_by_abundance():

                if ms_peak: continue
                #already assigned a molecular formula
               
                nominal_mz  = ms_peak.nominal_mz_exp

                #get mono isotopic peaks that was added a molecular formula obj
                #TODO update error variables

                possible_formulas_nominal = possible_formulas_dict.get(nominal_mz)
                
                if possible_formulas_nominal:

                    if is_adduct:

                            for m_formula in possible_formulas_nominal: m_formula.ion_type = Labels.adduct_ion
                        
                    ms_peak_indexes = SearchMolecularFormulaWorker().find_formulas(possible_formulas_nominal, min_abundance, mass_spectrum_obj, ms_peak)    

                    all_assigned_indexes.extend(ms_peak_indexes)
            
            
            #filter peaks by percentile threshold of found isotopologues 
            all_assigned_indexes = MolecularFormulaSearchFilters().filter_isotopologue(all_assigned_indexes, mass_spectrum_obj)

            #filter noise by kendrick density
            all_assigned_indexes = MolecularFormulaSearchFilters().filter_kendrick(all_assigned_indexes, mass_spectrum_obj)

            #filter per min peaks per mono isotopic class
            # this function should always be the last filter, 
            # thefore no need to return remaining indexes
            MolecularFormulaSearchFilters().check_min_peaks(all_assigned_indexes, mass_spectrum_obj)

        #error_average = self.mass_spectrum_obj.molecular_search_settings.mz_error_average
        
        kdm_base = self.mass_spectrum_obj.mspeaks_settings.kendrick_base
        
        self.mass_spectrum_obj.change_kendrick_base_all_mspeaks(kdm_base)

        ClusteringFilter().filter_kendrick(self.mass_spectrum_obj)

        min_abundance = self.mass_spectrum_obj.min_abundance

        for classe_tuple in assign_classes_order_tuples:

            classe_str  = classe_tuple[0]
            classe_dict = classe_tuple[1]
            
            is_adduct = check_adduct_class(classe_dict)
            set_min_max_dbe_by_oxygen(classe_dict)
            
            #if len(classe_dict.keys()) == 2:
            #    if classe_dict.get('S') == 1:
            #       continue
            # limits the dbe by the Ox class most abundant,
            # need to add other atoms contribution to be more accurate
            # but +-7 should be sufficient to cover the range 
            
            if self.mass_spectrum_obj.molecular_search_settings.isProtonated and not is_adduct:

                    print("Started molecular formula search for class %s, (de)protonated " % classe_str)

                    ion_type = Labels.protonated_de_ion

                    possible_formulas_dict = self.dict_molecular_lookup_table.get(ion_type).get(classe_str)
                    
                    if possible_formulas_dict:

                        run_search(possible_formulas_dict, self.mass_spectrum_obj, min_abundance, is_adduct=is_adduct)

            if self.mass_spectrum_obj.molecular_search_settings.isRadical and not is_adduct:

                    print("Started molecular formula search for class %s,  radical" % classe_str)

                    ion_type = Labels.radical_ion
                    
                    possible_formulas_dict = self.dict_molecular_lookup_table.get(ion_type).get(classe_str)
                    
                    if possible_formulas_dict:

                        run_search(possible_formulas_dict, self.mass_spectrum_obj, min_abundance, is_adduct=is_adduct)

            # looks for adduct, used_atom_valences should be 0 
            # this code does not support H exchance by halogen atoms
            if self.mass_spectrum_obj.molecular_search_settings.isAdduct and is_adduct:
                
                print("Started molecular formula search for class %s, adduct" % classe_str)
                
                ion_type = Labels.radical_ion
                
                possible_formulas_dict = self.dict_molecular_lookup_table.get(ion_type).get(classe_str)
                    
                if possible_formulas_dict:

                    run_search(possible_formulas_dict, self.mass_spectrum_obj, min_abundance, is_adduct=is_adduct)
        
        print("Finished molecular formula search")
    
    def get_dict_molecular_database(self, classes_str, nominal_mzs):
            
        dict_res = {}
        
        if self.mass_spectrum_obj.molecular_search_settings.isProtonated:
            
            ion_type = Labels.protonated_de_ion
            
            dict_res[ion_type] = self.sql_db.get_dict_entries(classes_str, ion_type, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings)
            
        if self.mass_spectrum_obj.molecular_search_settings.isRadical or self.mass_spectrum_obj.molecular_search_settings.isAdduct:

            ion_type = Labels.radical_ion

            dict_res[ion_type] = self.sql_db.get_dict_entries(classes_str, ion_type, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings)
    
        return dict_res
    
    def ox_classes_and_peaks_in_order_(self) -> dict:
        # order is only valid in python 3.4 and above
        # change to OrderedDict if your version is lower
        dict_ox_class_and_ms_peak = dict()
        
        for mspeak in self.mass_spectrum_obj.sort_by_abundance(reverse=True):
            
            #change this filter to cia filter, give more option here, confidence, number of isotopologue found etc

            ox_classe = mspeak.best_molecular_formula_candidate.class_label
            
            if ox_classe in dict_ox_class_and_ms_peak.keys():
                
                #get the most abundant of the same ox class
                if mspeak.abundance > dict_ox_class_and_ms_peak[ox_classe].abundance:

                    dict_ox_class_and_ms_peak[ox_classe] = (mspeak)
            else:
                    
                dict_ox_class_and_ms_peak[ox_classe] = (mspeak)
        
        return dict_ox_class_and_ms_peak

    def get_classes_in_order(self, dict_ox_class_and_ms_peak)-> [(str, dict)]: 
        ''' structure is 
            ('HC', {'HC': 1})'''
        
        usedAtoms = deepcopy(self.mass_spectrum_obj.molecular_search_settings.usedAtoms)
        
        usedAtoms.pop("C")
        usedAtoms.pop("H")
        usedAtoms.pop("O")

        min_n, max_n = usedAtoms.get('N')
        min_s, max_s = usedAtoms.get('S')
        min_p, max_p = usedAtoms.get('P')

        possible_n = [n for n in range(min_n, max_n + 1)]
        possible_s = [s for s in range(min_s, max_s + 1)]
        possible_p = [p for p in range(min_p, max_p + 1)]
        
        #used to enforce order for commum atoms 
        # and track the atom index in on the tuple in all_atoms_tuples
        atoms_in_order = ['N', 'S', 'P']
        
        #do number atoms prodcut and remove then from the usedAtoms dict
        all_atoms_tuples = product(possible_n, possible_s, possible_p)
        for atom in atoms_in_order:
            
            usedAtoms.pop(atom, None)
        
        #iterate over other atoms besides C,H, N, O, S and P
        
        for selected_atom_label, min_max_tuple in usedAtoms.items():
            
            min_x = min_max_tuple[0]
            max_x = min_max_tuple[1]

            possible_x = [x for x in range(min_x, max_x + 1)]
            all_atoms_tuples = product(all_atoms_tuples, possible_x)
            
            #merge tuples
            all_atoms_tuples = [all_atoms_combined[0] + (all_atoms_combined[1],) for all_atoms_combined in
                                all_atoms_tuples]
            
            #add atom label to the atoms_in_order list
            
            #important to index where the atom position is in on the tuple in all_atoms_tuples
            atoms_in_order.append(selected_atom_label)

        classes_strings_dict_tuples, hc_class = self.get_class_strings_dict(all_atoms_tuples, atoms_in_order)

        combined_classes = self.combine_ox_class_with_other(atoms_in_order, classes_strings_dict_tuples, dict_ox_class_and_ms_peak)
        
        combination_classes_ordered = self.sort_classes(atoms_in_order, combined_classes)
        
        oxygen_class_str_dict_tuple = [(ox_class, mspeak[0].class_dict) for ox_class, mspeak in dict_ox_class_and_ms_peak.items()] 

        ## add classes together and ignores classes selected from the main series
        for class_tuple in  combination_classes_ordered:
            if class_tuple not in oxygen_class_str_dict_tuple:
                oxygen_class_str_dict_tuple.append(class_tuple)
        
        return oxygen_class_str_dict_tuple

    @staticmethod
    def get_class_strings_dict(all_atoms_tuples, atoms_in_order) -> [(str, dict)]: 
        
        classe_list= []
        hc_class = []
        
        for all_atoms_tuple in all_atoms_tuples:
            
            classe_str = ''
            classe_dict = dict()
            
            for each_atoms_index, atoms_number in enumerate(all_atoms_tuple):
                
                if atoms_number != 0:
                    
                    classe_str = (classe_str + atoms_in_order[each_atoms_index] + str(atoms_number) + ' ')
                    
                    classe_dict[atoms_in_order[each_atoms_index]] = atoms_number

            classe_str = classe_str.strip()
            
            if len(classe_str) > 0:
            
                classe_list.append((classe_str,classe_dict))

            elif len(classe_str) == 0:

                hc_class.append(('HC', {'HC':1}))
        
        return classe_list, hc_class
    
    @staticmethod
    def combine_ox_class_with_other( atoms_in_order, classes_strings_dict_tuples, dict_ox_class_and_ms_peak) -> [dict]:
        
        #sort methods that uses the key of classes dictionary and the atoms_in_order as reference
        # c_tuple[1] = class_dict, because is one key:value map we loop through keys and get the first item only 
        # sort by len first then sort based on the atoms_in_order list
        atoms_in_order = Atoms.atoms_order

        Oxygen_mfs = dict_ox_class_and_ms_peak.values()
        
        
        #sort_method = lambda word: (len(word[0]), [atoms_in_order.index(atom) for atom in list( word[1].keys())])
        
        #print(classes_strings_dict_tuples)
        #classe_in_order = sorted(classes_strings_dict_tuples, key = sort_method)
        #print(classe_in_order)
        
        combination = []
        
        # _ ignoring the class_str
        for _ , other_classe_dict in classes_strings_dict_tuples:
          
           #combination.extend([[other_classe_str + ' ' + Oxygen_mf[0].class_label , {**other_classe_dict, **Oxygen_mf[0].class_dict}] for Oxygen_mf in Oxygen_mfs])
           combination.extend([{**other_classe_dict, **Oxygen_mf[0].class_dict} for Oxygen_mf in Oxygen_mfs])
 
        return combination
    
    @staticmethod
    def sort_classes( atoms_in_order, combination_tuples) -> [(str, dict)]: 
        
        join_list_of_list_classes = list()
        atoms_in_order =  ['N','S','P','O'] + atoms_in_order[3:]
        
        sort_method = lambda atoms_keys: [atoms_in_order.index(atoms_keys)] #(len(word[0]), print(word[1]))#[atoms_in_order.index(atom) for atom in list( word[1].keys())])
        for class_dict in combination_tuples:
            
            sorted_dict_keys = sorted(class_dict, key = sort_method)
            class_str = ' '.join([atom + str(class_dict[atom]) for atom in sorted_dict_keys])
            new_class_dict = { atom: class_dict[atom] for atom in sorted_dict_keys}
            join_list_of_list_classes.append((class_str, new_class_dict))
        
        return join_list_of_list_classes
 