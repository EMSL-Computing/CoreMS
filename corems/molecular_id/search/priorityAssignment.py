
import os,  sys
sys.path.append('.')
from copy import deepcopy
from threading import Thread
from itertools import product

from corems.encapsulation.constant import Labels
from corems.molecular_id.search.findOxigenPeaks import FindOxygenPeaks
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulaWorker, SearchMolecularFormulas
from corems.molecular_id.factory.MolecularLookupTable import MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL as molform_db
#from corems.molecular_id.factory.molecularMongo import MolForm_Mongo as molform_db

class OxigenPriorityAssignment(Thread):

    def __init__(self, mass_spectrum_obj):
        '''TODO:- add support for other atoms and adducts: Done
                - add dbe range on search runtime : Done
                - add docs
                - improve performace : Done 
        '''
        Thread.__init__(self)
        self.mass_spectrum_obj = mass_spectrum_obj
        #  initiated at create_molecular_database()
        self.dict_molecular_lookup_table = None
        self.assign_classes_order_str_dict_tuple_list = None
        self.dict_ox_class_and_ms_peak = None

    def run(self):
        
        # get oxigen classes dict and the associate mspeak class 
        # list_of_classes_min_max_dbe = self.class_and_dbes_in_ordem()
        # create database separated to give the user the chance to use mass spec filters
        
        self.create_data_base()

        if self.assign_classes_order_str_dict_tuple_list and self.dict_ox_class_and_ms_peak:
            
            self.run_worker_mass_spectrum(self.assign_classes_order_str_dict_tuple_list, self.dict_ox_class_and_ms_peak)
        
        else:
            raise RuntimeError('call create_data_base() first')
    
    def create_data_base(self):
        
        find_formula_thread = FindOxygenPeaks(self.mass_spectrum_obj)
        find_formula_thread.run()
        #mass spec obj indexes are set to interate over only the peaks with a molecular formula candidate
        find_formula_thread.set_mass_spec_indexes_by_found_peaks()
        
        self.dict_ox_class_and_ms_peak = self.ox_classes_and_dbes_in_ordem_()
        
        self.assign_classes_order_str_dict_tuple_list = self.get_classes_in_order(self.dict_ox_class_and_ms_peak)

        # reset back the massspec obj indexes since we already collected the mspeak indexes
        
        self.create_molecular_database()
        
    def run_worker_mass_spectrum(self, assign_classes_order_tuples, dict_ox_class_and_ms_peak):
        
        def check_adduct_class(classe_dict):
            return any([key in classe_dict.keys() for key in self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_neg])
        
        def check_min_peaks(ms_peak_indexes):
            
            if  self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter:

                if not len(ms_peak_indexes) >= self.mass_spectrum_obj.molecular_search_settings.min_peaks_per_class:
                    
                    for index in ms_peak_indexes: self.mass_spectrum_obj[index].clear_molecular_formulas()

        def set_min_max_dbe_by_oxygen(classe_dict):
            # calculates min and max DBE based on the oxigen number
            # ref :https://pubs.acs.org/doi/full/10.1021/ac200464q
            # if class does not has O it use the pha rule
            # ref : Vlad Lobodin manuscript to be include here
            '''
            atoms_exchanges = ['N']
            if 'O' in classe_dict.keys():
                
                oxigen_number = classe_dict.get("O")
                for atom in atoms_exchanges:
                    if atom in classe_dict.keys():
                        oxigen_number += classe_dict.get(atom)

                self.mass_spectrum_obj.molecular_search_settings.min_dbe = (oxigen_number/3) - 0.5 
                self.mass_spectrum_obj.molecular_search_settings.max_dbe = oxigen_number*3 + 0.5 + 2
            
            else:
            '''    
            self.mass_spectrum_obj.molecular_search_settings.use_pah_line_rule = True

        def run_search(possible_formulas_dict, mass_spectrum_obj, min_abundance, is_adduct=False):
            
            all_assigned_indexes = list()

            for ms_peak in mass_spectrum_obj.sort_by_abundance():

                #already assinged a molecular formula
               
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
            check_min_peaks(all_assigned_indexes)
        
        error_average = self.mass_spectrum_obj.molecular_search_settings.mz_error_average
        
        min_abundance = self.mass_spectrum_obj.min_abundance

        print('classes assignment order', assign_classes_order_tuples)
        
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
        
                    ion_type = Labels.protonated_de_ion

                    possible_formulas_dict = self.dict_molecular_lookup_table.get(ion_type).get(classe_str)
                    
                    if possible_formulas_dict:

                        run_search(possible_formulas_dict, self.mass_spectrum_obj, min_abundance, is_adduct=is_adduct)

            if self.mass_spectrum_obj.molecular_search_settings.isRadical and not is_adduct:
                
                    ion_type = Labels.radical_ion
                    
                    possible_formulas_dict = self.dict_molecular_lookup_table.get(ion_type).get(classe_str)
                    
                    if possible_formulas_dict:

                        run_search(possible_formulas_dict, self.mass_spectrum_obj, min_abundance, is_adduct=is_adduct)

            # looks for adduct, used_atom_valences should be 0 
            # this code does not support H exchance by halogen atoms
            if self.mass_spectrum_obj.molecular_search_settings.isAdduct and is_adduct:
                
                ion_type = Labels.radical_ion
                
                possible_formulas_dict = self.dict_molecular_lookup_table.get(ion_type).get(classe_str)
                    
                if possible_formulas_dict:

                    run_search(possible_formulas_dict, self.mass_spectrum_obj, min_abundance, is_adduct=is_adduct)

    def create_molecular_database(self):
        #number_of_process = multiprocessing.cpu_count()

        '''loading this on a shared memory would be better than having to serialize it for every process
            waiting for python 3.8 release'''
        
        min_o = min(self.mass_spectrum_obj, key=lambda msp: msp[0]['O'])[0]['O'] - 2
        
        max_o = max(self.mass_spectrum_obj, key=lambda msp: msp[0]['O'])[0]['O'] + 2

        #min_dbe = min(self.mass_spectrum_obj, key=lambda msp: msp[0].dbe)[0].dbe

        #max_dbe = max(self.mass_spectrum_obj, key=lambda msp: msp[0].dbe)[0].dbe

        #self.lookupTableSettings.use_pah_line_rule = False
        
        #self.lookupTableSettings.min_dbe = min_dbe/2#min_dbe - 7 if  (min_dbe - 7) > 0 else 0
        
        #self.lookupTableSettings.max_dbe = max_dbe * 2 #max_dbe + 7
        
        self.mass_spectrum_obj.reset_indexes()

        #initial_ox = deepcopy(self.mass_spectrum_obj.molecular_search_settings.usedAtoms)

        self.mass_spectrum_obj.molecular_search_settings.usedAtoms['O'] = (min_o, max_o)

        classes = MolecularCombinations().runworker(self.mass_spectrum_obj.molecular_search_settings)

        classes_str = [class_tuple[0] for class_tuple in classes]

        nominal_mzs = self.mass_spectrum_obj.nominal_mz_exp

        self.dict_molecular_lookup_table = self.get_dict_molecular_database(classes_str, nominal_mzs)
        
    def get_dict_molecular_database(self, classes_str, nominal_mzs):
            
        dict_res = {}
        
        #print (classes_str)
        if self.mass_spectrum_obj.molecular_search_settings.isProtonated:
            
            ion_type = Labels.protonated_de_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings)

        if self.mass_spectrum_obj.molecular_search_settings.isRadical:

            ion_type = Labels.radical_ion

            with molform_db() as sql_handle:

                dict_res[ion_type] = sql_handle.get_dict_entries(classes_str, ion_type, nominal_mzs, self.mass_spectrum_obj.molecular_search_settings)
        
        return dict_res

    def ox_classes_and_dbes_in_ordem_(self) -> dict:
        
        dict_ox_class_and_ms_peak = dict()
        
        for mspeak in sorted(self.mass_spectrum_obj, key=lambda msp: msp.abundance, reverse=True):
            
            ox_classe = mspeak.molecular_formula_lowest_error.class_label
            
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
        
        #used to inforce order for commum atoms 
        # and track the atom index in on the tuple in all_atomos_tuples
        atomos_in_ordem = ['N', 'S', 'P']
        
        #do number atoms prodcut and remove then from the usedAtoms dict
        all_atomos_tuples = product(possible_n, possible_s, possible_p)
        for atomo in atomos_in_ordem:
            
            usedAtoms.pop(atomo, None)
        
        #iterate over other atoms besides C,H, N, O, S and P
        
        for selected_atom_label, min_max_tuple in usedAtoms.items():
            
            min_x = min_max_tuple[0]
            max_x = min_max_tuple[1]
            possible_x = [x for x in range(min_x, max_x + 1)]
            all_atomos_tuples = product(all_atomos_tuples, possible_x)
            
            #merge tuples
            all_atomos_tuples = [all_atomos_combined[0] + (all_atomos_combined[1],) for all_atomos_combined in
                                all_atomos_tuples]
            
            #add atom label to the atomos_in_ordem list
            
            #important to index where the atom position is in on the tuple in all_atomos_tuples
            atomos_in_ordem.append(selected_atom_label)
       
        classes_strings_dict_tuples, hc_class = self.get_class_strings_dict(all_atomos_tuples, atomos_in_ordem)

        combined_classes = self.combine_ox_class_with_other(atomos_in_ordem, classes_strings_dict_tuples, dict_ox_class_and_ms_peak)
        
        combination_classes_ordered = self.sort_classes(atomos_in_ordem, combined_classes)
        
        oxygen_class_str_dict_tuple = [(oxclass, mspeak[0].class_dict) for oxclass, mspeak in dict_ox_class_and_ms_peak.items()] 

        return  oxygen_class_str_dict_tuple + combination_classes_ordered #+ classe_in_ordem + hc_class

    @staticmethod
    def get_class_strings_dict(all_atomos_tuples, atomos_in_ordem) -> [(str, dict)]: 
        
        classe_list= []
        hc_class = []
        
        for all_atomos_tuple in all_atomos_tuples:
            
            classe_str = ''
            classe_dict = dict()
            
            for each_atomos_index, atoms_number in enumerate(all_atomos_tuple):
                
                if atoms_number != 0:
                    
                    classe_str = (classe_str + atomos_in_ordem[each_atomos_index] + str(atoms_number) + ' ')
                    
                    classe_dict[atomos_in_ordem[each_atomos_index]] = atoms_number

            classe_str = classe_str.strip()
            
            if len(classe_str) > 0:
            
                classe_list.append((classe_str,classe_dict))

            elif len(classe_str) == 0:

                hc_class.append(('HC', {'HC':1}))
        
        return classe_list, hc_class
    
    @staticmethod
    def combine_ox_class_with_other( atomos_in_ordem, classes_strings_dict_tuples, dict_ox_class_and_ms_peak) -> [dict]:
        
        #sort methods that uses the key of classes dictonary and the atoms_in_ordem as referece
        # c_tuple[1] = class_dict, because is one key:value map we loop throught keys and get the first item only 
        # sort by len first then sort based on the atomos_in_ordem list
        oxigen_mfs = dict_ox_class_and_ms_peak.values()
        
        sort_method = lambda word: (len(word[0]), [atomos_in_ordem.index(atom) for atom in list( word[1].keys())])
        classe_in_ordem = sorted(classes_strings_dict_tuples, key = sort_method)
        combination = []
        
        # _ ignoring the class_str
        for _ , other_classe_dict in classe_in_ordem:
          
           #combination.extend([[other_classe_str + ' ' + oxigen_mf[0].class_label , {**other_classe_dict, **oxigen_mf[0].class_dict}] for oxigen_mf in oxigen_mfs])
           combination.extend([{**other_classe_dict, **oxigen_mf[0].class_dict} for oxigen_mf in oxigen_mfs])
 
        return combination
    
    @staticmethod
    def sort_classes( atomos_in_ordem, combination_tuples) -> [(str, dict)]: 
        
        join_list_of_list_classes = list()
        atomos_in_ordem =  ['N','S','P', 'O'] + atomos_in_ordem[3:]
        
        sort_method = lambda atoms_keys: [atomos_in_ordem.index(atoms_keys)] #(len(word[0]), print(word[1]))#[atomos_in_ordem.index(atom) for atom in list( word[1].keys())])
        for class_dict in combination_tuples:
            
            sorted_dict_keys = sorted(class_dict, key = sort_method)
            class_str = ' '.join([atom + str(class_dict[atom]) for atom in sorted_dict_keys])
            new_class_dict = { atom: class_dict[atom] for atom in sorted_dict_keys}
            join_list_of_list_classes.append((class_str, new_class_dict))
        
        return join_list_of_list_classes
        
'''   
if __name__ == "__main__":
    pass
    from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
    from corems.encapsulation.settings.molecular_id.MolecularIDSettings import self.mass_spectrum_obj.molecular_search_settings, MoleculaLookupTableSettings
    from corems.mass_spectrum.calc.CalibrationCalc import MZDomain_Calibration, FreqDomain_Calibration
    from matplotlib import pyplot, colors as mcolors

    def calibrate():
        
        self.mass_spectrum_obj.molecular_search_settings.error_method = 'average'
        self.mass_spectrum_obj.molecular_search_settings.min_mz_error = -2
        self.mass_spectrum_obj.molecular_search_settings.max_mz_error = 2
        self.mass_spectrum_obj.molecular_search_settings.mz_error_range = 1

        find_formula_thread = FindOxygenPeaks(mass_spectrum_obj, lookupTableSettings)
        find_formula_thread.run()
        mspeaks_results = find_formula_thread.get_list_found_peaks()
        
        calibrate = FreqDomain_Calibration(mass_spectrum_obj, mspeaks_results)
        calibrate.ledford_calibration()
        mass_spectrum_obj.clear_molecular_formulas()

    
    def assign_mf():
        
        self.mass_spectrum_obj.molecular_search_settings.error_method = 'symmetrical'
        self.mass_spectrum_obj.molecular_search_settings.min_mz_error = -1
        self.mass_spectrum_obj.molecular_search_settings.max_mz_error = 1
        self.mass_spectrum_obj.molecular_search_settings.mz_error_range = 2
        self.mass_spectrum_obj.molecular_search_settings.mz_error_average = 0
        self.mass_spectrum_obj.molecular_search_settings.min_abun_error = -30 # percentage
        self.mass_spectrum_obj.molecular_search_settings.max_abun_error = 70 # percentage
        self.mass_spectrum_obj.molecular_search_settings.isProtonated = True
        self.mass_spectrum_obj.molecular_search_settings.isRadical = False
        self.mass_spectrum_obj.molecular_search_settings.isAdduct = True

        assignOx = OxigenPriorityAssignment(mass_spectrum_obj, lookupTableSettings)
        assignOx.start()
        assignOx.join()

    def plot():
        colors = list(mcolors.XKCD_COLORS.keys())
        oxigens = range(6,21)
        for o in oxigens:
            #o_c = list()
            for mspeak in mass_spectrum_obj:
                if mspeak:
                    #molecular_formula = mspeak.molecular_formula_lowest_error
                    for molecular_formula in mspeak:
                        if molecular_formula['O'] == o:
                            if  not molecular_formula.is_isotopologue:
                                pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                                pyplot.plot(molecular_formula['C'], molecular_formula.dbe, "o",   color=colors[molecular_formula['O']])
                                pyplot.annotate(molecular_formula.class_label, (molecular_formula['C']+0.5, molecular_formula.dbe+0.5))

            pyplot.show()                
    
    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + os.path.normcase("ESI_NEG_SRFA.d/")
    
    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    lookupTableSettings = MoleculaLookupTableSettings()
    
    lookupTableSettings.usedAtoms['O'] = (1, 22)
    
    calibrate()

    assign_mf()

    plot()


'''
    