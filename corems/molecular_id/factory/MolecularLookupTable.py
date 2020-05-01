__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

from copy import deepcopy
import itertools
import multiprocessing
import pickle
import json
from corems.encapsulation.factory.processingSetting  import MolecularLookupDictSettings
from corems.encapsulation.constant import Labels
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula 
from corems.molecular_id.factory.molecularSQL import MolForm_SQL

class MolecularCombinations:
     
    '''
    runworker()
    Returns a dictionary of molecular formula obj inside a dict nominal masses and ion type
        {Labels.ion_type:{
            classes:{
                nominal_mass:{
                    [MolecularFormula obj,]
                }
            }
        }
    '''
    def __init__(self, sql_db = None):
            
        self.sql_db = sql_db    
    
    def check_database_get_class_list(self, molecular_search_settings):

        all_class_to_create = []
        
        classes_dict = self.get_classes_in_order(molecular_search_settings)
        
        class_str_set = set(classes_dict.keys())
        
        #print(classes_list)
        
        if molecular_search_settings.isProtonated:

            existing_classes = set(self.sql_db.get_all_classes(Labels.protonated_de_ion, molecular_search_settings))

            class_to_create = class_str_set - existing_classes

            for class_str in class_to_create:
                
                class_tuple =  (class_str, classes_dict.get(class_str)) 
                
                all_class_to_create.append((class_tuple, Labels.protonated_de_ion))

        if  molecular_search_settings.isRadical or molecular_search_settings.isAdduct:

            existing_classes = set(self.sql_db.get_all_classes(Labels.radical_ion, molecular_search_settings))

            class_to_create = class_str_set - existing_classes

            for class_str in class_to_create:
                
                class_tuple =  (class_str, classes_dict.get(class_str)) 
                
                all_class_to_create.append((class_tuple, Labels.radical_ion))
                
        #print (existing_classes)class_str_list

        return [(c_s, c_d) for c_s, c_d in classes_dict.items()], all_class_to_create       

    def runworker(self, molecular_search_settings) :
        
        if not self.sql_db:
            
            polarity = 1 if molecular_search_settings.ion_charge > 0 else - 1

            self.sql_db = MolForm_SQL(polarity, molecular_search_settings.url_database)
        
        from tqdm import tqdm
        
        print ("Querying database for existing classes")
        classes_list, class_to_create = self.check_database_get_class_list(molecular_search_settings)
        print ("Finished querying database for existing classes")
        print()

        if class_to_create:
            
            settings = MolecularLookupDictSettings()
            settings.usedAtoms = deepcopy(molecular_search_settings.usedAtoms)
            settings.ion_charge = molecular_search_settings.ion_charge
            settings.url_database = molecular_search_settings.url_database
            settings.db_jobs = molecular_search_settings.db_jobs
            c_h_combinations= self.get_c_h_combination(settings)
            
            number_of_process = int(settings.db_jobs)
            
            if settings.db_jobs > 1:
                
                print("Using %i logical CPUs for database entry generation"% number_of_process )
                #number_of_process = psutil.cpu_count(logical=False)
                print('creating database entry for %i classes' % len(class_to_create))
                print()

                worker_args = [(class_tuple, c_h_combinations, ion_type, settings) for class_tuple, ion_type in class_to_create]
                p = multiprocessing.Pool(number_of_process)
                for class_list in tqdm(p.imap_unordered(CombinationsWorker(), worker_args)):
                    # TODO this will slow down a bit, need to find a better way      
                    self.add_to_sql_session(class_list)
                    self.sql_db.commit()    
                    
                #p.map(, args)
                p.close()
                p.join()
            
            else:
                
                for class_tuple, ion_type in tqdm(class_to_create):

                    (class_tuple, c_h_combinations, ion_type, settings)
                
                    class_list = CombinationsWorker().get_combinations(class_tuple, c_h_combinations, ion_type, settings)
                    self.add_to_sql_session(class_list)
                    self.sql_db.commit()             
            
            
        return classes_list
   
    def add_to_sql_session(self, class_list):
        
        self.sql_db.add_all(class_list)

    def get_c_h_combination(self, settings):

        usedAtoms = settings.usedAtoms
        result = {}
        
        min_c, max_c = usedAtoms.get('C')
        min_h, max_h = usedAtoms.get('H')

        min_h_fix = self.get_fixed_initial_number_of_hydrogen(min_h, 'odd')
        possible_c = [c for c in range(min_c, max_c + 1)]

        possible_h = [h for h in range(min_h_fix, max_h + 2, 2)]

        list_products = [i for i in itertools.product(possible_c, possible_h)]

        result['odd'] = list_products
        
        min_h_fix = self.get_fixed_initial_number_of_hydrogen(min_h, 'even')
        possible_h = [h for h in range(min_h, max_h + 2, 2)]

        list_products = [i for i in itertools.product(possible_c, possible_h)]

        result['even'] = list_products

        return result

    def swap_class_order(self, class_first, class_second, new_list2):

        if class_first in new_list2:

            if class_second in new_list2:
                n_index, s_index = (
                    new_list2.index(class_first),
                    new_list2.index(class_second),
                )

                new_list2[n_index], new_list2[s_index] = (
                    new_list2[s_index],
                    new_list2[n_index],
                )

        return new_list2    
    
    
    def get_classes_in_order(self, molecular_search_settings):
        ''' structure is 
            ('HC', {'HC': 1})'''
        
        usedAtoms = deepcopy(molecular_search_settings.usedAtoms)
        
        usedAtoms.pop("C")
        usedAtoms.pop("H")

        min_n, max_n = usedAtoms.get('N')
        min_o, max_o = usedAtoms.get('O')
        min_s, max_s = usedAtoms.get('S')
        min_p, max_p = usedAtoms.get('P')

        possible_n = [n for n in range(min_n, max_n + 1)]
        possible_o = [o for o in range(min_o, max_o + 1)]
        possible_s = [s for s in range(min_s, max_s + 1)]
        possible_p = [p for p in range(min_p, max_p + 1)]
        
        atoms_in_order = ['N', 'O', 'S', 'P']

        classe_in_order = {}

        all_atoms_tuples = itertools.product(possible_n, possible_o,
                                            possible_s, possible_p)
        
        for atom in atoms_in_order:
            usedAtoms.pop(atom, None)
        
        for selected_atom, min_max_tuple in usedAtoms.items():
            
            min_x = min_max_tuple[0]
            max_x = min_max_tuple[1]
            

            possible_x = [x for x in range(min_x, max_x + 1)]

            all_atoms_tuples = itertools.product(all_atoms_tuples, possible_x)
            all_atoms_tuples = [all_atoms_combined[0] + (all_atoms_combined[1],) for all_atoms_combined in
                                all_atoms_tuples]
            atoms_in_order.append(selected_atom)
        
        for all_atoms_tuple in all_atoms_tuples:

            classe_str = ''
            classe_dict = {}
            
            for each_atoms_index, atom_number in enumerate(all_atoms_tuple):
                
                if atom_number != 0:
                    classe_str = (classe_str + atoms_in_order[each_atoms_index] + str(atom_number) +  ' ')
                    classe_dict[atoms_in_order[each_atoms_index]] = atom_number

            classe_str = classe_str.strip()
            
            if len(classe_str) > 0:
                
                classe_in_order[classe_str] =  classe_dict

            elif len(classe_str) == 0:

                classe_in_order['HC'] = {'HC': ''}
        
        classe_in_order_dict = self.sort_classes(atoms_in_order, classe_in_order)
        
        return classe_in_order_dict

    @staticmethod
    def sort_classes( atoms_in_order, combination_dict) -> [str]: 
        
        join_dict_classes = dict()
        atoms_in_order =  ['N','S','P','O'] + atoms_in_order[4:] + ['HC']
        
        sort_method = lambda atoms_keys: [atoms_in_order.index(atoms_keys)] #(len(word[0]), print(word[1]))#[atoms_in_order.index(atom) for atom in list( word[1].keys())])
        for class_str, class_dict in combination_dict.items():
            
            sorted_dict_keys = sorted(class_dict, key = sort_method)
            class_str = ' '.join([atom + str(class_dict[atom]) for atom in sorted_dict_keys])
            class_dict = { atom: class_dict[atom] for atom in sorted_dict_keys}
            join_dict_classes[class_str] =  class_dict
        
        return join_dict_classes

    @staticmethod
    def get_fixed_initial_number_of_hydrogen( min_h, odd_even):

        remaining_h = min_h % 2
        
        if odd_even == 'even':
            
            if remaining_h == 0: return remaining_h
            
            else: return remaining_h + 1    
        
        else:
            
            if remaining_h == 0: return remaining_h + 1
            
            else: return remaining_h    

class CombinationsWorker:

    # needs this warper to pass the class to multiprocessing
    
    def __call__(self, args):
    
        return self.get_combinations(*args)  # ,args[1]

    def get_combinations(self, classe_tuple,
                          c_h_combinations, ion_type, settings
                         ):

        min_dbe = settings.min_dbe

        max_dbe = settings.max_dbe

        ion_charge = settings.ion_charge
        
        min_mz = settings.min_mz
        
        max_mz = settings.max_mz
        
        isRadical = ion_type == Labels.radical_ion
        
        isProtonated = ion_type == Labels.protonated_de_ion

        #class_dict = classe_tuple[1]

        #print("isRadical", classe_tuple[0], isRadical) 
        
        #print("isProtonated", classe_tuple[0], isProtonated) 
        
        class_dict = classe_tuple[1]
        mf_data_list = []
        if isProtonated:

            ion_type = Labels.protonated_de_ion

            odd_or_even = self.get_h_odd_or_even(ion_type, class_dict, settings)

            carbon_hydrogen_combination = c_h_combinations.get(odd_or_even)

            list_mf = self.get_mol_formulas(carbon_hydrogen_combination, ion_type, classe_tuple, 
                                                    min_dbe, max_dbe,
                                                    min_mz, max_mz, ion_charge)
            
            mf_data_list.extend(list_mf)

        if isRadical:
           
            ion_type = Labels.radical_ion
            
            odd_or_even = self.get_h_odd_or_even(ion_type, class_dict, settings )

            carbon_hydrogen_combination = c_h_combinations.get(odd_or_even)

            list_mf = self.get_mol_formulas(carbon_hydrogen_combination, ion_type, classe_tuple, 
                                                    min_dbe, max_dbe,
                                                    min_mz, max_mz, ion_charge)

            mf_data_list.extend(list_mf)
        
        return  mf_data_list   
            
    @staticmethod
    def get_mol_formulas(carbon_hydrogen_combination,
                    ion_type,
                    classe_tuple,
                    min_dbe,
                    max_dbe,
                    min_mz,
                    max_mz, 
                    ion_charge,
                    ):
        
        class_dict = classe_tuple[1]
        class_str = classe_tuple[0]
        
        list_formulas = []
        for cada_possible in carbon_hydrogen_combination:
            
            c_number = cada_possible[0]
            h_number = cada_possible[1]
            
            formula_dict = {}
            for each_atom in class_dict.keys() :
                if each_atom != 'HC':
                    formula_dict[each_atom] = class_dict.get(each_atom)

            formula_dict['C'] = c_number
            formula_dict['H'] = h_number
            formula_dict[Labels.ion_type] = ion_type

            molecular_formula = MolecularFormula(formula_dict, ion_charge)
            
            mz = molecular_formula._calc_mz()
            dbe = molecular_formula._calc_dbe()
            
            formula_dict = molecular_formula.to_dict

            nominal_mass = molecular_formula.mz_nominal_calc

            if min_mz <= nominal_mass <= max_mz:
                
                if min_dbe <= dbe <= max_dbe:
                    
                    dict_results = {"mol_formula" : json.dumps(formula_dict),
                                    "mz" : mz,
                                    "nominal_mz" : nominal_mass,
                                    "ion_type" : ion_type,
                                    "ion_charge" : ion_charge,
                                    "classe" : class_str,
                                    "C" : molecular_formula.get('C'),
                                    "H" : molecular_formula.get('H'),
                                    "N" : molecular_formula.get('N'),
                                    "O" : molecular_formula.get('O'),
                                    "S" : molecular_formula.get('S'),
                                    "P" : molecular_formula.get('P'),
                                    "H_C" : molecular_formula.get('H')/molecular_formula.get('C'),
                                    "O_C" : molecular_formula.get('O')/molecular_formula.get('C'),
                                    "DBE" : dbe}
                                
                    yield dict_results
               
    def get_h_odd_or_even(self, ion_type, class_dict, molecular_search_settings):

        
        HAS_NITROGEN = 'N' in class_dict.keys()
        HAS_PHOSPHORUS = 'P' in class_dict.keys()

        number_of_halogen = self.get_total_halogen_atoms(class_dict, molecular_search_settings)

        if number_of_halogen > 0:

            TEM_HALOGEN = True

        else:

            TEM_HALOGEN = False

        if TEM_HALOGEN:

            remaining_halogen = number_of_halogen % 2

        else:

            remaining_halogen = 0

        if HAS_NITROGEN and HAS_PHOSPHORUS:

            number_of_n = class_dict.get('N') + class_dict.get('P')
            remaining_n = number_of_n % 2

        elif HAS_NITROGEN and not HAS_PHOSPHORUS:

            number_of_n = class_dict.get('N')
            remaining_n = number_of_n % 2

        elif HAS_PHOSPHORUS and not HAS_NITROGEN:

            number_of_n = class_dict.get('P')
            remaining_n = number_of_n % 2

        else:

            remaining_n = -1

        if ion_type == 'DE_OR_PROTONATED':
            # has nitrogen and is odd hydrogen, has to start as even
            if remaining_n > 0.0:
                if HAS_NITROGEN or HAS_PHOSPHORUS and remaining_n > 0:

                    if TEM_HALOGEN:
                        if remaining_halogen == 0:
                            return 'even'
                        else:
                            return 'odd'
                    else:

                        return 'even'

            if remaining_n == 0.0:

                if HAS_NITROGEN or HAS_PHOSPHORUS and remaining_n == 0:

                    if TEM_HALOGEN:
                        if remaining_halogen == 0:
                            return 'odd'
                        else:
                            return 'even'
                    else:
                        return 'odd'

            else:

                if TEM_HALOGEN:
                    if remaining_halogen == 0:
                        return 'odd'
                    else:
                        return 'even'
                else:
                    return 'odd'

        elif ion_type == 'RADICAL':
            
            # has nitrogen and is odd hydrogen, has to start as odd
            if remaining_n > 0.0:
                if HAS_NITROGEN or HAS_PHOSPHORUS:

                    if TEM_HALOGEN:
                        if remaining_halogen == 0:
                            return 'odd'
                        else:
                            return 'even'
                    else:
                        return 'odd'

            elif remaining_n == 0.0:

                if HAS_NITROGEN or HAS_PHOSPHORUS:

                    if TEM_HALOGEN:
                        if remaining_halogen == 0:
                            return 'even'
                        else:
                            return 'odd'
                    else:
                        return 'even'

            else:

                if TEM_HALOGEN:
                    if remaining_halogen == 0:
                        return 'even'
                    else:
                        return 'odd'
                else:
                    return 'even'

    def get_dbe_limits(self, classe_dict, use_pah_line_rule, formula_dict, min_dbe, max_dbe):

        sum_hetero_atoms = 0
        for i in classe_dict.keys():
            if i != Labels.ion_type:

                sum_hetero_atoms = sum_hetero_atoms + classe_dict.get(i)

        if not use_pah_line_rule:

            minDBE = min_dbe
            maxDBE = max_dbe

        else:

            minDBE = 0

            maxDBE = (int(formula_dict.get('C')) + sum_hetero_atoms) * 0.9

        return maxDBE, minDBE

    @staticmethod
    
    def get_total_halogen_atoms(class_dict, molecular_search_settings):

            atoms = ['F', 'Cl', 'Br', 'I']

            return 0

            total_number = 0
            
            for atom in atoms:

                TEM_HALOGEN = atom in class_dict.keys()

                if TEM_HALOGEN:

                    valencia = molecular_search_settings.used_atom_valences.get(atom)

                    if valencia != 0:
                        
                        total_number = total_number + class_dict.get(atom)
            
            return total_number    
