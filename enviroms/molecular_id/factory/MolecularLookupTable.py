__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

import itertools
import multiprocessing
from enviroms.molecular_id.factory.MolecularFormulaFactory import MolecularFormula
from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings
from enviroms.encapsulation.constant import Labels
from copy import deepcopy

class MolecularCombinations:
     
    '''
    runworker()
    Returns a dictionary of molecular formula objs inside a dict nominal masses and ion type
        {Labels.ion_type:{
            classes:{
                nominal_mass:{
                    [MolecularFormula objs,]
                }
            }
        }
    note:
    waiting for python 3.8 release (set 2019)  to have share memory between subprocesses;
    by then we need to pass the resulting dict and all settings into shared memory;
    the current serialization is adding 1.5s seconds for each ion type iteration
    '''

    def runworker(self,settings) :

        c_h_combinations= self.get_c_h_combination(settings)
        
        classes_list = self.get_classes_in_order(settings)
        
        number_of_process = int(multiprocessing.cpu_count())

        #number_of_process = psutil.cpu_count(logical=False)

        print('number_of_process', number_of_process)

        print('classes_list', len(classes_list))

        '''exited with code=0 in 6.9 seconds windows and 5.9 Linux with 12 logical CPUs'''
        p = multiprocessing.Pool(number_of_process)
        args = [(class_tuple, c_h_combinations, settings) for class_tuple in classes_list]
        results = p.map(CombinationsWorker(), args)
        p.close()
        p.join()
        '''
        args = [(class_tuple, isProtonated, isRadical, use_pah_line_rule, usedAtoms, 
                min_dbe, max_dbe, ion_charge, c_h_combinations) for class_tuple in classes_list]
        results = list()
        
        for arg in args:
            #exited with code=0 in 17.444 seconds
            results.append(CombinationsWorker().get_combinations(*arg))
        '''
        class_list_str= (x[0] for x in classes_list)
        
        return dict(zip(class_list_str, results))
        
    def get_c_h_combination(self, settings):

        # return dois dicionarios com produto das combinacooes de hidrogenio e carbono
        # para nitrogenio impar e par para radicais e protonados
        usedAtoms = settings.usedAtoms
        result = {}
        
        min_c, max_c = usedAtoms.get('C')
        min_h, max_h = usedAtoms.get('H')

        min_h_fix = self.get_fixed_initial_number_of_hidrogen(min_h, 'odd')
        possible_c = [c for c in range(min_c, max_c + 1)]

        possible_h = [h for h in range(min_h_fix, max_h + 2, 2)]

        list_products = [i for i in itertools.product(possible_c, possible_h)]

        result['odd'] = list_products
        
        min_h_fix = self.get_fixed_initial_number_of_hidrogen(min_h, 'even')
        possible_h = [h for h in range(min_h, max_h + 2, 2)]

        list_products = [i for i in itertools.product(possible_c, possible_h)]

        result['even'] = list_products

        return result

    def swap_class_order(self, class_frist, class_second, new_list2):

        if class_frist in new_list2:

            if class_second in new_list2:
                n_index, s_index = (
                    new_list2.index(class_frist),
                    new_list2.index(class_second),
                )

                new_list2[n_index], new_list2[s_index] = (
                    new_list2[s_index],
                    new_list2[n_index],
                )

        return new_list2    
    
    def get_classes_in_order(self, settings ):
        ''' structure is 
            ('HC', {'HC': 1})'''
        
        usedAtoms = deepcopy(settings.usedAtoms )
        
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
        
        atomos_in_ordem = ['N', 'O', 'S', 'P']

        classe_in_orderm = []

        all_atomos_tuples = itertools.product(possible_n, possible_o,
                                            possible_s, possible_p)
        
        for atomo in atomos_in_ordem:
            usedAtoms.pop(atomo, None)
        
        for selected_atomo, min_max_tuple in usedAtoms.items():
            
            min_x = min_max_tuple[0]
            max_x = min_max_tuple[1]
            # massa =  Recal_Assign_Settings.selection_of_atomos[selected_atomo][2]

            possible_x = [x for x in range(min_x, max_x + 1)]

            all_atomos_tuples = itertools.product(all_atomos_tuples, possible_x)
            all_atomos_tuples = [all_atomos_combined[0] + (all_atomos_combined[1],) for all_atomos_combined in
                                all_atomos_tuples]
            atomos_in_ordem.append(selected_atomo)
        
        for all_atomos_tuple in all_atomos_tuples:

            classe_str = ''
            classe_dict = {}
            
            for each_atomos_index, atom_number in enumerate(all_atomos_tuple):
                
                if atom_number != 0:
                    classe_str = (classe_str + atomos_in_ordem[each_atomos_index] + str(atom_number) +  ' ')
                    classe_dict[atomos_in_ordem[each_atomos_index]] = atom_number

            classe_str = classe_str.strip()
            
            if len(classe_str) > 0:
                
                classe_in_orderm.append((classe_str, classe_dict))

            elif len(classe_str) == 0:

                classe_in_orderm.append(('HC', {'HC': 1}))
        
        classe_in_orderm = self.sort_classes(atomos_in_ordem, classe_in_orderm)
        
        return classe_in_orderm

    @staticmethod
    def sort_classes( atomos_in_ordem, combination_tuples) -> [str]: 
        
        join_list_of_list_classes = list()
        atomos_in_ordem =  ['N','S','P', 'O'] + atomos_in_ordem[4:] + ['HC']
        
        sort_method = lambda atoms_keys: [atomos_in_ordem.index(atoms_keys)] #(len(word[0]), print(word[1]))#[atomos_in_ordem.index(atom) for atom in list( word[1].keys())])
        for class_tuple in combination_tuples:
            
            sorted_dict_keys = sorted(class_tuple[1], key = sort_method)
            class_str = ' '.join([atom + str(class_tuple[1][atom]) for atom in sorted_dict_keys])
            class_dict = { atom: class_tuple[1][atom] for atom in sorted_dict_keys}
            join_list_of_list_classes.append((class_str, class_dict))
        
        return join_list_of_list_classes

    
    def get_fixed_initial_number_of_hidrogen(self, min_h, odd_even):

        remaining_h = min_h % 2
        
        if odd_even == 'even':
            
            if remaining_h == 0: return remaining_h
            
            else: return remaining_h + 1    
        
        else:
            
            if remaining_h == 0: return remaining_h + 1
            
            else: return remaining_h    

class CombinationsWorker:

    # needs this wraper to pass the class to multiprocessing
    def __call__(self, args):

        return self.get_combinations(*args)  # ,args[1]

    def get_combinations(self, classe_tuple,
                          c_h_combinations,settings
                         ):

        usedAtoms = settings.usedAtoms
        
        min_dbe = settings.min_dbe

        max_dbe = settings.max_dbe

        use_pah_line_rule = settings.use_pah_line_rule

        isRadical = settings.isRadical

        isProtonated = settings.isProtonated

        ion_charge = settings.ion_charge
        
        min_mz = settings.min_mz
        
        max_mz = settings.max_mz
        
        hc_filter = settings.hc_filter
        
        oc_filter = settings.oc_filter
        
        dict_results = {}
        
        
        #class_str = classe_tuple[0]
        class_dict = classe_tuple[1]

        if isProtonated:

            ion_type = 'DE_OR_PROTONATED'

            par_ou_impar = self.get_h_impar_ou_par(ion_type, class_dict)

            carbon_hidrogen_combination = c_h_combinations.get(par_ou_impar)

            dict_results_per_ion_type = self.get_theoretical_formulas_per_class(carbon_hidrogen_combination, ion_type, class_dict, 
                                                    use_pah_line_rule, usedAtoms, min_dbe, max_dbe,
                                                    min_mz, max_mz, hc_filter,oc_filter, ion_charge)
            dict_results[ion_type] = dict_results_per_ion_type
        
        if isRadical:

            ion_type = 'RADICAL'

            par_ou_impar = self.get_h_impar_ou_par(ion_type, class_dict)

            carbon_hidrogen_combination = c_h_combinations.get(par_ou_impar)

            dict_results_per_ion_type = self.get_theoretical_formulas_per_class(carbon_hidrogen_combination, ion_type, class_dict, 
                                                    use_pah_line_rule, usedAtoms, min_dbe, max_dbe, 
                                                    min_mz, max_mz, hc_filter,oc_filter, ion_charge)
            dict_results[ion_type] = dict_results_per_ion_type
        
        return dict_results

    def get_theoretical_formulas_per_class(self,
                                           carbon_hidrogen_combination,
                                           ion_type,
                                           class_dict,
                                           use_pah_line_rule,
                                           usedAtoms,
                                           min_dbe,
                                           max_dbe,
                                           min_mz,
                                           max_mz, 
                                           hc_filter, 
                                           oc_filter,
                                           ion_charge,
                                           ):
       
        dict_results = {}
        for cada_possible in carbon_hidrogen_combination:

            c_number = cada_possible[0]
            h_number = cada_possible[1]
            o_number = class_dict.get('O')
            continuar = True

            if float(h_number / c_number) >= hc_filter:
                
                if o_number:

                    if float(o_number) / c_number > oc_filter:
                        continuar = False

                if continuar:

                    formula_dict = {}
                    for each_atom in class_dict.keys() :
                        if each_atom != 'HC':
                            formula_dict[each_atom] = class_dict.get(each_atom)

                    formula_dict['C'] = c_number
                    formula_dict['H'] = h_number
                    formula_dict[Labels.ion_type] = ion_type

                    molecular_formula = MolecularFormula(formula_dict, ion_charge)
                    DBE = molecular_formula.dbe
                    nominal_mass = molecular_formula.mz_nominal_theo

                    # one second overhead to create and serialize the molecular formula object
                    #DBE = self.get_DBE(formula_dict, 1)
                    #nominal_mass = int(self.getMass(formula_dict, ion_charge))

                    if min_mz <= nominal_mass <= max_mz:
                        maxDBE, minDBE = self.get_dbe_limits(
                            class_dict,
                            use_pah_line_rule,
                            formula_dict,
                            min_dbe,
                            max_dbe,
                        )

                        if minDBE <= DBE <= maxDBE:
                            
                            if nominal_mass in dict_results.keys():

                                dict_results[nominal_mass].append(molecular_formula)

                            else:

                                dict_results[nominal_mass] = [molecular_formula]
            
        return dict_results
    
    def get_h_impar_ou_par(self, ion_type, class_dict):

        TEM_NITROGENIO = 'N' in class_dict.keys()
        TEM_PHOSPOROUS = 'P' in class_dict.keys()

        number_of_halogen = self.get_total_halogen_atoms(class_dict)

        if number_of_halogen > 0:

            TEM_HALOGEN = True

        else:

            TEM_HALOGEN = False

        if TEM_HALOGEN:

            remaining_halogen = number_of_halogen % 2

        else:

            remaining_halogen = 0

        if TEM_NITROGENIO and TEM_PHOSPOROUS:

            number_of_n = class_dict.get('N') + class_dict.get('P')
            remaining_n = number_of_n % 2

        elif TEM_NITROGENIO and not TEM_PHOSPOROUS:

            number_of_n = class_dict.get('N')
            remaining_n = number_of_n % 2

        elif TEM_PHOSPOROUS and not TEM_NITROGENIO:

            number_of_n = class_dict.get('P')
            remaining_n = number_of_n % 2

        else:

            remaining_n = -1

        if ion_type == 'DE_OR_PROTONATED':
            # tem nitrogenio e eh impar hidrogenio tem que comecar com par
            if remaining_n > 0.0:
                if TEM_NITROGENIO or TEM_PHOSPOROUS and remaining_n > 0:

                    if TEM_HALOGEN:
                        if remaining_halogen == 0:
                            return 'even'
                        else:
                            return 'odd'
                    else:

                        return 'even'

            if remaining_n == 0.0:

                if TEM_NITROGENIO or TEM_PHOSPOROUS and remaining_n == 0:

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
            # tem nitrogenio e eh impar hidrogenio tem que comecar com impar

            if remaining_n > 0.0:
                if TEM_NITROGENIO or TEM_PHOSPOROUS:

                    if TEM_HALOGEN:
                        if remaining_halogen == 0:
                            return 'odd'
                        else:
                            return 'even'
                    else:
                        return 'odd'

            elif remaining_n == 0.0:

                if TEM_NITROGENIO or TEM_PHOSPOROUS:

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

        sum_hetero_atomos = 0
        for i in classe_dict.keys():
            if i != Labels.ion_type:

                sum_hetero_atomos = sum_hetero_atomos + classe_dict.get(i)

        # sum_hetero_atomos = sum(classe_dict.values())

        if not use_pah_line_rule:

            minDBE = min_dbe
            maxDBE = max_dbe

        else:

            minDBE = 0

            maxDBE = (int(formula_dict.get('C')) + sum_hetero_atomos) * 0.9

        return maxDBE, minDBE

    @staticmethod
    
    def get_total_halogen_atoms(class_dict):

            atoms = ['F', 'Cl', 'Br', 'I']

            total_number = 0
            '''
            for atom in atoms:

                TEM_HALOGEN = in class_dict.keys(atom)

                if TEM_HALOGEN:

                    valencia = Recal_Assign_Settings.selection_of_atomos.get(atom)[3]

                    if valencia == 1:
                        total_number = total_number + class_dict.get(atom)
            '''
            return total_number    
