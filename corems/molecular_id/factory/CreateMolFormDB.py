__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

from copy import deepcopy
import itertools
import multiprocessing
import pickle
import json
import cProfile
import io
import pstats
import contextlib
import time

from tqdm import tqdm

from sqlalchemy.orm import scoped_session
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, func

from corems.encapsulation.factory.processingSetting  import MolecularLookupDictSettings
from corems.encapsulation.constant import Atoms
from corems.molecular_id.factory.molformSQL import CarbonHydrogen, HeteroAtoms, MolecularFormulaLink
from corems.encapsulation.factory.parameters import MSParameters
from corems import timeit
from corems.molecular_id.factory.molformSQL import Base

@contextlib.contextmanager
def profiled():
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    ps.print_stats()
    # uncomment this to see who's calling what
    # ps.print_callers()
    print(s.getvalue())

class NewMolecularCombinations:
     
    def __init__(self, sql_db = None, url='sqlite:///'):

        self.engine = create_engine(url, echo = False)
        session_factory = sessionmaker(bind=self.engine )
        Session = scoped_session(session_factory)
        self.session = session_factory()
        
        Base.metadata.create_all(self.engine)
    
    def cProfile_worker(self, args):
        
        cProfile.runctx('self.get_mol_formulas(*args)', globals(), locals(), 'mf_database_cprofile.prof')

    def check_database_get_class_list(self, molecular_search_settings):
        
        all_class_to_create = []
        
        classes_dict = self.get_classes_in_order(molecular_search_settings)
        
        class_str_set = set(classes_dict.keys())
        
        existing_classes_objs = self.session.query(HeteroAtoms).distinct().all()
        
        existing_classes_str = set([classe.name for classe in existing_classes_objs])

        self.len_existing_classes = len(existing_classes_str)

        class_to_create = class_str_set - existing_classes_str
        
        for class_str in class_to_create:
            
            class_tuple =  (class_str, classes_dict.get(class_str)) 
            
            all_class_to_create.append(class_tuple)

        #print (existing_classes)class_str_list
        
        return [(c_s, c_d) for c_s, c_d in classes_dict.items()], all_class_to_create       
    
    def get_carbonsHydrogens(self, settings, odd_even):
        
        operator = '==' if odd_even == 'even' else '!=' 

        usedAtoms = settings.usedAtoms
        user_min_c, user_max_c = usedAtoms.get('C')
        user_min_h, user_max_h = usedAtoms.get('H')

        carbonHydrogenObjs =  eval("self.session.query(CarbonHydrogen).filter(" 
                                       "CarbonHydrogen.C >= user_min_c,"
                                        "CarbonHydrogen.H >= user_min_h,"
                                        "CarbonHydrogen.C <= user_max_c,"
                                        "CarbonHydrogen.H <= user_max_h,"
                                        "CarbonHydrogen.C % 2" + operator+ "odd_even).all()")
        
        return carbonHydrogenObjs                                        
    
    def add_carbonsHydrogens(self, settings):

        usedAtoms = settings.usedAtoms

        self.session.query(CarbonHydrogen).distinct().all()

        user_min_c, user_max_c = usedAtoms.get('C')
        user_min_h, user_max_h = usedAtoms.get('H')

        query_obj = self.session.query(func.max(CarbonHydrogen.C).label("max_c"), 
                        func.min(CarbonHydrogen.C).label("min_c"),
                        func.max(CarbonHydrogen.H).label("max_h"),
                        func.min(CarbonHydrogen.C).label("min_h"),
                        )

        database = query_obj.first()
        
        if database.max_c == user_max_c and database.min_c == user_min_c and database.max_h == user_max_h and database.min_h == user_min_h:   
            #yeah we are good, 
            pass

        else:
            
            databaseCarbon = set(self.session.query(CarbonHydrogen.C).all())
            databaseHydrogen = set(self.session.query(CarbonHydrogen.H).all())
            
            userCarbon = set(range(user_min_c, user_max_c + 1))
            userHydrogen = set(range(user_min_h, user_max_h + 1))
            
            carbonCreate = databaseCarbon ^ userCarbon
            hydrogenCreate = databaseHydrogen ^ userHydrogen

            carbon_hydrogen_objs = [CarbonHydrogen(C=i[0],H=i[1]) for i in itertools.product(carbonCreate, hydrogenCreate)]

            self.session.add_all(carbon_hydrogen_objs)  

            

    def runworker(self, molecular_search_settings):
        
        classes_list, class_to_create = self.check_database_get_class_list(molecular_search_settings)
        
        if class_to_create:
            
            settings = MolecularLookupDictSettings()
            settings.usedAtoms = deepcopy(molecular_search_settings.usedAtoms)
            settings.url_database = molecular_search_settings.url_database
            settings.db_jobs = molecular_search_settings.db_jobs
            
            print("OK2")
            self.add_carbonsHydrogens(settings)
            print("OK")
            self.odd_ch = self.get_carbonsHydrogens(settings,'odd')
            self.even_ch = self.get_carbonsHydrogens(settings, 'even')
            print("OK3")
            for class_index, class_tuple in tqdm(enumerate(class_to_create)):
                
                self.populate_combinations(class_tuple, settings, class_index)
            
        return classes_list
    
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

    def calc_mz(self, datadict, class_mass=0):
        
        mass = class_mass
        
        for each_atom in datadict.keys():
            
            if each_atom != 'HC':    
                
                mass = mass + Atoms.atomic_masses[each_atom]  *  datadict.get(each_atom)
            
        return mass 
        
    def calc_dbe(self, datadict, init_dbe=0, final_calc=True):
            
            for atom in datadict.keys():
                  
                n_atom = int(datadict.get(atom))
                
                clean_atom = ''.join([i for i in atom if not i.isdigit()]) 
                
                valencia = MSParameters.molecular_search.used_atom_valences.get(clean_atom)
                
                if valencia and valencia > 0:
                    #print atom, valencia, n_atom, init_dbe
                    init_dbe = init_dbe + (n_atom * (valencia - 2))
                else:
                    continue
            
            if final_calc:
                return 1 + (0.5 * init_dbe)
            else:    
                return init_dbe

    def populate_combinations(self, classe_tuple, settings, class_index):
        
        ion_charge =  0
        
        class_str = classe_tuple[0]
        class_dict = classe_tuple[1]
        
        odd_or_even = self.get_h_odd_or_even(class_dict, settings)
        
        #get_mol_formulas
        #self.cProfile_worker((carbon_hydrogen_objs, classe_tuple, settings))
        self.get_mol_formulas(odd_or_even, classe_tuple, settings, class_index)
        #self.session.commit()    

    def get_or_add(self, SomeClass, kw):
            
        obj = self.session.query(SomeClass).filter_by(**kw).first()
        if not obj:
            obj = SomeClass(**kw)
        return obj
    
    def get_mol_formulas(self, odd_even_tag, classe_tuple, settings, class_index):
        
        class_str = classe_tuple[0]
        class_dict = classe_tuple[1]
        
        #heteroAtom_obj = HeteroAtoms(name=class_str)
        #self.session.add(heteroAtom_obj)
        #if not heteroAtom_obj.id:
            
        #    heteroAtom_obj.id = self.len_existing_classes + (class_index+1) 
        
        #results = list()
        
        if 'HC' in class_dict:
            del class_dict['HC']
        
        class_mass = self.calc_mz(class_dict)
        class_dbe = self.calc_dbe(class_dict)
        
        sum_init_creation = 0
        sum_obj_creation = 0
        
        #self.session.query(CarbonHydrogen).filter(CarbonHydrogen.C%2 == 0).all()
        
        
        self.odd_ch[0].C
        
        #carbonHydrogen_objs = self.odd_ch if odd_even_tag == 'even' else self.even_ch 
        '''
        #carbonHydrogen_objs = self.objs_carbonHydrogen[odd_even_tag]
        #carbonHydrogen_dict = self.dict_carbonHydrogen[odd_even_tag]
        
        for index, carbonHydrogen_obj in enumerate(carbonHydrogen_objs):
            
            timestart = time.time()
            #C, H = carbonHydrogen_obj.C, carbonHydrogen_obj.H,
            #{'C': C, 'H': H}
            dict_ch = carbonHydrogen_dict[index]
            timefinish = time.time()
            sum_init_creation = sum_init_creation + (timefinish - timestart)

            mass = self.calc_mz(dict_ch, class_mass)
            
            DBE = self.calc_dbe(dict_ch, class_dbe, final_calc=True)
            
            if settings.min_mz <= mass <= settings.max_mz:
                
                if settings.min_dbe <= DBE <= settings.max_dbe:
                    
                    #mf = MolecularFormula(mass=mass, DBE = DBE)
                    #mf.carbon_hydrogen = carbonHydrogen_obj
                    timestart = time.time()
                    
                    molecularFormula = MolecularFormulaLink(heteroAtoms=heteroAtom_obj, 
                                        carbonHydrogen=carbonHydrogen_obj, mass=mass, DBE=DBE)
                    
                    results.append(molecularFormula)
                    
                    sum_obj_creation = sum_obj_creation + (time.time() - timestart)
                    #heteroAtom_obj.carbonHydrogen.append(carbonHydrogen_obj)
                    #self.sqldb.session.add(mf)
                    #ha = HeteroAtoms(class_label='N5')
                    #mf = MolecularFormula(mass=203.1, DBE = 11)
                    #ch = CarbonHydrogen(C=10, H=27)

                    #mf.carbon_hydrogen = ch
                    #ha.carbon_hydrogen.append(mf)
        
        self.session.add_all(results)
        self.session.commit()
        '''
        
    def get_h_odd_or_even(self, class_dict, molecular_search_settings):

        
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

    @staticmethod
    def get_total_halogen_atoms(class_dict, molecular_search_settings):

            atoms = ['F', 'Cl', 'Br', 'I']

            total_number = 0
            
            for atom in atoms:

                if atom in class_dict.keys():

                    total_number = total_number + class_dict.get(atom)
            
            return total_number    
