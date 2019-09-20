__author__ = "Yuri E. Corilo"
__date__ = "Set 06, 2019"

import re
from threading import Thread

import pandas as pd

from enviroms.encapsulation.Constants import Atoms

class MolecularLookUpDictExport(Thread):
    '''
    TODO: add MSPeak indexes: done
    
    '''
    def __init__(self, out_file_path, dict_molecular_lookup_dict, output_type='pandas'):
        
        '''
        output_type:str 
            'excel', 'csv' or 'pandas' 
        '''
        Thread.__init__(self)

        self.file_location = out_file_path

        #'excel', 'csv' or 'pandas'
        self.output_type = output_type
        
        self.dict_molecular_lookup_dict  = dict_molecular_lookup_dict

        # collect all assigned atoms and order them accordingly to the Atoms.atoms_order list 
        self.atomos_order_list = self.get_all_used_atoms_in_ordem()

        #column labels in order
        self.columns = ['Heteroatom Class',
                        'DBE', 
                        'H/C' ,
                        'O/C' ,
                        'Ion Type' ,
                        'Is Isotopologue',
                        ]
    
    @property
    def output_type(self):
        return self._output_type
    
    @output_type.setter
    def output_type(self, output_type):
        output_types = ['excel', 'csv', 'pandas']
        if output_type in output_types:
            self._output_type = output_type
        else:
            raise TypeError('Supported types are "excel", "csv" or "pandas", %s entered' % output_type)    

    def run(self):
        
        dict_data_list = self.get_list_dict_data()
        
        #add atoms labels to the columns 
        self.columns.extend(self.atomos_order_list)

        if self.output_type == 'excel':
            self.to_excel(dict_data_list)
        elif self.output_type == 'csv': 
            self.to_csv(dict_data_list)
        elif self.output_type == 'pandas':
            self.to_pandas(dict_data_list)
        else:
            raise ValueError("Unkown output type: %s; it can be 'excel', 'csv' or 'pandas'" %self.output_type)
    
    def get_pandas_df(self):
        
        self.columns.extend(self.atomos_order_list)
        dict_data_list = self.get_list_dict_data()
        df = pd.DataFrame(dict_data_list, columns=self.columns)
        df.name =  self.output_type
        return df

    def to_pandas(self, dict_data_list):
        
        df = pd.DataFrame(dict_data_list, columns=self.columns)
        
        df.to_pickle(self.file_location + '.pkl',  index=False)

    def to_excel(self, dict_data_list):
        
        df = pd.DataFrame(dict_data_list, columns=self.columns)
        
        df.to_excel(self.file_location + '.xlsx',  index=False)

    def to_csv(self, dict_data_list):
        
        import csv
        try:
            with open(self.file_location, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.columns)
                writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)
        except IOError as ioerror:
            print(ioerror) 
    
    def get_all_used_atoms_in_ordem(self):
        
        atomos_in_order = Atoms.atoms_order
        all_used_atoms = set()
        for classe in self.dict_molecular_lookup_dict.keys():
             if classe != 'HC1':
                atoms = classe.split(' ')
                for atom in atoms:
                    #removes atoms qnt but retain atomic mass number
                    clean_atom = re.sub('([0-9]+)$','',atom)    
                    all_used_atoms.add(clean_atom)

        sort_method = lambda atom: [atomos_in_order.index(atom)]

        return sorted(all_used_atoms, key = sort_method)

    def get_list_dict_data(self, include_no_match=False, include_isotopolgues=True, 
                            isotopologue_inline=False, no_match_inline=False):
        
        dict_data_list = []

        def add_isotopologue_data():
            
            return {'Heteroatom Class' : m_formula.class_label,
                    'DBE' :  m_formula.dbe,
                    'H/C' :  m_formula.H_C,
                    'O/C' :  m_formula.O_C,
                    'Ion Type' : m_formula.ion_type.lower(),
                    'Is Isotopologue' : int(m_formula.is_isotopologue),
                    'Probability Ratio' : m_formula.prob_ratio
                    }    
            
        def add_mfobj_dict_data():
            
            return {'Heteroatom Class' : m_formula.class_label,
                    'DBE' :  m_formula.dbe,
                    'H/C' :  m_formula.H_C,
                    'O/C' :  m_formula.O_C,
                    'Ion Type' : m_formula.ion_type.lower(),
                    'Is Isotopologue' : int(m_formula.is_isotopologue),
                    }

        for dict_ion_type in self.dict_molecular_lookup_dict.values():
            
            for dict_nominal_mass in dict_ion_type.values():
                
                for m_formula in dict_nominal_mass.values():

                    if m_formula.is_isotopologue:
                        dict_result = add_isotopologue_data()
                    else:
                        dict_result = add_mfobj_dict_data()     
                    
                    add_mfobj_dict_data()
                    
                    formula_dict = m_formula.to_dict

                    for atom in self.atomos_order_list:
                        if atom in formula_dict.keys():
                            dict_result[atom] =  formula_dict.get(atom)
                    
                    dict_data_list.append(dict_result)

        return dict_data_list           