__author__ = "Yuri E. Corilo"
__date__ = "Set 06, 2019"

from enviroms.encapsulation.constant import Atoms
from threading import Thread
import pandas as pd

class MassSpecExport(Thread):
    '''
    TODO: add MSPeak indexes: done
    
    '''
    def __init__(self, out_file_path, mass_spectrum, output_type='excel'):
        
        '''
        output_type:str 
            'excel', 'csv' or 'pandas' 
        '''
        Thread.__init__(self)

        self.output_file = out_file_path

        #'excel', 'csv' or 'pandas'
        self.output_type = output_type
        
        self.mass_spectrum  = mass_spectrum

        # collect all assigned atoms and order them accordingly to the Atoms.atoms_order list 
        self.atomos_order_list = self.get_all_used_atoms_in_ordem()
        
        self._init_columns()

        
    def _init_columns(self):
        
        #column labels in order
        self.columns = [    'Index',
                            'Measured m/z',
                            'Calibrated m/z',
                            'Calculated m/z',
                            'Measured Abundance', 
                            'Measured Resolving Power',
                            'Signal/Noise',
                            'Mass Error (ppm)',
                            'DBE',
                            'H/C', 
                            'O/C',  
                            'Heteroatom Class',
                            'Ion Type',
                            'Is Isotopologue',
                            #'Aromaticity Index',
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
        df.name =  self.output_file
        return df

    def to_pandas(self, dict_data_list):
        
        df = pd.DataFrame(dict_data_list, columns=self.columns)
        
        df.to_pickle(self.output_file + '.pkl')

    def to_excel(self, dict_data_list):
        
        df = pd.DataFrame(dict_data_list, columns=self.columns)
        
        df.to_excel(self.output_file + '.xlsx')

    def to_csv(self, dict_data_list):
        
        import csv
        try:
            with open(self.output_file +'.csv', 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.columns)
                writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)
        except IOError as ioerror:
            print(ioerror) 
    
    def get_all_used_atoms_in_ordem(self):
        
        atomos_in_order = Atoms.atoms_order
        all_used_atoms = set()
        for ms_peak in self.mass_spectrum:
             for m_formula in ms_peak:
                 if ms_peak:
                     for atom in m_formula.atoms:
                        all_used_atoms.add(atom)

        sort_method = lambda atom: [atomos_in_order.index(atom)]

        return sorted(all_used_atoms, key = sort_method)

    def get_list_dict_data(self, include_no_match=True, include_isotopolgues=True, 
                            isotopologue_inline=False, no_match_inline=False):
        
        dict_data_list = []

        def add_no_match_dict_data():
            
            dict_result = { 'Index': index,
                            'Measured m/z':  ms_peak.mz_exp,
                            'Calibrated m/z': ms_peak.mz_exp,
                            'Measured Abundance': ms_peak.abundance,
                            'Measured Resolving Power': ms_peak.resolving_power,
                            'Signal/Noise':  ms_peak.signal_to_noise}
                     
            dict_data_list.append(dict_result)

        def add_match_dict_data():
            
            formula_dict = m_formula.to_dict
            dict_result = { 'Index': index,
                            'Measured m/z':  ms_peak.mz_exp,
                            'Calibrated m/z': ms_peak.mz_exp,
                            'Calculated m/z': m_formula.mz_theor,
                            'Measured Abundance': ms_peak.abundance,
                            'Measured Resolving Power': ms_peak.resolving_power,
                            'Signal/Noise':  ms_peak.signal_to_noise,
                            'Mass Error (ppm)':m_formula._calc_assigment_mass_error(ms_peak.mz_exp),
                            'DBE' :  m_formula.dbe,
                            'Heteroatom Class' : m_formula.class_label,
                            'H/C' :  m_formula.H_C,
                            'O/C' :  m_formula.O_C,
                            'Ion Type' : m_formula.ion_type.lower(),
                            'Is Isotopologue' : int(m_formula.is_isotopologue),
                            }

            for atom in self.atomos_order_list:
                if atom in formula_dict.keys():
                    dict_result[atom] =  formula_dict.get(atom)
            
            dict_data_list.append(dict_result)
        
        for index, ms_peak in enumerate(self.mass_spectrum.sort_by_mz()):
            
            #check if there is a molecular formula candidate for the msPeak
            if ms_peak:
                #m_formula = ms_peak.molecular_formula_lowest_error
                for m_formula in ms_peak:
                    
                    if m_formula.is_isotopologue: # isotopologues inline
                        if include_isotopolgues and isotopologue_inline: add_match_dict_data()
                    else: add_match_dict_data() #add monoisotopic peak
                           
            else:
                # include not_match
                if include_no_match and no_match_inline:
                    add_no_match_dict_data()
        
        if include_isotopolgues and not isotopologue_inline:
            for index, ms_peak in enumerate(self.mass_spectrum.sort_by_mz()):
                for m_formula in ms_peak:
                    if m_formula.is_isotopologue:
                        add_match_dict_data()
                        
        if include_no_match and not no_match_inline:
            for index, ms_peak in enumerate(self.mass_spectrum.sort_by_mz()):
                if not ms_peak:
                    add_no_match_dict_data()
        return dict_data_list            



