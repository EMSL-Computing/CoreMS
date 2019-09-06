__author__ = "Yuri E. Corilo"
__date__ = "Set 06, 2019"

from enviroms.encapsulation.Constants import Atoms
from threading import Thread

class MassSpectoCSV(Thread):
    '''
    classdocs
    '''
    def __init__(self, out_file_path, mass_spectrum):
        
        '''
        Constructor
        '''
        Thread.__init__(self)

        self.file_location = out_file_path
        
        self.mass_spectrum  = mass_spectrum

        self.atomos_order_list = self.get_all_used_atoms_in_ordem()

        self.csv_columns = ['Measured m/z',
                            'Calibrated m/z',
                            'Calculated m/z',
                            'Measured Abundance', 
                            'Measured Resolving Power',
                            'Signal/Noise',
                            'Mass Error (ppm)',
                            'DBE',
                            'H/C', 
                            'O/C',  
                            #'Aromaticity Index',
                            'Heteroatom class']
    def run(self):
        
        dict_data_list = self.get_list_dict_data()
        
        self.csv_columns.extend(self.atomos_order_list)

        self.write_to_file(dict_data_list)

    def write_to_file(self, dict_data_list):
        
        import csv
        try:
            with open(self.file_location, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.csv_columns)
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

    def get_list_dict_data(self, include_no_match=False, include_isotopolgues=True, ):
        
        dict_data_list = []
        for ms_peak in self.mass_spectrum.sort_by_mz():
            
            #check if there is a molecular formula candidate for the msPeak
            if ms_peak:
                #m_formula = ms_peak.molecular_formula_lowest_error
                for m_formula in ms_peak:
                    
                    formula_dict = m_formula.to_dict
                    dict_result = { 'Measured m/z':  ms_peak.mz_exp,
                                    'Calibrated m/z': ms_peak.mz_exp,
                                    'Calculated m/z': m_formula.mz_theor,
                                    'Measured Abundance': ms_peak.abundance,
                                    'Measured Resolving Power': ms_peak.resolving_power,
                                    'Signal/Noise':  ms_peak.signal_to_noise,
                                    'Mass Error (ppm)':m_formula._calc_assigment_mass_error(ms_peak.mz_exp),
                                    'DBE' :  m_formula.dbe,
                                    'H/C' :  m_formula.H_C,
                                    'O/C' :  m_formula.O_C,
                                    #'Aromaticity Index': m_formula.ai,
                                    'Heteroatom class' : m_formula.class_label
                    }
                    
                    for atom in self.atomos_order_list:
                        if atom in formula_dict.keys():
                            dict_result[atom] =  formula_dict.get(atom)

                    if m_formula.is_isotopologue: #skip or not the isotopologues
                        if include_isotopolgues : dict_data_list.append(dict_result)
                    else: dict_data_list.append(dict_result)
                            
            else:
               # include not_match
               if include_no_match:
                     dict_result = { 'Measured m/z':  ms_peak.exp_mz,
                                    'Calibrated m/z': ms_peak.exp_mz,
                                    'Calculated m/z': m_formula.mz_theor,
                                    'Measured Abundance': ms_peak.abundance,
                                    'Measured Resolving Power': ms_peak.rp,
                                    'Signal/Noise':  ms_peak.stn
                     }
                     
                     dict_data_list.append(dict_result)
        
        return dict_data_list            



