__author__ = "Yuri E. Corilo"
__date__ = "Set 06, 2019"

from threading import Thread

from corems.encapsulation.constant import Atoms
from corems.encapsulation.settings.io import settings_parsers
from numpy import string_, array, NaN
from pandas import DataFrame


class MassSpecExport(Thread):
    '''
    TODO: add MSPeak indexes: done
    '''
    def __init__(self, out_file_path, mass_spectrum, output_type='excel'):
        '''
        output_type:str
            'excel', 'csv', 'hdf5' or 'pandas'
        '''
        Thread.__init__(self)

        self.output_file = out_file_path

        # 'excel', 'csv' or 'pandas'
        self.output_type = output_type

        self.mass_spectrum = mass_spectrum

        # collect all assigned atoms and order them accordingly to the Atoms.atoms_order list
        self.atomos_order_list = self.get_all_used_atoms_in_ordem()

        self._init_columns()

    def _init_columns(self):

        # column labels in order
        self.columns_label = ['Index',
                        'm/z',
                        'Calibrated m/z',
                        'Calculated m/z',
                        'Abundance',
                        'Resolving Power',
                        'S/N',
                        'Ion Charge',
                        'Mass Error (ppm)',
                        'DBE',
                        'H/C',
                        'O/C',
                        'Heteroatom Class',
                        'Ion Type',
                        'Is Isotopologue',
                        # 'Aromaticity Index',
                        ]

    @property
    def output_type(self):
        return self._output_type

    @output_type.setter
    def output_type(self, output_type):
        output_types = ['excel', 'csv', 'pandas', 'hdf5']
        if output_type in output_types:
            self._output_type = output_type
        else:
            raise TypeError(
                'Supported types are "excel", "csv" or "pandas", %s entered' % output_type)

    def save(self):
        '''wrapper to run in a separated thread'''

        if self.output_type == 'excel':
            self.to_excel()
        elif self.output_type == 'csv':
            self.to_csv()
        elif self.output_type == 'pandas':
            self.to_pandas()
        elif self.output_type == 'hdf5':
            self.to_pandas()    
        else:
            raise ValueError(
                "Unkown output type: %s; it can be 'excel', 'csv' or 'pandas'" % self.output_type)

    def run(self):
        ''' run is called when the Tread start
            call exportMS.start() '''
        self.save()

    def get_pandas_df(self):

        columns = self.columns_label + self.atomos_order_list
        dict_data_list = self.get_list_dict_data()
        df = DataFrame(dict_data_list, columns=columns)
        df.name = self.output_file
        return df

    def to_pandas(self):
        
        columns = self.columns_label + self.atomos_order_list

        dict_data_list = self.get_list_dict_data()

        df = DataFrame(dict_data_list, columns=columns)

        df.to_pickle(self.output_file + '.pkl')

    def to_excel(self):

        columns = self.columns_label + self.atomos_order_list

        dict_data_list = self.get_list_dict_data()

        df = DataFrame(dict_data_list, columns=columns)

        df.to_excel(self.output_file + '.xlsx')

    def to_csv(self):
        
        columns = self.columns_label + self.atomos_order_list

        dict_data_list = self.get_list_dict_data()

        import csv
        try:
            with open(self.output_file + '.csv', 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)
        except IOError as ioerror:
            print(ioerror)

    
    def to_hdf(self, ):
        
        import h5py
        import json

        list_results = self.list_dict_to_list()
        
        dict_ms_attrs = self.none_to_nan(self.get_mass_spec_attrs())
        
        setting_dicts = settings_parsers.get_dict_data()

        with h5py.File(self.output_file + '.hdf5', 'w') as f:
        
            dset = f.create_dataset(str(self.mass_spectrum.scan_number), data=list_results)
            dset.attrs['ColumnsLabels'] = json.dumps(self.columns_label + self.atomos_order_list)
            dset.attrs['MassSpecAttrs'] = json.dumps(dict_ms_attrs)
            
            for setting_label, setting_dict in setting_dicts.items(): 
                if setting_label != 'Transient':
                    dset.attrs[setting_label+'Setting'] =  json.dumps(self.none_to_nan(setting_dict))

            #dset.attrs.update(dict_ms_attrs)
          
    def none_to_nan(self, dict_data):
        
        for key, values in dict_data.items():
            if not values: dict_data[key] = NaN
        
        return dict_data   

    def get_mass_spec_attrs(self):

        dict_ms_attrs = {}
        dict_ms_attrs['polarity'] =     self.mass_spectrum.polarity
        dict_ms_attrs['rt'] =     self.mass_spectrum.rt
        dict_ms_attrs['mobility_scan'] =     self.mass_spectrum.mobility_scan
        dict_ms_attrs['mobility_rt'] =     self.mass_spectrum.mobility_rt
        dict_ms_attrs['Aterm'] =  self.mass_spectrum.Aterm
        dict_ms_attrs['Bterm'] =  self.mass_spectrum.Bterm
        dict_ms_attrs['Cterm'] =  self.mass_spectrum.Cterm
        dict_ms_attrs['baselise_noise'] =  self.mass_spectrum.baselise_noise
        dict_ms_attrs['baselise_noise_std'] =  self.mass_spectrum.baselise_noise_std
        
        
        return dict_ms_attrs


    def get_all_used_atoms_in_ordem(self):

        atomos_in_order = Atoms.atoms_order
        all_used_atoms = set()
        for ms_peak in self.mass_spectrum:
            for m_formula in ms_peak:
                if ms_peak:
                    for atom in m_formula.atoms:
                        all_used_atoms.add(atom)

        def sort_method(atom):
            return [atomos_in_order.index(atom)]

        return sorted(all_used_atoms, key=sort_method)

    def list_dict_to_list(self):
        
        column_labels = self.columns_label + self.atomos_order_list

        dict_list = self.get_list_dict_data()
        
        all_lines = []
        for dict_res in dict_list:
            
            result_line = [NaN] * len(column_labels)
            
            for label, value in dict_res.items():
                
                label_index = column_labels.index(label)
                result_line[label_index] =  value   
            
            all_lines.append(result_line)
        return  all_lines       

    
    def get_list_dict_data(self, include_no_match=True, include_isotopolgues=True,
                           isotopologue_inline=False, no_match_inline=False):

        dict_data_list = []

        def add_no_match_dict_data():

            dict_result = {'Index': index,
                           'm/z':  ms_peak.mz_exp,
                           'Calibrated m/z': ms_peak.mz_exp,
                           'Abundance': ms_peak.abundance,
                           'Resolving Power': ms_peak.resolving_power,
                           'S/N':  ms_peak.signal_to_noise,
                           'Ion Charge': ms_peak.ion_charge,
                           }

            dict_data_list.append(dict_result)

        def add_match_dict_data():

            formula_dict = m_formula.to_dict
            dict_result = {'Index': index,
                           'm/z':  ms_peak.mz_exp,
                           'Calibrated m/z': ms_peak.mz_exp,
                           'Calculated m/z': m_formula.mz_theor,
                           'Abundance': ms_peak.abundance,
                           'Resolving Power': ms_peak.resolving_power,
                           'S/N':  ms_peak.signal_to_noise,
                           'Ion Charge': ms_peak.ion_charge,
                           'Mass Error (ppm)': m_formula._calc_assigment_mass_error(ms_peak.mz_exp),
                           'DBE':  m_formula.dbe,
                           'Heteroatom Class': m_formula.class_label.encode('utf-8'),
                           'H/C':  m_formula.H_C,
                           'O/C':  m_formula.O_C,
                           'Ion Type': m_formula.ion_type.lower().encode('utf-8'),
                           'Is Isotopologue': int(m_formula.is_isotopologue),
                           }

            for atom in self.atomos_order_list:
                if atom in formula_dict.keys():
                    dict_result[atom] = formula_dict.get(atom)

            dict_data_list.append(dict_result)

        for index, ms_peak in enumerate(self.mass_spectrum.sort_by_mz()):

            # check if there is a molecular formula candidate for the msPeak
            if ms_peak:
                #m_formula = ms_peak.molecular_formula_lowest_error
                for m_formula in ms_peak:

                    if m_formula.is_isotopologue:  # isotopologues inline
                        if include_isotopolgues and isotopologue_inline:
                            add_match_dict_data()
                    else:
                        add_match_dict_data()  # add monoisotopic peak

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
