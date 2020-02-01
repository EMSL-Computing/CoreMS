__author__ = "Yuri E. Corilo"
__date__ = "Nov 11, 2019"

from threading import Thread
from pathlib import Path

from numpy import string_, array, NaN, empty
from pandas import DataFrame

from corems.encapsulation.constant import Atoms
from corems.encapsulation.constant import Labels
from corems.encapsulation.settings.io import settings_parsers
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq

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

        self.output_file = Path(out_file_path)

        # 'excel', 'csv' or 'pandas'
        self.output_type = output_type

        self.mass_spectrum = mass_spectrum

        # collect all assigned atoms and order them accordingly to the Atoms.atoms_order list
        self.atoms_order_list = self.get_all_used_atoms_in_ordem(self.mass_spectrum)

        self._init_columns()

    def _init_columns(self):

        # column labels in order
        self.columns_label = ['Index',
                        'm/z',
                        'Calibrated m/z',
                        'Calculated m/z',
                        'Peak Height',
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
            self.to_hdf()
        else:
            raise ValueError(
                "Unkown output type: %s; it can be 'excel', 'csv' or 'pandas'" % self.output_type)

    def run(self):
        ''' run is called when the Tread start
            call exportMS.start() '''
        self.save()

    def get_pandas_df(self):

        columns = self.columns_label + self.get_all_used_atoms_in_ordem(self.mass_spectrum)
        dict_data_list = self.get_list_dict_data(self.mass_spectrum)
        df = DataFrame(dict_data_list, columns=columns)
        df.name = self.output_file
        return df

    def write_settings(self, output_path, mass_spectrum):
        
        import json
        
        dict_setting = settings_parsers.get_dict_data_ms(mass_spectrum)

        dict_setting['MassSpecAttrs'] = self.get_mass_spec_attrs(mass_spectrum)
        dict_setting['analyzer'] = mass_spectrum.analyzer
        dict_setting['instrument_label'] = mass_spectrum.instrument_label
        dict_setting['sample_name'] = mass_spectrum.sample_name

        with open(output_path.with_suffix('.json'), 'w', encoding='utf8', ) as outfile:

            output = json.dumps(dict_setting, sort_keys=True, indent=4, separators=(',', ': '))
            outfile.write(output)
    
    def to_pandas(self):
        
        columns = self.columns_label + self.get_all_used_atoms_in_ordem(self.mass_spectrum)

        dict_data_list = self.get_list_dict_data(self.mass_spectrum)

        df = DataFrame(dict_data_list, columns=columns)

        df.to_pickle(self.output_file.with_suffix('.pkl'))
        
        self.write_settings(self.output_file, self.mass_spectrum)
    
    def to_excel(self):

        columns = self.columns_label + self.get_all_used_atoms_in_ordem(self.mass_spectrum)

        dict_data_list = self.get_list_dict_data(self.mass_spectrum)

        df = DataFrame(dict_data_list, columns=columns)

        df.to_excel(self.output_file.with_suffix('.xlsx'))

        self.write_settings(self.output_file, self.mass_spectrum)

    def to_csv(self):
        
        columns = self.columns_label + self.get_all_used_atoms_in_ordem(self.mass_spectrum)

        dict_data_list = self.get_list_dict_data(self.mass_spectrum)

        import csv
        try:
            with open(self.output_file.with_suffix('.csv'), 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)
            
            self.write_settings(self.output_file, self.mass_spectrum)
        
        except IOError as ioerror:
            print(ioerror)
    
   
    def to_hdf(self):
        
        import h5py
        import json
        from datetime import datetime, timezone

        with h5py.File(self.output_file.with_suffix('.hdf5'), 'a') as hdf_handle:
            
            list_results = self.list_dict_to_list(self.mass_spectrum)
            
            dict_ms_attrs = self.get_mass_spec_attrs(self.mass_spectrum)
            
            setting_dicts = settings_parsers.get_dict_data_ms(self.mass_spectrum)

            columns_labels = json.dumps(self.columns_label + self.get_all_used_atoms_in_ordem(self.mass_spectrum))

            if not hdf_handle.attrs.get('date_utc'):

                timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
                hdf_handle.attrs['date_utc'] = timenow
                hdf_handle.attrs['file_name'] = self.mass_spectrum.filename.name
                hdf_handle.attrs['data_structure'] = 'mass_spectrum'
                hdf_handle.attrs['analyzer'] = self.mass_spectrum.analyzer
                hdf_handle.attrs['instrument_label'] = self.mass_spectrum.instrument_label
                hdf_handle.attrs['sample_name'] = self.mass_spectrum.sample_name

            group_key = str(self.mass_spectrum.scan_number)
            
            if not group_key in hdf_handle.keys():
                
                scan_group = hdf_handle.create_group(str(self.mass_spectrum.scan_number))

                if list(self.mass_spectrum.abundance_profile):
                
                    mz_abun_array = empty(shape=(2,len(self.mass_spectrum.abundance_profile)))
                    
                    mz_abun_array[0] = self.mass_spectrum.abundance_profile
                    mz_abun_array[1] = self.mass_spectrum.mz_exp_profile
                    
                    raw_ms_dataset = scan_group.create_dataset('raw_ms', data=mz_abun_array, dtype="f8")

                else:
                    #create empy dataset for missing raw data
                    raw_ms_dataset = scan_group.create_dataset('raw_ms', dtype="f8")

                raw_ms_dataset.attrs['MassSpecAttrs'] = self.to_json(dict_ms_attrs)
                
                if isinstance(self.mass_spectrum, MassSpecfromFreq):
                    raw_ms_dataset.attrs['TransientSetting'] = self.to_json(setting_dicts.get('TransientSetting'))

            else:
                
                scan_group = hdf_handle.get(group_key)
    
            # if there is not processed data len = 0, otherwise len() will return next index
            index_processed_data = str(len(scan_group.keys()))
            
            print('index_processed_data', index_processed_data)
            
            timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
            
            processed_dset = scan_group.create_dataset(index_processed_data, data=list_results)
        
            processed_dset.attrs['date_utc'] = timenow

            processed_dset.attrs['ColumnsLabels'] = columns_labels
            
            processed_dset.attrs['MoleculaSearchSetting'] = self.to_json(setting_dicts.get('MoleculaSearch'))
            
            processed_dset.attrs['MassSpecPeakSetting'] = self.to_json(setting_dicts.get('MassSpecPeak'))

            processed_dset.attrs['MassSpectrumSetting'] = self.to_json(setting_dicts.get('MassSpectrum'))


    @staticmethod     
    def to_json(dict_data, nan=False):
        
        import json
        #for key, values in dict_data.items():
        #    if not values: dict_data[key] = NaN
        
        #output = json.dumps(dict_data, sort_keys=True, indent=4, separators=(',', ': '))
        return json.dumps(dict_data)

    def get_mass_spec_attrs(self, mass_spectrum):

        dict_ms_attrs = {}
        dict_ms_attrs['polarity'] = mass_spectrum.polarity
        dict_ms_attrs['rt'] =     mass_spectrum.rt
        dict_ms_attrs['tic'] =  mass_spectrum.tic
        dict_ms_attrs['mobility_scan'] =     mass_spectrum.mobility_scan
        dict_ms_attrs['mobility_rt'] =     mass_spectrum.mobility_rt
        dict_ms_attrs['Aterm'] =  mass_spectrum.Aterm
        dict_ms_attrs['Bterm'] =  mass_spectrum.Bterm
        dict_ms_attrs['Cterm'] =  mass_spectrum.Cterm
        dict_ms_attrs['baselise_noise'] =  mass_spectrum.baselise_noise
        dict_ms_attrs['baselise_noise_std'] =  mass_spectrum.baselise_noise_std

        return dict_ms_attrs


    def get_all_used_atoms_in_ordem(self, mass_spectrum):

        atoms_in_order = Atoms.atoms_order
        all_used_atoms = set()
        if mass_spectrum:
            for ms_peak in mass_spectrum:
                if ms_peak:
                    for m_formula in ms_peak:
                        for atom in m_formula.atoms:
                            all_used_atoms.add(atom)

        def sort_method(atom):
            return [atoms_in_order.index(atom)]

        return sorted(all_used_atoms, key=sort_method)

    def list_dict_to_list(self, mass_spectrum):
        
        column_labels = self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum)

        dict_list = self.get_list_dict_data(mass_spectrum)
        
        all_lines = []
        for dict_res in dict_list:
            
            result_line = [NaN] * len(column_labels)
            
            for label, value in dict_res.items():
                
                label_index = column_labels.index(label)
                result_line[label_index] =  value   
            
            all_lines.append(result_line)
        return  all_lines       

    
    def get_list_dict_data(self, mass_spectrum, include_no_match=True, include_isotopolgues=True,
                           isotopologue_inline=False, no_match_inline=False):

        dict_data_list = []

        def add_no_match_dict_data():

            dict_result = {'Index': index,
                           'm/z':  ms_peak._mz_exp,
                           'Calibrated m/z': ms_peak.mz_exp,
                           'Peak Height': ms_peak.abundance,
                           'Resolving Power': ms_peak.resolving_power,
                           'S/N':  ms_peak.signal_to_noise,
                           'Ion Charge': ms_peak.ion_charge,
                           'Heteroatom Class' : Labels.unassigned
                           }

            dict_data_list.append(dict_result)

        def add_match_dict_data():

            formula_dict = m_formula.to_dict
            
            dict_result = {'Index': index,
                           'm/z':  ms_peak._mz_exp,
                           'Calibrated m/z': ms_peak.mz_exp,
                           'Calculated m/z': m_formula.mz_theor,
                           'Peak Height': ms_peak.abundance,
                           'Resolving Power': ms_peak.resolving_power,
                           'S/N':  ms_peak.signal_to_noise,
                           'Ion Charge': ms_peak.ion_charge,
                           'Mass Error (ppm)': m_formula.mz_error,
                           'DBE':  m_formula.dbe,
                           'Heteroatom Class': str(m_formula.class_label.encode('utf-8')),
                           'H/C':  m_formula.H_C,
                           'O/C':  m_formula.O_C,
                           'Ion Type': str(m_formula.ion_type.lower().encode('utf-8')),
                           'Is Isotopologue': int(m_formula.is_isotopologue),
                           }

            for atom in self.atoms_order_list:
                if atom in formula_dict.keys():
                    dict_result[atom] = formula_dict.get(atom)

            dict_data_list.append(dict_result)

        for index, ms_peak in enumerate(mass_spectrum):

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
            for index, ms_peak in enumerate(mass_spectrum):
                for m_formula in ms_peak:
                    if m_formula.is_isotopologue:
                        add_match_dict_data()

        if include_no_match and not no_match_inline:
            for index, ms_peak in enumerate(mass_spectrum):
                if not ms_peak:
                    add_no_match_dict_data()
        
        return dict_data_list
