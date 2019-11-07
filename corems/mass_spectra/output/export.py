__author__ = "Yuri E. Corilo"
__date__ = "Nov 06, 2019"


from pathlib import Path

from numpy import   NaN
from pandas import DataFrame

from corems.mass_spectrum.output.export import MassSpecExport
from corems.encapsulation.constant import Atoms
from corems.encapsulation.settings.io import settings_parsers



class MassSpectraExport(MassSpecExport):
    '''
    TODO: add MSPeak indexes: done
    '''
    def __init__(self, out_file_path, mass_spectra, output_type='excel'):
        '''
        output_type:str
            'excel', 'csv', 'hdf5' or 'pandas'
        '''
        
        self.output_file = Path(out_file_path).name

        self.dir_loc = Path(out_file_path / ".corems")
        
        self.dir_loc.mkdir(exist_ok=True)

        # 'excel', 'csv' or 'pandas'
        self._output_type = output_type

        self.mass_spectra = mass_spectra

        self._init_columns()

        
    def get_pandas_df(self):

        list_df = []
        
        for mass_spectrum in self.mass_spectra:

            columns = self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum)
            
            dict_data_list = self.get_list_dict_data(mass_spectrum)
            
            df = DataFrame(dict_data_list, columns=columns)
            
            scan_number = mass_spectrum.scan_number
            
            df.name = str(self.output_file) + '_' + str(scan_number)
            
            list_df.append(df)
        
        return list_df

    def write_settings(self, output_path, mass_spectrum):
        
        import json
        
        dict_setting = settings_parsers.get_dict_data()
        del dict_setting['DataInput']
        dict_setting['MassSpecAttrs'] = self.none_to_nan_and_json(self.get_mass_spec_attrs(mass_spectrum))
        
        with open(output_path.with_suffix('.json'), 'w', encoding='utf8', ) as outfile:

            output = json.dumps(dict_setting, sort_keys=True, indent=4, separators=(',', ': '))
            outfile.write(output)

    def to_pandas(self):
        
        for mass_spectrum in self.mass_spectra:

            columns = self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum)

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number

            out_filename = Path("%s_scan%s%s" % (self.output_file, str(scan_number), '.pkl'))
            
            df.to_pickle(self.dir_loc / out_filename)

            self.write_settings(self.dir_loc / out_filename.with_suffix(''), mass_spectrum)

    def to_excel(self):

        for mass_spectrum in self.mass_spectra:
            
            columns = self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum)

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number
            
            out_filename = Path("%s_scan%s%s" % (self.output_file, str(scan_number), '.xlsx'))
            
            df.to_excel(self.dir_loc / out_filename)

            self.write_settings(self.dir_loc / out_filename.with_suffix(''), mass_spectrum)
            

    def to_csv(self):
        
        import csv

        for mass_spectrum in self.mass_spectra:
            
            columns = self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum)

            scan_number = mass_spectrum.scan_number

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            out_filename = Path("%s_scan%s%s" % (self.output_file, str(scan_number), '.csv'))
            
            
            with open(self.dir_loc / out_filename, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)

            self.write_settings(self.dir_loc / out_filename.with_suffix(''), mass_spectrum)                    
    
    def to_hdf(self, ):
        
        import h5py
        import json

        with h5py.File(self.output_file + '.hdf5', 'w') as hdf_handle:
        
            setting_dicts = settings_parsers.get_dict_data()
            
            self.write_metadata_hdf(setting_dicts, hdf_handle)

            for mass_spectrum in self.mass_spectra:
                
                list_results = self.list_dict_to_list()
                
                dset = hdf_handle.create_dataset(str(self.mass_spectrum.scan_number), data=list_results)
                
                columns_labels = json.dumps(self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum))

                dict_ms_attrs = self.none_to_nan_and_json(self.get_mass_spec_attrs(mass_spectrum))
                
                moleculaSearchSetting = self.none_to_nan_and_json(setting_dicts.get('MoleculaSearch'))
                
                massSpecPeakSetting  = self.none_to_nan_and_json(setting_dicts.get('MassSpecPeak'))

                dset.attrs['ColumnsLabels'] = columns_labels
                
                dset.attrs['MassSpecAttrs'] = dict_ms_attrs
                
                #needs to store the molecular search in a dict by scan number in the search class then it can be accessed here
                dset.attrs['MoleculaSearchSetting'] = moleculaSearchSetting
               
                dset.attrs['MassSpecPeakSetting'] = massSpecPeakSetting
            
            #dset.attrs.update(dict_ms_attrs)
    
    @staticmethod     
    def none_to_nan_and_json( dict_data):
        import json
        for key, values in dict_data.items():
            if not values: dict_data[key] = NaN
        
        return json.dumps(dict_data)

    @staticmethod
    def get_mass_spec_attrs( mass_spectrum):

        dict_ms_attrs = {}
        dict_ms_attrs['polarity'] = mass_spectrum.polarity
        dict_ms_attrs['rt'] = mass_spectrum.rt
        dict_ms_attrs['mobility_scan'] = mass_spectrum.mobility_scan
        dict_ms_attrs['mobility_rt'] = mass_spectrum.mobility_rt
        dict_ms_attrs['Aterm'] = mass_spectrum.Aterm
        dict_ms_attrs['Bterm'] = mass_spectrum.Bterm
        dict_ms_attrs['Cterm'] = mass_spectrum.Cterm
        dict_ms_attrs['baselise_noise'] = mass_spectrum.baselise_noise
        dict_ms_attrs['baselise_noise_std'] = mass_spectrum.baselise_noise_std

        return dict_ms_attrs


    def get_all_used_atoms_in_ordem(self, mass_spectrum):

        def sort_method(atom):
            return [Atoms.atoms_order.index(atom)]
            
        all_used_atoms = set()
        for ms_peak in mass_spectrum:
            for m_formula in ms_peak:
                if ms_peak:
                    for atom in m_formula.atoms:
                        all_used_atoms.add(atom)

        return sorted(all_used_atoms, key=sort_method)

    def list_of_dict_data_to_list_of_table_data(self, mass_spectrum):
        
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

            for atom in self.get_all_used_atoms_in_ordem(mass_spectrum):
                if atom in formula_dict.keys():
                    dict_result[atom] = formula_dict.get(atom)

            dict_data_list.append(dict_result)

        for index, ms_peak in enumerate(mass_spectrum.sort_by_mz()):

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
            for index, ms_peak in enumerate(mass_spectrum.sort_by_mz()):
                for m_formula in ms_peak:
                    if m_formula.is_isotopologue:
                        add_match_dict_data()

        if include_no_match and not no_match_inline:
            for index, ms_peak in enumerate(mass_spectrum.sort_by_mz()):
                if not ms_peak:
                    add_no_match_dict_data()
        
        return dict_data_list
