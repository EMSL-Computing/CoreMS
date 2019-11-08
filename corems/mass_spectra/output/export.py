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
        output_type: str
            'excel', 'csv', 'hdf5' or 'pandas'
        '''
        
        self.output_file = Path(out_file_path).name

        self.dir_loc = Path(out_file_path + ".corems")
        
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
                
                list_results = self.list_dict_to_list(mass_spectrum)
                
                dset = hdf_handle.create_dataset(str(mass_spectrum.scan_number), data=list_results)
                
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
    
   