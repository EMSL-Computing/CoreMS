__author__ = "Yuri E. Corilo"
__date__ = "Nov 06, 2019"


from pathlib import Path

from numpy import  NaN, concatenate
from pandas import DataFrame

from corems.mass_spectrum.output.export import MassSpecExport
from corems.encapsulation.constant import Atoms
from corems.encapsulation.settings.io import settings_parsers
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq

class MassSpectraExport(MassSpecExport):
    '''
    
    '''
    def __init__(self, out_file_path, mass_spectra, output_type='excel'):
        '''
        output_type: str
            'excel', 'csv', 'hdf5' or 'pandas'
        '''
        
        self.output_file = Path(out_file_path)

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
    
    def get_mass_spectra_attrs(self, mass_spectra):

        dict_ms_attrs = {}
        dict_ms_attrs['analyzer'] = self.mass_spectra.analyzer
        dict_ms_attrs['instrument_label'] = self.mass_spectra.instrument_label
        dict_ms_attrs['sample_name'] = self.mass_spectra.sample_name
        
        import json  
        return json.dumps(dict_ms_attrs)

    def to_hdf(self):
        
        import h5py
        import json
        from numpy import array
        from datetime import datetime, timezone

        with h5py.File(self.output_file.with_suffix('.hdf5'), 'a') as hdf_handle:
            
            if not hdf_handle.attrs.get('date_utc'):

                    timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
                    hdf_handle.attrs['date_utc'] = timenow
                    hdf_handle.attrs['filename'] = self.mass_spectra.file_location.name
                    hdf_handle.attrs['data_structure'] = 'mass_spectra'
                    hdf_handle.attrs['analyzer'] = self.mass_spectra.analyzer
                    hdf_handle.attrs['instrument_label'] = self.mass_spectra.instrument_label
                    hdf_handle.attrs['sample_name'] = self.mass_spectra.sample_name

            for mass_spectrum in self.mass_spectra:

                list_results = self.list_dict_to_list(mass_spectrum)
                
                dict_ms_attrs = self.get_mass_spec_attrs(mass_spectrum)
                
                setting_dicts = settings_parsers.get_dict_data_ms(mass_spectrum)

                columns_labels = json.dumps(self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum))

                group_key = str(mass_spectrum.scan_number)

                if not group_key in hdf_handle.keys():
                
                    scan_group = hdf_handle.create_group(str(mass_spectrum.scan_number))

                    if list(mass_spectrum.abundance_profile):
                    
                       mz_abun_array = concatenate((mass_spectrum.abundance_profile, mass_spectrum.mz_exp_profile), axis=0)
                        
                       
                       raw_ms_dataset = scan_group.create_dataset('raw_ms', data=mz_abun_array, dtype="f8")

                    else:
                        #create empy dataset for missing raw data
                        raw_ms_dataset = scan_group.create_dataset('raw_ms', dtype="f8")

                    raw_ms_dataset.attrs['MassSpecAttrs'] = self.to_json(dict_ms_attrs)
                    
                    if isinstance(mass_spectrum, MassSpecfromFreq):
                        raw_ms_dataset.attrs['TransientSetting'] = self.to_json(setting_dicts.get('TransientSetting'))

                else:
                    
                    scan_group = hdf_handle.get(group_key)
        
                    # if there is not processed data len = 0, otherwise len() will return next index
                index_processed_data = str(len(scan_group.keys()))
                
                timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
                
                processed_dset = scan_group.create_dataset(index_processed_data, data=list_results)
            
                processed_dset.attrs['date_utc'] = timenow

                processed_dset.attrs['ColumnsLabels'] = columns_labels
                
                processed_dset.attrs['MoleculaSearchSetting'] = self.to_json(setting_dicts.get('MoleculaSearch'))
                
                processed_dset.attrs['MassSpecPeakSetting'] = self.to_json(setting_dicts.get('MassSpecPeak'))

                processed_dset.attrs['MassSpectrumSetting'] = self.to_json(setting_dicts.get('MassSpectrum'))