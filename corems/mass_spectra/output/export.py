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
    
    def get_mass_spectra_attrs(self, mass_spectra):

        dict_ms_attrs = {}
        dict_ms_attrs['polarity'] = mass_spectra.polarity
        dict_ms_attrs['rt'] =     mass_spectra.rt
        dict_ms_attrs['tic'] =  mass_spectra.tic
        dict_ms_attrs['mobility_scan'] =     mass_spectra.mobility_scan
        dict_ms_attrs['mobility_rt'] =     mass_spectra.mobility_rt
        dict_ms_attrs['Aterm'] =  mass_spectra.Aterm
        dict_ms_attrs['Bterm'] =  mass_spectra.Bterm
        dict_ms_attrs['Cterm'] =  mass_spectra.Cterm
        dict_ms_attrs['baselise_noise'] =  mass_spectra.baselise_noise
        dict_ms_attrs['baselise_noise_std'] =  mass_spectra.baselise_noise_std

    def to_hdf(self, analyzer, instrument_label):
        
        import h5py
        import json
        from numpy import array
        from datetime import datetime, timezone

        with h5py.File(self.output_file.with_suffix('.hdf5'), 'w') as hdf_handle:
            
            if not hdf_handle['date']:

                    timenow = datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z")
                    hdf_handle['date_utc'] = timenow
                    hdf_handle['filename'] = self.mass_spectra.filename.name
                    hdf_handle['analyzer'] = analyzer
                    hdf_handle['instrument_label'] = instrument_label
                    hdf_handle['data_structure'] = 'mass_spectra'

            for mass_spectrum in self.mass_spectra:

                list_results = self.list_dict_to_list(mass_spectrum)
                
                dict_ms_attrs = self.none_to_nan_and_json(self.get_mass_spec_attrs(mass_spectrum))
                
                setting_dicts = settings_parsers.get_dict_data()

                columns_labels = json.dumps(self.columns_label + self.get_all_used_atoms_in_ordem(mass_spectrum))

                scan_group = hdf_handle.create_group(str(mass_spectrum.scan_number))
    
                if mass_spectrum.abundance_profile:
                    
                    mz_abun_array = array(mass_spectrum.abundance_profile, mass_spectrum.mz_profile)
                    
                    raw_ms_dataset = scan_group.create_dataset('raw_ms', data=mz_abun_array, dtype="f8")
                
                else:
                    
                    #create empy dataset for missing raw data
                    raw_ms_dataset = scan_group.create_dataset('raw_ms', dtype="f8")

                raw_ms_dataset.attrs['MassSpecAttrs'] = dict_ms_attrs
                raw_ms_dataset.attrs['TransientSetting'] = self.none_to_nan_and_json(setting_dicts.get('TransientSetting'))

                # if there is not processed data len = 0, otherwise len() will return next index
                index_processed_data = str(len(scan_group.keys()))

                processed_dset = scan_group.create_dataset(index_processed_data, data=list_results,  dtype="S")
            
                processed_dset['date_utc'] = timenow

                processed_dset.attrs['ColumnsLabels'] = columns_labels
                
                processed_dset.attrs['MoleculaSearchSetting'] = self.none_to_nan_and_json(setting_dicts.get('MoleculaSearch'))
                
                processed_dset.attrs['MassSpecPeakSetting'] = self.none_to_nan_and_json(setting_dicts.get('MassSpecPeak'))

                processed_dset.attrs['MassSpectrumSetting'] = self.none_to_nan_and_json(setting_dicts.get('MassSpectrum'))

        
   