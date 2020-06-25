__author__ = "Yuri E. Corilo"
__date__ = "Nov 06, 2019"


from pathlib import Path

from numpy import  NaN, concatenate
from pandas import DataFrame, ExcelWriter, read_excel
from openpyxl import load_workbook

from corems.mass_spectrum.output.export import HighResMassSpecExport
from corems.encapsulation.constant import Atoms
from corems.encapsulation.output import parameter_to_dict
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq

class LowResGCMSExport():
    
    def __init__(self, out_file_path, gcms):
        '''
        output_type: str
            'excel', 'csv', 'hdf5' or 'pandas'
        '''
        
        self.output_file = Path(out_file_path)

        self.gcms = gcms

        self._init_columns()
        
    def _init_columns(self):

        columns =  ['Sample name', 'Retention Time', 'Retention Time Ref', 'Peak Height',
                'Peak Area', 'Retention index', 'Retention index Ref','Retention Index Score',
                'Similarity Score',
                'Spectral Similarity Score',
                'Compound Name']
        
        if self.gcms.molecular_search_settings.exploratory_mode:
                
                columns.extend(['Weighted Cosine Correlation',
                                'Cosine Correlation',
                                'Stein Scott Similarity',
                                'Pearson Correlation',
                                'Spearman Correlation',
                                'Kendall Tau Correlation', 
                                'Euclidean Distance',
                                'Manhattan Distance',
                                'Jaccard Similarity',
                                'DWT Correlation',
                                'DFT Correlation' ])
                                                
        
        return columns        
    
    def get_pandas_df(self, highest_score=True):

        columns = self._init_columns()
        
        dict_data_list = self.get_list_dict_data(self.gcms, highest_score=highest_score)
        
        df = DataFrame(dict_data_list, columns=columns)
        
        df.name = self.gcms.sample_name

        return df

    def get_json(self, nan=False, highest_score=True):
        
        import json

        dict_data_list = self.get_list_dict_data(self.gcms, highest_score=highest_score)
        
        return json.dumps(dict_data_list)

    def to_pandas(self, highest_score=True):
        
        columns = self._init_columns() 
        
        dict_data_list = self.get_list_dict_data(self.gcms, highest_score=highest_score)

        df = DataFrame(dict_data_list, columns=columns)

        df.to_pickle(self.output_file.with_suffix('.pkl'))
        
        self.write_settings(self.output_file, self.gcms)
               
    def to_excel(self, write_mode='a', highest_score=True):

        out_put_path = self.output_file.with_suffix('.xlsx')

        columns = self._init_columns() 
        
        dict_data_list = self.get_list_dict_data(self.gcms, highest_score=highest_score)

        df = DataFrame(dict_data_list, columns=columns)

        if write_mode == 'a' and out_put_path.exists():
            
            writer = ExcelWriter(out_put_path, engine='openpyxl')
            # try to open an existing workbook
            writer.book = load_workbook(out_put_path)
            # copy existing sheets
            writer.sheets = dict((ws.title, ws) for ws in writer.book.worksheets)
            # read existing file
            reader = read_excel(out_put_path)
            # write out the new sheet
            df.to_excel(writer,index=False,header=False,startrow=len(reader)+1)

            writer.close()
        else:
        
            df.to_excel(self.output_file.with_suffix('.xlsx'), index=False, engine='openpyxl')

        self.write_settings(self.output_file, self.gcms)

    def to_csv(self, write_mode='a', highest_score=True) :
        
        import csv
        
        columns = self._init_columns() 
        
        dict_data_list = self.get_list_dict_data(self.gcms, highest_score=highest_score)

        out_put_path = self.output_file.with_suffix('.csv')

        write_header = not out_put_path.exists()
        
        try:
            with open(out_put_path, write_mode, newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                if write_header:
                    writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)
            
            self.write_settings(self.output_file, self.gcms)
        
        except IOError as ioerror:
            print(ioerror)                 
    
    
    def write_settings(self, output_path, gcms):
        
        import json
        
        dict_setting = parameter_to_dict.get_dict_data_gcms(gcms)
        dict_setting['calibration_rt_ri_pairs_ref'] = gcms.ri_pairs_ref
        dict_setting['analyzer'] = gcms.analyzer
        dict_setting['instrument_label'] = gcms.instrument_label
        dict_setting['sample_name'] = [gcms.sample_name]

        with open(output_path.with_suffix('.json'), 'w', encoding='utf8', ) as outfile:

            output = json.dumps(dict_setting, sort_keys=True, indent=4, separators=(',', ': '))
            outfile.write(output)

    def get_list_dict_data(self, gcms, include_no_match=True, no_match_inline=False, highest_score=False) :

        dict_data_list = []

        def add_match_dict_data():

            out_dict = {'Sample name': gcms.sample_name,
                        'Retention Time': gc_peak.rt,
                        'Retention Time Ref': compound_obj.rt,
                        'Peak Height': gc_peak.tic,
                        'Peak Area': gc_peak.area,
                        'Retention index': gc_peak.ri,
                        'Retention index Ref':  compound_obj.ri,
                        'Retention Index Score': compound_obj.ri_score,
                        'Spectral Similarity Score': compound_obj.spectral_similarity_score,
                        'Similarity Score': compound_obj.similarity_score,
                        'Compound Name' : compound_obj.name
            }    

            if self.gcms.molecular_search_settings.exploratory_mode:
                out_dict.update({
                    'Weighted Cosine Correlation': compound_obj.spectral_similarity_scores.get("weighted_cosine_correlation"), 
                    'Cosine Correlation': compound_obj.spectral_similarity_scores.get("cosine_correlation"), 
                    'Stein Scott Similarity': compound_obj.spectral_similarity_scores.get("stein_scott_similarity"), 
                    'Pearson Correlation': compound_obj.spectral_similarity_scores.get("pearson_correlation"), 
                    'Spearman Correlation': compound_obj.spectral_similarity_scores.get("spearman_correlation"), 
                    'Kendall Tau Correlation': compound_obj.spectral_similarity_scores.get("kendall_tau_correlation"), 
                    'Euclidean Distance': compound_obj.spectral_similarity_scores.get("euclidean_distance"), 
                    'Manhattan Distance': compound_obj.spectral_similarity_scores.get("manhattan_distance"),
                    'Jaccard Similarity': compound_obj.spectral_similarity_scores.get("jaccard_distance"),
                    'DFT Correlation': compound_obj.spectral_similarity_scores.get("dft_correlation"),
                    'DWT Correlation': compound_obj.spectral_similarity_scores.get("dwt_correlation"),
                })
        
            dict_data_list.append(out_dict)
        
        def add_no_match_dict_data():

            dict_data_list.append( {'Sample name': gcms.sample_name,
                           'Retention Time': gc_peak.rt,
                           'Peak Height': gc_peak.tic,
                           'Peak Area': gc_peak.area,
                           'Retention index': gc_peak.ri,
                           } )

           
        for gc_peak in gcms:

            # check if there is a compound candidate 
            if gc_peak:
                
                if highest_score:
                    compound_obj = gc_peak.highest_score_compound
                    add_match_dict_data()
                    
                else:
                    for compound_obj in gc_peak:
                        add_match_dict_data()  # add monoisotopic peak

            else:
                # include not_match
                if include_no_match and no_match_inline:
                    add_no_match_dict_data()

        if include_no_match and not no_match_inline:
            for gc_peak in gcms:
                if not gc_peak:
                    add_no_match_dict_data()
        
        return dict_data_list        

    
class HighResMassSpectraExport(HighResMassSpecExport):
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

            columns = self.columns_label + self.get_all_used_atoms_in_order(mass_spectrum)
            
            dict_data_list = self.get_list_dict_data(mass_spectrum)
            
            df = DataFrame(dict_data_list, columns=columns)
            
            scan_number = mass_spectrum.scan_number
            
            df.name = str(self.output_file) + '_' + str(scan_number)
            
            list_df.append(df)
        
        return list_df

    def to_pandas(self):
        
        for mass_spectrum in self.mass_spectra:

            columns = self.columns_label + self.get_all_used_atoms_in_order(mass_spectrum)

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number

            out_filename = Path("%s_scan%s%s" % (self.output_file, str(scan_number), '.pkl'))
            
            df.to_pickle(self.dir_loc / out_filename)

            self.write_settings(self.dir_loc / out_filename.with_suffix(''), mass_spectrum)

    def to_excel(self):

        for mass_spectrum in self.mass_spectra:
            
            columns = self.columns_label + self.get_all_used_atoms_in_order(mass_spectrum)

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number
            
            out_filename = Path("%s_scan%s%s" % (self.output_file, str(scan_number), '.xlsx'))
            
            df.to_excel(self.dir_loc / out_filename)

            self.write_settings(self.dir_loc / out_filename.with_suffix(''), mass_spectrum)
            

    def to_csv(self):
        
        import csv

        for mass_spectrum in self.mass_spectra:
            
            columns = self.columns_label + self.get_all_used_atoms_in_order(mass_spectrum)

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

                list_results = self.list_dict_to_list(mass_spectrum, is_hdf5=True)
                
                dict_ms_attrs = self.get_mass_spec_attrs(mass_spectrum)
                
                setting_dicts = parameter_to_dict.get_dict_data_ms(mass_spectrum)

                columns_labels = json.dumps(self.columns_label + self.get_all_used_atoms_in_order(mass_spectrum))

                group_key = str(mass_spectrum.scan_number)

                if not group_key in hdf_handle.keys():
                
                    scan_group = hdf_handle.create_group(str(mass_spectrum.scan_number))

                    if list(mass_spectrum.abundance_profile):
                    
                       mz_abun_array = concatenate((mass_spectrum.abundance_profile, mass_spectrum.mz_exp_profile), axis=0)
                        
                       
                       raw_ms_dataset = scan_group.create_dataset('raw_ms', data=mz_abun_array, dtype="f8")

                    else:
                        #create empy dataset for missing raw data
                        raw_ms_dataset = scan_group.create_dataset('raw_ms', dtype="f8")

                    raw_ms_dataset.attrs['MassSpecAttrs'] = json.dumps(dict_ms_attrs)
                    
                    if isinstance(mass_spectrum, MassSpecfromFreq):
                        raw_ms_dataset.attrs['TransientSetting'] = json.dumps(setting_dicts.get('TransientSetting'))

                else:
                    
                    scan_group = hdf_handle.get(group_key)
        
                    # if there is not processed data len = 0, otherwise len() will return next index
                index_processed_data = str(len(scan_group.keys()))
                
                timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
                
                processed_dset = scan_group.create_dataset(index_processed_data, data=list_results)
            
                processed_dset.attrs['date_utc'] = timenow

                processed_dset.attrs['ColumnsLabels'] = columns_labels
                
                processed_dset.attrs['MoleculaSearchSetting'] = json.dumps(setting_dicts.get('MoleculaSearch'))
                
                processed_dset.attrs['MassSpecPeakSetting'] = json.dumps(setting_dicts.get('MassSpecPeak'))

                processed_dset.attrs['MassSpectrumSetting'] = json.dumps(setting_dicts.get('MassSpectrum'))