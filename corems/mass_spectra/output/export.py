__author__ = "Yuri E. Corilo"
__date__ = "Nov 06, 2019"


import csv
from pathlib import Path

from numpy import  NaN, concatenate
from openpyxl import load_workbook
from pandas import DataFrame, ExcelWriter, read_excel

from corems.mass_spectrum.output.export import HighResMassSpecExport
from corems.molecular_id.calc.SpectralSimilarity import methods_name
from corems.encapsulation.constant import Atoms
from corems.encapsulation.output import parameter_to_dict
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq
from corems import __version__, corems_md5
import uuid

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

        columns =  ['Sample name', 'Peak Index',  'Retention Time', 'Retention Time Ref', 'Peak Height',
                'Peak Area', 'Retention index', 'Retention index Ref','Retention Index Score',
                'Similarity Score',
                'Spectral Similarity Score',
                'Compound Name']
        
        if self.gcms.molecular_search_settings.exploratory_mode:
                
                columns.extend(['Weighted Cosine Correlation',
                                'Cosine Correlation',
                                'Stein Scott Similarity',
                                'Stein Scott Similarity Nist',
                                'Pearson Correlation',
                                'Spearman Correlation',
                                'Kendall Tau Correlation', 
                                'Euclidean Distance',
                                'Manhattan Distance',
                                'Jaccard Distance',
                                'DWT Correlation',
                                'DFT Correlation' ])

                columns.extend(list(methods_name.values()))
                                                
        
        return columns        
    
    def get_pandas_df(self, id_label="corems:"):

        columns = self._init_columns()

        dict_data_list = self.get_list_dict_data(self.gcms)

        df = DataFrame(dict_data_list, columns=columns)

        df.name = self.gcms.sample_name

        return df

    def get_json(self, nan=False, id_label="corems:"):

        import json

        dict_data_list = self.get_list_dict_data(self.gcms)

        return json.dumps(dict_data_list, sort_keys=False, indent=4, separators=(',', ': '))

    def to_pandas(self, id_label="corems:"):

        columns = self._init_columns() 

        dict_data_list = self.get_list_dict_data(self.gcms)

        df = DataFrame(dict_data_list, columns=columns)

        df.to_pickle(self.output_file.with_suffix('.pkl'))

        self.write_settings(self.output_file, self.gcms, id_label="corems:")

    def to_excel(self, write_mode='w', id_label="corems:"):

        out_put_path = self.output_file.with_suffix('.xlsx')

        columns = self._init_columns()

        dict_data_list = self.get_list_dict_data(self.gcms)

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
            df.to_excel(writer, index=False, header=False, startrow=len(reader) + 1)

            writer.close()
        else:

            df.to_excel(self.output_file.with_suffix('.xlsx'), index=False, engine='openpyxl')

        self.write_settings(self.output_file, self.gcms, id_label="corems:")

    def to_csv(self, separate_output=False, id_label="corems:"):

        if separate_output:
            # set write mode to write
            # this mode will overwrite the file without warning
            write_mode = 'w'
        else:
            # set write mode to append
            write_mode = 'a'
        
        columns = self._init_columns() 
        
        dict_data_list = self.get_list_dict_data(self.gcms)

        out_put_path = self.output_file.with_suffix('.csv')

        write_header = not out_put_path.exists()
        
        try:
            with open(out_put_path, write_mode, newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                if write_header:
                    writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)
            
            self.write_settings(self.output_file, self.gcms, id_label=id_label)
        
        except IOError as ioerror:
            print(ioerror)                 
    
    def to_hdf(self, id_label="corems:"):
        
        # save sample at a time
        def add_compound(gc_peak, compound_obj):

            compound_group = peak_group.create_group(compound_obj.name.replace('/', ''))
            compound_group.attrs["retention_time"] = compound_obj.rt
            compound_group.attrs["retention_index"] = compound_obj.ri
            compound_group.attrs["retention_index_score"] = compound_obj.ri_score
            compound_group.attrs["spectral_similarity_score"] = compound_obj.spectral_similarity_score
            compound_group.attrs["similarity_score"] = compound_obj.similarity_score
            
            compond_mz = compound_group.create_dataset('mz', data=np.array(compound_obj.mz), dtype="f8")  
            compond_abundance = compound_group.create_dataset('abundance', data=np.array(compound_obj.abundance), dtype="f8")

            if self.gcms.molecular_search_settings.exploratory_mode:
                
                compound_group.attrs['Spectral Similarities'] =  json.dumps(compound_obj.spectral_similarity_scores, sort_keys=False, indent=4, separators=(',', ': '))
                
        import h5py
        import json
        import numpy as np
        from datetime import datetime, timezone

        output_path = self.output_file.with_suffix('.hdf5')

        with h5py.File(output_path, 'w') as hdf_handle:

            timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
            hdf_handle.attrs['time_stamp'] = timenow
            hdf_handle.attrs['data_structure'] = 'gcms'
            hdf_handle.attrs['analyzer'] = self.gcms.analyzer
            hdf_handle.attrs['instrument_label'] = self.gcms.instrument_label

            hdf_handle.attrs['sample_id'] = "self.gcms.id"
            hdf_handle.attrs['sample_name'] = self.gcms.sample_name
            hdf_handle.attrs['input_data'] = str(self.gcms.file_location)
            hdf_handle.attrs['output_data'] = str(output_path)
            hdf_handle.attrs['output_data_id'] = id_label + uuid.uuid4().hex
            hdf_handle.attrs['corems_version'] = __version__

            hdf_handle.attrs["Stats"] = json.dumps(self.get_data_stats(self.gcms), sort_keys=False, indent=4, separators=(',', ': '))
            hdf_handle.attrs["Calibration"] = json.dumps(self.get_calibration_stats(self.gcms, id_label), sort_keys=False, indent=4, separators=(',', ': '))
            hdf_handle.attrs["Blank"] = json.dumps(self.get_blank_stats(self.gcms), sort_keys=False, indent=4, separators=(',', ': '))

            corems_dict_setting = parameter_to_dict.get_dict_data_gcms(self.gcms)
            hdf_handle.attrs["CoreMSParameters"] = json.dumps(corems_dict_setting, sort_keys=False, indent=4, separators=(',', ': '))
            
            scans_dataset = hdf_handle.create_dataset('scans', data=np.array(self.gcms.scans_number), dtype="f8")                
            rt_dataset = hdf_handle.create_dataset('rt', data=np.array(self.gcms.retention_time), dtype="f8")                
            tic_dataset = hdf_handle.create_dataset('tic', data=np.array(self.gcms.tic), dtype="f8")            
            processed_tic_dataset = hdf_handle.create_dataset('processed_tic', data=np.array(self.gcms.processed_tic), dtype="f8")

            output_score_method = self.gcms.molecular_search_settings.output_score_method

            for gc_peak in self.gcms:

                # print(gc_peak.rt)
                # print(gc_peak.tic)
                
                # check if there is a compound candidate 
                peak_group = hdf_handle.create_group(str(gc_peak.rt))
                peak_group.attrs["deconvolution"] = int(self.gcms.chromatogram_settings.use_deconvolution)
                
                peak_group.attrs["start_index"] = gc_peak.start_index
                peak_group.attrs["index"] = gc_peak.index
                peak_group.attrs["final_index"] = gc_peak.final_index

                peak_group.attrs["retention_index"] = gc_peak.ri
                peak_group.attrs["retention_time"] = gc_peak.rt
                peak_group.attrs["area"] = gc_peak.area
                
                mz = peak_group.create_dataset('mz', data=np.array(gc_peak.mass_spectrum.mz_exp), dtype="f8")
                abundance = peak_group.create_dataset('abundance', data=np.array(gc_peak.mass_spectrum.abundance), dtype="f8")

                if gc_peak:

                    if output_score_method == 'highest_sim_score':
                        compound_obj = gc_peak.highest_score_compound
                        add_compound(gc_peak, compound_obj)

                    elif output_score_method == 'highest_ss':
                        compound_obj = gc_peak.highest_ss_compound
                        add_compound(gc_peak, compound_obj)

                    else:

                        for compound_obj in gc_peak:
                            add_compound(gc_peak, compound_obj)
                    
    def get_data_stats(self, gcms):
        
        matched_peaks = gcms.matched_peaks
        no_matched_peaks = gcms.no_matched_peaks
        unique_metabolites = gcms.unique_metabolites
        
        peak_matchs_above_0p85 = 0
        unique_peak_match_above_0p85 = 0
        for match_peak in matched_peaks:
            gc_peak_above_85 = 0
            matches_above_85 = list(filter(lambda m: m.similarity_score >= 0.85, match_peak))
            if matches_above_85:
                peak_matchs_above_0p85 +=1
            if len(matches_above_85) == 1:
                unique_peak_match_above_0p85 += 1

        data_stats = {}
        data_stats['average_signal_noise'] = "ni"
        data_stats['chromatogram_dynamic_range'] = gcms.dynamic_range
        data_stats['total_number_peaks'] = len(gcms)
        data_stats['total_peaks_matched'] = len(matched_peaks)
        data_stats['total_peaks_without_matches'] = len(no_matched_peaks)
        data_stats['total_matches_above_similarity_score_0.85'] = peak_matchs_above_0p85
        data_stats['single_matches_above_similarity_score_0.85'] = unique_peak_match_above_0p85
        data_stats['unique_metabolites'] = len(unique_metabolites)

        return data_stats

    def get_calibration_stats(self, gcms, id_label):
        
        calibration_parameters = {}

        calibration_parameters['calibration_rt_ri_pairs_ref'] = gcms.ri_pairs_ref
        calibration_parameters['file_path'] = str(gcms.cal_file_path)
        calibration_parameters['file_id'] = id_label + corems_md5(gcms.cal_file_path)
        calibration_parameters['sample_name'] = str(gcms.cal_file_path.stem)
        calibration_parameters['calibration_method'] = ""

        return calibration_parameters

    def get_blank_stats(self, gcms):
        
        blank_parameters = {}

        blank_parameters['sample_name'] = "ni"
        blank_parameters['blank_id'] = "ni"
        blank_parameters['file_path'] = "ni"
        blank_parameters['file_id'] = "ni"
        blank_parameters['common_features_to_blank'] = "ni"

        return blank_parameters

    def get_parameters_json(self, gcms, id_label, output_path):

        import json

        if isinstance(output_path, str):
            output_path = Path(output_path)

        paramaters_path = output_path.with_suffix('.json')

        if paramaters_path.exists():
            with paramaters_path.open() as current_param:
                output_parameters_dict = json.load(current_param)

        else:
            output_parameters_dict = {}
      
            output_parameters_dict['sample_name'] = []
            output_parameters_dict['sample_id'] = []
            output_parameters_dict['input_data'] = []
            output_parameters_dict['input_data_id'] = []

        blank_parameters = {}
        
        output_parameters_dict['analyzer'] = gcms.analyzer
        output_parameters_dict['instrument_label'] = gcms.instrument_label
        
        output_parameters_dict['sample_name'].append(gcms.sample_name)
        output_parameters_dict['sample_id'].append("sample_id")
        output_parameters_dict['input_data'].append(str(gcms.file_location))
        output_parameters_dict['input_data_id'].append(id_label + corems_md5(gcms.file_location))
        
        output_parameters_dict['output_data'] = str(output_path)
        output_parameters_dict['output_data_id'] = id_label + uuid.uuid4().hex
        output_parameters_dict['corems_version'] = __version__
        
        output_parameters_dict["Stats"] = self.get_data_stats(gcms)
        output_parameters_dict["Calibration"] = self.get_calibration_stats(gcms, id_label)
        output_parameters_dict["Blank"] = self.get_blank_stats(gcms)

        corems_dict_setting = parameter_to_dict.get_dict_data_gcms(gcms)
        output_parameters_dict["CoreMSParameters"] = corems_dict_setting

        output = json.dumps(output_parameters_dict, sort_keys=False, indent=4, separators=(',', ': '))

        return output

    def write_settings(self, output_path, gcms, id_label="corems:"):

        output = self.get_parameters_json(gcms, id_label, output_path)

        with open(output_path.with_suffix('.json'), 'w', encoding='utf8', ) as outfile:

            outfile.write(output)

    def get_list_dict_data(self, gcms, include_no_match=True, no_match_inline=False) :

        output_score_method = gcms.molecular_search_settings.output_score_method

        dict_data_list = []

        def add_match_dict_data():

            out_dict = {'Sample name': gcms.sample_name,
                        'Peak Index': gcpeak_index, 
                        'Retention Time': gc_peak.rt,
                        'Retention Time Ref': compound_obj.rt,
                        'Peak Height': gc_peak.tic,
                        'Peak Area': gc_peak.area,
                        'Retention index': gc_peak.ri,
                        'Retention index Ref': compound_obj.ri,
                        'Retention Index Score': compound_obj.ri_score,
                        'Spectral Similarity Score': compound_obj.spectral_similarity_score,
                        'Similarity Score': compound_obj.similarity_score,
                        'Compound Name': compound_obj.name
                        }

            if self.gcms.molecular_search_settings.exploratory_mode:

                out_dict.update({
                    'Weighted Cosine Correlation': compound_obj.spectral_similarity_scores.get("weighted_cosine_correlation"), 
                    'Cosine Correlation': compound_obj.spectral_similarity_scores.get("cosine_correlation"), 
                    'Stein Scott Similarity': compound_obj.spectral_similarity_scores.get("stein_scott_similarity"), 
                    'Stein Scott Similarity Nist': compound_obj.spectral_similarity_scores.get("stein_scott_similarity_nist"),                     
                    'Pearson Correlation': compound_obj.spectral_similarity_scores.get("pearson_correlation"), 
                    'Spearman Correlation': compound_obj.spectral_similarity_scores.get("spearman_correlation"), 
                    'Kendall Tau Correlation': compound_obj.spectral_similarity_scores.get("kendall_tau_correlation"), 
                    'DFT Correlation': compound_obj.spectral_similarity_scores.get("dft_correlation"),
                    'DWT Correlation': compound_obj.spectral_similarity_scores.get("dwt_correlation"),
                    'Euclidean Distance': compound_obj.spectral_similarity_scores.get("euclidean_distance"), 
                    'Manhattan Distance': compound_obj.spectral_similarity_scores.get("manhattan_distance"),
                    'Jaccard Distance': compound_obj.spectral_similarity_scores.get("jaccard_distance"),
                })
                
                for method in methods_name:
                    
                    out_dict[methods_name.get(method)] = compound_obj.spectral_similarity_scores.get(method)

            dict_data_list.append(out_dict)
        
        def add_no_match_dict_data():

            dict_data_list.append({'Sample name': gcms.sample_name,
                                    'Peak Index': gcpeak_index,
                                    'Retention Time': gc_peak.rt,
                                    'Peak Height': gc_peak.tic,
                                    'Peak Area': gc_peak.area,
                                    'Retention index': gc_peak.ri,
                                  })

        for gcpeak_index, gc_peak in enumerate(gcms.sorted_gcpeaks):

            # check if there is a compound candidate 
            if gc_peak:
                
                if output_score_method == 'highest_sim_score':
                    compound_obj = gc_peak.highest_score_compound
                    add_match_dict_data()

                elif output_score_method == 'highest_ss':
                    
                    compound_obj = gc_peak.highest_ss_compound
                    add_match_dict_data()
                    
                else:
                    for compound_obj in gc_peak:
                        add_match_dict_data()  # add monoisotopic peak

            else:
                # include not_match
                if include_no_match and no_match_inline:
                    add_no_match_dict_data()

        if include_no_match and not no_match_inline:
            for gcpeak_index , gc_peak in enumerate(gcms.sorted_gcpeaks):
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
        return json.dumps(dict_ms_attrs, sort_keys=False, indent=4, separators=(',', ': '))

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

                    raw_ms_dataset.attrs['MassSpecAttrs'] = json.dumps(dict_ms_attrs, sort_keys=False, indent=4, separators=(',', ': '))
                    
                    if isinstance(mass_spectrum, MassSpecfromFreq):
                        raw_ms_dataset.attrs['TransientSetting'] = json.dumps(setting_dicts.get('TransientSetting'), sort_keys=False, indent=4, separators=(',', ': '))

                else:
                    
                    scan_group = hdf_handle.get(group_key)
        
                    # if there is not processed data len = 0, otherwise len() will return next index
                index_processed_data = str(len(scan_group.keys()))
                
                timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
                
                processed_dset = scan_group.create_dataset(index_processed_data, data=list_results)
            
                processed_dset.attrs['date_utc'] = timenow

                processed_dset.attrs['ColumnsLabels'] = columns_labels
                
                processed_dset.attrs['MoleculaSearchSetting'] = json.dumps(setting_dicts.get('MoleculaSearch'), sort_keys=False, indent=4, separators=(',', ': '))
                
                processed_dset.attrs['MassSpecPeakSetting'] = json.dumps(setting_dicts.get('MassSpecPeak'), sort_keys=False, indent=4, separators=(',', ': '))

                processed_dset.attrs['MassSpectrumSetting'] = json.dumps(setting_dicts.get('MassSpectrum'), sort_keys=False, indent=4, separators=(',', ': '))