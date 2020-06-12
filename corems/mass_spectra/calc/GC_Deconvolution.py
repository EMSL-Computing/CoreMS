

from matplotlib import pyplot as plt
import numpy as np

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroidLowRes
from corems.chroma_peak.factory.ChromaPeakClasses import GCPeak
from corems.mass_spectra.calc import SignalProcessing as sp

class MassDeconvolution:

    def run_deconvolution(self):

        maximum_tic = max(self._processed_tic)

        eic_dict = self.ion_extracted_chroma(self._ms)
        peaks_entity_data = self.find_peaks_entity(eic_dict)
        
        # select model peaks, create Mass Spectrum objs, GCPeak objs, store results in GC_Class gcpeaks obj 
        self.deconvolution(peaks_entity_data, maximum_tic)

    def centroid_detector(self, tic, rt):
        
        # need to change the parameter to accommodate EIC peak picking 
        # needs a better algorithm to detect start and end of a peak
        
        noise_std = self.chromatogram_settings.std_noise_threshold

        method = self.chromatogram_settings.noise_threshold_method
        
        # peak picking
        min_height = self.chromatogram_settings.peak_height_min_percent
        min_datapoints = self.chromatogram_settings.min_peak_datapoints
        
        # baseline detection
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent
        max_height = self.chromatogram_settings.peak_height_max_percent
        
        peak_indexes_generator = sp.peak_detector_generator(tic, noise_std, method, rt, max_height, min_height, max_prominence, min_datapoints)

        return peak_indexes_generator

    def ion_extracted_chroma(self, mass_spectra_obj):
        
        eic_dict = {} 
        
        for scan_number, ms_obj in mass_spectra_obj.items():

            mz_list = ms_obj.mz_exp
            abundance_list = ms_obj.abundance
            # add list of scan numbers
            for index, mz in enumerate(mz_list):

                # dict of mz and tuple (mass spectrum abundances index, and scan number)
                if not mz in eic_dict.keys(): 
                   
                    eic_dict[mz] = [ [abundance_list[index]], [ms_obj.rt] ]
                
                else:
                    
                    eic_dict[mz][0].append(ms_obj.abundance[index])
                    eic_dict[mz][1].append(ms_obj.rt)
        
        return eic_dict           

    def find_peaks_entity(self, eic_dict):
        
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent
        max_height = self.chromatogram_settings.peak_height_max_percent

        peaks_entity_data = {}

        for mz, eic_scan_index_rt in eic_dict.items():
            
            eic = eic_scan_index_rt[0]
            rt_list = eic_scan_index_rt[1]
                                    
            if len(eic) > 10:

                smooth_eic = self.smooth_tic(eic)

                sp.find_minima_derivative(rt_list, smooth_eic,  max_height, max_prominence, correct_baseline=False)

                include_indexes = list(self.centroid_detector(smooth_eic, rt_list))
                
                for initial_scan, apex_scan, final_scan in include_indexes:

                        rt_corrected_therm = self.quadratic_interpolation(rt_list, smooth_eic, apex_scan)
                        
                        ref_apex_rt = round(rt_list[apex_scan] + rt_corrected_therm,4)
                        apex_rt = rt_list[apex_scan]
                        apex_abundance = smooth_eic[apex_scan]
                        
                        #maximum_tic = apex_abundance if apex_abundance > maximum_tic else maximum_tic
                        
                        for scan_index in range(initial_scan, final_scan):
                            
                            peak_rt = rt_list[scan_index]
                            peak_abundance = smooth_eic[scan_index]

                            dict_data = {peak_rt: { 'mz': [mz] , 
                                                    'abundance':[peak_abundance], 
                                                    'scan_number': [scan_index] },
                                        "ref_apex_rt": ref_apex_rt
                                        }

                            if not apex_rt in peaks_entity_data.keys():
                                
                                peaks_entity_data[apex_rt] = dict_data
                            
                            else:
                                
                                if not peak_rt in peaks_entity_data[apex_rt].keys():
                                    
                                    peaks_entity_data[apex_rt][peak_rt] = dict_data.get(peak_rt)

                                else:    
                                    
                                    existing_data = peaks_entity_data[apex_rt].get(peak_rt)
                                    
                                    existing_data['mz'].append(mz)
                                    existing_data['abundance'].append(peak_abundance)
                                    existing_data['scan_number'].append(scan_index)    
        
        return peaks_entity_data                            
        
    def deconvolution(self, peaks_entity_data, maximum_tic):
        
        i = 0
        tic_list = []
        rt_list = []

        for apex_rt, datadict in sorted(peaks_entity_data.items()):
            
            if apex_rt in datadict.keys():
                
                apex_data = datadict[apex_rt] 
                
                ref_apex_rt = datadict["ref_apex_rt"]

                tic = sum(apex_data.get('abundance'))
                
                norm_smooth_tic = (tic/ maximum_tic)*100

                if norm_smooth_tic > self.chromatogram_settings.peak_height_min_percent and len(apex_data['mz']) > 3:
                    
                    scan_index = apex_data['scan_number'][0]
                    
                    mz_list, abundance_list = zip(*sorted(zip(apex_data['mz'], apex_data['abundance'])))

                    data_dict = {Labels.mz: mz_list, Labels.abundance: abundance_list}
                    
                    d_params = default_parameters(self._ms[scan_index]._filename)
                    
                    d_params["rt"] = apex_rt

                    d_params["scan_number"] = scan_index

                    d_params['label'] = Labels.gcms_centroid

                    d_params["polarity"] = self._ms[scan_index].polarity

                    d_params['analyzer'] = self._ms[scan_index].analyzer

                    d_params['instrument_label'] = self._ms[scan_index].instrument_label
                    
                    d_params["filename_path"] = self._ms[scan_index].instrument_label

                    ms = MassSpecCentroidLowRes(data_dict, d_params )
                    
                    #needs to define peak start and end, passing just minus and plus one from apex pos for now
                    gc_peak =  GCPeak(ms, (i-1,i, i+1))
                    i += 1
                    self.gcpeaks.append(gc_peak)

                    tic_list.append(tic)
                    rt_list.append(ref_apex_rt)

                    peak_rt = []
                    peak_tic = []
                    
                    for rt, each_datadict in datadict.items():
                        
                        if rt != "ref_apex_rt":
                            peak_rt.append(rt)
                            peak_tic.append(sum(each_datadict["abundance"]))
                    
                    peak_rt, peak_tic = zip(*sorted(zip(peak_rt, peak_tic)))

                    #ax = plt.gca()

                    #markerline_a, stemlines_a, baseline_a  = ax.stem(data[0], data[1], linefmt='-',  markerfmt=" ", use_line_collection =True, label=rt)
                    
                    #plt.setp(markerline_a, 'color', c, 'linewidth', 2)
                    #plt.setp(stemlines_a, 'color', c, 'linewidth', 2)
                    #plt.setp(baseline_a, 'color', c, 'linewidth', 2)

                    #ax.set_xlabel("$\t{m/z}$", fontsize=12)
                    #ax.set_ylabel('Abundance', fontsize=12)
                    #ax.tick_params(axis='both', which='major', labelsize=12)

                    #ax.axes.spines['top'].set_visible(False)
                    #ax.axes.spines['right'].set_visible(False)

                    #ax.get_yaxis().set_visible(False)
                    #ax.spines['left'].set_visible(False)
                    plt.plot(peak_rt, peak_tic)
                    #plt.legend()
                    #plt.show()
                    #plt.close()
        #self.rt_clustering(rt_list, tic_list)
        plt.plot(self.retention_time, self._processed_tic, c='black')
        plt.plot(rt_list, tic_list, c='black', marker= '^', linewidth=0)
        plt.show()
        #print(i)

    
    def quadratic_interpolation(self, rt_list, tic_list, apex_index):
        
        rt_list = np.array(rt_list)
        tic_list = np.array(tic_list)
        three_highest_i = [i for i in range(apex_index-1, apex_index+2)]
        
        z = np.poly1d(np.polyfit(rt_list[three_highest_i], tic_list[three_highest_i], 2))
        a = z[2]
        b = z[1]

        corrected_apex_rt = -b/(2*a)
        initial_rt = rt_list[apex_index]
        
        return initial_rt- corrected_apex_rt
