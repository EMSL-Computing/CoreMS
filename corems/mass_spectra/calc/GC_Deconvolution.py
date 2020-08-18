

from matplotlib import pyplot as plt
import numpy as np

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroidLowRes
from corems.chroma_peak.factory.ChromaPeakClasses import GCPeak
from corems.mass_spectra.calc import SignalProcessing as sp

class MassDeconvolution:

    def run_deconvolution(self):

        eic_dict = self.ion_extracted_chroma(self._ms)
        
        peaks_entity_data = self.find_peaks_entity(eic_dict)
        
        # select model peaks, create Mass Spectrum objs, GCPeak objs, store results in GC_Class gcpeaks obj 
        self.deconvolution(peaks_entity_data)

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

    def hc(self, X, Y, max_rt_distance=0.025):
        
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.cluster.hierarchy import fcluster
        from matplotlib import pyplot as plt

        Z = linkage(np.reshape(X, (len(X), 1)), method  = "ward")
        #Z = linkage(X, method  = "ward")
        
        max_d = max_rt_distance
        distance_clusters = fcluster(Z, max_d, criterion='distance')
        # print("distance")
        # print(distance_clusters)

        #inconsistency_cluster = fcluster(Z, 2, depth=2)
        #max_cluster = fcluster(Z, 2, criterion='maxclust')

        grouped_rt = {}
        for index_obj, group in enumerate(distance_clusters):
           
            if not group in grouped_rt.keys():
                grouped_rt[group] = [X[index_obj]]
            else:
                grouped_rt[group].append(X[index_obj])
        
        return grouped_rt      

        #plt.figure(figsize=(10, 8))
        # plt.scatter(X, Y, c=distance_clusters, cmap='prism')  # plot points with cluster dependent colors
        # plt.show()
        #labelList = range(int(min(X)), int(max(X)))

        # plt.figure(figsize=(10, 7))
        # dendrogram(Z,
        #            orientation='top',
        #            distance_sort='descending',
        #            show_leaf_counts=True)
        # plt.show()
        # print(Z)

    
    def find_peaks_entity(self, eic_dict):
        
        ''' combine eic with mathing rt apexes''' 
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent
        
        max_height = self.chromatogram_settings.peak_height_max_percent

        signal_threshold = self.chromatogram_settings.eic_signal_threshold

        correct_baseline = False
        peaks_entity_data = {}

        max_eic = 0
        for mz, eic_scan_index_rt in eic_dict.items():
            
            ind_max_eic = max(eic_scan_index_rt[0])
            max_eic = ind_max_eic if ind_max_eic > max_eic else max_eic
        
        for mz, eic_scan_index_rt in eic_dict.items():
            
            eic = eic_scan_index_rt[0]
            rt_list = eic_scan_index_rt[1]

            if len(eic) > 10:

                smooth_eic = self.smooth_tic(eic)

                include_indexes = sp.find_minima_derivative(rt_list, smooth_eic,  max_height, max_prominence, max_eic, 
                                                            signal_threshold=signal_threshold,  correct_baseline=correct_baseline)

                #include_indexes = list(self.centroid_detector(smooth_eic, rt_list))
                
                for initial_scan, apex_scan, final_scan in include_indexes:

                        rt_corrected_therm =  0 #self.quadratic_interpolation(rt_list, smooth_eic, apex_scan)
                        
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
        
    def smooth_signal(self, signal):
            
        implemented_smooth_method = self.chromatogram_settings.implemented_smooth_method
        
        pol_order = self.chromatogram_settings.savgol_pol_order

        window_len = self.chromatogram_settings.smooth_window

        window = self.chromatogram_settings.smooth_method

        return sp.smooth_signal(signal, window_len, window, pol_order, implemented_smooth_method)

    def deconvolution(self, peaks_entity_data):
        
        i = 0
        tic_list = []
        rt_list = []

        domain = self.retention_time
        signal = self._processed_tic
        max_height = self.chromatogram_settings.peak_height_max_percent
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent
        
        max_signal = max(signal)
        
        signal_threshold = self.chromatogram_settings.peak_height_min_percent
        correct_baseline = False
        
        max_rt_distance = self.chromatogram_settings.max_rt_distance

        include_indexes = sp.find_minima_derivative(domain, signal,  max_height, max_prominence, max_signal, 
                                                    signal_threshold=signal_threshold, correct_baseline=correct_baseline, plot_res=False)
        
        ''' deconvolution window is defined by the TIC peak region'''
        all_apexes_rt = np.array(list(peaks_entity_data.keys()))

        for indexes_tuple in include_indexes:
            
            start_rt = self.retention_time[indexes_tuple[0]]
            apex_rt = self.retention_time[indexes_tuple[1]]
            final_rt = self.retention_time[indexes_tuple[2]]

            #plt.plot(self.retention_time[indexes_tuple[1]], signal[indexes_tuple[1]], marker= '^', c='red',  linewidth=0)
            #plt.plot(self.retention_time[indexes_tuple[0]], signal[indexes_tuple[0]], marker= '^', c='blue',  linewidth=0)
            #plt.plot(self.retention_time[indexes_tuple[2]], signal[indexes_tuple[2]], marker= '^', c='blue',  linewidth=0)
            
            ''' find all features within TIC peak window'''
            peak_features_indexes = np.where((all_apexes_rt >= start_rt) & (all_apexes_rt <= final_rt))[0]
            peak_features_rts = all_apexes_rt[peak_features_indexes]

            filtered_features_rt = []
            filtered_features_abundance = []
            
            for each_apex_rt in peak_features_rts:
                
                apex_data = peaks_entity_data.get(each_apex_rt).get(each_apex_rt)

                peak_features_tic = sum(peaks_entity_data.get(each_apex_rt).get(each_apex_rt).get('abundance'))
            
                norm_smooth_tic = (peak_features_tic/ max_signal)*100

                ''' TODO: 
                    Improve Peak Filtering

                    Calculate peaks sharpness here and filter it out (Amax - An /n)?
                    Peak Fit and Calculate Peak Gaussian Similarity?
                    Currentely using flat % tic relative abundance threshold and min 3 m/z per mass spectrum
                '''
                if norm_smooth_tic > signal_threshold and len(apex_data['mz']) > 1:
                       
                       #print(len(apex_data['mz']))
                       filtered_features_rt.append(each_apex_rt)
                       filtered_features_abundance.append(peak_features_tic)
            
            if len(filtered_features_rt) > 1: 
                ''' more than one peak feature identified inside a TIC peak  '''
                # plt.plot(self.retention_time[indexes_tuple[0]:indexes_tuple[2]], signal[indexes_tuple[0]:indexes_tuple[2]], c='black')

                #print(filtered_features_rt)
                grouped_rt = self.hc(filtered_features_rt, filtered_features_abundance, max_rt_distance=max_rt_distance)
                #print(grouped_rt)
                
                for group, apex_rt_list in grouped_rt.items():
                    ''' each group is a peak feature defined by the hc algorithm
                        summing the apex data only, for now
                    '''
                    group_datadict = {}
                    group_datadict['ref_apex_rt'] = []

                    for each_group_apex_rt in apex_rt_list:
                        
                        datadict = peaks_entity_data.get(each_group_apex_rt)

                        for rt, each_datadict in datadict.items():
                            
                            if rt == "ref_apex_rt":
                                
                                group_datadict['ref_apex_rt'].append(each_datadict)
                            
                            else:
                                
                                if rt in group_datadict.keys():
                                    
                                    mz_list = each_datadict.get("mz")
                                    abundance_list = each_datadict.get("abundance")
                                    
                                    each_mz_abun = dict(zip(mz_list, abundance_list)) 

                                    for index_mz, mz in enumerate(group_datadict[rt].get("mz")):
                                        if mz in each_mz_abun.keys():
                                            
                                            each_mz_abun[mz] = each_mz_abun[mz] + group_datadict[rt].get("abundance")[index_mz]
                                        
                                        else:    
                                            
                                            each_mz_abun[mz] = group_datadict[rt].get("abundance")[index_mz]
                                    
                                    group_datadict[rt] = { 'mz': list(each_mz_abun.keys()) , 
                                                        'abundance': list(each_mz_abun.values()), 
                                                        'scan_number': each_datadict.get('scan_number') }
                                                        

                                else:
                                    
                                    group_datadict[rt] = each_datadict
                            
                    peak_rt = []
                    peak_tic = []

                    #print(group_datadict.get('ref_apex_rt'))
                    for rt, each_datadict in group_datadict.items():
                        if rt != "ref_apex_rt":
                            peak_rt.append(rt)
                            peak_tic.append(sum(each_datadict["abundance"]))
                    
                    peak_rt, peak_tic = zip(*sorted(zip(peak_rt, peak_tic)))
                    
                    smoothed_tic = self.smooth_signal(peak_tic)
                    
                    include_indexes = sp.find_minima_derivative(peak_rt, smoothed_tic,  max_height, max_prominence, max_signal, 
                                                                            signal_threshold=signal_threshold,  correct_baseline=False, plot_res=False)
                    
                    #include_indexes = sp.find_minima_derivative(domain, signal,  max_height, max_prominence, max_signal, signal_threshold=signal_threshold, correct_baseline=correct_baseline)

                    include_indexes = list(include_indexes)
                    #print("start include_indexes")
                    #print(list(include_indexes))
                    #print("end include_indexes")
                    #print(include_indexes)
                    if include_indexes:

                        if len(include_indexes) > 1:
                            
                            for new_apex_index in include_indexes:
                                
                                if start_rt <= peak_rt[new_apex_index[1]] <= final_rt:
                                    
                                    #pass
                                    plt.plot(peak_rt[new_apex_index[1]], smoothed_tic[new_apex_index[1]], c='red', marker= '^', linewidth=0)  
                                    plt.plot(peak_rt[new_apex_index[0]:new_apex_index[2]], smoothed_tic[new_apex_index[0]:new_apex_index[2]])

                                else:
                                    
                                    print('new apex outside deconvolution window')    
                        
                        else:
                            
                            new_apex_index = include_indexes[0]

                            plt.plot(peak_rt[new_apex_index[1]], smoothed_tic[new_apex_index[1]], c='red', marker= '^', linewidth=0)  
                            plt.plot(peak_rt[new_apex_index[0]:new_apex_index[2]], smoothed_tic[new_apex_index[0]:new_apex_index[2]])
                    
                    #rt_list.append(peak_rt[include_indexes])
                    #tic_list.append(sum(smoothed_tic))

            elif len(filtered_features_rt) == 1:
                ''' only one peak feature inside deconvolution window '''
                
                each_apex_rt = filtered_features_rt[0]

                datadict = peaks_entity_data.get(each_apex_rt)

                apex_data = datadict.get(each_apex_rt)

                tic = sum(apex_data.get('abundance'))

                scan_index = apex_data['scan_number'][0]

                ref_apex_rt = datadict["ref_apex_rt"]
                
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

                ''' Only stores the apex data 
                    TODO Store all ms on the GCPeak obj
                    It Might be necessary to to create a new obj or change workflow with out deconvolution
                    Indexes on original objs are related to the TIC 
                '''
                
                ms = MassSpecCentroidLowRes(data_dict, d_params )
                
                '''TODO fake indexes, area calculation will fail 
                    workaround is to calculate the area here and store the result after the gcpeak obj is created
                '''
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
                    
                #peak_rt, peak_tic = zip(*sorted(zip(peak_rt, peak_tic)))

                #plt.plot(peak_rt, self.smooth_signal(peak_tic))
                

            else:
                
                #print('no data after filter')
                pass
        '''        
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
        
        '''

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
