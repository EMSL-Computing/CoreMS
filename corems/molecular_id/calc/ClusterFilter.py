
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import MeanShift, estimate_bandwidth



import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

class ClusteringFilter():
    
    def get_mass_error_matrix_data(self, ms_peaks):
        
        mass_list = list()
        error_list = list()
        list_indexes_mass_spec = []
        
        for index, mspeak in enumerate(ms_peaks):

            if mspeak.is_assigned:
                    
                #print(mspeak.mz_exp, len(mspeak))
                for mformula in mspeak:
                    mass_list.append(mspeak.mz_exp)
                    error_list.append(mformula.mz_error)
                    list_indexes_mass_spec.append(index)
        
        kendrick_dict = {'mass': mass_list, 'error': error_list}  
        df = pd.DataFrame(kendrick_dict) 
        matrix_data = df.values.astype("float32", copy = False)
        return matrix_data, list_indexes_mass_spec

    def get_kendrick_matrix_data(self, mass_spectrum):
        
        km = mass_spectrum.kendrick_mass
        kdm = mass_spectrum.kmd
        kendrick_dict = {'km': km, 'kdm': kdm}  
        df = pd.DataFrame(kendrick_dict) 
        matrix_data = df.values.astype("float32", copy = False)
        return matrix_data

    def filter_kendrick(self, mass_spectrum):
        
        matrix_data = self.get_kendrick_matrix_data(mass_spectrum)

        stdscaler = StandardScaler().fit(matrix_data)
        
        matrix_data_scaled = stdscaler.transform(matrix_data)

        clusters = DBSCAN(eps = .75, min_samples=50).fit_predict(matrix_data_scaled)
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(clusters)) - (1 if -1 in clusters else 0)
        n_noise_ = list(clusters).count(-1)
        
        indexes = []
        for i in range(len(clusters)):
            if clusters[i] == -1:
                indexes.append(i)
        
        print('Estimated number of clusters: %d' % n_clusters_)
        print('Estimated number of noise points: %d' % n_noise_)
        print()
        mass_spectrum.filter_by_index(indexes)
        #from matplotlib import pyplot as plt
        #plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="jet")
        #plt.xlabel("km")
        #plt.ylabel("kdm")
        #plt.show()
        #plt.close()

    def filter_kendrick_by_index(self, ms_peak_indexes, mass_spectrum_obj):
        
        min_samples = mass_spectrum_obj.molecular_search_settings.min_peaks_per_class

        kendrick_dict = {'km': list(), 'kmd': list()}  

        if len(ms_peak_indexes) <= 1: return []
        
        for index, _ in ms_peak_indexes:
           kendrick_dict["km"].append(mass_spectrum_obj[index].kendrick_mass)
           kendrick_dict["kmd"].append(mass_spectrum_obj[index].kmd)
           
        # check min data points otherwise StandardScaler().fit(0 will fail
        
        df = pd.DataFrame(kendrick_dict) 
        matrix_data = df.values.astype("float32", copy = False)

        stdscaler = StandardScaler().fit(matrix_data)
        matrix_data_scaled = stdscaler.transform(matrix_data)

        clusters = DBSCAN(eps = .8, min_samples=min_samples).fit_predict(matrix_data_scaled)
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(clusters)) - (1 if -1 in clusters else 0)
        n_noise_ = list(clusters).count(-1)
        
        print('Estimated number of clusters: %d' % n_clusters_)
        print('Estimated number of noise points: %d' % n_noise_)
        print()

        noise_idx = []
        
        other_peaks_idx = []

        for i in range(len(clusters)):
            
            if clusters[i] == -1:
                noise_idx.append(ms_peak_indexes[i])
            
            else:
                other_peaks_idx.append(ms_peak_indexes[i])    

        #mfs = [mass_spectrum_obj[index].best_molecular_formula_candidate.string for index in other_peaks_idx]
        
        #mfs_noise = [mass_spectrum_obj[index].best_molecular_formula_candidate.string for index in noise_idx]
        
        #print(mfs)
        #print(mfs_noise)

        #from matplotlib import pyplot as plt
        #plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="jet")
        #plt.xlabel("km")
        #plt.ylabel("kdm")
        #plt.show()
        #plt.close()
        
        return noise_idx      

    def remove_assignment_by_mass_error(self, mass_spectrum):
        
        #data need to be binned by mz unit or more to be able to use clustering
        
        matrix_data, list_indexes_mass_spec = self.get_mass_error_matrix_data(mass_spectrum)

        stdscaler = StandardScaler().fit(matrix_data)
        
        matrix_data_scaled = stdscaler.transform(matrix_data)
        
        #bandwidth = estimate_bandwidth(matrix_data_scaled, quantile=0.3, n_samples=int(len(ms_peaks)/3))

        #clusters = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit_predict(matrix_data_scaled)
        
        #eps and min_samp need to be optimized by precision and number of mspeaks
        clusters = DBSCAN(eps = .15).fit_predict(matrix_data_scaled)
        
        indexes = []
        
        #from matplotlib import pyplot as plt
        #plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="plasma")
        #plt.xlabel("km")
        #plt.ylabel("kdm")
        #plt.show()
        #plt.close()

        for i in range(len(clusters)):
            if clusters[i] == -1:
                indexes.append(list_indexes_mass_spec[i])
        
        mass_spectrum.remove_assignment_by_index(indexes)    