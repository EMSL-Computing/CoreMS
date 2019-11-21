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
                    error_list.append(mformula._calc_assigment_mass_error(mspeak.mz_exp))
                    list_indexes_mass_spec.append(index)
        
        dict = {'mass': mass_list, 'error': error_list}  
        df = pd.DataFrame(dict) 
        matrix_data = df.values.astype("float32", copy = False)
        return matrix_data, list_indexes_mass_spec

    def get_kendrick_matrix_data(self, mass_spectrum):
        km = mass_spectrum.kendrick_mass
        kdm = mass_spectrum.kmd
        dict = {'km': km, 'kdm': kdm}  
        df = pd.DataFrame(dict) 
        matrix_data = df.values.astype("float32", copy = False)
        return matrix_data

    def filter_kendrick(self, mass_spectrum):
        
        matrix_data = self.get_kendrick_matrix_data(mass_spectrum)

        stscaler = StandardScaler().fit(matrix_data)
        
        matrix_data_scaled = stscaler.transform(matrix_data)

        clusters = DBSCAN(eps = .15, min_samples=15).fit_predict(matrix_data_scaled)
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(clusters)) - (1 if -1 in clusters else 0)
        n_noise_ = list(clusters).count(-1)
        
        indexes = []
        for i in range(len(clusters)):
            if clusters[i] == -1:
                indexes.append(i)
        
        print('Estimated number of clusters: %d' % n_clusters_)
        print('Estimated number of noise points: %d' % n_noise_)
        mass_spectrum.filter_by_index(indexes)
        from matplotlib import pyplot as plt
        plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="spectral")
        plt.xlabel("km")
        plt.ylabel("kdm")
        plt.show()
        plt.close()
       

    def remove_assigment_by_mass_error(self, mass_spectrum):
        
        #data need to be binned by mz unit or more to be able to use clustering
        
        matrix_data, list_indexes_mass_spec = self.get_mass_error_matrix_data(mass_spectrum)

        stscaler = StandardScaler().fit(matrix_data)
        
        matrix_data_scaled = stscaler.transform(matrix_data)
        
        #bandwidth = estimate_bandwidth(matrix_data_scaled, quantile=0.3, n_samples=int(len(ms_peaks)/3))

        #clusters = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit_predict(matrix_data_scaled)
        
        #eps and min_samp need to be optmized by precision and number of mspeaks
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
        
        mass_spectrum.remove_assigment_by_index(indexes)    