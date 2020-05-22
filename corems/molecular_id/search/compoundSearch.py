
from threading import Thread
from pathlib import Path

from pandas import DataFrame

from numpy import power, dot, absolute, subtract, intersect1d, where, average
from numpy.linalg import norm
from scipy.spatial.distance import cosine, jaccard, euclidean, cityblock
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics.pairwise import cosine_similarity
from math import exp
from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite
from scipy import fft
from scipy.fftpack import rfft

class LowResMassSpectralMatch(Thread):

    def __init__(self, gcms_obj, sql_obj=None, calibration=False):
        
        '''TODO:
        '''
        Thread.__init__(self)
        
        self.gcms_obj = gcms_obj

        #  initiated at create_molecular_database()
        #self.dict_molecular_lookup_table = None
        self.calibration = calibration
        # reading local file for now, 
        if not sql_obj:
            self.sql_obj = EI_LowRes_SQLite(url=self.gcms_obj.molecular_search_settings.url_database)
        else:
            self.sql_obj = sql_obj
    
    def metabolite_detector_score(self, gc_peak, ref_obj):

        spectral_similarity_scores = {}
        spectral_similarity_scores["cosine_correlation"] = self.cosine_correlation(gc_peak.mass_spectrum, ref_obj)
        
        if self.gcms_obj.molecular_search_settings.exploratory_mode:
            
            spectral_similarity_scores["weighted_cosine_correlation"] = self.weighted_cosine_correlation(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["stein_scott_similarity"] = self.stein_scott(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["pearson_correlation"] = self.pearson_correlation(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["spearman_correlation"] = self.spearman_correlation(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["kendall_tau_correlation"] = self.kendall_tau(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["euclidean_distance"] = self.euclidean_distance(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["manhattan_distance"] = self.manhattan_distance(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["jaccard_distance"] = self.jaccard_distance(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["dft_correlation"] = self.dft_correlation(gc_peak.mass_spectrum, ref_obj)
            spectral_similarity_scores["dwt_correlation"] = self.dwt_correlation(gc_peak.mass_spectrum, ref_obj)

        #print(ref_obj.get('ri'), gc_peak.ri, self.gcms_obj.molecular_search_settings.ri_window)

        ri_score = exp( -1*(power((gc_peak.ri - ref_obj.get('ri')), 2 )  / (2 * power(self.gcms_obj.molecular_search_settings.ri_std, 2)) ))

        similarity_score = ((spectral_similarity_scores.get("cosine_correlation")**2) * (ri_score))**(1/3)

        return spectral_similarity_scores, ri_score, similarity_score
        
    def weighted_cosine_correlation(self, mass_spec, ref_obj):
        
        a , b = 0.5, 1.3    
        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # weight exp data 

        exp_abun = list(ms_mz_abun_dict.values())
        exp_mz = list(ms_mz_abun_dict.keys())

        xc = power(exp_abun, a) *  power(exp_mz, b) 
        
        # track back to individual mz
        weighted_exp_dict = dict(zip(ms_mz_abun_dict.keys(), xc))

        # weight ref data
        yc = power(ref_obj.get("abundance"), a) *  power(ref_obj.get("mz"), b) 
        
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), yc))

        #parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([weighted_exp_dict, ref_mz_abun_dict])

        # fill missing mz with weight {abun**a}{m/z**b} to 0
        df.fillna(0, inplace=True)
        
        #calculate cosine correlation, 
        x = df.T[0]
        y = df.T[1]

        #correlation = (1 - cosine(x, y))
        
        correlation = dot(x, y)/(norm(x)*norm(y))

        return correlation

    def cosine_correlation(self, mass_spec, ref_obj):

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        #print(ref_mz_abun_dict)
        #print(ms_mz_abun_dict)

        #parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        #calculate cosine correlation, 
        x = df.T[0]
        y = df.T[1]

        correlation = (1 - cosine(x, y))
        
        #correlation = dot(x, y)/(norm(x)*norm(y))

        return correlation

    def stein_scott(self, mass_spec, ref_obj):

        ms_mz_abun_dict = mass_spec.mz_abun_dict
 
        exp_abun = list(ms_mz_abun_dict.values())
        exp_mz = list(ms_mz_abun_dict.keys())
 
        ref_abun = ref_obj.get("abundance")
        ref_mz = ref_obj.get("mz")
 
        # important: I assume ref_mz and ref_abun are in order, and one-to-one; this needs to be be verified
        ref_mz_abun_dict = dict(zip(ref_mz, ref_abun))
 
        # filter out the mass values that have zero intensities in exp_mz
        exp_mz_filtered = set([k for k in exp_mz if ms_mz_abun_dict[k] != 0])
 
        # filter out the mass values that have zero intensities in ref_mz
        ref_mz_filtered = set([k for k in ref_mz if ref_mz_abun_dict[k] != 0])
 
        # find the intersection/common mass values of both ref and exp, and sort them
        common_mz_values = sorted(list(exp_mz_filtered.intersection(ref_mz_filtered)))

        # find the number of common mass values (after filtering 0s)
        n_x_y = len(common_mz_values)
 
        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in exp_abun)
 
        s_r_x_y = 0
        for i in range(1,n_x_y):
            current_value = common_mz_values[i]
            previous_value = common_mz_values[i-1]
 
            y_i = ref_mz_abun_dict[current_value]
            y_i_minus1 = ref_mz_abun_dict[previous_value]
            x_i = ms_mz_abun_dict[current_value]
            x_i_minus1 = ms_mz_abun_dict[previous_value]
 
            temp_computation = (y_i/y_i_minus1) * (x_i_minus1/x_i)
            n = 0
            if temp_computation <= 1:
                n = -1
            else:
                n = 1
            s_r_x_y = s_r_x_y + temp_computation ** n
 
        # finish the calculation of S_R(X,Y)
        s_r_x_y = s_r_x_y / n_x_y
 
        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation(mass_spec, ref_obj)
 
        # final step
        s_ss_x_y = (n_x * s_wc_x_y + n_x_y * s_r_x_y) / (n_x + n_x_y)

        return s_ss_x_y

    def pearson_correlation(self, mass_spec, ref_obj):

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Pearson correlation
        correlation = pearsonr(df.T[0], df.T[1])

        return correlation[0]

    def spearman_correlation(self, mass_spec, ref_obj):

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Spearman correlation
        ### TODO - Check axis
        correlation = spearmanr(df.T[0], df.T[1], axis=0)

        return correlation[0]

    def kendall_tau(self, mass_spec, ref_obj):

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Kendall's tau
        correlation = kendalltau(df.T[0], df.T[1])

        return correlation[0]

    def euclidean_distance(self, mass_spec, ref_obj):

        def euclidean_distance_manual(qlist,rlist):

            T1=sum(subtract(qlist,rlist)**2)

            T2=sum((qlist)**2)

            return (1+T1/T2)**(-1)
            
        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Pearson correlation
        #correlation = euclidean(df.T[0], df.T[1])

        correlation = euclidean_distance_manual(df.T[0], df.T[1])
        return correlation

    def manhattan_distance(self, mass_spec, ref_obj):
        
        def mann_distance_manual(qlist,rlist):
        
            T1=sum(absolute(subtract(qlist,rlist)))
            T2=sum(qlist)
            return (1+T1/T2)**(-1)

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate manhattan correlation
        #correlation = cityblock(df.T[0], df.T[1])
        correlation = mann_distance_manual(df.T[0], df.T[1])
        
        return correlation

    def dwt_correlation(self, mass_spec, ref_obj):

        return 0
        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofilland tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate DFT correlation
        #
        # return correlation

    def dft_correlation(self, mass_spec, ref_obj):

        ms_mz_abun_dict = mass_spec.mz_abun_dict
 
        exp_abun = list(ms_mz_abun_dict.values())
        exp_mz = list(ms_mz_abun_dict.keys())
 
        ref_abun = ref_obj.get("abundance")
        ref_mz = ref_obj.get("mz")
 
        # important: I assume ref_mz and ref_abun are in order, and one-to-one; this needs to be be verified
        ref_mz_abun_dict = dict(zip(ref_mz, ref_abun))
 
        # filter out the mass values that have zero intensities in exp_mz
        exp_mz_filtered = set([k for k in exp_mz if ms_mz_abun_dict[k] != 0])
 
        # filter out the mass values that have zero intensities in ref_mz
        ref_mz_filtered = set([k for k in ref_mz if ref_mz_abun_dict[k] != 0])
 
        # find the intersection/common mass values of both ref and exp, and sort them
        
        common_mz_values = exp_mz_filtered.intersection(ref_mz_filtered)
        
        # find the number of common mass values (after filtering 0s)
        n_x_y = len(common_mz_values)
 
        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in exp_abun)

        #parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])
 
        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        #calculate cosine correlation, 
        x = df.T[0]
        y = df.T[1]
 
        # get the Fourier transform of x and y
        x_dft = rfft(x)
        y_dft = rfft(y)
 
        s_dft_xy = dot(x_dft, y_dft)/(norm(x_dft)*norm(y_dft))
 
        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation(mass_spec, ref_obj)
 
        # final step
        s_dft = (n_x * s_wc_x_y + n_x_y * s_dft_xy) / (n_x + n_x_y)

        return s_dft

    def jaccard_distance(self, mass_spec, ref_obj):

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofilland tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate jaccard correlation
        correlation = jaccard(df.T[0], df.T[1])

        return correlation
    #@timeit
    def run(self):
        # TODO select the best gcms peak    
        import tqdm
        
        if not self.gcms_obj:
            
            self.gcms_obj.process_chromatogram()
        list_cpu = []
        
        for gc_peak in tqdm.tqdm(self.gcms_obj):
            
            if not self.calibration:
                
                window = self.gcms_obj.molecular_search_settings.ri_search_range

                ri = gc_peak.ri
                
                min_mat_ri = (ri-window, ri+window)    
                
                ref_objs = self.sql_obj.query_min_max_ri(min_mat_ri)
                
            else:

                compound_names = self.gcms_obj.molecular_search_settings.ri_calibration_compound_names

                window = self.gcms_obj.molecular_search_settings.rt_search_range

                rt = gc_peak.rt

                min_mat_rt = (rt-window, rt+window)    
                
                ref_objs = self.sql_obj.query_names_and_rt(min_mat_rt, compound_names)
                
            for ref_obj in ref_objs:
                # uses spectral similarly and uses a threshold to only select peaks with high data correlation
                if self.calibration:
                    
                    spectral_similarity_scores = {}
                    spectral_similarity_scores["cosine_correlation"] = self.cosine_correlation(gc_peak.mass_spectrum, ref_obj)

                    #print(w_correlation_value,correlation_value )
                    if spectral_similarity_scores["cosine_correlation"] >= self.gcms_obj.molecular_search_settings.correlation_threshold:
                    
                        gc_peak.add_compound(ref_obj, spectral_similarity_scores)

                # use score, usually a combination of Retention index and Spectral Similarity
                # Threshold is implemented by not necessarily used
                else:

                    # m/q developed methods will be implemented here
                    spectral_similarity_scores, ri_score, similarity_score = self.metabolite_detector_score(gc_peak, ref_obj)   
                    
                    #TODO need to add similarity score option in the parameters encapsulation class
                    if spectral_similarity_scores.get("cosine_correlation") >= self.gcms_obj.molecular_search_settings.score_threshold:
                    
                        gc_peak.add_compound(ref_obj, spectral_similarity_scores, ri_score, similarity_score)
                
        self.sql_obj.session.close()
        self.sql_obj.engine.dispose()