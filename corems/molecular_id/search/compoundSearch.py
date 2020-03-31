
from threading import Thread
from pathlib import Path

from pandas import DataFrame

from numpy import power, dot, absolute, subtract, sum
from numpy.linalg import norm
from scipy.spatial.distance import cosine, jaccard, euclidean, cityblock
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics.pairwise import cosine_similarity
from math import exp

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.encapsulation.settings.processingSetting import CompoundSearchSettings, GasChromatographSetting

class LowResMassSpectralMatch(Thread):

    def __init__(self, gcms_obj, ref_lib_path, calibration=False):
        
        '''TODO:
        '''
        Thread.__init__(self)
        
        self.gcms_obj = gcms_obj

        #  initiated at create_molecular_database()
        #self.dict_molecular_lookup_table = None
        self.calibration = calibration
        # reading local file for now, 
        self.sqlLite_obj = ReadNistMSI(ref_lib_path).get_sqlLite_obj()

    def metabolite_detector_score(self, gc_peak, ref_obj):

        spectral_similarity_score = self.cosine_correlation(gc_peak.mass_spectrum, ref_obj)

        #print(ref_obj.get('ri'), gc_peak.ri, CompoundSearchSettings.ri_window)

        ri_score = exp(-((gc_peak.ri - ref_obj.get('ri'))**2 ) / (2 * CompoundSearchSettings.ri_window**2))

        similarity_score = ((spectral_similarity_score**2) * (ri_score))**(1/3)

        return spectral_similarity_score, ri_score, similarity_score
        
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

        #correlation = (1 - cosine(x, y))
        
        correlation = dot(x, y)/(norm(x)*norm(y))

        return correlation

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

        return correlation

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

        return correlation

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

        return correlation


    def euclidean_distance(self, mass_spec, ref_obj):

        def euclidean_distance_manual(qlist,rlist):

            T1=sum(subtract(qlist,rlist)**2)

            T2=sum((qlist)**2)

            T3=(1+T1/T2)**(-1)
            
        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Pearson correlation
        correlation = euclidean(df.T[0], df.T[1])

        return correlation

    def manhattan_distance(self, mass_spec, ref_obj):
        
        def mann_distance_manual(qlist,rlist):
        
            T1=sum(absolute(subtract(qlist,rlist)))
            T2=sum(qlist)
            T3=(1+T1/T2)**(-1)

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Pearson correlation
        correlation = cityblock(df.T[0], df.T[1])

        return correlation

    def jaccard_distance(self, mass_spec, ref_obj):

        # create dict['mz'] = abundance, for experimental data
        ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

        # parse to dataframe, easier to zerofilland tranpose
        df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        # calculate Pearson correlation
        correlation = jaccard(df.T[0], df.T[1])

        return correlation

    def run(self):
        
        import tqdm

        if not self.gcms_obj:
            
            self.gcms_obj.process_chromatogram()

        for gc_peak in tqdm.tqdm(self.gcms_obj):
            
            if not self.calibration:
                
                window = CompoundSearchSettings.ri_search_range

                ri = gc_peak.ri

                min_mat_ri = (ri-window, ri+window)    
                
                ref_objs = self.sqlLite_obj.query_min_max_ri(min_mat_ri)
                
            else:

                window = CompoundSearchSettings.rt_search_range

                rt = gc_peak.rt

                min_mat_rt = (rt-window, rt+window)    
                
                ref_objs = self.sqlLite_obj.query_min_max_rt(min_mat_rt)
                
            for ref_obj in ref_objs:
                # uses spectral similarly and uses a threshold to only select peaks with high data correlation
                if self.calibration:
                    
                    correlation_value = self.cosine_correlation(gc_peak.mass_spectrum, ref_obj)

                    #print(w_correlation_value,correlation_value )
                    if correlation_value >= CompoundSearchSettings.correlation_threshold:
                    
                        gc_peak.add_compound(ref_obj, correlation_value)

                # use score, usually a combination of Retention index and Spectral Similarity
                # Threshold is implemented by not necessarily used
                else:

                    # TODO: add other scoring methods
                    # m/q developed methods will be implemented here
                    spectral_similarity_score, ri_score, similarity_score = self.metabolite_detector_score(gc_peak, ref_obj)    
                    
                    if similarity_score >= CompoundSearchSettings.score_threshold:
                    
                        gc_peak.add_compound(ref_obj, spectral_similarity_score, ri_score, similarity_score)
                
        
        self.sqlLite_obj.session.close()
        self.sqlLite_obj.engine.dispose()