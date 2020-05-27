
from threading import Thread
from pathlib import Path

from pandas import DataFrame
from math import exp
from numpy import power
from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite
from corems.molecular_id.calc.SpectralSimilarity import *

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
        spectral_similarity_scores["cosine_correlation"] = cosine_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
        
        if self.gcms_obj.molecular_search_settings.exploratory_mode:
            
            spectral_similarity_scores["weighted_cosine_correlation"] = weighted_cosine_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["stein_scott_similarity"] = stein_scott(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["pearson_correlation"] = pearson_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["spearman_correlation"] = spearman_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["kendall_tau_correlation"] = kendall_tau(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["euclidean_distance"] = euclidean_distance(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["manhattan_distance"] = manhattan_distance(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["jaccard_distance"] = jaccard_distance(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["dft_correlation"] = dft_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)
            spectral_similarity_scores["dwt_correlation"] = dwt_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)

        #print(ref_obj.get('ri'), gc_peak.ri, self.gcms_obj.molecular_search_settings.ri_window)

        ri_score = exp( -1*(power((gc_peak.ri - ref_obj.get('ri')), 2 )  / (2 * power(self.gcms_obj.molecular_search_settings.ri_std, 2)) ))

        similarity_score = ((spectral_similarity_scores.get("cosine_correlation")**2) * (ri_score))**(1/3)

        return spectral_similarity_scores, ri_score, similarity_score
        
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
                    spectral_similarity_scores["cosine_correlation"] = cosine_correlation(gc_peak.mass_spectrum.mz_abun_dict, ref_obj)

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