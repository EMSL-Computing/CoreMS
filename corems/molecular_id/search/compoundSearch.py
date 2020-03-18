
from threading import Thread
from pathlib import Path

from pandas import DataFrame
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import cosine
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

        cosine_correlation = self.cosine_correlation(gc_peak.mass_spectrum, ref_obj)

        #print(ref_obj.get('ri'), gc_peak.ri, CompoundSearchSettings.ri_window)

        ri_score = exp(-((gc_peak.ri - ref_obj.get('ri'))**2 ) / (2 * CompoundSearchSettings.ri_window))

        score = (cosine_correlation * (ri_score**2))**(1/3)

        return score
        
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
        correlation = (1 - cosine(df.T[0], df.T[1]))

        return correlation

    def run(self):
        
        if not self.gcms_obj:
            
            self.gcms_obj.process_chromatogram()

        for gc_peak in self.gcms_obj:
            
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

                    if correlation_value >= CompoundSearchSettings.correlation_threshold:
                    
                        gc_peak.add_compound(ref_obj, correlation_value)

                # use score, usually a combination of Retention index and Spectral Similarity
                # Threshold is implemented by not necessarily used
                else:

                    # TODO: add other scoring methods
                    # m/q developed methods will be implemented here
                    score_value = self.metabolite_detector_score(gc_peak, ref_obj)    

                    if score_value >= CompoundSearchSettings.score_threshold:
                    
                        gc_peak.add_compound(ref_obj, score_value)
                
        
        self.sqlLite_obj.session.close()
        self.sqlLite_obj.engine.dispose()