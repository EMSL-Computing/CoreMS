
from threading import Thread
from pathlib import Path

from pandas import DataFrame
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import cosine

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.encapsulation.settings.processingSetting import CompoundSearchSettings, GasChromatographSetting

class LowResMassSpectralMatch(Thread):

    def __init__(self, gcms_obj, ref_lib_path):
        
        '''TODO:
        '''
        Thread.__init__(self)
        
        self.gcms_obj = gcms_obj

        #  initiated at create_molecular_database()
        #self.dict_molecular_lookup_table = None
        
        # reading local file for now, 
        self.sqlLite_obj = ReadNistMSI(ref_lib_path).get_sqlLite_obj()

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
        
        window = CompoundSearchSettings.rt_search_range

        if not self.gcms_obj:
            self.gcms_obj.process_chromatogram()

        for gc_peak in self.gcms_obj:
            
            rt = gc_peak.mass_spectrum.rt

            min_mat_rt = (rt-window, rt+window)    
            
            ref_objs = self.sqlLite_obj.query_min_max_rt(min_mat_rt)
            
            for ref_obj in ref_objs:
            
                correlation_value = self.cosine_correlation(gc_peak.mass_spectrum, ref_obj)
                
                if correlation_value >= CompoundSearchSettings.similarity_threshold:
                    gc_peak.add_compound(ref_obj, correlation_value)
        
        self.sqlLite_obj.session.close()
        self.sqlLite_obj.engine.dispose()