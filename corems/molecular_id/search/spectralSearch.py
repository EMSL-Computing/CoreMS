
from threading import Thread
from pathlib import Path

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.encapsulation.settings.processingSetting import GasChromatographSetting

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

    def run(self):
        
        window = GasChromatographSetting.rt_window

        if not self.gcms_obj:
            self.gcms_obj.process_chromatogram()

        for gc_peak in self.gcms_obj:
            
            rt = gc_peak.mass_spectrum.rt

            min_mat_rt = (rt-window, rt+window)    
            
            ref_objs = self.sqlLite_obj.query_min_max_rt(min_mat_rt)
            
            for ref_obj in ref_objs:
                
                print(gc_peak.mass_spectrum.rt, gc_peak.mass_spectrum.tic, ref_obj.get("rt"))

            #ms.cosine_similarity(ref_objs.pair)

