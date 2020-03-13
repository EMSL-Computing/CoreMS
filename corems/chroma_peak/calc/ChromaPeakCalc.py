
from numpy import array, polyfit, poly1d, where

from corems.encapsulation.settings.processingSetting import CompoundSearchSettings

__author__ = "Yuri E. Corilo"
__date__ = "March 11, 2020"

class GCPeakCalculation(object):

    '''
    classdocs
    '''
    def linear_rt(self,  left_ri,left_rt,right_rt ):

        return left_ri + (CompoundSearchSettings.ri_spacing * (self.rt - left_rt)/(right_rt-left_rt))

    def calc_ri(self, dict_ref):
        # find left and right peaks, gets retention time for both, e calculates retention index 
        #dict_ref has to be comprehensive enough to cover all rt range or it will fail
        
        rts = array(list(dict_ref.keys())) 
        
        right_peak_index = where(rts >= self.rt)[0][0]
        left_peak_index = where(rts < self.rt)[0][-1]

        left_rt  = rts[left_peak_index]
        left_ri = dict_ref[left_rt]

        right_rt = rts[right_peak_index]
            
        self._ri = self.linear_rt(left_ri,left_rt,right_rt)
