
from numpy import array, polyfit, poly1d

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

        rts = array(list(dict_ref.keys())) 
        ris = array(list(dict_ref.values())) 
        
        z = polyfit(ris,rts, 1)
        poli = poly1d(z)

        diff = rts - self.rt
        
        i = 0
        
        while diff[i] < 0:
            
            if i == len(diff)-1:
                break 
            if diff[i] > 0:
                break
            else:
                i +=1
        #print(i, len(rts))
        
        if i == 0:
            
            left_ri = dict_ref[rts[0]] - 200 
            left_rt = (poli(left_ri ))
            right_rt = rts[i]

        elif i ==len(rts):

            
            left_rt  = rts[i-1]  
            left_ri = dict_ref[left_rt]    
            right_rt = (poli(dict_ref[rts[-1]] + 200 ))
        
        else:    
            
            left_rt  = rts[i-1]
            left_ri = dict_ref[left_rt]

            right_rt = rts[i]
            #right_ri = dict_ref[right_rt]
        
        self._ri = self.linear_rt(left_ri,left_rt,right_rt)

        #print((left_rt, self.rt, right_rt, left_ri, self._ri, right_ri))