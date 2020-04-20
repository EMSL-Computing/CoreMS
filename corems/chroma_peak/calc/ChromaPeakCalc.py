
from numpy import array, polyfit, poly1d, where, trapz

from bisect import bisect_left

__author__ = "Yuri E. Corilo"
__date__ = "March 11, 2020"

class GCPeakCalculation(object):

    '''
    classdocs
    '''
    def calc_area(self, tic, dx):
        
        yy = tic[self.start_index:self.final_index]
        
        self._area = trapz(yy, dx = dx)

    def linear_ri(self, right_ri, left_ri,left_rt,right_rt ):

        #ri = ((right_ri -  left_ri) * (self.rt - left_rt)/ (right_rt - left_rt)) + left_ri   

        return left_ri + ((right_ri -  left_ri) * (self.rt - left_rt)/(right_rt-left_rt))

    def calc_ri(self, rt_ri_pairs):
        current_rt = self.rt

        rts = [rt_ri[0] for rt_ri in rt_ri_pairs]
        index = bisect_left(rts, current_rt)
        
        if index >= len(rt_ri_pairs):
            index -= 1
        
        current_ref = rt_ri_pairs[index]
        
        if current_rt == current_ref[0]:
            #print(current_rt, current_ref)
            self._ri = current_ref[1]
            return 1
        
        else:
            if index == 0:
                index += 1
            
            left_rt = rt_ri_pairs[index-1][0]
            left_ri = rt_ri_pairs[index-1][1]

            right_rt = rt_ri_pairs[index][0]
            right_ri =rt_ri_pairs[index][1]

            self._ri = self.linear_ri(right_ri, left_ri,left_rt,right_rt)

            #print(rt_ri_pairs[index-1], current_rt, rt_ri_pairs[index], left_ri, self._ri, right_ri)

            return 1