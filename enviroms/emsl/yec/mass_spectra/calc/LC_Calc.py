'''
Created on Jun 14, 2019

@author: eber373
'''

class LC_Calculations(object):
    
    '''
    classdocs
    '''
    
    def find_nearest_scan(self, rt):

        from numpy import abs as absolute
        from numpy import array

        array_rt = array(self.retention_time_list)

        scan_index = (absolute(array_rt - rt)).argmin()

        real_scan = self.scans_number_list[scan_index]

        return real_scan + 1
    