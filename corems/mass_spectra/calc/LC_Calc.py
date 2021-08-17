'''
Created on Jun 14, 2019

@author: eber373
'''

__author__ = "Yuri E. Corilo"
__date__ = "Aug 10, 2021"

from corems.mass_spectra.calc import SignalProcessing as sp

class LC_Calculations:
    
    '''
    classdocs
    '''
    def smooth_tic(self, tic):
            
        implemented_smooth_method = self.chromatogram_settings.implemented_smooth_method
        
        pol_order = self.chromatogram_settings.savgol_pol_order

        window_len = self.chromatogram_settings.smooth_window

        window = self.chromatogram_settings.smooth_method

        return sp.smooth_signal(tic, window_len, window, pol_order, implemented_smooth_method)

    def centroid_detector(self, tic, rt):
        
        noise_std = self.chromatogram_settings.std_noise_threshold

        method = self.chromatogram_settings.noise_threshold_method
        
        #peak picking
        min_height = self.chromatogram_settings.peak_height_min_percent 
        min_datapoints = self.chromatogram_settings.min_peak_datapoints   
        
        # baseline detection
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent 
        max_height = self.chromatogram_settings.peak_height_max_percent 
        
        peak_indexes_generator = sp.peak_detector_generator(tic, noise_std, method, rt, max_height, min_height, max_prominence, min_datapoints)

        return peak_indexes_generator

    def find_nearest_scan(self, rt):

        from numpy import abs as absolute
        from numpy import array

        array_rt = array(self.retention_time)

        scan_index = (absolute(array_rt - rt)).argmin()

        real_scan = self.scans_number[scan_index]

        return real_scan + 1
    