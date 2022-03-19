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
    @staticmethod
    def get_max_eic(eic_data: dict):
        
        max_eic = 0
        for eic_data in eic_data.values():

            ind_max_eic = max(eic_data.get('EIC'))
            max_eic = ind_max_eic if ind_max_eic > max_eic else max_eic

        return max_eic
    
    def smooth_tic(self, tic):
            
        implemented_smooth_method = self.chromatogram_settings.implemented_smooth_method
        
        pol_order = self.chromatogram_settings.savgol_pol_order

        window_len = self.chromatogram_settings.smooth_window

        window = self.chromatogram_settings.smooth_method

        return sp.smooth_signal(tic, window_len, window, pol_order, implemented_smooth_method)

    def eic_centroid_detector(self, rt, eic, max_eic):
        
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent

        max_height = self.chromatogram_settings.peak_height_max_percent

        signal_threshold = self.chromatogram_settings.eic_signal_threshold

        min_peak_datapoints = self.chromatogram_settings.min_peak_datapoints

        peak_derivative_threshold = self.chromatogram_settings.peak_derivative_threshold

        correct_baseline = False

        include_indexes = sp.peak_picking_first_derivative(rt, eic, max_height, max_prominence, max_eic,
                                                           min_peak_datapoints, peak_derivative_threshold,
                                                           signal_threshold=signal_threshold,
                                                           correct_baseline=correct_baseline)

        return include_indexes

    def centroid_detector(self, rt, tic):
        
        #noise_std = self.chromatogram_settings.std_noise_threshold

        #method = self.chromatogram_settings.noise_threshold_method
        #peak picking
        #min_height = self.chromatogram_settings.peak_height_min_percent 
        min_peak_datapoints = self.chromatogram_settings.min_peak_datapoints   
        
        peak_derivative_threshold = self.chromatogram_settings.peak_derivative_threshold
        signal_threshold = self.chromatogram_settings.peak_height_min_percent

        # baseline detection
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent 
        max_height = self.chromatogram_settings.peak_height_max_percent 
        
        correct_baseline = False
        
        #peak_indexes_generator = sp.peak_detector_generator(tic, noise_std, method, rt, max_height, min_height, max_prominence, min_datapoints)

        include_indexes = sp.peak_picking_first_derivative(rt, tic, max_height, max_prominence, max(tic),
                                                           min_peak_datapoints, peak_derivative_threshold,
                                                           signal_threshold=signal_threshold,
                                                           correct_baseline=correct_baseline)

        return include_indexes

    def find_nearest_scan(self, rt):

        from numpy import abs as absolute
        from numpy import array

        array_rt = array(self.retention_time)

        scan_index = (absolute(array_rt - rt)).argmin()

        real_scan = self.scans_number[scan_index]

        return real_scan + 1
    