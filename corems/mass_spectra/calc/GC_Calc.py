__author__ = "Yuri E. Corilo"
__date__ = "Feb 13, 2020"

from scipy.signal import savgol_filter
from numpy import hstack, inf, isnan, where, hanning, convolve, blackman, bartlett, ones, r_, sum
from corems.encapsulation.settings.processingSetting import GasChromatographSetting

class GC_Calculations:
    
    '''
    classdocs
    '''
    
    def smooth(self, x):
        
        """smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        output:
            the smoothed signal
            
        see also: 
        
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.savgol_filter
    
        """

        window_len = GasChromatographSetting.smooth_window

        window = GasChromatographSetting.smooth_method

        pol_order = GasChromatographSetting.savgol_pol_order
        
        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")

        if x.size < window_len:
            raise ValueError("Input array needs to be bigger than window size")

        if window_len < 3:
            return x

        if not window in GasChromatographSetting.implemented_smooth_method:
            raise ValueError("Window method should be 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

        s = r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
        
        if window == 'savgol':
            
            return savgol_filter(self.tic, window_len, pol_order)
        
        elif window == 'flat':  # moving average
           
            w = ones(window_len, 'd')
        
        else:
            
            w = eval( window + '(window_len)')

            y = convolve(w / w.sum(), s, mode='valid')

        print( [(window_len / 2 - 1), -(window_len / 2)] )

        return y[int(window_len / 2 - 1):int(-window_len / 2)]
        
    def smooth_tic(self,):
            
            return self.smooth(self.tic)
            '''
            window = GasChromatographSetting.smooth_window
            
            method = GasChromatographSetting.smooth_method

            pol_order = GasChromatographSetting.savgol_pol_order
            
            if GasChromatographSetting.smooth_method == "savgol":
                
                return savgol_filter(self.tic, window, pol_order)
            
            else:
                
                w = eval(method+'(window)')
                
                y = convolve(w / w.sum(), self.tic, mode='valid')

                return y[(window/2 - 1):-(window/2)]
            '''
    
    def peaks_detector(self, mass, abund):

        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            
            return mass[indexes], abund[indexes]
    
    def baseline_detector(self, mass, abund):

        abund = -abund

        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            
            return mass[indexes], abund[indexes]

    def find_nearest_scan(self, rt):

        from numpy import abs as absolute
        from numpy import array

        array_rt = array(self.retention_time)

        scan_index = (absolute(array_rt - rt)).argmin()

        real_scan = self.scans_number[scan_index]

        return real_scan + 1
