__author__ = "Yuri E. Corilo"
__date__ = "Feb 13, 2020"

from scipy.signal import savgol_filter, boxcar

from numpy import hstack, inf, isnan, where, hanning, convolve, blackman, bartlett, ones, r_, sum, empty, nan, array, nan_to_num

from pandas import Series
from corems.encapsulation.settings.processingSetting import GasChromatographSetting

class GC_Calculations:
    
    '''
    classdocs
    '''
    
    def smooth(self, x, window_len, window):
        
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
        print(window)
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
            
            pol_order = GasChromatographSetting.savgol_pol_order
            
            return savgol_filter(x, window_len, pol_order)
        
        
        elif window == 'boxcar':  # moving average
            
            print('yes')
            
            w = boxcar(window_len)
            
            y = convolve(w, s, mode='valid')

            
        elif window == 'flat':  # moving average
            
            w = ones(window_len, 'd')

            y = convolve(w / w.sum(), s, mode='valid')
            
        else:
            
            w = eval( window + '(window_len)')

            y = convolve(w / w.sum(), s, mode='valid')

        return y[int(window_len / 2 - 1):int(-window_len / 2)]
        
    def smooth_tic(self, tic):
            
            window_len = GasChromatographSetting.smooth_window

            window = GasChromatographSetting.smooth_method

            return self.smooth(tic, window_len, window)
    
    def peak_detector(self, tic, max_tic):

        def find_minima(index, right=True):
            
            j = index
            tic_len = len(tic)

            if right: minima = tic[j] > tic[j+1]
            else: minima = tic[j] > tic[j-1]

            while minima:
                
                if j == 1 or j == tic_len -2:
                    break
                
                if right: 
                    j += 1

                    minima = tic[j] >= tic[j+1]

                else: 
                    j -= 1
                    minima = tic[j] >= tic[j-1]
                        

            if right: return j
            else: return j

        dy = tic[1:] - tic[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(tic))[0]
        
        if indices_nan.size:
            
            tic[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        for index in indexes:
            
            start_index = find_minima(index, right=False)
            final_index = find_minima(index)

            yield (start_index, index, final_index)

    def centroid_detector(self, tic, max_tic):

        peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_tic )*100

        for start_index, index, final_index in self.peak_detector(tic, max_tic):

             #abundance min threshold        
            if final_index-start_index > GasChromatographSetting.min_peak_datapoints:

                if (tic[index]/max_tic) * 100 > GasChromatographSetting.peak_height_min_percent:
                #if self.retention_time[final_index]-self.retention_time[start_index] < GasChromatographSetting.max_peak_width:
                    
                    #calculates prominence and filter  
                    if  min( peak_height_diff(index,start_index), peak_height_diff(index,final_index) )> GasChromatographSetting.peak_min_prominence_percent :   
                        
                        yield (start_index, index, final_index)

    def minima_detector(self, tic, max_tic):

        peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_tic )*100

        for start_index, index, final_index in self.peak_detector(tic, max_tic):

            #abundance max threshold    
            if (tic[index]/max_tic) * 100 < GasChromatographSetting.peak_height_max_percent:

                    #calculates prominence and filter   
                    if  min(peak_height_diff(index,start_index), peak_height_diff(index,final_index) )< GasChromatographSetting.peak_max_prominence_percent :   
                        
                        yield (start_index, final_index)    

    def baseline_detector(self, tic):
            
        #maximum_abundance = max(tic)
        
        #dy = tic[1:] - tic[:-1]
        
        #indices_nan = where(isnan(tic))[0]
        
        #if indices_nan.size:
            
        #    tic[indices_nan] = inf
        #    dy[where(isnan(dy))[0]] = inf
        
        #indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        indexes = list(i for i in self.minima_detector(tic, max(tic)))
        
        tic = -tic
        
        baseline = empty(len(tic))

        baseline.fill(nan)

        baseline[indexes] = tic[indexes]

        s = Series(baseline).interpolate(method='nearest')
        
        s = nan_to_num(self.smooth(s, 10, 'bartlett'))

        return s
  
    def remove_outliers(self, data):
        
        from numpy import percentile
        q25, q75 = percentile(data, 25), percentile(data, 75)
        iqr = q75 - q25
        print('Percentiles: 25th=%.3f, 75th=%.3f, IQR=%.3f' % (q25, q75, iqr))
        # calculate the outlier cutoff
        cut_off = iqr * 1.5
        lower, upper = q25 - cut_off, q75 + cut_off
        # identify outliers
        outliers = [x for x in data if x < lower or x > upper]
        print('Identified outliers: %d' % len(outliers))
        # remove outliers
        nanfilled_outliers = Series([x if lower <= x <= upper else nan for x in data])

        return nanfilled_outliers
        
   
