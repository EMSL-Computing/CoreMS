

import numpy as np

from pandas import Series, DataFrame

from scipy.signal import savgol_filter, boxcar
from scipy import interpolate


def peak_detector(tic, max_tic):

        dy = derivate(tic)

        indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]
        
        for index in indexes:
            
            start_index = find_minima(index, tic, right=False)
            final_index = find_minima(index, tic)

            yield (start_index, index, final_index)

def find_minima(index, tic, right=True):
            
    j = index
    apex_abundance = tic[index]
    tic_len = len(tic)

    if right: minima = tic[j] >= tic[j+1] 
    else: minima = tic[j] >= tic[j-1] 

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

def derivate(data_array):
        dy = data_array[1:] - data_array[:-1]
    
        '''replaces nan for infinity'''
        indices_nan = np.where(np.isnan(data_array))[0]
        
        if indices_nan.size:
            
            data_array[indices_nan] = np.inf
            dy[where(isnan(dy))[0]] = np.inf
        
        return dy

def minima_detector(tic, max_tic, peak_height_max_percent, peak_max_prominence_percent):

        peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_tic )*100

        for start_index, index, final_index in peak_detector(tic, max_tic):

            #abundance max threshold    
            if (tic[index]/max_tic) * 100 < peak_height_max_percent:

                    #calculates prominence and filter   
                    if  peak_height_diff(index,start_index) and peak_height_diff(index,final_index) < peak_max_prominence_percent :   
                        
                            yield from (start_index, final_index)

def baseline_detector( tic, rt,  peak_height_max_percent, peak_max_prominence_percent, do_interpolation = True):
        
    rt = np.array(rt)

    max_tic =  max(tic)
    
    indexes = sorted(list(set(i for i in minima_detector(tic, max_tic, peak_height_max_percent, peak_max_prominence_percent))))
    
    y = -tic
    
    x1 = rt[indexes]
    
    y1 = y[indexes]

    if len(x1) <= 5:
        return tic

    if not do_interpolation:
        
        y0 = np.zeros(tic.shape)
        y0[indexes] = y[indexes]
        
        return y0
    
    else:    
        
        f1 = interpolate.interp1d(x1, y1, kind='quadratic',fill_value="extrapolate")
        
        ynew1 = f1(list(rt))
        
        # from matplotlib import pyplot as plt   
        # if self.deconv_rt_list and  self.deconv_mz == 51:
            
        #   plt.plot(rt, tic-(-1* ynew1), color='green')

        # plt.plot(rt, -1* ynew1, c='black')
        
        # s = self.smooth(s, 10, 'blackman')

        # plt.plot(self.retention_time, -s)
        
        # plt.show()
        
        return -1* ynew1

def peak_detector_generator(tic, stds, method, rt, max_height, min_height, max_prominence, min_datapoints):
        
    max_tic = max(tic)

    if method=='manual_relative_abundance':
        
        tic = tic - baseline_detector(tic, rt, max_height, max_prominence)
        
        norm_tic = (tic/max_tic)*100
        
        remove_indexes = np.where(norm_tic < min_height)[0]

        #if self.deconv_rt_list and  self.deconv_mz == 51:
        #    plt.plot(self.deconv_rt_list, tic, label=self.deconv_mz)
        
    elif method == 'auto_relative_abundance':
        
        tic = tic - baseline_detector(tic, rt, max_height, max_prominence)

        baseline = baseline_detector(tic, rt, max_height, max_prominence)

        peak_detect_threshold = ((np.nanmean(baseline) + (stds * np.std(baseline))))

        remove_indexes = np.where(tic < peak_detect_threshold)[0]

    elif method == 'second_derivative':
        
        remove_indexes = second_derivative_threshold(tic, stds, rt, max_height, max_prominence)

    else:

        NotImplemented(method)
    
    peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_tic )*100
    
    dy = derivate(tic)

    include_indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]

    final_indexes = sorted(set(include_indexes)-set(remove_indexes))

    #from matplotlib import pyplot as plt   
    
    #plt.plot(self.retention_time, tic, color='black')
    #plt.scatter(self.retention_time[remove_indexes], tic[remove_indexes], color='red')
    #plt.scatter(self.retention_time[include_indexes], tic[include_indexes], color='blue')
    #plt.scatter(self.retention_time[final_indexes], tic[final_indexes], color='blue')
    
    #plt.show()

    for index in final_indexes:
            
        start_index = find_minima(index, tic, right=False)
        final_index = find_minima(index, tic)
        
        if final_index-start_index > min_datapoints:

            #if min( peak_height_diff(index,start_index), peak_height_diff(index,final_index) )> self.chromatogram_settings.peak_min_prominence_percent :   
                
                yield (start_index, index, final_index)

def smooth_signal(x, window_len, window, pol_order, implemented_smooth_method):
    
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
    x= np.array(x)

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input array needs to be bigger than window size")

    #if window_len < 3:
    #    return x

    if not window in implemented_smooth_method:
        raise ValueError("Window method should be 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'savgol':
        
        return savgol_filter(x, window_len, pol_order)


    elif window == 'boxcar':  # moving average
        
        w = boxcar(window_len)
        
        y = np.convolve(w, s, mode='valid')

        
    elif window == 'flat':  # moving average
        
        w = np.ones(window_len, 'd')

        y = np.convolve(w / w.sum(), s, mode='valid')
        
    else:
        
        w = eval( window + '(window_len)')

        y = np.convolve(w / w.sum(), s, mode='valid')

    return y[int(window_len / 2 - 1):int(-window_len / 2)]

def second_derivative_threshold(tic, stds, rt, peak_height_max_percent, peak_max_prominence_percent):
          
    dy = derivate(tic)
    
    dydy = derivate(dy)
    dydy = np.hstack((dydy, 0))
    dydy = np.hstack((0, dydy))

    baseline = baseline_detector(dydy, rt, peak_height_max_percent, peak_max_prominence_percent, do_interpolation=False)
    
    threshold_median = np.median(baseline) - (stds * np.std(baseline))
    
    remove_indexes = np.where(dydy > threshold_median)[0]
    
    return remove_indexes
        