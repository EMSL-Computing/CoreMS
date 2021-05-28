

import numpy as np

from pandas import Series, DataFrame

from scipy.signal import savgol_filter, boxcar
from scipy import interpolate
from matplotlib import pyplot as plt
from numpy import abs
from numpy import array, poly1d, polyfit

def peak_detector(tic, max_tic):

    dy = derivate(tic)

    indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]
    
    for index in indexes:
        
        start_index = find_minima(index, tic, right=False)
        final_index = find_minima(index, tic)

        yield (start_index, index, final_index)

def find_nearest_scan(data, nodes):

    array_data = array(nodes)

    scan_index = (abs(array_data - data)).argmin()

    return nodes[scan_index]

def peak_picking_first_derivative(domain, signal, max_height, max_prominence, max_signal, min_peak_datapoints,
                           signal_threshold=0.2, correct_baseline=True, plot_res=False):
    
    if correct_baseline:
        signal = signal - baseline_detector(signal, domain, max_height, max_prominence)

    domain = np.array(domain)
    signal = np.array(signal)

    # dy_signal = derivate(signal)

    # dydy = derivate(dy_signal)
    
    # where returns a tuple of indexes and data type and we only need  the indexes
    # right_indexes = np.where((np.hstack( (0, dy_signal)) < neg_dy_threshold))[0]
    # left_indexes = np.where((np.hstack( (0, dy_signal)) > pos_dy_threshold))[0]
  
    # dy_left_signal = derivate(left_indexes)  
    # dy_left_signal_zero_filled = np.hstack( (dy_left_signal,0))

    # start_peak = []

    # needs a more efficient way of doing that 
    # for index in range(1,len(dy_left_signal_zero_filled)-1):
    #    delta = dy_left_signal_zero_filled[index]
    #    delta_previous = dy_left_signal_zero_filled[index - 1]
    #    delta_next = dy_left_signal_zero_filled[index + 1]
    #    if delta == 1 and delta_previous != 1:
    #        start_peak.append(left_indexes[index]-1)
    # start_peak = array(start_peak)

    dy = derivate(signal)
    apex_indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]
    # min_indexes = np.where((np.hstack((dy, 0)) > 0) & (np.hstack((0, dy)) < 0))[0]

    start_peak = []
    end_peak = []

    pos_dy_threshold = max(dy) * 0.0005
    neg_dy_threshold = min(dy) * 0.0005
    len_dy = len(dy)
    # take apex_index and move left to find start
    for index in apex_indexes:
        # catch for starting position
        if index == 0:
            index_start = index
        else:
            index_start = index - 1

        # catch for ending position
        if index == (len(dy) - 1):
            index_end = index
        else:
            index_end = index + 1

        # while dy[index_start-1] > 0 and index_start != 0:
        while dy[index_start - 1] > pos_dy_threshold and index_start > 0:
            index_start = index_start - 1
        start_peak.append(index_start)

        # while dy[index_end] < 0 and index_end != (len(dy) - 1):
        while dy[index_end] < neg_dy_threshold and index_end != (len_dy - 1):
            index_end = index_end + 1
        end_peak.append(index_end)

    start_peak = array(start_peak)
    end_peak = array(end_peak)

    # left_index = []
    # right_index = []

    last_closest_right = 0

    for apex_index in apex_indexes:

        index_gt_apex = np.where(end_peak >= apex_index)[0]
        index_lt_apex = np.where(start_peak <= apex_index)[0]
        # index_gt_apex = np.where(min_indexes >= apex_index)[0]
        # index_lt_apex = np.where(min_indexes <= apex_index)[0]

        if not index_gt_apex.size == 0 and not index_lt_apex.size == 0:

            closest_right = find_nearest_scan(apex_index, end_peak[index_gt_apex])
            closest_left = find_nearest_scan(apex_index,  start_peak[index_lt_apex])
            # closest_right = find_nearest_scan(apex_index, min_indexes[index_gt_apex])
            # closest_left = find_nearest_scan(apex_index, min_indexes[index_lt_apex])

            x = [closest_left, closest_right]
            y = [signal[closest_left], signal[closest_right]]

            pol = poly1d(polyfit(x, y, 1))

            corrected_peak_height = signal[apex_index] - pol(apex_index)

            if (corrected_peak_height / max_signal) * 100 > signal_threshold:

                if (closest_right - closest_left) >= min_peak_datapoints:

                    if plot_res:

                        # plt.plot(domain[closest_left: closest_right+1], dydy[closest_left:closest_right+1], c='black')
                        plt.plot(domain[closest_left: closest_right+1], dy[closest_left:closest_right+1], c='red')
                        plt.plot(domain[[apex_index, apex_index]], [signal[apex_index], pol(apex_index)], c='red')
                        plt.plot(domain[start_peak], signal[start_peak], c='blue', linewidth='0', marker="^")
                        # plt.plot(domain[[closest_left,apex_index,closest_right]], signal[[closest_left,apex_index,closest_right]], c='blue', linewidth='0', marker="s")    
                        # plt.plot(domain[[closest_left,apex_index,closest_right]], pol([closest_left, apex_index, closest_right]), c='red')
                        plt.plot(domain[closest_left: closest_right+1], signal[closest_left:closest_right+1], c='black')
                        plt.title(str((corrected_peak_height/max_signal)*100))
                        
                        # plt.show()

                    # if not closest_right <= apex_index:
                    #right_index.append(closest_right)
                    #print (closest_right, apex_index)

                    #left_index.append(closest_left)
                    #print (closest_left, apex_index)
                    #plt.show()
                    #if closest_left < last_closest_right:
                    #    closest_left = last_closest_right
                        
                    yield (closest_left, apex_index, closest_right)
                    last_closest_right = closest_right
    
    #plt.plot(domain, dydy, c='black')
    #plt.plot(domain, np.hstack((0, dy)), c='green')
    #plt.plot(domain, [0 for i in range(len(domain))], c='blue')
    #plt.plot(domain[right_index], signal[right_index], c='blue', linewidth='0', marker="s")
    #plt.plot(domain[left_index], signal[left_index], c='yellow', linewidth='0', marker="^")
    #plt.plot(domain[apex_indexes], signal[apex_indexes], c='red', linewidth='0', marker="^")
    #plt.plot(domain[right_indexes], signal[right_indexes], c='green', linewidth='0', marker="^")
    #plt.show()

def find_minima(index, tic, right=True):
            
    j = index
    #apex_abundance = tic[index]
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

    data_array = np.array(data_array)

    dy = data_array[1:] - data_array[:-1]

    '''replaces nan for infinity'''
    indices_nan = np.where(np.isnan(data_array))[0]

    if indices_nan.size:

        data_array[indices_nan] = np.inf
        dy[np.where(np.isnan(dy))[0]] = np.inf

    return dy

def minima_detector(tic, max_tic, peak_height_max_percent, peak_max_prominence_percent):

    peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_tic )*100

    for start_index, index, final_index in peak_detector(tic, max_tic):

        # abundance max threshold    
        if (tic[index] / max_tic) * 100 < peak_height_max_percent:

            # calculates prominence and filter   
            if peak_height_diff(index, start_index) and peak_height_diff(index, final_index) < peak_max_prominence_percent:
                
                    yield from (start_index, final_index)

def baseline_detector(tic, rt, peak_height_max_percent, peak_max_prominence_percent, do_interpolation=True):

    rt = np.array(rt)

    max_tic = max(tic)

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

        return -1 * ynew1

def peak_detector_generator(tic, stds, method, rt, max_height, min_height, max_prominence, min_datapoints):

    max_tic = max(tic)

    if method == 'manual_relative_abundance':

        tic = tic - baseline_detector(tic, rt, max_height, max_prominence)

        norm_tic = (tic / max_tic) * 100

        remove_indexes = np.where(norm_tic < min_height)[0]

        # if self.deconv_rt_list and  self.deconv_mz == 51:
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
        