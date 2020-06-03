__author__ = "Yuri E. Corilo"
__date__ = "Feb 13, 2020"

from scipy.signal import savgol_filter, boxcar
from numpy import hstack, inf, isnan, where, hanning, convolve, blackman, bartlett, ones, r_, sum, empty, nan, array, nan_to_num, std, nanmean, median, average, zeros
from pandas import Series


class GC_Calculations:
    
    '''
    classdocs
    '''
   
    def lowres_deconvolution(self):
        from matplotlib import pyplot as plt

        deconvolve_dict = {} 
        maximum_tic = max(self._processed_tic)
        
        for scan_number, ms_obj in self._ms.items():

            mz_list = ms_obj.mz_exp
            abundance_list = ms_obj.abundance
            #add list of scan numbers
            for index, mz in enumerate(mz_list):

                #dict of mz and tuple (mass spectrum abundances index, and scan number )
                if not mz in deconvolve_dict.keys(): 
                   
                    deconvolve_dict[mz] = [ [abundance_list[index]], [scan_number], [index], [ms_obj.rt] ]
                
                else:
                    
                    deconvolve_dict[mz][0].append(ms_obj.abundance[index])
                    deconvolve_dict[mz][1].append(scan_number)
                    deconvolve_dict[mz][2].append(index)
                    deconvolve_dict[mz][3].append(ms_obj.rt)

        dict_deconvolution_data = {}

        for mz, eic_scan_index_rt in deconvolve_dict.items():

            eic = eic_scan_index_rt[0]
            rt_list = eic_scan_index_rt[3]
            
            self.deconv_max_tic = self._ms[index]._processed_tic
            self.deconv_rt_list = rt_list
            self.deconv_mz = mz

            if len(eic) > 5:

                #print(len(eic), len(rt_list))    
                
                smooth_eic = self.smooth_tic(eic)

                include_indexes = list(self.centroid_detector(smooth_eic, rt_list))
                
                norm_smooth_eic = (smooth_eic/maximum_tic)*100
                
                remove_indexes = where(smooth_eic < self.chromatogram_settings.peak_height_min_abun)[0]
                final_indexes = sorted(set(include_indexes)-set(remove_indexes))
                #print(len(include_indexes), len(remove_indexes), len(final_indexes))
                
                rt_apex = [rt_list[i[1]] for i in final_indexes]
                abundance_apex = [smooth_eic[i[1]] for i in final_indexes]

                for peak_index, deconvoled_rt in enumerate(rt_apex):
                    
                    deconvoled_abundance = abundance_apex[peak_index]
                    
                    if not deconvoled_rt in dict_deconvolution_data.keys():
                        dict_deconvolution_data[deconvoled_rt] = [ [mz], [deconvoled_abundance], [scan_number] ]
                    else:
                        
                        dict_deconvolution_data[deconvoled_rt][0].append(mz)
                        dict_deconvolution_data[deconvoled_rt][1].append(deconvoled_abundance)  
                        dict_deconvolution_data[deconvoled_rt][2].append(scan_number)    

                #if rt_apex and mz == 51:
                    
                    #plt.plot(rt_list, smooth_eic, label=mz)
                    #plt.plot(rt_apex, abundance_apex, linewidth= 0, c='black', marker= '^')
                    #plt.legend()

        print(len(dict_deconvolution_data.keys()))
        c = 'green'
        
        i = 0
        for rt, data in dict_deconvolution_data.items():
            
            norm_smooth_tic = (sum(data[1])/maximum_tic)*100

            if norm_smooth_tic > self.chromatogram_settings.peak_height_min_abun and len(data[1]) > 3:
                i += 1
                ax = plt.gca()

                markerline_a, stemlines_a, baseline_a  = ax.stem(data[0], data[1], linefmt='-',  markerfmt=" ", use_line_collection =True, label=rt)
                
                plt.setp(markerline_a, 'color', c, 'linewidth', 2)
                plt.setp(stemlines_a, 'color', c, 'linewidth', 2)
                plt.setp(baseline_a, 'color', c, 'linewidth', 2)

                ax.set_xlabel("$\t{m/z}$", fontsize=12)
                ax.set_ylabel('Abundance', fontsize=12)
                ax.tick_params(axis='both', which='major', labelsize=12)

                ax.axes.spines['top'].set_visible(False)
                ax.axes.spines['right'].set_visible(False)

                ax.get_yaxis().set_visible(False)
                ax.spines['left'].set_visible(False)
                #plt.plot(rt_apex, abundance_apex, linewidth= 0, c='black', marker= '^')
                plt.legend()
                plt.show()
            #plt.close()
        print(i)

    def calibrate_ri(self, ref_dict):
        
        if not self:
            
            self.process_chromatogram()

        for gcms_peak in self:
        
            gcms_peak.calc_ri(ref_dict)

        self.ri_pairs_ref = ref_dict
        
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
        x= array(x)

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")

        if x.size < window_len:
            raise ValueError("Input array needs to be bigger than window size")

        #if window_len < 3:
        #    return x

        if not window in self.chromatogram_settings.implemented_smooth_method:
            raise ValueError("Window method should be 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

        s = r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
        
        if window == 'savgol':
            
            pol_order = self.chromatogram_settings.savgol_pol_order
            
            return savgol_filter(x, window_len, pol_order)
        
        
        elif window == 'boxcar':  # moving average
            
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
            
            window_len = self.chromatogram_settings.smooth_window

            window = self.chromatogram_settings.smooth_method

            return self.smooth(tic, window_len, window)
    
    @staticmethod
    def derivate(data_array):
            dy = data_array[1:] - data_array[:-1]
        
            '''replaces nan for infinity'''
            indices_nan = where(isnan(data_array))[0]
            
            if indices_nan.size:
                
                data_array[indices_nan] = inf
                dy[where(isnan(dy))[0]] = inf
            
            return dy
    
    @staticmethod
    def find_minima(index, tic, right=True):
            
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

    def second_derivative_threshold(self, tic, stds, rt):
          
        dy = self.derivate(tic)
        
        dydy = self.derivate(dy)
        dydy = hstack((dydy, 0))
        dydy = hstack((0, dydy))

        baseline = self.baseline_detector(dydy, rt, do_interpolation=False)
        
        threshold_median = median(baseline) - (stds * std(baseline))
        
        #threshold_average = average(baseline) - (stds * std(baseline))
        
        remove_indexes = where(dydy > threshold_median)[0]
        
        #sorted(set(include_indexes)-set(remove_indexes))
        #from matplotlib import pyplot as plt   
        
        #plt.plot(self.retention_time, dydy, color='red')
        #plt.plot(self.retention_time, [threshold_median]*len(tic), color='red')
        
        #plt.scatter(self.retention_time[remove_indexes], dydy[remove_indexes], color='green')
        
        return remove_indexes
        #plt.plot(self.retention_time, hstack((dy, 0)), color='green', label="First Derivative")
        #ax = plt.gca()

        #ax.plot(self.retention_time, dydy, color='red', label="Second Derivative")
        
        #ax.plot(self.retention_time, [median(baseline)]*len(threshold_average), color='blue', label="Baseline")

        #ax.plot(self.retention_time, threshold_median, color='black', label="Auto Second Derivative Threshold")
        
        #ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')
        #plt.legend()
        #plt.xlim(27.4, 27.7)
        #plt.ylim(-50000, 300000)
        
        
    def peak_detector(self, tic, max_tic):

        dy = self.derivate(tic)

        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        for index in indexes:
            
            start_index = self.find_minima(index, tic, right=False)
            final_index = self.find_minima(index, tic)

            yield (start_index, index, final_index)

    def peak_detector_generator(self, tic, stds, method, rt):
        
        if method=='manual_relative_abundance':
            
            tic = tic - self.baseline_detector(tic, rt)
            
            if not self.deconv_max_tic:
                max_height = max(tic)
            else:
                max_height = self.deconv_max_tic

            norm_tic = (tic/max_height)*100
            
            remove_indexes = where(norm_tic < self.chromatogram_settings.peak_height_min_abun)[0]

            #if self.deconv_rt_list and  self.deconv_mz == 51:
            #    plt.plot(self.deconv_rt_list, tic, label=self.deconv_mz)

            #plt.plot(self.retention_time, hstack((dy, 0)), color='green', label="First Derivative")
            
        elif method == 'auto_relative_abundance':
            
            tic = tic - self.baseline_detector(tic, rt)

            if not self.deconv_max_tic:
                max_height = max(tic)
            else:
                max_height = self.deconv_max_tic

            baseline = self.baseline_detector(tic, rt)

            peak_detect_threshold = ((nanmean(baseline) + (stds * std(baseline))))

            remove_indexes = where(tic < peak_detect_threshold)[0]

        elif method == 'second_derivative':
            
            if not self.deconv_max_tic:
                max_height = max(tic)
            else:
                max_height = self.deconv_max_tic

            remove_indexes = self.second_derivative_threshold(tic, stds, rt)

        else:

            NotImplemented(method)
        
        peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_height )*100
        
        dy = self.derivate(tic)

        include_indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        final_indexes = sorted(set(include_indexes)-set(remove_indexes))

        #from matplotlib import pyplot as plt   
        
        #plt.plot(self.retention_time, tic, color='black')
        #plt.scatter(self.retention_time[remove_indexes], tic[remove_indexes], color='red')
        #plt.scatter(self.retention_time[include_indexes], tic[include_indexes], color='blue')
        #plt.scatter(self.retention_time[final_indexes], tic[final_indexes], color='blue')
        
        #plt.show()

        for index in final_indexes:
                
            start_index = self.find_minima(index, tic, right=False)
            final_index = self.find_minima(index, tic)
            
            if final_index-start_index > self.chromatogram_settings.min_peak_datapoints:

               #if min( peak_height_diff(index,start_index), peak_height_diff(index,final_index) )> self.chromatogram_settings.peak_min_prominence_percent :   
                    
                    yield (start_index, index, final_index)

    def centroid_detector(self, tic, rt):
        
        stds = self.chromatogram_settings.std_noise_threshold

        method = self.chromatogram_settings.noise_threshold_method
            
        peak_indexes_generator = self.peak_detector_generator(tic, stds, method, rt)

        return peak_indexes_generator

    def minima_detector(self, tic, max_tic):

        peak_height_diff = lambda hi, li : ((tic[hi] - tic[li]) / max_tic )*100

        for start_index, index, final_index in self.peak_detector(tic, max_tic):

            #abundance max threshold    
            if (tic[index]/max_tic) * 100 < self.chromatogram_settings.peak_height_max_percent:

                    #calculates prominence and filter   
                    if  peak_height_diff(index,start_index) and peak_height_diff(index,final_index) < self.chromatogram_settings.peak_max_prominence_percent :   
                        
                            yield from (start_index, final_index)

    def baseline_detector(self, tic, rt, do_interpolation=True):
        
        #from matplotlib import pyplot as plt   

        from scipy import interpolate
        
        rt = array(rt)

        if not self.deconv_max_tic:
           max_tic =  max(tic)
        else:
            max_tic = self.deconv_max_tic

        indexes = sorted(list(set(i for i in self.minima_detector(tic, max_tic))))
        
        y = -tic
        
        x1 = rt[indexes]
        
        y1 = y[indexes]

        if len(x1) <= 5:
            return tic

        if not do_interpolation:
            
            y0 = zeros(tic.shape)
            y0[indexes] = y[indexes]
            
            return y0
        
        else:    
            
            f1 = interpolate.interp1d(x1, y1, kind='quadratic',fill_value="extrapolate")
            
            ynew1 = f1(list(rt))
            
            #if self.deconv_rt_list and  self.deconv_mz == 51:
                
            #    plt.plot(rt, tic-(-1* ynew1), color='green')

            #plt.plot(rt, -1* ynew1, c='black')
            
            #s = self.smooth(s, 10, 'blackman')

            #plt.plot(self.retention_time, -s)
            
            #plt.show()
            
            return -1* ynew1
  
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
        
   
