import time
from typing import Tuple

from numpy import where, average, std, isnan, inf, hstack, median, argmax, percentile, log10, histogram, nan
#from scipy.signal import argrelmax
from corems import chunks
import warnings

#from matplotlib import pyplot
__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"


class NoiseThresholdCalc:
    """Class for noise threshold calculation.

    Parameters
    ----------
    mass_spectrum : MassSpectrum
        The mass spectrum object.
    settings : MSParameters
        The mass spectrum parameters object.
    is_centroid : bool
        Flag indicating whether the mass spectrum is centroid or profile.
    baseline_noise : float
        The baseline noise.
    baseline_noise_std : float
        The baseline noise standard deviation.
    max_signal_to_noise : float
        The maximum signal to noise.
    max_abundance : float
        The maximum abundance.
    abundance : np.array
        The abundance array.
    abundance_profile : np.array
        The abundance profile array.
    mz_exp : np.array
        The experimental m/z array.
    mz_exp_profile : np.array
        The experimental m/z profile array.

    Attributes
    ----------
    None

    Methods
    -------
    * get_noise_threshold(). Get the noise threshold.    
    * cut_mz_domain_noise(). Cut the m/z domain to the noise threshold regions.  
    * get_noise_average(ymincentroid). 
        Get the average noise and standard deviation.   
    * get_abundance_minima_centroid(abun_cut)
        Get the abundance minima for centroid data.   
    * run_log_noise_threshold_calc(). 
        Run the log noise threshold calculation.  
    * run_noise_threshold_calc(). 
        Run the noise threshold calculation.  
    """


    def get_noise_threshold(self) -> Tuple[Tuple[float, float], Tuple[float,float ]]:
        """ Get the noise threshold.

        Returns
        -------
        Tuple[Tuple[float, float], Tuple[float, float]]
            A tuple containing the m/z and abundance noise thresholds.
            (min_mz, max_mz), (noise_threshold, noise_threshold)
        """
       
        if self.is_centroid:

            x = min(self.mz_exp), max((self.mz_exp))
            
            if self.settings.noise_threshold_method == 'minima':
                
                abundance_threshold = self.baseline_noise + (self.settings.noise_threshold_min_std * self.baseline_noise_std)
                y = (abundance_threshold, abundance_threshold)

            elif self.settings.noise_threshold_method == 'signal_noise':

                normalized_threshold = (self.max_abundance * self.settings.noise_threshold_min_s2n )/self.max_signal_to_noise
                y = (normalized_threshold, normalized_threshold)
            
            elif self.settings.noise_threshold_method == "relative_abundance":

                normalized_threshold = (max(self.abundance)/100)*self.settings.noise_threshold_min_relative_abundance
                y = (normalized_threshold, normalized_threshold)    

            elif self.settings.noise_threshold_method == "absolute_abundance":

                normalized_threshold = self.abundance*self.settings.noise_threshold_absolute_abundance
                y = (normalized_threshold, normalized_threshold)
            #log noise method not tested for centroid data
            else:
                    raise  Exception("%s method was not implemented, please refer to corems.mass_spectrum.calc.NoiseCalc Class" % self.settings.noise_threshold_method)
                
            return x, y    

        else:

            if self.baseline_noise and self.baseline_noise_std:
                
                x = (self.mz_exp_profile.min(), self.mz_exp_profile.max())
                y = (self.baseline_noise_std, self.baseline_noise_std)
                
                if self.settings.noise_threshold_method == 'minima':
                
                    #print(self.settings.noise_threshold_min_std)
                    abundance_threshold = self.baseline_noise + (self.settings.noise_threshold_min_std * self.baseline_noise_std)
                    
                    y = (abundance_threshold, abundance_threshold)

                elif self.settings.noise_threshold_method == 'signal_noise':

                    max_sn = self.abundance_profile.max()/self.baseline_noise_std

                    normalized_threshold = (self.abundance_profile.max() * self.settings.noise_threshold_min_s2n )/max_sn
                    y = (normalized_threshold, normalized_threshold)

                elif self.settings.noise_threshold_method == "relative_abundance":

                    normalized_threshold = (self.abundance_profile.max()/100)*self.settings.noise_threshold_min_relative_abundance
                    y = (normalized_threshold, normalized_threshold)

                elif self.settings.noise_threshold_method == "absolute_abundance":

                    normalized_threshold = self.settings.noise_threshold_absolute_abundance
                    y = (normalized_threshold, normalized_threshold)

                elif self.settings.noise_threshold_method == "log":
                    normalized_threshold = self.settings.noise_threshold_log_nsigma * self.baseline_noise_std
                    y = (normalized_threshold, normalized_threshold)

                else:
                    raise  Exception("%s method was not implemented, \
                        please refer to corems.mass_spectrum.calc.NoiseCalc Class" % self.settings.noise_threshold_method)
                
                return x, y
            
            else:
                
                warnings.warn(
                    "Noise Baseline and Noise std not specified,\
                    defaulting to 0,0 run process_mass_spec() ?"
                )    
                return (0,0) , (0,0)

    def cut_mz_domain_noise(self):
        """Cut the m/z domain to the noise threshold regions.

        Returns
        -------
        Tuple[np.array, np.array]
            A tuple containing the m/z and abundance arrays of the truncated spectrum region.
        """
        min_mz_whole_ms = self.mz_exp_profile.min()
        max_mz_whole_ms = self.mz_exp_profile.max()

        if self.settings.noise_threshold_method == 'minima':
            
            # this calculation is taking too long (about 2 seconds)
            number_average_molecular_weight = self.weight_average_molecular_weight(
                profile=True)
           
            # +-200 is a guess for testing only, it needs adjustment for each type of analysis
            # need to check min mz here or it will break
            min_mz_noise = number_average_molecular_weight - 100
            # need to check max mz here or it will break
            max_mz_noise = number_average_molecular_weight + 100

        else:

            min_mz_noise = self.settings.noise_min_mz
            max_mz_noise = self.settings.noise_max_mz

        if min_mz_noise < min_mz_whole_ms:
            min_mz_noise = min_mz_whole_ms

        if max_mz_noise > max_mz_whole_ms:
            max_mz_noise = max_mz_whole_ms

        #print(min_mz_noise, max_mz_noise)
        low_mz_index = (where(self.mz_exp_profile >= min_mz_noise)[0][0])
        #print(self.mz_exp_profile[low_mz_index])
        # low_mz_index = (argmax(self.mz_exp_profile <= min_mz_noise))
        
        high_mz_index = (where(self.mz_exp_profile <= max_mz_noise)[-1][-1])
        
        #high_mz_index = (argmax(self.mz_exp_profile <= max_mz_noise))
        
        if high_mz_index > low_mz_index:
            # pyplot.plot(self.mz_exp_profile[low_mz_index:high_mz_index], self.abundance_profile[low_mz_index:high_mz_index])
            # pyplot.show()
            return self.mz_exp_profile[high_mz_index:low_mz_index], self.abundance_profile[low_mz_index:high_mz_index]
        else:
            # pyplot.plot(self.mz_exp_profile[high_mz_index:low_mz_index], self.abundance_profile[high_mz_index:low_mz_index])
            # pyplot.show()
            return self.mz_exp_profile[high_mz_index:low_mz_index], self.abundance_profile[high_mz_index:low_mz_index]
      

    def get_noise_average(self, ymincentroid):
        """ Get the average noise and standard deviation.

        Parameters
        ----------
        ymincentroid : np.array
            The ymincentroid array.
        
        Returns
        -------
        Tuple[float, float]
            A tuple containing the average noise and standard deviation.
            
        """
        # assumes noise to be gaussian and estimate noise level by 
        # calculating the valley. 
        
        auto = True if self.settings.noise_threshold_method == 'minima' else False

        average_noise = median((ymincentroid))*2 if auto else median(ymincentroid)
        
        s_deviation = ymincentroid.std()*3 if auto else ymincentroid.std()
            
        return average_noise, s_deviation

    def get_abundance_minima_centroid(self, abun_cut):
        """Get the abundance minima for centroid data.

        Parameters
        ----------
        abun_cut : np.array
            The abundance cut array.

        Returns
        -------
        np.array
            The abundance minima array.
        """ 
        maximum = self.abundance_profile.max()
        threshold_min = (maximum * 1.00)

        y = -abun_cut

        dy = y[1:] - y[:-1]
        '''replaces NaN for Infinity'''
        indices_nan = where(isnan(y))[0]
        
        if indices_nan.size:

            y[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf

        
        indices = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        if indices.size and threshold_min is not None:
            indices = indices[abun_cut[indices] <= threshold_min]
  
        return abun_cut[indices]

    def run_log_noise_threshold_calc(self):
        """ Run the log noise threshold calculation.


        Returns
        -------
        Tuple[float, float]
            A tuple containing the average noise and standard deviation.
            
        Notes 
        --------
        Method for estimating the noise based on decimal log of all the data point

        Idea is that you calculate a histogram of of the log10(abundance) values. 
        The maximum of the histogram == the standard deviation of the noise. 


        For aFT data it is a gaussian distribution of noise - not implemented here!
        For mFT data it is a Rayleigh distribution, and the value is actually 10^(abu_max)*0.463.


        See the publication cited above for the derivation of this. 

        References
        --------
        1. dx.doi.org/10.1021/ac403278t | Anal. Chem. 2014, 86, 3308−3316

        """

        if self.is_centroid:
            raise  Exception("log noise Not tested for centroid data")
        else:
            # cut the spectrum to ROI
            mz_cut, abundance_cut = self.cut_mz_domain_noise()
            # If there are 0 values, the log will fail
            # But we may have negative values for aFT data, so we check if 0 exists
            # Need to make a copy of the abundance cut values so we dont overwrite it....
            tmp_abundance = abundance_cut.copy()
            if 0 in tmp_abundance:
                tmp_abundance[tmp_abundance==0] = nan
                tmp_abundance = tmp_abundance[~isnan(tmp_abundance)]
                # It seems there are edge cases of sparse but high S/N data where the wrong values may be determined. 
                # Hard to generalise - needs more investigation.

            # calculate a histogram of the log10 of the abundance data
            hist_values = histogram(log10(tmp_abundance),bins=self.settings.noise_threshold_log_nsigma_bins) 
            #find the apex of this histogram
            maxvalidx = where(hist_values[0] == max(hist_values[0]))
            # get the value of this apex (note - still in log10 units)
            log_sigma = hist_values[1][maxvalidx]
            # If the histogram had more than one maximum frequency bin, we need to reduce that to one entry
            if len(log_sigma)>1:
                log_sigma = average(log_sigma)
            ## To do : check if aFT or mFT and adjust method
            noise_mid = 10**log_sigma
            noise_1std = noise_mid*self.settings.noise_threshold_log_nsigma_corr_factor #for mFT 0.463
            return float(noise_mid), float(noise_1std)

    def run_noise_threshold_calc(self):
        """ Runs noise threshold calculation (not log based method)
        
        Returns
        -------
        Tuple[float, float]
            A tuple containing the average noise and standard deviation.

        """
        if self.is_centroid:
            # calculates noise_baseline and noise_std
            # needed to run auto noise threshold mode
            # it is not used for signal to noise nor 
            # relative abudance methods
            abundances_chunks = chunks(self.abundance, 50)
            each_min_abund = [min(x) for x in abundances_chunks]

            return average(each_min_abund), std(each_min_abund)
        
        else:

            mz_cut, abundance_cut = self.cut_mz_domain_noise()
            
            if self.settings.noise_threshold_method == 'minima':

                yminima = self.get_abundance_minima_centroid(abundance_cut)
                
                return self.get_noise_average(yminima)

            else:
                
                # pyplot.show()
                return self.get_noise_average(abundance_cut)
