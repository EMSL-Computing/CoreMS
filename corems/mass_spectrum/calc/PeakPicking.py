'''
@author: Yuri E. Corilo
@date: Jun 27, 2019
'''

import math
from numpy import hstack, inf, isnan, poly1d, polyfit, where
from corems.encapsulation.constant import Labels

class PeakPicking:

    def cut_mz_domain_peak_picking(self):

        max_picking_mz = self.settings.max_picking_mz
        min_picking_mz = self.settings.min_picking_mz

        min_final = where(self.mz_exp_profile > min_picking_mz)[-1][-1]
        min_comeco = where(self.mz_exp_profile > min_picking_mz)[0][0]

        mz_domain_X_low_cutoff, mz_domain_low_Y_cutoff, = self.mz_exp_profile[min_comeco: min_final], self.abundance_profile[min_comeco: min_final]

        max_final = where(self.mz_exp_profile < max_picking_mz)[-1][-1]
        max_comeco = where(self.mz_exp_profile < max_picking_mz)[0][0]

        if self.has_frequency:

            if self.freq_exp_profile.any():

                freq_domain_low_Y_cutoff = self.freq_exp_profile[min_comeco:min_final]

                return mz_domain_X_low_cutoff[max_comeco:max_final], mz_domain_low_Y_cutoff[max_comeco:max_final], freq_domain_low_Y_cutoff[max_comeco:max_final]

        else:

            return mz_domain_X_low_cutoff[max_comeco:max_final], mz_domain_low_Y_cutoff[max_comeco:max_final], None

    def do_peak_picking(self):

        mz, abudance, freq = self.cut_mz_domain_peak_picking()

        self.mz_exp_profile = mz
        self.abundance_profile = abudance
        self.freq_exp_profile = freq
        
        if self.label == Labels.bruker_frequency or self.label == Labels.midas_frequency:

            self.calc_centroid(mz, abudance, freq)

        elif self.label == Labels.thermo_profile:
            self.calc_centroid(mz, abudance, self.freq_exp_profile)

        elif self.label == Labels.bruker_profile:
            self.calc_centroid(mz, abudance, self.freq_exp_profile)

        elif self.label == Labels.booster_profile:
            self.calc_centroid(mz, abudance, self.freq_exp_profile)

        elif self.label == Labels.simulated_profile:
            self.calc_centroid(mz, abudance, self.freq_exp_profile)

        else: raise Exception("Unknow mass spectrum type", self.label)

    def find_minima(self, apex_index, abundance, len_abundance, right=True):
            
            j = apex_index
            
            if right: minima = abundance[j] > abundance[j+1]
            else: minima = abundance[j] > abundance[j-1]

            while minima:
                
                if j == 1 or j == len_abundance -2:
                    break
                
                if right: 
                    j += 1

                    minima = abundance[j] >= abundance[j+1]

                else: 
                    j -= 1
                    minima = abundance[j] >= abundance[j-1]
            
            if right: return j
            else: return j

    def calculate_resolving_power(self, intes, massa, current_index):
            
            '''this is a conservative calculation of resolving power,
               the peak need to be resolved at least at the half-maximum magnitude,
               otherwise, the combined full width at half maximum is used to calculate resolving power'''

            peak_height = intes[current_index]
            target_peak_height = peak_height/2

            peak_height_minus = peak_height
            peak_height_plus = peak_height

            index_minus = current_index
            while peak_height_minus  >= target_peak_height:

                index_minus = index_minus -1
                peak_height_minus = intes[index_minus]
                #print "massa", "," , "intes", "," , massa[index_minus], "," , peak_height_minus
            x = [ massa[index_minus],  massa[index_minus+1]]
            y = [ intes[index_minus],  intes[index_minus+1]]
            coefficients = polyfit(x, y, 1)

            a = coefficients[0]
            b = coefficients[1]

            y_intercept =  intes[index_minus] + ((intes[index_minus+1] - intes[index_minus])/2)
            massa1 = (y_intercept -b)/a

            index_plus = current_index
            while peak_height_plus  >= target_peak_height:

                index_plus = index_plus + 1
                peak_height_plus = intes[index_plus]
                #print "massa", "," , "intes", "," , massa[index_plus], "," , peak_height_plus

            x = [massa[index_plus],  massa[index_plus - 1]]
            y = [intes[index_plus],  intes[index_plus - 1]]

            coefficients = polyfit(x, y, 1)
            a = coefficients[0]
            b = coefficients[1]

            y_intercept =  intes[index_plus - 1] + ((intes[index_plus] - intes[index_plus - 1])/2)
            massa2 = (y_intercept -b)/a

            if massa1 > massa2:

                resolvingpower =  massa[current_index]/(massa1-massa2)

            else:

                resolvingpower =  massa[current_index]/(massa2-massa1)

            return resolvingpower

    def cal_minima(self, mass, abun):

        abun = -abun

        dy = abun[1:] - abun[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abun))[0]
        
        if indices_nan.size:
            
            abun[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            
            return mass[indexes], abun[indexes]
    
    def calc_centroid(self, mass, abund, freq):
        #TODO: remove peaks that minimum is one data point from the maximum
        # to remove artifacts 

        len_abundance = len(abund)
        
        max_abundance = max(abund)
        
        peak_height_diff = lambda hi, li : ((abund[hi] - abund[li]) / max_abundance )*100

        abundance_threshold, factor = self.get_threshold(abund)
        #print(abundance_threshold, factor)
        # find indices of all peaks
        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        # noise threshold
        if indexes.size and abundance_threshold is not None:
            indexes = indexes[abund[indexes]/factor >= abundance_threshold]
            
        for current_index in indexes: 
            
            if self.label == Labels.simulated_profile: 

                mz_exp_centroid, intes_centr, peak_indexes = self.use_the_max(mass, abund, current_index, len_abundance, peak_height_diff)
                if mz_exp_centroid:
                    
                    peak_resolving_power = self.calculate_resolving_power( abund, mass, current_index)
                    s2n = intes_centr/self.baselise_noise_std
                    freq_centr = None
                    self.add_mspeak(self.polarity, mz_exp_centroid, abund[current_index] , peak_resolving_power, s2n, peak_indexes, exp_freq=freq_centr, ms_parent=self)
            
            else:
            
                mz_exp_centroid, freq_centr, intes_centr, peak_indexes = self.find_apex_fit_quadratic(mass, abund, freq, current_index, len_abundance, peak_height_diff)
                if mz_exp_centroid:
                    
                    peak_resolving_power = self.calculate_resolving_power( abund, mass, current_index)
                    s2n = intes_centr/self.baselise_noise_std
                    self.add_mspeak(self.polarity, mz_exp_centroid, abund[current_index] , peak_resolving_power, s2n, peak_indexes, exp_freq=freq_centr, ms_parent=self)
            
        
    def get_threshold(self, intes):
        
        threshold_method = self.settings.threshold_method

        if threshold_method == 'auto':
            
            #print(self.settings.noise_threshold_std)
            abundance_threshold = self.baselise_noise + (self.settings.noise_threshold_std * self.baselise_noise_std)
            factor = 1

        elif threshold_method == 'signal_noise':

            abundance_threshold = self.settings.s2n_threshold
            factor = self.baselise_noise_std

        elif threshold_method == "relative_abundance":

            abundance_threshold = self.settings.relative_abundance_threshold
            factor = intes.max()/100

        else:
            raise  Exception("%s method was not implemented, please refer to corems.mass_spectrum.calc.NoiseCalc Class" % threshold_method)
        
        return abundance_threshold, factor
        
    def find_apex_fit_quadratic(self, mass, abund, freq, current_index, len_abundance, peak_height_diff):
        
        # calc prominence
        peak_indexes = self.check_prominence(abund, current_index, len_abundance, peak_height_diff )
        
        if not peak_indexes:        
            
            return None, None, None, None           
        
        else:    
            
            # fit parabola to three most abundant datapoints
            list_mass = [mass[current_index - 1], mass[current_index], mass[current_index +1]]
            list_y = [abund[current_index - 1],abund[current_index], abund[current_index +1]]
            
            z = poly1d(polyfit(list_mass, list_y, 2))
            a = z[2]
            b = z[1]

            calculated = -b/(2*a)
            
            if math.isnan(calculated):
                
                mz_exp_centroid = list_mass[1]

            elif calculated < 1 or int(calculated) != int(list_mass[1]):

                mz_exp_centroid = list_mass[1]
            
            else:
                
                mz_exp_centroid = calculated 
            
            if self.label == Labels.bruker_frequency or self.label == Labels.midas_frequency:
                
                # fit parabola to three most abundant frequency datapoints
                list_freq = [freq[current_index - 1], freq[current_index], freq[current_index +1]]
                z = poly1d(polyfit(list_freq, list_y, 2))
                a = z[2]
                b = z[1]

                calculated_freq = -b/(2*a)

                if calculated_freq < 1 or int(calculated_freq) != freq[current_index]:
                    freq_centr = list_freq[1]

                else:
                    freq_centr = calculated_freq
            
            else:
                    freq_centr = None
                    
            return mz_exp_centroid, freq_centr, abund[current_index], peak_indexes
    
    def check_prominence(self, abun, current_index, len_abundance, peak_height_diff ):

        final_index = self.find_minima(current_index, abun, len_abundance, right=True)
            
        start_index = self.find_minima(current_index, abun, len_abundance, right=False)
            
        peak_indexes = (start_index, current_index, final_index)

        if min( peak_height_diff(current_index,start_index), peak_height_diff(current_index,final_index) ) >  self.mspeaks_settings.peak_min_prominence_percent :   
            
            return peak_indexes
        
        else:
            
            return False

    def use_the_max(self, mass, abund, current_index, len_abundance, peak_height_diff):

        peak_indexes = self.check_prominence(abund, current_index, len_abundance, peak_height_diff )
        
        if not peak_indexes:        

            return None, None, None
        
        else:    
            
            return mass[current_index], abund[current_index], peak_indexes


                    