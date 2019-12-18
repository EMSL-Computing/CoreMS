'''
@author: Yuri E. Corilo
@date: Jun 27, 2019
'''


from numpy import hstack, inf, isnan, poly1d, polyfit, where

from corems.encapsulation.settings.input.ProcessingSetting import MassSpectrumSetting
from corems.encapsulation.constant import Labels

class PeakPicking(object):

    def cut_mz_domain_peak_picking(self):

        max_picking_mz = MassSpectrumSetting.max_picking_mz
        min_picking_mz =  MassSpectrumSetting.min_picking_mz
        
        final =  where(self.mz_exp_profile  > min_picking_mz)[-1][-1]
        comeco =  where(self.mz_exp_profile  > min_picking_mz)[0][0]

        mz_domain_X_low_cutoff, mz_domain_low_Y_cutoff,  = self.mz_exp_profile [comeco:final], self.abundance_profile[comeco:final]

        final =  where(self.mz_exp_profile < max_picking_mz)[-1][-1]
        comeco =  where(self.mz_exp_profile < max_picking_mz)[0][0]

        if self.freq_exp_profile.any():

            freq_domain_Y_cutoff  = self.freq_exp_profile[comeco:final]

            return mz_domain_X_low_cutoff[comeco:final], mz_domain_low_Y_cutoff[comeco:final], freq_domain_Y_cutoff[comeco:final]

        else:

            return mz_domain_X_low_cutoff[comeco:final], mz_domain_low_Y_cutoff[comeco:final], None

    def do_peak_picking(self):

            if self.label == Labels.bruker_frequency or self.label == Labels.midas_frequency:
                
                mz, abudance, freq = self.cut_mz_domain_peak_picking()
                self.calc_centroid(mz, abudance, freq)
            
            elif self.label == Labels.thermo_profile:
                self.calc_centroid(self.mz_exp_profile, self.abundance_profile, self.freq_exp_profile)
            
            elif self.label == Labels.bruker_profile:
                self.calc_centroid(self.mz_exp_profile, self.abundance_profile, self.freq_exp_profile)
            
            elif self.label == Labels.booster_profile:
                self.calc_centroid(self.mz_exp_profile, self.abundance_profile, self.freq_exp_profile)

            elif self.label == Labels.simulated_profile:
                self.calc_centroid(self.mz_exp_profile, self.abundance_profile, self.freq_exp_profile)
            
           
            else: raise Exception("Unknow mass spectrum type", self.label)
            
           
    
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

            y_intercep =  intes[index_minus] + ((intes[index_minus+1] - intes[index_minus])/2)
            massa1 = (y_intercep -b)/a

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

            y_intercep =  intes[index_plus - 1] + ((intes[index_plus] - intes[index_plus - 1])/2)
            massa2 = (y_intercep -b)/a

            if massa1 > massa2:

                resolvingpower =  massa[current_index]/(massa1-massa2)

            else:

                resolvingpower =  massa[current_index]/(massa2-massa1)

            return resolvingpower

    def cal_minima(self, mass, abund):

        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            
            return mass[indexes], abund[indexes]
    
    def calc_centroid(self, mass, abund, freq):
        #TODO: remove peaks that minimum is one data point from the maximum
        # to remove artifacts 

        abundance_threshold, factor = self.get_threshold(abund)
        # find indices of all peaks
        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        if indexes.size and abundance_threshold is not None:
            indexes = indexes[abund[indexes]/factor >= abundance_threshold]
        
        for current_index in indexes: 
            
            if self.label == Labels.bruker_frequency:                                                               
                
                mz_exp_centroid, freq_centr, intes_centr = self.find_apex_fit_quadratic(mass, abund, freq, current_index)
                peak_resolving_power = self.calculate_resolving_power( abund, mass, current_index)
                s2n = intes_centr/self.baselise_noise_std
            
            elif self.label == Labels.bruker_profile or self.label == Labels.booster_profile or self.label.thermo_profile: 

                peak_resolving_power = self.calculate_resolving_power( abund, mass, current_index)
                mz_exp_centroid, freq_centr, intes_centr = self.find_apex_fit_quadratic(mass, abund, freq, current_index)
                s2n = intes_centr/self.baselise_noise_std

            elif self.label == Labels.simulated_profile: 

                peak_resolving_power = self.calculate_resolving_power( abund, mass, current_index)
                mz_exp_centroid, intes_centr = self.use_the_max(mass, abund, current_index)
                s2n = intes_centr/self.baselise_noise_std
                freq_centr = None

            else: raise Exception("Label '%s' not recognized inside : %s" % (self.label, self.__str__()))
            
            self.add_mspeak(self.polarity, mz_exp_centroid, abund[current_index] , peak_resolving_power, s2n, current_index, exp_freq=freq_centr)
        
    def get_threshold(self, intes):
        
        threshold_method = MassSpectrumSetting.threshold_method

        if threshold_method == 'auto':
            
            #print(MassSpectrumSetting.noise_threshold_stds)
            abundance_threshold = self.baselise_noise + (MassSpectrumSetting.noise_threshold_stds * self.baselise_noise_std)
            factor = 1

        elif threshold_method == 'signal_noise':

            abundance_threshold = MassSpectrumSetting.s2n_threshold
            factor = self.baselise_noise_std

        elif threshold_method == "relative_abundance":

            abundance_threshold = MassSpectrumSetting.relative_abundance_threshold
            factor = intes.max()/100

        else:
            raise  Exception("%s method was not implemented, please refer to corems.mass_spectrum.calc.NoiseCalc Class" % threshold_method)
        
        return abundance_threshold, factor
        
    def find_apex_fit_quadratic(self, mass, abund, freq, current_index):
        
        list_mass = [mass[current_index - 1], mass[current_index], mass[current_index +1]]
        list_y = [abund[current_index - 1],abund[current_index], abund[current_index +1]]
        
        z = poly1d(polyfit(list_mass, list_y, 2))
        a = z[2]
        b = z[1]

        calculated = -b/(2*a)
        
        if calculated < 1 or int(calculated) != int(list_mass[1]):

            mz_exp_centroid = list_mass[1]
        
        else:
            
            mz_exp_centroid = calculated 
        
        if self.label == Labels.bruker_frequency or self.label == Labels.midas_frequency:
            
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
                
        return mz_exp_centroid, freq_centr, abund[current_index]
    
    def use_the_max(self, mass, abund, current_index):
        
            return mass[current_index], abund[current_index]

    def old_calc_centroid(self, massa, intes, freq_exp): #pragma: no cover

        #this function is too slow, may need slice and apply multi processing,
        abundance_threshold, factor = self.get_threshold(intes)
        
        do_freq = freq_exp.any()
        
        for x in range(len(intes)-1):
            if (intes[x]/factor) > abundance_threshold:
                
                if  intes[x] > intes[x +1] and intes[x] > intes[x - 1]:# and (intes[x]/factor) > abundance_threshold:#and :
                    intes_centr = intes[x]
                    mz_exp_centroid = massa[x]
                    
                    if do_freq:
                        freq_centr = freq_exp[x]
                    else:
                        freq_centr = None

                    peak_resolving_power = self.calculate_resolving_power(intes, massa, x)

                    #parms ion_charge, mz_exp, abundance, resolving_power, signal_to_noise, massspec_index,
                    self.add_mspeak(self.polarity, mz_exp_centroid, intes_centr, peak_resolving_power, intes_centr/self.baselise_noise_std, x, exp_freq=freq_centr)

    def calc_min(self, massa, intes, freq_exp): #pragma: no cover

        #this function is too slow, may need slice and apply multi processing,
        
        min_mz = []
        min_abu = list()
        for x in range(len(intes)-1):
                
            if  intes[x] < intes[x +1] and intes[x] < intes[x - 1]:# and (intes[x]/factor) > abundance_threshold:#and :
                intes_centr = intes[x]
                mz_exp_centroid = massa[x]
                
                min_mz.append(mz_exp_centroid)
                min_abu.append(intes_centr)
        
        return min_mz, min_abu 

                    