'''
Created on Jun 27, 2019

@author: eber373
'''
from enviroms.emsl.yec.encapsulation.settings.ProcessingSetting import MassSpectrumSetting
from numpy import where,  isnan, inf, hstack, average, std

class NoiseThreshouldCalc():

        def cut_mz_domain_noise(self, auto):

            if auto:

                number_average_molecular_weight = self.weight_average_molecular_weight(profile=True)

                #+-200 is a guess for testing only, it needs adjustment for each type of analysis
                #need to check min mz here or it will break
                min_mz_noise = number_average_molecular_weight - 100
                #need to check max mz here or it will break
                max_mz_noise = number_average_molecular_weight + 100

                min_mz_whole_ms = self.exp_mz.min()
                max_mz_whole_ms = self.exp_mz.max()

                if min_mz_noise < min_mz_whole_ms:
                    min_mz_noise = min_mz_whole_ms

                if max_mz_noise < max_mz_whole_ms:
                    max_mz_noise = max_mz_whole_ms


            else:

                min_mz_noise = MassSpectrumSetting.min_noise_mz
                max_mz_noise =  MassSpectrumSetting.max_noise_mz


            final = where(self.exp_mz > min_mz_noise)[-1][-1]
            comeco = where(self.abundance > min_mz_noise)[0][0]

            mz_domain_low_Y_cutoff = self.abundance[comeco:final]

            final = where(self.exp_mz  < max_mz_noise)[-1][-1]
            comeco = where(self.exp_mz  < max_mz_noise)[0][0]

            return  mz_domain_low_Y_cutoff[comeco:final]

        def get_noise_average(self, ymincentroid):

            average_noise = average(ymincentroid)*2
            s_desviation = std(ymincentroid)*2
            print( "Baseline noise level is %.2f, and the standard deviation is: %.2f" % (average_noise, s_desviation))
            return average_noise, s_desviation
        
        def get_abundance_minima_centroide(self, intes):
        
            maximum = intes.max()
            
            thresould_min = (maximum * 5) / 100
            
            y = -intes
            
            dy = y[1:] - y[:-1]
            
            '''replaces NaN for Infinity'''
            indices_nan = where(isnan(y))[0]
            
            if indices_nan.size:
                
                y[indices_nan] = inf
                dy[where(isnan(dy))[0]] = inf
            
            indices = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
            
            if indices.size and thresould_min is not None:
                indices = indices[y[indices] <= thresould_min]
            
            return intes[indices]
        
        def run_noise_threshould_calc(self, auto):

            Y_cut = self.cut_mz_domain_noise(auto)

            if auto:

                yminima = self.get_abundance_minima_centroide(Y_cut)
                return self.get_noise_average(yminima)

            else:

                return self.get_noise_average(Y_cut)
