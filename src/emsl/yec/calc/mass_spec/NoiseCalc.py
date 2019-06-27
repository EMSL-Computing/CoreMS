'''
Created on Jun 27, 2019

@author: eber373
'''
from numpy import where, average, std
from emsl.yec.encapsulation.settings.ProcessingSetting import MassSpectrumSetting

class NoiseThreshouldCalc():        
        
        def cut_mz_domain(self, auto):

            if auto: 
                
                number_average_molecular_weight = self.weight_average_molecular_weight()
                 
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
                
                min_mz_noise = MassSpectrumSetting.max_noise_mz
                max_mz_noise =  MassSpectrumSetting.min_noise_mz
            
            
            final = where(self.exp_mz > min_mz_noise)[-1][-1]
            comeco = where(self.magnitude > min_mz_noise)[0][0]

            mz_domain_low_Y_cutoff = self.magnitude[comeco:final]

            final = where(self.exp_mz  < max_mz_noise)[-1][-1]
            comeco = where(self.exp_mz  < max_mz_noise)[0][0]
            
            return  mz_domain_low_Y_cutoff[comeco:final]
        
        def get_noise_average(self, ymincentroid):
    
            average_noise = average(ymincentroid)*2
            s_desviation = 2*std(ymincentroid)
            print (average_noise, s_desviation)
            return average_noise, s_desviation
    
        def get_magnitude_minima_centroide(self, intes):
    
            maximum = max(intes)
    
            thresould_min = (maximum * 5) / 100
            
            intes_min =  []
            
            x = 0
            while (len(intes) - 1) > x:
                if intes[x - 1] > intes[x] < intes[x + 1] and intes[x] < thresould_min:
                    intes_min.append(intes[x]) # cria lista de intensidade centroide
                    x = x + 1
                x = x + 1
            
            del intes
            
            return intes_min
        
        def run_noise_threshould_calc(self, auto):
            
            Y_cut = self.cut_mz_domain(auto)
        
            if auto:
                
                yminima = self.get_magnitude_minima_centroide(Y_cut)
                
                return self.get_noise_average(yminima)
            
            else: 
                
                return self.get_noise_average(Y_cut) 