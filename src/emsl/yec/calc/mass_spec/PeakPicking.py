'''
@author: Yuri E. Corilo
@date: Jun 27, 2019
'''

import threading
import time

from numpy import where, polyfit, poly1d

from emsl.yec.encapsulation.settings.ProcessingSetting import MassSpectrumSetting


class PeakPicking():
    
    def cut_mz_domain_peak_picking(self):
        
        max_picking_mz = MassSpectrumSetting.max_picking_mz
        min_picking_mz =  MassSpectrumSetting.min_picking_mz
        
        final =  where(self.exp_mz  > min_picking_mz)[-1][-1]
        comeco =  where(self.exp_mz  > min_picking_mz)[0][0]
        
        mz_domain_X_low_cutoff, mz_domain_low_Y_cutoff,  = self.exp_mz [comeco:final], self.abundance[comeco:final]
    
        final =  where(self.exp_mz < max_picking_mz)[-1][-1]
        comeco =  where(self.exp_mz < max_picking_mz)[0][0]
           
        if self.frequency_domain.any():
            
            freq_domain_Y_cutoff  = self.frequency_domain[comeco:final] 
            
            return mz_domain_X_low_cutoff[comeco:final], mz_domain_low_Y_cutoff[comeco:final], freq_domain_Y_cutoff[comeco:final]
        
        else:
            
            return mz_domain_X_low_cutoff[comeco:final], mz_domain_low_Y_cutoff[comeco:final], None
    
    def do_peak_picking(self):
            
            
            mz, abudance, freq_exp = self.cut_mz_domain_peak_picking()
            print("Done selecting m/z range")
            
            x = threading.Thread(target=self.calc_centroid, args=(mz, abudance, freq_exp))
            x.start()
            x.join()
            #while x.is_alive():
                # progress bar here so it wont freeze application
                #print( len(self.mspeaks))
            #self.calc_centroid(mz, abudance, freq_exp, abudance_thresould)
            
            print("Done picking peaks")
            
            '''
            print "HELL YES"
            max_peaks_x_y = SpectrumProcess().peakdetect_mob(ymaxcentroid, x_axis = xmaxcentroid, lookahead = 2)
            
            max_peaks_x, max_peaks_y = max_peaks_x_y[0], max_peaks_x_y[1]
            
            del max_peaks_x_y
            '''
    def calculate_resolving_power(self, intes, massa, current_index):
            
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
        
    def calc_centroid(self, massa, intes, freq_exp):    
        
        #this function is too slow, may need slice and apply multi processing, 
        threshold_method = MassSpectrumSetting.threshold_method 
            
        if threshold_method == 'auto':
            print(MassSpectrumSetting.noise_threshold_stds)
            abudance_thresould = self.baselise_noise + (MassSpectrumSetting.noise_threshold_stds * self.baselise_noise_std)
            factor = 1
            
        elif threshold_method == 'signal_noise':
            
            abudance_thresould = MassSpectrumSetting.s2n_threshold
            factor = self.baselise_noise_std
        
        elif threshold_method == "relative_abudance":
        
            abudance_thresould = MassSpectrumSetting.relative_abundace_threshold
            factor = intes.max()/100
        
        else: 
            raise  Exception("%s method was not implemented, please refer to emsl.yec.calc.mass_spec.NoiseCalc Class" % threshold_method)        
        
        time0 = time.time()    
        
        x = 0
        while (len(intes)-1) > x:
                
                        if  intes[x] > intes[x +1] and intes[x] > intes[x - 1] and (intes[x]/factor) > abudance_thresould:#and :
                           
                            z = poly1d(polyfit([massa[x - 1], massa[x], massa[x +1]], [intes[x - 1],intes[x], intes[x +1]], 2))
                            a = z[2]
                            b = z[1]
                            
                            calculated = -b/(2*a)
                            if calculated < 1 or int(calculated) != int(massa[x]):
                                
                                intes_centr = intes[x] 
                                exp_mz_centroid = massa[x]
                            else:   
                                #use the tallest point 
                                intes_centr = intes[x] 
                                exp_mz_centroid = calculated # cria lista de intensidade centroide
                            
                            if freq_exp.any():
                                
                                z = poly1d(polyfit([freq_exp[x - 1],freq_exp[x],freq_exp[x +1]], [intes[x - 1],intes[x],intes[x +1]], 2))
                                a = z[2]
                                b = z[1]
                                
                                calculated_freq = -b/(2*a)
                                
                                if calculated_freq < 1 or int(calculated_freq) != freq_exp[x]:
                                    freq_centr = freq_exp[x]
                                
                                else:    
                                    freq_centr = calculated_freq # cria lista de intensidade centroide
                            else:
                                freq_centr = None
                            
                            peak_resolving_power = self.calculate_resolving_power(intes, massa, x)
                            
                            #parms ion_charge, exp_mz, abundance, resolving_power, signal_to_noise, massspec_index,
                            self.add_mspeak(self.polarity, exp_mz_centroid, intes_centr, peak_resolving_power, intes_centr/self.baselise_noise_std, x, exp_freq=freq_centr)
                            
                            x = x + 1
                        x = x + 1
                    #x = x + 1    
        print( round(time.time() - time0, 2), "seconds to find the peaks")            