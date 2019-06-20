'''
Created on Jun 19, 2019

@author: eber373
'''

import gc

from numpy import hamming, hanning, blackman, zeros, fft, sqrt, arange, where, \
    power


class TransientCalculations(object):
    
    '''
    classdocs
    '''
    
    def cal_transient_time(self):
        
        return (1 / self.bandwidth) * ((self.number_data_points) / 2)
        
    def zero_fill(self, number_fills, transient ):
        
        zeros_filled_transient = zeros(len(transient)*(number_fills+1))
                
        zeros_filled_transient[0:len(transient)] = transient    
        
        del transient
        
        gc.collect()
        
        return  zeros_filled_transient 
    
    def truncation(self, transient, number_of_truncations):
        
        data_count = len(transient)
            
        for _ in range(number_of_truncations):
        
            data_count = data_count / 2
            
        time_domain_truncated = transient[0:data_count]
        
        del transient
                
        gc.collect()  
        
        return time_domain_truncated
    
    def apodization(self, transient, apodi_method):
        
        length = self.number_data_points
            
        if apodi_method == "Hamming":
                H_function = hamming(length)
        elif apodi_method == "Hanning":
                H_function = hanning(length)
        elif apodi_method == "Blackman":
                H_function = blackman(length)
            
        S_x = transient * H_function
        
        del transient
        gc.collect()  
            
        return S_x
    
    def calculate_frequency_domain(self, number_data_points):
        
        qntpoints = arange(0,(number_data_points))
        
        factor_distancy = (self.bandwidth)/(number_data_points)  
                
        frequency_domain = qntpoints * factor_distancy
        
        del qntpoints   
        del factor_distancy
        gc.collect()  
        
        return frequency_domain  
    
    def cut_freq_domain(self, freqdomain_X, freqdomain_Y):
      
        if self.exc_low_freq > self.exc_high_freq:
            
            final =  where(freqdomain_X > self.exc_high_freq)[-1][-1]
            comeco =  where(freqdomain_X > self.exc_high_freq())[0][0]
        
        else:
            
            final =  where(freqdomain_X > self.exc_low_freq)[-1][-1]
            comeco =  where(freqdomain_X > self.exc_low_freq)[0][0]
            
        
        return freqdomain_X[comeco:final], freqdomain_Y[comeco:final]
        del freqdomain_X, freqdomain_Y
        gc.collect()
    
    def perform_magniture_mode_ft(self, transient, number_of_zero_fills):
        
        A = fft.fft(transient)
        
        factor = int(number_of_zero_fills+1)
        Max_index = int(len(A)/factor)
        
        A = A[0:Max_index]
        
        datapoints = len(A)
        
        freqdomain_X = self.calculate_frequency_domain(datapoints)
        
        magnitude_Y = sqrt((A.real**2) + (A.imag**2))
        
        freqdomain_X_cut, magnitude_Y_cut = self.cut_freq_domain(freqdomain_X, magnitude_Y)  
        
        del transient 
        #del freqdomain_X
        #del magnitude_Y
        gc.collect()
        
        return freqdomain_X_cut, magnitude_Y_cut
    
    def correct_dc_offset(self):
        
        pass
    
    def f_to_mz(self):
        
        print( self.Atherm, self.BTherm, self.frequency_domain)
        m_z = (self.Atherm/ self.frequency_domain ) + (self.BTherm / power(self.frequency_domain, 2))
       
        return m_z 
    
    