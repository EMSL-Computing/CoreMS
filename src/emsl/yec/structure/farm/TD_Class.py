'''
Created on Jun 12, 2019

TIME DOMAIN CLASS
'''
from os.path import basename, dirname

from numpy import linspace

from emsl.yec.calc.TD_Calc import TransientCalculations

import matplotlib.pyplot as plt
        
__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

fig = plt.figure()

fig.patch.set_facecolor(None)
        
fig.patch.set_alpha(0)


class Transient(TransientCalculations):
    
    def __init__(self, data, d_params):
        
        self._transient_data = data
        
        self.__set__parameters__objects(d_params)
        
        self.__set__transient__time()
        
                
    def __set__parameters__objects(self, d_params):
        
        self.full_filename_path  = d_params.get("filename_path")
        
        self.calibration_terms = (d_params.get("Aterm"), d_params.get("Bterm"), d_params.get("Cterm"))
         
        self._exc_high_freq = d_params.get("exc_high_freq")
        
        self._exc_low_freq = d_params.get("exc_low_freq")
        
        self.bandwidth = d_params.get("bandwidth")
        
        self.number_data_points = d_params.get("number_data_points")
        
        self.polarity = d_params.get("polarity")
        
        self.location = 230
        
    def __set__transient__time(self):
        
        ### needs __set__parameters__ 
        self.transient_time = self.cal_transient_time()
    
    ''''move it to main code and store frequency domain at MS class???'''
    def get_frequency_domain(self, apodization_method, number_of_truncations, number_of_zero_fills ):
        
        if number_of_truncations > 0:
            
            new_time_domain = self.truncation(self._transient_data, number_of_truncations)
            
        else: 
            
            new_time_domain = self._transient_data
            
        if apodization_method != "None":
                
            new_time_domain = self.apodization(new_time_domain, apodization_method)  
        
        #self.plot_transient(self.transient_data)
        
        #self.plot_transient(new_time_domain)
        
        time_domain_y_zero_filled = self.zero_fill(number_of_zero_fills, new_time_domain)     
        
        #self.transient_time = self.transient_time*(number_of_zero_fills+1)
        
        #self.plot_transient(time_domain_y_zero_filled)
        
        frequency_domain, magnitude = self.perform_magniture_mode_ft(time_domain_y_zero_filled, number_of_zero_fills) 
        #creat MS Object here? 
        return frequency_domain, magnitude 
    
    def generate_mass_spec(self, apodization_method, number_of_truncations, number_of_zero_fills):
        
        from emsl.yec.structure.farm.MS_Class import MassSpec
        
        frequency_domain, magnitude = self.get_frequency_domain(apodization_method, number_of_truncations, number_of_zero_fills)
        
        return MassSpec(self.calibration_terms, self.polarity, frequency_domain, magnitude, )
        
    @property
    def get_filename(self):
        
        return basename(self.full_filename_path)
    
    @property
    def get_dir_location(self):
        
        return dirname(self.full_filename_path)
    
    @property
    def get_A_therm(self):
        
        return self.calibration_terms[0]
    
    @property
    def get_B_therm(self):
        
        return self.calibration_terms[1]
    
    @property
    def get_C_therm(self):
        
        return self.calibration_terms[2]
    
    def plot_transient(self, transient_data):
       
        
        self.location +=1
        print( self.location)
        time_axis = linspace(0, self.transient_time, num=len(transient_data))
        plt.subplot(self.location)
        plt.plot(time_axis, transient_data, color='green')
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        #plt.show()    
    
    '''move it to mass_spec_object'''
    def plot_frequency_domain(self):
        
        self.location +=1
        plt.subplot(self.location)
        plt.plot(self.frequency_domain, self.magnitude, color='green')
        plt.xlabel("Hz")
        plt.ylabel("Magnitude")
        #plt.show()    
    
    '''move to ms_class'''
    def plot_mz_domain(self):
        
        self.location +=1
        mz = self.f_to_mz()
        plt.subplot(self.location)
        plt.plot(mz, self.magnitude, color='green')
        plt.xlabel("m/z")
        plt.ylabel("Magnitude")
        plt.show()        
        
            