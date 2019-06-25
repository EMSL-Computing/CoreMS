'''
Created on Jun 12, 2019

TIME DOMAIN CLASS
'''
from os.path import basename, dirname

from numpy import linspace

from emsl.yec.calc.TransientCalc import TransientCalculations
from emsl.yec.structure.factory.MassSpectrumClasses import MassSpecBase
import matplotlib.pyplot as plt


__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

fig = plt.figure()

fig.patch.set_facecolor(None)
        
fig.patch.set_alpha(0)


class Transient(TransientCalculations):
    
    def __init__(self, data, d_params):
        
        self._transient_data = data
        
        self.d_params = d_params
        
        self.frequency_domain = None 
        
        self.magnitude = None
        
        self.__set__parameters__objects(d_params)
        
        self.__set__transient__time()
        
                
    def __set__parameters__objects(self, d_params):
        
        self._full_filename_path  = d_params.get("filename_path")
        
        self.calibration_terms = (d_params.get("Aterm"), d_params.get("Bterm"), d_params.get("Cterm"))
         
        self._exc_high_freq = d_params.get("exc_high_freq")
        
        self._exc_low_freq = d_params.get("exc_low_freq")
        
        self.bandwidth = d_params.get("bandwidth")
        
        self.number_data_points = d_params.get("number_data_points")
        
        self.polarity = d_params.get("polarity")
        
        self.location = 220
        
        self.apodization_method = None
        
        self.number_of_truncations = None
        
        self.number_of_zero_fills = None
        
        
    def __set__transient__time(self):
        
        ### needs __set__parameters__ 
        self.transient_time = self.cal_transient_time()
    
    
    def set_processing_parameter(self, apodization_method, number_of_truncations, number_of_zero_fills):
        
        self.apodization_method = apodization_method
        
        self.number_of_truncations = number_of_truncations
        
        self.number_of_zero_fills = number_of_zero_fills
        
        
    def get_frequency_domain(self,  plot_result=True):
        
        if not self.apodization_method and self.number_of_truncations and self.number_of_zero_fills:
            
            raise Exception("you need to call set_processing parameter before fFt")
        
        else:
            
            if self.number_of_truncations > 0:
                
                new_time_domain = self.truncation(self._transient_data, self.number_of_truncations)
                
            else: 
                
                new_time_domain = self._transient_data
                
            if self.apodization_method != None:
                    
                new_time_domain = self.apodization(new_time_domain, self.apodization_method)  
            
            if plot_result:
                
                self.plot_transient(self._transient_data)
            
                self.plot_transient(new_time_domain)
            
            time_domain_y_zero_filled = self.zero_fill(self.number_of_zero_fills, new_time_domain)     
            
            self.transient_time = self.transient_time*(self.number_of_zero_fills+1)
            
            if plot_result:
            
                self.plot_transient(time_domain_y_zero_filled)
            
            return self.perform_magniture_mode_ft(time_domain_y_zero_filled, self.number_of_zero_fills) 
        #return frequency_domain, magnitude 
    
    def generate_mass_spec(self, plot_result=True):
        
       
        frequency_domain, magnitude = self.get_frequency_domain()
        
        if plot_result:
            
            self._plot_frequency_domain(frequency_domain, magnitude)
        
        return MassSpecBase(frequency_domain, magnitude , self.d_params)
        
    
    @property
    def filename(self):
        
        return basename(self._full_filename_path)
    
    @property
    def dir_location(self):
        
        return dirname(self._full_filename_path)
    
    @property
    def A_therm(self):
        
        return self.calibration_terms[0]
    
    @property
    def B_therm(self):
        
        return self.calibration_terms[1]
    
    @property
    def C_therm(self):
        
        return self.calibration_terms[2]
    
        
    def plot_transient(self, transient_data):
       
        #self.location +=1
        #print( self.location)
        time_axis = linspace(0, self.transient_time, num=len(transient_data))
        #plt.subplot(self.location)
        plt.plot(time_axis, transient_data, color='green')
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        plt.show()    
    
    def plot_zerofilled_transient(self):
       
        new_time_domain = self.apodization(self.transient_time, self.apodization_method)
        time_domain_y_zero_filled = self.zero_fill(self.number_of_zero_fills, new_time_domain)     
        self.transient_time = self.transient_time*(self.number_of_zero_fills+1)  
        time_axis = linspace(0, time_domain_y_zero_filled, num=len(time_domain_y_zero_filled))
        #plt.subplot(self.location)
        plt.plot(time_axis, new_time_domain, color='green')
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        plt.show()   
        
    def plot_apodized_transient(self):
       
        #self.location +=1
        #print( self.location)
        new_time_domain = self.apodization(self.transient_time, self.apodization_method)  
        time_axis = linspace(0, new_time_domain, num=len(new_time_domain))
        #plt.subplot(self.location)
        plt.plot(time_axis, new_time_domain, color='green')
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        plt.show()    
    
    
    def plot_frequency_domain(self):
        
        #self.location +=1
        #plt.subplot(self.location)
        frequency_domain, magnitude = self.get_frequency_domain(plot_result=False)
        plt.plot(frequency_domain, magnitude, color='green')
        plt.xlabel("Hz")
        plt.ylabel("Magnitude")
        plt.show()    
    
    
    def _plot_frequency_domain(self,frequency_domain, magnitude):
        
        #self.location +=1
        #plt.subplot(self.location)
        plt.plot(frequency_domain, magnitude, color='green')
        plt.xlabel("Hz")
        plt.ylabel("Magnitude")
        plt.show()    
    
          
        
            