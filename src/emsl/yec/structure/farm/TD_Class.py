'''
Created on Jun 12, 2019

TIME DOMAIN CLASS
'''
from os.path import basename, dirname

from numpy import linspace

from emsl.yec.calc.TD_Calc import TransientCalculations


__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

class Transient(TransientCalculations):
    
    def __init__(self, data, d_params):
        
        self.transient_data = data
        
        self.frequency_domain = None
        
        self.magnitude = None
        
        self.__set__parameters__objects(d_params)
        
        self.__set__transient__time()
        
        
    def __set__parameters__objects(self, d_params):
        
        self.full_filename_path  = d_params.get("filename_path")
        
        self.Atherm = d_params.get("Atherm")
        
        self.BTherm = d_params.get("Btherm")
        print('BBB', self.BTherm)
        self.CTherm = d_params.get("CTherm")
        
        self.exc_high_freq = d_params.get("exc_high_freq")
        
        self.exc_low_freq = d_params.get("exc_low_freq")
        
        self.bandwidth = d_params.get("bandwidth")
        
        self.number_data_points = d_params.get("number_data_points")
        
    def __set__transient__time(self):
        
        ### needs __set__parameters__ 
        self.transient_time = self.cal_transient_time()
    
    ''''move it to main code'''
    def set_frequency_domain(self):
        
        apodization_method = 'Hanning'
        number_of_truncations = 0
        number_of_zero_fills = 1
        
        if number_of_truncations > 0:
            
            new_time_domain = self.truncation(self.transient_data, number_of_truncations)
            
        else: 
            
            new_time_domain = self.transient_data
            
        if apodization_method != "None":
                
            new_time_domain = self.apodization(new_time_domain, apodization_method)  
        
        self.plot_transient(self.transient_data)
        
        self.plot_transient(new_time_domain)
        
        time_domain_y = self.zero_fill(number_of_zero_fills, new_time_domain)     
        
        
        self.plot_transient(time_domain_y)
        
        self.frequency_domain, self.magnitude = self.perform_magniture_mode_ft(time_domain_y, number_of_zero_fills) 
    
    @property
    def get_filename(self):
        
        return basename(self.full_filename_path)
    
    @property
    def get_dir_location(self):
        
        return dirname(self.full_filename_path)
    
    def plot_transient(self, transient_data):
        
        import matplotlib.pyplot as plt
        time_axis = linspace(0, self.transient_time, num=len(transient_data))
        plt.plot(time_axis, transient_data)
        plt.show()    
    
    '''move it to mass_spec_object'''
    def plot_frequency_domain(self):
        
        import matplotlib.pyplot as plt
        plt.plot(self.frequency_domain, self.magnitude)
        plt.show()    
    
    '''move to ms_class'''
    def plot_mz_domain(self):
        mz = self.f_to_mz()
        import matplotlib.pyplot as plt
        plt.plot(mz, self.magnitude)
        plt.show()        
        
            