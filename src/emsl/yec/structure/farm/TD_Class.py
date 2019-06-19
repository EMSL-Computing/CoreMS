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
        
        self.__set__parameters__objects(d_params)
        
        self.__set__transient__time()
        
        
        
    def __set__parameters__objects(self, d_params):
        
        self.full_filename_path  = d_params.get("filename_path")
        
        self.Atherm = d_params.get("Atherm")
        
        self.BTherm = d_params.get("Btherm")
        
        self.CTherm = d_params.get("CTherm")
        
        self.exc_high_freq = d_params.get("exc_high_freq")
        
        self.exc_low_freq = d_params.get("exc_low_freq")
        
        self.bandwidth = d_params.get("bandwidth")
        
        self.number_data_points = d_params.get("number_data_points")
        
    def __set__transient__time(self):
        
        ### needs __set__parameters__ 
        self.transient_time = self.cal_transient_time()
    
    @property
    def get_filename(self):
        
        return basename(self.full_filename_path)
    
    @property
    def get_dir_location(self):
        
        return dirname(self.full_filename_path)
    
    def plot_transient(self):
        
        import matplotlib.pyplot as plt
        time_axis = linspace(0, self.transient_time, num=self.number_data_points)
        plt.plot(time_axis, self.transient_data)
        plt.show()    
    