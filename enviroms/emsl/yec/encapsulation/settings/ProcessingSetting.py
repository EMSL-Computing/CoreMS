__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


class TransientSetting():
    
    implemented_apodization_function = ["Hamming", "Hanning", "Blackman"]
    
    _apodization_method = 'Hanning'
    _number_of_truncations = 0
    _number_of_zero_fills = 1 
    
    @property
    def apodization_method(self): return self._apodization_method
    
    @property
    def number_of_zero_fills(self): return self._number_of_zero_fills    
    
    @property
    def number_of_truncations(self):  return self._number_of_truncations
    
    @apodization_method.setter
    def apodization_method(self, apodization_method):
    
        if apodization_method in self.implemented_apodization_function: 
            return self._apodization_method 
        else: 
            raise Exception("%s function is not implemented, please refer to TransientCalc Class" % apodization_method) 
        
    @number_of_truncations.setter
    def number_of_truncations(self, number_of_truncations): return self.check_negative_number(number_of_truncations)
    
    @number_of_zero_fills.setter
    def number_of_zero_fills(self, number_of_zero_fills): return self.check_negative_number(number_of_zero_fills)
    
    
    def check_negative_number(self, number):
        if number >= 0: 
            return number
        else: 
            raise Exception("Can not be negative") 
        
            
class MassSpectrumSetting():
    
    threshold_method = "auto"
    implemented_noise_threshold_methods = {"auto", "signal_noise", "relative_abudance"}
    noise_threshold_stds = 6.0
    s2n_threshold = 4
    relative_abundace_threshold = 5# from 1-100
    
    min_noise_mz = 100.0
    max_noise_mz = 200.0    
    
    min_picking_mz = 100.0
    max_picking_mz = 1000.0
    
class MassSpecPeakSetting():
    
    '''needs to clear previous results from peak_picking'''
    kendrick_base =  {"C": 1, "H": 2}   
  