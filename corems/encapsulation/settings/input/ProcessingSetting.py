__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

class TransientSetting:
    
    implemented_apodization_function = ("Hamming", "Hanning", "Blackman")
    apodization_method = "Hanning"
    number_of_truncations = 0
    number_of_zero_fills = 1
            
class MassSpectrumSetting:
    
    threshold_method = "auto"
    implemented_noise_threshold_methods = ("auto", "signal_noise", "relative_abundance")
    noise_threshold_std = 6
    s2n_threshold = 4
    relative_abundance_threshold = 6 # from 1-100
    
    min_noise_mz = 100.0
    max_noise_mz = 1200.0  
    
    min_picking_mz = 100.0
    max_picking_mz = 1200.0
    
class MassSpecPeakSetting:
    
    '''needs to clear previous results from peak_picking'''
    #kendrick_base =  {"C": 1, "H":2}   
    kendrick_base =  {"C": 1, "H":2, "O":1}   
