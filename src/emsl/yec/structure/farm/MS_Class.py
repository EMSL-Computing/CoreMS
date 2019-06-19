'''
Created on Jun 12, 2019
'''

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from emsl.yec.calc.MassSpecCalc import MassSpecCalculations

class MassSpec(MassSpecCalculations):
    '''
    classdocs
    '''

    #change for panda object input
    def __init__(self, scan_number,polarity, noise_std, f_rt, ): #input_dataframe
                 
                #  polarity, a_exp_mz, a_exp_freqs, a_magnitude, scan_number, rt, 
                # noise_std=None, l_base_noise=None, l_signal_to_noises=None, l_charge=None, l_resolution=None ):
        
        self.ion_charge = polarity
        self.scan_number = scan_number
        self.retention_time = f_rt
        self.noise_std = noise_std
        
        
        self.molecular_formulas = None 
        