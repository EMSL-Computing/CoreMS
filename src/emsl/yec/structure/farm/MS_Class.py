'''
Created on Jun 12, 2019
'''

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from emsl.yec.calc.MS_Calc import MassSpecCalculations

class MassSpec(MassSpecCalculations):
    '''
    classdocs
    '''

    #change for panda object input
    def __init__(self, calibration_terms,  polarity, magnitude, frequency=None, mz=None, scan_number=1): #input_dataframe
                 
                #  polarity, a_exp_mz, a_exp_freqs, a_magnitude, scan_number, rt, 
                # noise_std=None, l_base_noise=None, l_signal_to_noises=None, l_charge=None, l_resolution=None ):
        
        self.calibration_terms = calibration_terms #(A, B, C) 
        self.frequency_domain = frequency
        self.magnitude = magnitude
        self.mz = mz
        self.polarity = polarity
        self.scan_number = scan_number
        self.peakObjs = None 
        
        if not mz:
            self.mz = self.f_to_mz()
    
    @property
    def get_A_therm(self):
        
        return self.calibration_terms[0]
    
    @property
    def get_B_therm(self):
        
        return self.calibration_terms[1]
    
    @property
    def get_C_therm(self):
        
        return self.calibration_terms[2]        