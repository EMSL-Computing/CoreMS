'''
Created on Jun 12, 2019
'''
from os.path import basename, dirname

from emsl.yec.calc.MassSpectrumCalc import MassSpecCalculations
from emsl.yec.structure.factory.MSPeakClasses import MSPeak
import matplotlib.pyplot as plt


__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
    
class MassSpecBase(MassSpecCalculations):
    '''
    classdocs
    '''
    def __init__(self,exp_mz, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        self.magnitude = magnitude
        self._exp_mz = exp_mz
        self.mspeaks = None
        
        self._set_parameters_objects(d_params)
        
        #for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value) 
        
    def _set_parameters_objects(self, d_params):
        
        self._full_filename_path  = d_params.get("filename_path")
        
        self._calibration_terms = (d_params.get("Aterm"), d_params.get("Bterm"), d_params.get("Cterm"))
        
        self.polarity = d_params.get("polarity")
        
        self.scan_number = d_params.get("scan_number")
        
        self.rt = d_params.get("rt")
        
        self.location = 220 
        
    @property
    def exp_mz(self): return self._exp_mz
    
    @property
    def Aterm(self):
        
        return self._calibration_terms[0]
    
    @property
    def Bterm(self):
        
        return self._calibration_terms[1]
    
    @property
    def Cterm(self):
        
        return self._calibration_terms[2]        
    
    @property
    def filename(self):
        
        return basename(self.full_filename_path)
    
    @property
    def dir_location(self):
        
        return dirname(self.full_filename_path)
    '''move to ms_class'''
    
    def set_mspeaks(self):
        
        #pass
        self.mspeaks = [MSPeak(1, 200.1, 100, 1e6, 4, 0)] 
        #self.mspeaks = [MSPeak(self.ion_charge, exp_mass, exp_freq, relative_abundance) for exp_mass, exp_freq, relative_abundance in zip(exp_masses,exp_freqs, relative_abundances)]
            
    def assign_molecular_formulas(self):
        
        '''call assigment algorithms here'''
        
        formula_dict = {'C':10,'H':20, 'O':10, "IonType": 'closed_shell' }
         
        self.mspeaks[0].molecular_formula = formula_dict
        
    def plot_mz_domain_profile(self):
        
        #self.location +=1
        #plt.subplot(self.location)
        #print(self.exp_mz)
        #print(self.magnitude)
        
        plt.plot(self.exp_mz, self.magnitude, color='green')
        plt.xlabel("m/z")
        plt.ylabel("Magnitude")
        plt.show()          
        
class MassSpecfromFreq(MassSpecBase):
       
    def __init__(self,frequency_domain, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        super().__init__(None, magnitude, d_params)
        
        self.frequency_domain = frequency_domain
        
        self._set_mz_domain()
        
        self.mspeaks = None
        
        """implement here, code in MassSpecCalc"""
        #self.calc_noise_threshould
        #self.peak_picking
        #self.calc_resolving_power  
   
    def _set_mz_domain(self):
            
        self._exp_mz = self.f_to_mz()
           
class MassSpecCentroid(MassSpecBase):
    
    '''
    classdocs
    '''
    
    def __init__(self, exp_mz, magnitude, dataframe, d_params, **kwargs): 
                 
        '''
        Constructor
        dataframe will contain l_base_noise, l_signal_to_noise, l_charge, l_resolving_power
        needs to simulate peak shape and pass as exp_mz and magnitude.
        '''
        super().__init__(exp_mz, magnitude, d_params)
        
        #self.__process__from__centroid(l_base_noise, l_signal_to_noise, l_charge, l_resolving_power)
        
    def __process__from__centroid(self, l_base_noise, l_signal_to_noise, l_charge, l_resolving_power):
        
        if l_base_noise.any():
            self.set_base_noise(l_base_noise)
        
        if l_signal_to_noise.any():
            self.set_noise(l_signal_to_noise)
        
        if l_charge.any():
            self.set_charge(l_charge)
        
        if l_resolving_power.any():
            self.set_resolving_power(l_resolving_power)  
        
class MassSpecProfile(MassSpecBase):
    
    '''
    classdocs
    '''
    
    def __init__(self, exp_mz, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        super().__init__(exp_mz, magnitude, d_params)
        
        self.calc_noise_threshould
        self.peak_picking
        self.calc_resolving_power  
        
        #for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value)            