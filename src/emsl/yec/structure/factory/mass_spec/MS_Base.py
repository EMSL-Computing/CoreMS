'''
Created on Jun 12, 2019
'''
from os.path import basename, dirname

from emsl.yec.structure.factory.MSPeakClasses import MSPeak
import matplotlib.pyplot as plt
from emsl.yec.calc.mass_spec.MassSpectrumCalc import MassSpecCalc


__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"
    
class MassSpecBase(MassSpecCalc):
    '''
    classdocs
    '''
    def __init__(self,exp_mz, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        self._magnitude = magnitude
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
    def magnitude(self): return self._magnitude
    
    @property
    def baselise_noise(self): return self._baselise_noise
    
    @property
    def baselise_noise_std(self): return self._baselise_noise_std
    
    @property
    def Aterm(self): return self._calibration_terms[0]
    
    @property
    def Bterm(self): return self._calibration_terms[1]
    
    @property
    def Cterm(self): return self._calibration_terms[2]        
    
    @property
    def filename(self): return basename(self.full_filename_path)
    
    @property
    def dir_location(self): return dirname(self.full_filename_path)
    
    def set_noise_treshould(self, auto=True):
        
        self._baselise_noise, self._baselise_noise_std  = self.run_noise_threshould_calc(auto)
        
    def set_mspeaks(self):
        
        #pass
        self.mspeaks = [MSPeak(1, 200.1, 100, 1e6, 4, 0)] 
        #self.mspeaks = [MSPeak(self.ion_charge, exp_mass, exp_freq, relative_abundance) for exp_mass, exp_freq, relative_abundance in zip(exp_masses,exp_freqs, relative_abundances)]
            
    def assign_molecular_formulas(self):
        
        '''call assigment algorithms here'''
        
        formula_dict = {'C':10,'H':20, 'O':10, "IonType": 'closed_shell' }
         
        self.mspeaks[0].molecular_formula = formula_dict
        
    def plot_mz_domain_profile_and_noise_threshold(self):
        
        if self.baselise_noise and self.baselise_noise:
            x = (self.exp_mz.min(), self.exp_mz.max())
            y = (self.baselise_noise, self.baselise_noise) 
        
            threshold = (self.baselise_noise + (3*self.baselise_noise_std))
            plt.plot(self.exp_mz, self.magnitude, color='green')
            plt.plot(x, (threshold,threshold), color='yellow')
            plt.plot(x, y, color='red')
            plt.xlabel("m/z")
            plt.ylabel("Magnitude")
            plt.show()      
        else: 
            
            raise Exception("Calculate noise threshold first")
        
    def plot_mz_domain_profile(self):
        
        #self.location +=1
        #plt.subplot(self.location)
        #print(self.exp_mz)
        #print(self.magnitude)
        
        plt.plot(self.exp_mz, self.magnitude, color='green')
        plt.xlabel("m/z")
        plt.ylabel("Magnitude")
        plt.show()          
        

           

        
 