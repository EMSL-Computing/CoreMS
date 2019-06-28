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
    Mass Spectrum class object with common features and functions
    '''
    def __init__(self,exp_mz, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        self._magnitude = magnitude
        self._exp_mz = exp_mz
        self.mspeaks = []
        
        self._set_parameters_objects(d_params)
        
        #for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value) 
        
    def _set_parameters_objects(self, d_params):
        
        self._calibration_terms = (d_params.get("Aterm"), d_params.get("Bterm"), d_params.get("Cterm"))
        
        self.polarity = int(d_params.get("polarity"))
        
        self.scan_number = d_params.get("scan_number")
        
        self.rt = d_params.get("rt")
        
        self._filename = d_params.get("filename")
        
        self._dir_location = d_params.get("dir_location")
         
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
    def filename(self): return self._filename
    
    @property
    def dir_location(self): return self._dir_location
    
    def find_peaks(self):
        
        self.do_peak_picking()
        
        print( "A total of %i peaks were found" % len(self.mspeaks))
    
    def set_noise_treshould(self, auto=True):
        
        self._baselise_noise, self._baselise_noise_std  = self.run_noise_threshould_calc(auto)
        
    def add_mspeak(self, ion_charge, exp_mz, magnitude, resolving_power, signal_to_noise, massspec_index, exp_freq=None):
        #parms ion_charge, exp_mz, magnitude, resolving_power, signal_to_noise, massspec_index, 
        #print( self.mspeaks)
        
        self.mspeaks.append(MSPeak(ion_charge, exp_mz, magnitude, resolving_power, signal_to_noise, massspec_index, exp_freq=exp_freq)) 
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
        
        plt.plot(self.exp_mz, self.magnitude, color='green')
        plt.xlabel("m/z")
        plt.ylabel("Magnitude")
        plt.show()          
   