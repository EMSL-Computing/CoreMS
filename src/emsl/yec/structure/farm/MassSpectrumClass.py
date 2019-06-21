'''
Created on Jun 12, 2019
'''
from os.path import basename, dirname

from emsl.yec.calc.MassSpectrumCalc import MassSpecCalculations
from emsl.yec.structure.farm.MolecularFormulaClass import MSPeak
import matplotlib.pyplot as plt


__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class MassSpec(MassSpecCalculations):
    '''
    classdocs
    '''
    def __init__(self,frequency_domain, magnitude, d_params, **kwargs): 
                 
        '''
        Constructor
        '''
        self.frequency_domain = frequency_domain
        self.magnitude = magnitude
        self.mz_domain = None
        self.molecular_formulas = None
        
        self.__set__parameters__objects(d_params)
        
        self.__set_mz_domain()
        #for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value) 
        
    def __set__parameters__objects(self, d_params):
        
        self._full_filename_path  = d_params.get("filename_path")
        
        self.calibration_terms = (d_params.get("Aterm"), d_params.get("Bterm"), d_params.get("Cterm"))
        
        self.polarity = d_params.get("polarity")
        
        self.scan_number = d_params.get("scan_number")
        
        self.rt = d_params.get("rt")
        
        self.location = 220 
    
    def __set_mz_domain(self):
            
        self.exp_mz = self.f_to_mz()
    
    def set_creat_molecular_formulas(self):
        
        pass
        #self.molecular_formulas = [MSPeak(self.ion_charge, exp_mass, exp_freq, relative_abundance) for exp_mass, exp_freq, relative_abundance in zip(exp_masses,exp_freqs, relative_abundances)]
    
    
    def __process__from__centroid(self, l_base_noise, l_signal_to_noise, l_charge, l_resolving_power):
        
        if l_base_noise.any():
            self.set_base_noise(l_base_noise)
        
        if l_signal_to_noise.any():
            self.set_noise(l_signal_to_noise)
        
        if l_charge.any():
            self.set_charge(l_charge)
        
        if l_resolving_power.any():
            self.set_resolving_power(l_resolving_power)  

    
    @property
    def A_therm(self):
        
        return self.calibration_terms[0]
    
    @property
    def B_therm(self):
        
        return self.calibration_terms[1]
    
    @property
    def C_therm(self):
        
        return self.calibration_terms[2]        
    
    @property
    def filename(self):
        
        return basename(self.full_filename_path)
    
    @property
    def dir_location(self):
        
        return dirname(self.full_filename_path)
    '''move to ms_class'''
    
    def plot_mz_domain_profile(self):
        
        #self.location +=1
        #plt.subplot(self.location)
        plt.plot(self.mz_domain, self.magnitude, color='green')
        plt.xlabel("m/z")
        plt.ylabel("Magnitude")
        plt.show()          