

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from copy import deepcopy
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from numpy import nan
from corems.ms_peak.calc.MSPeakCalc import MSPeakCalculation
from corems.mass_spectra.calc import SignalProcessing as sp
from numpy import NaN, power

class _MSPeak(MSPeakCalculation):
    '''
    classdocs
    '''
    def __init__(self, ion_charge, mz_exp, abundance, resolving_power, 
                 signal_to_noise, indexes, index, ms_parent=None, exp_freq=None):

        self._ms_parent = ms_parent
        # needed to create the object
        self.ion_charge = int(ion_charge)
        self._mz_exp = float(mz_exp)
        self.mass = float(mz_exp) / float(ion_charge)
        self.abundance = float(abundance)
        self.resolving_power = float(resolving_power)
        self.signal_to_noise = float(signal_to_noise)
        # profile indexes
        self.start_scan = int(indexes[0]) 
        self.apex_scan = int(indexes[1])
        self.final_scan = int(indexes[2]) 
        
        #mass spec obj index
        self.index = int(index)
        # parent mass spectrum obj instance

        # updated after mass error prediction'
        self.predicted_std = None
        # updated after calibration'
        self.mz_cal = None
        # updated individual calculation'
        self.baseline_noise = None

        if exp_freq:
            self.freq_exp = float(exp_freq)

        if self._ms_parent is not None:
            kendrick_dict_base = self._ms_parent.mspeaks_settings.kendrick_base
        else:
            kendrick_dict_base = {'C':1, 'H':2}
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)

        'updated after molecular formula ID'

        self.molecular_formulas = []
        self._confidence_score = None
        # placeholder for found isotopologues index 
        self.isotopologue_indexes = []
        # placeholder for found isotopologues molecular formula obj
        self.found_isotopologues = {}

        # Label for what type of peak it is - real signal, noise, sinc wiggle, magnetron or harmonic peak, etc. 
        self.peak_type = None

    def __len__(self):

        return len(self.molecular_formulas)

    def __setitem__(self, position, molecular_formula_obj):

        self.molecular_formulas[position] = molecular_formula_obj

    def __getitem__(self, position) -> MolecularFormula:

        return self.molecular_formulas[position]

    def change_kendrick_base(self, kendrick_dict_base):
        '''kendrick_dict_base = {"C": 1, "H": 2}'''
        self._kdm, self._kendrick_mass, self._nominal_km = self._calc_kdm(
            kendrick_dict_base)

    def add_molecular_formula(self, molecular_formula_obj):

        # freeze state
        molecular_formula_obj._mspeak_parent = self

        # new_mol_formula = deepcopy(molecular_formula_obj)
        # add link mass spectrum obj instance

        # new_mol_formula.mspeak_parent = self

        self.molecular_formulas.append(molecular_formula_obj)

        return molecular_formula_obj

    def remove_molecular_formula(self, mf_obj):

        self.molecular_formulas.remove(mf_obj)

    def clear_molecular_formulas(self):

        self.molecular_formulas = []

    @property
    def mz_exp(self):
        if self.mz_cal:
            return self.mz_cal
        else:
            return self._mz_exp

    @mz_exp.setter
    def mz_exp(self, mz_exp):
        self._mz_exp = mz_exp

    @property
    def area(self): 
        if self._ms_parent.is_centroid:
            return nan
        else:    
            return self.calc_area()

    @property
    def nominal_mz_exp(self): return int(self.mz_exp)

    @property
    def kmd(self): return self._kdm

    @property
    def kendrick_mass(self): return self._kendrick_mass

    @property
    def knm(self): return self._nominal_km

    @property
    def is_assigned(self) -> bool:

        return bool(self.molecular_formulas)

    def plot_simulation(self, sim_type="lorentz", ax=None, color="green",
                        oversample_multiplier=1, delta_rp = 0, mz_overlay=1):


        if self._ms_parent:

            import matplotlib.pyplot as plt

            x, y = eval("self."+sim_type+"(oversample_multiplier="+str(oversample_multiplier)+", delta_rp="+str(delta_rp)+", mz_overlay="+str(mz_overlay)+")")

            if ax is None:
                    ax = plt.gca()
            ax.plot(x, y, color=color, label="Simulation")
            ax.set(xlabel='m/z', ylabel='abundance')

            plt.legend()
            return ax

    def plot(self, ax=None, color="black", derivative=True, deriv_color='red'): #pragma: no cover
        
        if self._ms_parent:

            import matplotlib.pyplot as plt

            if ax is None:
                ax = plt.gca()
            x = self._ms_parent.mz_exp_profile[self.start_scan: self.final_scan]
            y =  self._ms_parent.abundance_profile[self.start_scan: self.final_scan]

            ax.plot(x, y, color=color, label="Data")
            ax.set(xlabel='m/z', ylabel='abundance')
            if derivative and not self._ms_parent.is_centroid:
                dy = sp.derivate(self._ms_parent.abundance_profile[self.start_index: self.final_index+1])
                ax.plot(x, dy, c=deriv_color)
            else:
                ax.plot((self.mz_exp, self.mz_exp),(0, self.abundance), color=color, label="Data")
                            
            #plt.legend()

            return ax

        else:
            print("Centroid Peak Object")

    @property
    def best_molecular_formula_candidate(self):
        
        if self._ms_parent.molecular_search_settings.score_method == "N_S_P_lowest_error":
            return self.cia_score_N_S_P_error()
        
        elif self._ms_parent.molecular_search_settings.score_method == "S_P_lowest_error":
            return self.cia_score_S_P_error()

        elif self._ms_parent.molecular_search_settings.score_method == "lowest_error":
            return self.molecular_formula_lowest_error()    
        
        elif self._ms_parent.molecular_search_settings.score_method == "air_filter_error":
            return self.molecular_formula_air_filter()    

        elif self._ms_parent.molecular_search_settings.score_method == "water_filter_error":
            return self.molecular_formula_water_filter()    

        elif self._ms_parent.molecular_search_settings.score_method == "earth_filter_error":
            return self.molecular_formula_earth_filter()   

        elif self._ms_parent.molecular_search_settings.score_method == "prob_score":
            return self.molecular_formula_highest_prob_score()
        else:
            
            raise TypeError("Unknown score method selected: % s, \
                            Please check score_method at \
                            encapsulation.settings.molecular_id.MolecularIDSettings.MolecularFormulaSearchSettings", 
                            self._ms_parent.parameters.molecular_search.score_method)    

class ICRMassPeak(_MSPeak):

    def __init__(self, *args, ms_parent=None, exp_freq=None):

        super().__init__(*args,exp_freq=exp_freq, ms_parent=ms_parent)

    def resolving_power_calc(self, B, T):
        
        '''
        low pressure limits, 
        T: float 
            transient time
        B: float
            Magnetic Filed Strength (Tesla)    
        
        reference
        Marshall et al. (Mass Spectrom Rev. 1998 Jan-Feb;17(1):1-35.)
        DOI: 10.1002/(SICI)1098-2787(1998)17:1<1::AID-MAS1>3.0.CO;2-K
        '''
        return (1.274e7 * self.ion_charge * B * T)/ (self.mz_exp*self.ion_charge)

    def set_calc_resolving_power(self, B, T):

        self.resolving_power = self.resolving_power_calc(B, T) 

    def _mz_to_f_bruker(self):
        '''
        Convert a peak m/z value to frequency
        '''
        '''
        # Currently Broken - Not sure why
        if self.mz_cal:
            mz_val = self.mz_cal
        else:
            mz_val = self.mz_exp
        Aterm, Bterm, Cterm = self._ms_parent.Aterm, self._ms_parent.Bterm, self._ms_parent.Cterm
        # Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        
        if Cterm == 0:
            
            if Bterm == 0:
                #uncalibrated data
                freq_domain = Aterm / mz_val
                
            else:
                
                freq_domain = (Aterm / (mz_val)) - Bterm

        # @will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:

            freq_domain = (Aterm / mz_val) + (Bterm / power(mz_val, 2)) + Cterm

        return freq_domain
        '''
        print('Function not confirmed to work, disabled.')
        return None





    
        
class TOFMassPeak(_MSPeak):

    def __init__(self, *args, exp_freq=None):

        super().__init__(*args,exp_freq=exp_freq)

    def set_calc_resolving_power(self):
        return 0

class OrbiMassPeak(_MSPeak):

    def __init__(self, *args, exp_freq=None):

        super().__init__(*args,exp_freq=exp_freq)

    def set_calc_resolving_power(self):
        return 0       

