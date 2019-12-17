import sys
sys.path.append(".")
from pathlib import Path

from numpy import linspace
from matplotlib import pyplot
from numpy import hstack, inf, isnan, poly1d, polyfit, where

from corems.transient.input.BrukerSolarix import ReadBrukerSolarix

def calc_maximum(mass, abund):

        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            
            return mass[indexes], abund[indexes]

def calc_minimum(mass, abund):

        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) > 0) & (hstack((0, dy)) < 0))[0]

        if indexes.size:
            
            return mass[indexes], abund[indexes]


file_path = Path.cwd() / "ESI_NEG_SRFA.d"

reader_obj = ReadBrukerSolarix(file_path)

transient_obj = reader_obj.get_transient()

mass_spectrum_obj = transient_obj.get_mass_spectrum(plot_result=False)

for peak_obj_idx, peak_obj in enumerate(mass_spectrum_obj):
    
    rp = peak_obj.resolving_power
    
    if  peak_obj_idx != 0 or peak_obj_idx != len(mass_spectrum_obj)-1 :

        sim_mz, sim_abun = peak_obj.lorentz_pdf()
        
        #pyplot.plot(sim_mz,sim_abun)
        
        #pyplot.show()    

        next_peak_obj = mass_spectrum_obj[peak_obj_idx + 1]
        next_sim_mz, next_sim_abun = next_peak_obj.lorentz_pdf()
        
        previous_peak_obj = mass_spectrum_obj[peak_obj_idx - 1]
        previous_sim_mz, previous_sim_abun = previous_peak_obj.lorentz_pdf()
        
        summed_peaks_abun = sim_abun + next_sim_abun + previous_sim_abun # fix

        min_mz = min(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
        max_mz = max(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
        
        if peak_obj.nominal_mz_exp == next_peak_obj.nominal_mz_exp and peak_obj.nominal_mz_exp == previous_peak_obj.nominal_mz_exp:
            
            summed_mz_domain = linspace(min_mz, max_mz, 10000)

            mz_centroid, abund_centroid = calc_maximum(summed_mz_domain,summed_peaks_abun)    

            mz_min_valley, abund_min_valley = calc_minimum(summed_mz_domain,summed_peaks_abun)    

            print(previous_peak_obj.mz_exp, peak_obj.mz_exp, next_peak_obj.mz_exp,)
            print(previous_peak_obj.abundance, peak_obj.abundance, next_peak_obj.abundance)
            
            pyplot.plot(mz_centroid, abund_centroid, 'o')

            pyplot.plot(mz_min_valley, abund_min_valley, 'o', c='g')


            pyplot.plot(summed_mz_domain,summed_peaks_abun)
            
            pyplot.show()  


