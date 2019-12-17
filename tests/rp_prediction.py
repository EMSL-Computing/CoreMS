import sys
sys.path.append(".")

from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
from pathlib import Path
from numpy import linspace
from matplotlib import pyplot

file_path = Path.cwd() / "ESI_NEG_SRFA.d"

reader_obj = ReadBrukerSolarix(file_path)

transient_obj = reader_obj.get_transient()

mass_spectrum_obj = transient_obj.get_mass_spectrum(plot_result=False).sort_by_mz()

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
            
            print(previous_peak_obj.mz_exp, peak_obj.mz_exp, next_peak_obj.mz_exp,)
            print(previous_peak_obj.abundance, peak_obj.abundance, next_peak_obj.abundance)
            
            summed_mz_domain = linspace(min_mz, max_mz, 10000)

            pyplot.plot(summed_mz_domain,summed_peaks_abun)
            
            pyplot.show()    

