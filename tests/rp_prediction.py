import sys
sys.path.append(".")
from pathlib import Path

from numpy import linspace
from matplotlib import pyplot
from numpy import hstack, inf, isnan, where
from scipy.stats import norm, cauchy

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

        return mass[indexes], abund[indexes]

max_mz = 1200

file_path = Path.cwd() / "ESI_NEG_SRFA.d"

reader_obj = ReadBrukerSolarix(file_path)

transient_obj = reader_obj.get_transient()

mass_spectrum_obj = transient_obj.get_mass_spectrum(plot_result=False)

for peak_obj_idx, peak_obj in enumerate(mass_spectrum_obj):
    
    if  peak_obj_idx != 0 or peak_obj_idx != len(mass_spectrum_obj)-1 :

        next_peak_obj = mass_spectrum_obj[peak_obj_idx + 1]
        
        previous_peak_obj = mass_spectrum_obj[peak_obj_idx - 1]
        
        if peak_obj.nominal_mz_exp < max_mz and  peak_obj.nominal_mz_exp == next_peak_obj.nominal_mz_exp and peak_obj.nominal_mz_exp == previous_peak_obj.nominal_mz_exp:
            
            sim_mz, sim_abun = peak_obj.lorentz_pdf()
            next_sim_mz, next_sim_abun = next_peak_obj.lorentz_pdf()
            previous_sim_mz, previous_sim_abun = previous_peak_obj.lorentz_pdf()
        
            summed_peaks_abun = (sim_abun + next_sim_abun + previous_sim_abun) # fix
            summed_peaks_abun = summed_peaks_abun/(max(summed_peaks_abun))

            min_mz = min(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
            max_mz = max(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
            
            summed_mz_domain = linspace(min_mz, max_mz, 10000)

            mz_centroid, abund_centroid = calc_maximum(summed_mz_domain,summed_peaks_abun)    

            mz_min_valley, abund_min_valley = calc_minimum(summed_mz_domain,summed_peaks_abun)    

            print(previous_peak_obj.mz_exp, peak_obj.mz_exp, next_peak_obj.mz_exp,)
            print(previous_peak_obj.abundance, peak_obj.abundance, next_peak_obj.abundance)
           
            delta_rp = 10000
            
            if len(mz_centroid) == 3:
                
                while  abund_min_valley[0] > 0.01 or abund_min_valley[1] > 0.01:
                    
                    previous_sim_mz, previous_sim_abun = previous_peak_obj.lorentz_pdf(delta_rp=delta_rp)
                    sim_mz, sim_abun = peak_obj.lorentz_pdf(delta_rp=delta_rp)
                    next_sim_mz, next_sim_abun = next_peak_obj.lorentz_pdf(delta_rp=delta_rp)
                    
                    summed_peaks_abun = (sim_abun + next_sim_abun + previous_sim_abun) # fix
                    summed_peaks_abun = summed_peaks_abun/(max(summed_peaks_abun))

                    min_mz = min(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
                    max_mz = max(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
                    
                    summed_mz_domain = linspace(min_mz, max_mz, 10000)

                    mz_centroid, abund_centroid = calc_maximum(summed_mz_domain,summed_peaks_abun)    

                    mz_min_valley, abund_min_valley = calc_minimum(summed_mz_domain,summed_peaks_abun)  

                    if len(abund_min_valley) == 0:
                        break
                    
                    delta_rp += 10000

                    print (abund_min_valley)

                    #pyplot.plot(mz_centroid, abund_centroid, 'o')

                    #pyplot.plot(mz_min_valley, abund_min_valley, 'o', c='g')

                    #pyplot.plot(summed_mz_domain,summed_peaks_abun)
                    
                    #pyplot.show()  
            
            '''
            delta_rp = 1000
            while abund_min_valley[1] > 1:
                
                previous_sim_mz, previous_sim_abun = previous_peak_obj.lorentz_pdf()
                sim_mz, sim_abun = peak_obj.lorentz_pdf(delta_rp=delta_rp)
                next_sim_mz, next_sim_abun = next_peak_obj.lorentz_pdf(delta_rp=delta_rp)
                
                summed_peaks_abun = (sim_abun + next_sim_abun + previous_sim_abun) # fix
                summed_peaks_abun = summed_peaks_abun/(max(summed_peaks_abun))

                min_mz = min(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
                max_mz = max(list(sim_mz) + list(previous_sim_mz) + list(next_sim_mz))
                
                summed_mz_domain = linspace(min_mz, max_mz, 10000)

                mz_centroid, abund_centroid = calc_maximum(summed_mz_domain,summed_peaks_abun)    

                mz_min_valley, abund_min_valley = calc_minimum(summed_mz_domain,summed_peaks_abun)  

                delta_rp += 100

                print (abund_min_valley)
            '''
            
            pyplot.plot(mz_centroid, abund_centroid, 'o')

            pyplot.plot(mz_min_valley, abund_min_valley, 'o', c='g')

            pyplot.plot(summed_mz_domain,summed_peaks_abun)
            
            pyplot.show()  


