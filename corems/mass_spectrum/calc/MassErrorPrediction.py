__author__ = 'Yuri E. Corilo'
__date__ = "03/31/2020"

from threading import Thread
from pandas import DataFrame
from numpy import hstack, inf, isnan, where, array
from tqdm import tqdm

class MassErrorPrediction(Thread):
    
    def __init__(self, mass_spectrum, mz_overlay=10, rp_increments=10000, base_line_target=0.01, max_interation=1000, interpolation='linear'):
        
        Thread.__init__(self)
        
        self.mass_spectrum_obj = mass_spectrum

        self.mz_overlay = mz_overlay

        self.rp_increments = 10000

        self.base_line_target = 0.01 

        self.max_interation = 1000

        self.df = None

        self.interpolation = interpolation
    
    def run(self):
            
        self.df = self.calc_error_dist()

    def get_results(self):

        if not self.df:
            self.run()

        return self.df

    def calc_error_dist(self):
        
        results_list = []
        
        indexes_without_results = list(range(len(self.mass_spectrum_obj)))
        # loop trough mass spectrum

        for peak_obj_idx, peak_obj in enumerate(tqdm(self.mass_spectrum_obj)):
            
            # access ms peaks triplets ( peak_obj_idx -1, peak_obj_idx, and peak_obj_idx + 1)
            # check lower and upper boundaries to not excesses mass spectrum range
            
            if  peak_obj_idx != 0 and peak_obj_idx != len(self.mass_spectrum_obj)-1:
                
                # current peak_obj initialted in the loop expression
                # geting the peak on the left (previous_peak_obj) and the one in the right position (next_peak_obj)
                next_peak_obj = self.mass_spectrum_obj[peak_obj_idx + 1]
                previous_peak_obj = self.mass_spectrum_obj[peak_obj_idx - 1]
                
                # check mz range defined in max_mz variable and check if peaks have same nominal mz
                # keeping same mz for better plotting representation only, remove it for production
                if  peak_obj.nominal_mz_exp == next_peak_obj.nominal_mz_exp and peak_obj.nominal_mz_exp == previous_peak_obj.nominal_mz_exp:
                    
                    #simulate peak shape
                    sim_mz, sim_abun = peak_obj.gaussian(mz_overlay=self.mz_overlay)
                    #update_plot(sim_mz,sim_abun, 0.5)
                    
                    #simulate peak shape
                    next_sim_mz, next_sim_abun = next_peak_obj.gaussian(mz_overlay=self.mz_overlay)
                    #update_plot(next_sim_mz, next_sim_abun, 0.5)
                    
                    
                    #simulate peak shape
                    previous_sim_mz, previous_sim_abun = previous_peak_obj.gaussian(mz_overlay=self.mz_overlay)
                    #update_plot(previous_sim_mz,  previous_sim_abun, 0.5)
                    
                    sim_mz_domain,  summed_peaks_abun = self.sum_data( ((previous_sim_mz,previous_sim_abun),  (sim_mz,sim_abun), (next_sim_mz, next_sim_abun)) )
                    #update_plot(sim_mz_domain,summed_peaks_abun, 0.5)
                    
                    #sum simulated abundances 
                    #summed_peaks_abun = (sim_abun + next_sim_abun + previous_sim_abun) 
                    
                    #normalize abundances to 0-1
                    #summed_peaks_abun = summed_peaks_abun/(max(summed_peaks_abun))

                    #find appexes location (mz) and magnitude
                    mz_centroid, abund_centroid = self.find_peak_apex(sim_mz_domain,summed_peaks_abun)    

                    #find valley location (mz_min_valley) and magnitude (abund_min_valley)
                    mz_min_valley, abund_min_valley = self.find_peak_valley(sim_mz_domain, summed_peaks_abun)  

                    # clear delta_rp (global implementation) and store choose resolving power increments   
                    delta_rp = self.rp_increments
                    
                    # used to limited number of iterations
                    i = 0
                    j = 0
                    
                    # TODO: fit peak shape and decide best fit #gaussian, lorentz and voigt 
                    #plot_triplets(mz_centroid,abund_centroid, mz_min_valley, abund_min_valley, sim_mz_domain, summed_peaks_abun )
                    if len(mz_centroid) == 2 :
                            
                        while len(mz_centroid) < 3 and i <= self.max_interation:
                            
                            previous_sim_mz, previous_sim_abun = previous_peak_obj.gaussian(delta_rp=delta_rp, mz_overlay=self.mz_overlay)
                            
                            sim_mz, sim_abun = peak_obj.gaussian(delta_rp=delta_rp, mz_overlay=self.mz_overlay)
                            
                            next_sim_mz, next_sim_abun = next_peak_obj.gaussian(delta_rp=delta_rp, mz_overlay=self.mz_overlay)

                            sim_mz_domain,  summed_peaks_abun = self.sum_data( ((previous_sim_mz,previous_sim_abun),  (sim_mz,sim_abun), (next_sim_mz, next_sim_abun)) )
                            
                            #update_plot(sim_mz_domain,  summed_peaks_abun, 0.01)

                            mz_centroid, abund_centroid = self.find_peak_apex(sim_mz_domain,summed_peaks_abun)    

                            delta_rp += self.rp_increments
                            
                            i += 1

                        mz_min_valley, abund_min_valley = self.find_peak_valley(sim_mz_domain, summed_peaks_abun)      

                    if len(mz_centroid) == 3 and len(abund_min_valley) == 2:
                        # increase all three peak resolving power until both valley magnitude is bellow the defined target
                        # calculate peak shapes with the needed resolving power to have a baseline resolution for all peaks
                        # calculate mass difference (ppm) between original centroid and the new simulated peak. 
                        
                        while  abund_min_valley[0] > self.base_line_target or abund_min_valley[1] > self.base_line_target and j <= self.max_interation:
                            
                            previous_sim_mz, previous_sim_abun = previous_peak_obj.gaussian(delta_rp=delta_rp, mz_overlay=self.mz_overlay)
                            
                            sim_mz, sim_abun = peak_obj.gaussian(delta_rp=delta_rp, mz_overlay=self.mz_overlay)
                            
                            next_sim_mz, next_sim_abun = next_peak_obj.gaussian(delta_rp=delta_rp, mz_overlay=self.mz_overlay)

                            sim_mz_domain,  summed_peaks_abun = self.sum_data( ((previous_sim_mz,previous_sim_abun),  (sim_mz,sim_abun), (next_sim_mz, next_sim_abun)) )
                            
                            #update_plot(sim_mz_domain,  summed_peaks_abun, 0.001)
                            
                            #summed_peaks_abun = (sim_abun + next_sim_abun + previous_sim_abun) 
                            
                            
                            #find appexes location (mz) and magnitude
                            mz_centroid, abund_centroid = self.find_peak_apex(sim_mz_domain,summed_peaks_abun)    
                            
                            #find valley location (mz_min_valley) and magnitude (abund_min_valley)
                            summed_peaks_abun = summed_peaks_abun/(summed_peaks_abun.max())
                            mz_min_valley, abund_min_valley = self.find_peak_valley(sim_mz_domain, summed_peaks_abun)  

                            if len(abund_min_valley) != 2:
                                break
                            
                            delta_rp += self.rp_increments
                            j += 1
                            
                            #plot_triplets(mz_centroid,abund_centroid, mz_min_valley, abund_min_valley, sim_mz_domain, summed_peaks_abun )
                        
                        
                        #plot_triplets(mz_centroid,abund_centroid, mz_min_valley, abund_min_valley, sim_mz_domain, summed_peaks_abun )

                        mass_shift_ppp = self.calc_error(mz_centroid[1], peak_obj.mz_exp, 1000000)
                        #delta_mz = mz_centroid[1] - peak_obj.mz_exp
                        height_shift_per = self.calc_error(abund_centroid[1], peak_obj.abundance, 100)
                        #excitation_amplitude = str(mass_spectrum_obj.filename.stem).split("ex")[1].split("pc")[0]
                        #ion_time = str(mass_spectrum_obj.filename.stem).split("0pt")[1].split("s")[0]
                        peak_obj.predicted_std = mass_shift_ppp
                        
                        results_list.append( {
                        "ms_index_position" : peak_obj_idx,
                        "predicted_std": mass_shift_ppp,
                        "mz_exp": peak_obj.mz_exp,
                        "nominal_mz_exp": peak_obj.nominal_mz_exp,
                        "predicted_mz": mz_centroid[1],
                        "s2n" : peak_obj.signal_to_noise,
                        "peak_height" : peak_obj.abundance,
                        "predicted_peak_height" : abund_centroid[1],
                        "peak_height_error" : height_shift_per,
                        "resolving_power" : peak_obj.resolving_power,
                        #"excitation_amplitude" : excitation_amplitude,
                        #"ion_time" : ion_time
                        })
                        
                        indexes_without_results.remove(peak_obj_idx)
                    #elif len(mz_centroid) == 3 and len(abund_min_valley) != 2:

        for peak_obj_idx in indexes_without_results:

            results_list.append( {
            "ms_index_position" : peak_obj_idx,
            "mz_exp": self.mass_spectrum_obj[peak_obj_idx].mz_exp,
            "nominal_mz_exp": self.mass_spectrum_obj[peak_obj_idx].nominal_mz_exp,
            "s2n" : self.mass_spectrum_obj[peak_obj_idx].signal_to_noise,
            "peak_height" : self.mass_spectrum_obj[peak_obj_idx].abundance,
            "resolving_power" : self.mass_spectrum_obj[peak_obj_idx].resolving_power,
            #"excitation_amplitude" : excitation_amplitude,
            #"ion_time" : ion_time
            } )

        df = DataFrame(results_list).sort_values("mz_exp")
        
        df.interpolate(method ='linear', limit_direction ='backward',  inplace=True)
        df.interpolate(method ='linear', limit_direction ='forward',  inplace=True)

        #TODO improve interpolation for missing data
        #f1 = interpolate.interp1d(x1, y1, kind='quadratic',fill_value="extrapolate")

        
        for peak_obj_idx in indexes_without_results:

            predicted_std = df.loc[peak_obj_idx].predicted_std
            
            self.mass_spectrum_obj[peak_obj_idx].predicted_std = predicted_std

        return df

    def sum_data(self, tuple_mz_abun_list):

        all_mz = {}

        for mz_list, abun_list in tuple_mz_abun_list:
            
            for index, mz in enumerate(mz_list):

                abundance = abun_list[index]

                if mz in all_mz:
                    all_mz[mz] = all_mz[mz] + abundance    
                else: 
                    all_mz[mz] = abundance
        
        mz_all = []
        abun_all = []

        for mz in sorted (all_mz) : 
            mz_all.append(mz)
            abun_all.append(all_mz[mz])

        return array(mz_all), array(abun_all)    

    def calc_error(self, mass_ref, mass_sim, factor):

        return (mass_sim-mass_ref/mass_ref)*factor

    def find_peak_apex(self, mass, abund):

        dy = abund[1:] - abund[:-1]

        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]

        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf

        indexes = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indexes.size:
            
            return mass[indexes], abund[indexes]

    def find_peak_valley(self, mz, abund):

        dy = abund[1:] - abund[:-1]
        
        '''replaces nan for infinity'''
        indices_nan = where(isnan(abund))[0]
        
        if indices_nan.size:
            
            abund[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf
        
        indexes = where((hstack((dy, 0)) > 0) & (hstack((0, dy)) < 0))[0]

        return mz[indexes], abund[indexes]    