import time

from numpy import where, average, std, isnan, inf, hstack, median, argmax, percentile
from corems import chunks
import warnings

#from matplotlib import pyplot
__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"

class NoiseThresholdCalc:

    def get_noise_threshold(self) -> ( (float, float), (float,float) ):
        ''' return two tuples (min_mz, max_mz) , (noise_threshold, noise_threshold)'''
        
        if self.is_centroid:

            x = min(self.mz_exp), max((self.mz_exp))
            
            if self.settings.threshold_method == 'auto':
                
                abundance_threshold = self.baselise_noise + (self.settings.noise_threshold_std * self.baselise_noise_std)
                y = (abundance_threshold, abundance_threshold)

            elif self.settings.threshold_method == 'signal_noise':

                normalized_threshold = (self.max_abundance * self.settings.s2n_threshold )/self.max_signal_to_noise
                y = (normalized_threshold, normalized_threshold)
            
            elif self.settings.threshold_method == "relative_abundance":

                normalized_threshold = (max(self.abundance)/100)*self.settings.relative_abundance_threshold
                y = (normalized_threshold, normalized_threshold)    
            
            else:
                    raise  Exception("%s method was not implemented, please refer to corems.mass_spectrum.calc.NoiseCalc Class" % self.settings.threshold_method)
                
            return x, y    

        else:

            if self.baselise_noise and self.baselise_noise_std:
                
                x = (self.mz_exp_profile.min(), self.mz_exp_profile.max())
                y = (self.baselise_noise_std, self.baselise_noise_std)
                
                if self.settings.threshold_method == 'auto':
                
                    #print(self.settings.noise_threshold_std)
                    abundance_threshold = self.baselise_noise + (self.settings.noise_threshold_std * self.baselise_noise_std)
                    
                    y = (abundance_threshold, abundance_threshold)

                elif self.settings.threshold_method == 'signal_noise':

                    max_sn = self.abundance_profile.max()/self.baselise_noise_std

                    normalized_threshold = (self.abundance_profile.max() * self.settings.s2n_threshold )/max_sn
                    y = (normalized_threshold, normalized_threshold)

                elif self.settings.threshold_method == "relative_abundance":

                    normalized_threshold = (self.abundance_profile.max()/100)*self.settings.relative_abundance_threshold
                    y = (normalized_threshold, normalized_threshold)
                    

                else:
                    raise  Exception("%s method was not implemented, \
                        please refer to corems.mass_spectrum.calc.NoiseCalc Class" % self.settings.threshold_method)
                
                return x, y
            
            else:
                
                warnings.warn(
                    "Noise Baseline and Noise std not specified,\
                    defaulting to 0,0 run process_mass_spec() ?"
                )    
                return (0,0) , (0,0)

    def cut_mz_domain_noise(self, auto):
        
        min_mz_whole_ms = self.mz_exp_profile.min()
        max_mz_whole_ms = self.mz_exp_profile.max()

        if auto:
            
            # this calculation is taking too long (about 2 seconds)
            number_average_molecular_weight = self.weight_average_molecular_weight(
                profile=True)
           
            # +-200 is a guess for testing only, it needs adjustment for each type of analysis
            # need to check min mz here or it will break
            min_mz_noise = number_average_molecular_weight - 100
            # need to check max mz here or it will break
            max_mz_noise = number_average_molecular_weight + 100

        else:

            min_mz_noise = self.settings.min_noise_mz
            max_mz_noise = self.settings.max_noise_mz

        if min_mz_noise < min_mz_whole_ms:
            min_mz_noise = min_mz_whole_ms

        if max_mz_noise > max_mz_whole_ms:
            max_mz_noise = max_mz_whole_ms

        #print(min_mz_noise, max_mz_noise)
        low_mz_index = (where(self.mz_exp_profile >= min_mz_noise)[0][0])
        #print(self.mz_exp_profile[low_mz_index])
        # low_mz_index = (argmax(self.mz_exp_profile <= min_mz_noise))
        
        high_mz_index = (where(self.mz_exp_profile <= max_mz_noise)[-1][-1])
        
        #high_mz_index = (argmax(self.mz_exp_profile <= max_mz_noise))
        
        if high_mz_index > low_mz_index:
            # pyplot.plot(self.mz_exp_profile[low_mz_index:high_mz_index], self.abundance_profile[low_mz_index:high_mz_index])
            # pyplot.show()
            return self.mz_exp_profile[high_mz_index:low_mz_index], self.abundance_profile[low_mz_index:high_mz_index]
        else:
            # pyplot.plot(self.mz_exp_profile[high_mz_index:low_mz_index], self.abundance_profile[high_mz_index:low_mz_index])
            # pyplot.show()
            return self.mz_exp_profile[high_mz_index:low_mz_index], self.abundance_profile[high_mz_index:low_mz_index]


    def from_posterior(self, param, samples):

        import pymc3 as pm
        import numpy as np
        import theano.tensor as tt
        from theano import as_op
        from scipy.stats import gaussian_kde
        
        smin, smax = np.min(samples), np.max(samples)
        width = smax - smin
        x = np.linspace(smin, smax, 100)
        y = gaussian_kde(samples)(x)
        
        # what was never sampled should have a small probability but not 0,
        # so we'll extend the domain and use linear approximation of density on it
        x = np.concatenate([[x[0] - 3 * width], x, [x[-1] + 3 * width]])
        y = np.concatenate([[0], y, [0]])
        
        return pm.distributions.Interpolated(param, x, y)

    def error_model_from_trace(self, trace, ymincentroid):
         
        import pymc3 as pm
        #from pymc3 import traceplot, plot_posterior
        
        with pm.Model() as model2:
            
            sd = self.from_posterior('sd', trace['sd'])
            y = pm.HalfNormal('y', sd=sd, observed=ymincentroid)
            start = pm.find_MAP()
            step = pm.NUTS() # Hamiltonian MCMC with No U-Turn Sampler
            trace = pm.sample(1000, step, start, random_seed=123, progressbar=True, tune=1000)
            pm.summary(trace)
            #plot_posterior(trace)
            #traceplot(trace)    
            return pm.summary(trace)['mean'].values[0] 

    def simple_model_error_dist(self,  ymincentroid):
        
        import pymc3 as pm
        # from pymc3 import traceplot, plot_posterior
        #import seaborn as sns
        #f, ax = pyplot.subplots(figsize=(6, 6))
        #sns.distplot(ymincentroid)
        #sns.kdeplot(ymincentroid, ax=ax, shade=True, color="g")
        #sns.rugplot(ymincentroid, color="black", ax=ax)
        #ax.set(xlabel= "Peak Minima Magnitude", ylabel= "Density")
        #pyplot.show()

        with pm.Model() as model:
            
            #mu = pm.Uniform('mu', lower=-1, upper=1)
            lower = ymincentroid.min()
            upper = ymincentroid.max()
            
            sd = pm.Uniform('sd', lower=lower , upper=upper)
            
            y = pm.HalfNormal('y', sd=sd, observed=ymincentroid)
            
            start = pm.find_MAP()
            step = pm.NUTS() # Hamiltonian MCMC with No U-Turn Sampler
            trace = pm.sample(1000, step, start, random_seed=123, progressbar=True, tune=1000)
            
            return pm.summary(trace)['mean'].values[0] 
            

    def get_noise_average(self, ymincentroid, auto=True, bayes=False):
        # assumes noise to be gaussian and estimate noise level by 
        # calculating the valley. If bayes is enable it will 
        # model the valley distributuion as half-Normal and estimate the std
        
        average_noise = median((ymincentroid))*2 if auto else median(ymincentroid)
        
        if bayes:
            
            s_deviation = self.simple_model_error_dist(ymincentroid)
        
        else:
            
            s_deviation = ymincentroid.std()*3 if auto else ymincentroid.std()
            
        return average_noise, s_deviation

    def get_abundance_minima_centroid(self, mz_cut, abun_cut):

        maximum = self.abundance_profile.max()
        threshold_min = (maximum * 1.00)

        y = -abun_cut

        dy = y[1:] - y[:-1]
        '''replaces NaN for Infinity'''
        indices_nan = where(isnan(y))[0]
        
        if indices_nan.size:

            y[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf

        
        indices = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]
        
        if indices.size and threshold_min is not None:
            indices = indices[abun_cut[indices] <= threshold_min]

        # pyplot.plot(mz_cut[indices], abun_cut[indices], linewidth=0, marker='o', color='red', markersize=2)
        # pyplot.show()

        # pyplot.hist(x=abun_cut[indices], bins='auto', color='#0504aa',
                            # alpha=0.7, rwidth=0.85)
        
        # cutoff = percentile(abun_cut[indices], 0.98)
        
        # print(cutoff)
        
        # pyplot.show()    
        return abun_cut[indices]


    def run_noise_threshold_calc(self, auto, bayes=False):
        
        if self.is_centroid:
            # calculates noise_baseline and noise_std
            # needed to run auto noise threshold mode
            # it is not used for signal to noise nor 
            # relative abudance methods
            abundances_chunks = chunks(self.abundance, 50)
            each_min_abund = [min(x) for x in abundances_chunks]

            return average(each_min_abund), std(each_min_abund)
        
        else:

            mz_cut, abundance_cut = self.cut_mz_domain_noise(auto)
            
            if auto:

                yminima = self.get_abundance_minima_centroid(mz_cut, abundance_cut)
                
                return self.get_noise_average(yminima, auto=auto, bayes=bayes)

            else:
                
                # pyplot.show()
                return self.get_noise_average(abundance_cut,auto=auto, bayes=bayes)
