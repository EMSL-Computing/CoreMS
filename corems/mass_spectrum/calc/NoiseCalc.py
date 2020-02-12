import time

from numpy import where, average, std, isnan, inf, hstack

from corems.encapsulation.settings.processingSetting import MassSpectrumSetting

__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"


class NoiseThresholdCalc:

    def cut_mz_domain_noise(self, auto):
        
        if auto:
            
            # this calculation is taking too long (about 2 seconds)
            number_average_molecular_weight = self.weight_average_molecular_weight(
                profile=True)
           
            # +-200 is a guess for testing only, it needs adjustment for each type of analysis
            # need to check min mz here or it will break
            min_mz_noise = number_average_molecular_weight - 100
            # need to check max mz here or it will break
            max_mz_noise = number_average_molecular_weight + 100

            min_mz_whole_ms = self.mz_exp_profile.min()
            max_mz_whole_ms = self.mz_exp_profile.max()

            if min_mz_noise < min_mz_whole_ms:
                min_mz_noise = min_mz_whole_ms

            if max_mz_noise < max_mz_whole_ms:
                max_mz_noise = max_mz_whole_ms

        else:

            min_mz_noise = MassSpectrumSetting.min_noise_mz
            max_mz_noise = MassSpectrumSetting.max_noise_mz
            
        final = where(self.mz_exp_profile > min_mz_noise)[-1][-1]
        comeco = where(self.mz_exp_profile > min_mz_noise)[0][0]

        mz_domain_low_Y_cutoff = self.abundance_profile[comeco:final]

        final = where(self.mz_exp_profile < max_mz_noise)[-1][-1]
        comeco = where(self.mz_exp_profile < max_mz_noise)[0][0]

        return mz_domain_low_Y_cutoff[comeco:final]

    def simple_model_error_dist(self,  ymincentroid):
        
        import pymc3 as pm

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
            
            print(pm.summary(trace))

            return pm.summary(trace)['mean'].values[0] 
            

    def get_noise_average(self, ymincentroid, bayes=False):
        #assumes noise to be gaussian and estimate noise level by calculating the valley 
        # if bayes is enable it will model the valley distributuion as half-Normal and estimate the std
        
        average_noise = (ymincentroid*2).mean()
        
        if bayes:
            
            s_deviation = self.simple_model_error_dist(ymincentroid)
        
        else:
            
            s_deviation = ymincentroid.std() * 2
        
        return average_noise, s_deviation

    def get_abundance_minima_centroid(self, intes):

        maximum = intes.max()

        threshold_min = (maximum * 0.05)

        y = -intes

        dy = y[1:] - y[:-1]

        '''replaces NaN for Infinity'''
        indices_nan = where(isnan(y))[0]

        if indices_nan.size:

            y[indices_nan] = inf
            dy[where(isnan(dy))[0]] = inf

        indices = where((hstack((dy, 0)) < 0) & (hstack((0, dy)) > 0))[0]

        if indices.size and threshold_min is not None:
            indices = indices[intes[indices] <= threshold_min]

        return intes[indices]


    def run_noise_threshold_calc(self, auto, bayes=False):

        Y_cut = self.cut_mz_domain_noise(auto)
        
        if auto:

            yminima = self.get_abundance_minima_centroid(Y_cut)
            
            return self.get_noise_average(yminima, bayes=bayes)

        else:

            return self.get_noise_average(Y_cut, bayes=bayes)
