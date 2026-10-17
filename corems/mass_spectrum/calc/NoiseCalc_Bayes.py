"""
This code is for Bayesian estimation of the noise levels.
It is it not implemented or used in the current code base.
The packages it uses are not part of the requirements.
If you want to use it, you will need to install them manually.
"""

from corems.mass_spectrum.calc.NoiseCalc import NoiseThresholdCalc


class BayesNoiseCalc(NoiseThresholdCalc):
    def from_posterior(self, param, samples):
        """
        # Legacy code for Bayesian efforts - not used.
        pymc3 is not installed by default,
            if have plans to use it manual installation of pymc3
            package before using this method is needed
        """

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
        """
        # Legacy code for Bayesian efforts - not used.
        pymc3 is not installed by default,
            if have plans to use it manual installation of pymc3
            package before using this method is needed
        """
        import pymc3 as pm
        # from pymc3 import traceplot, plot_posterior

        with pm.Model() as model2:
            sd = self.from_posterior("sd", trace["sd"])
            y = pm.HalfNormal("y", sd=sd, observed=ymincentroid)
            start = pm.find_MAP()
            step = pm.NUTS()  # Hamiltonian MCMC with No U-Turn Sampler
            trace = pm.sample(
                1000, step, start, random_seed=123, progressbar=True, tune=1000
            )
            pm.summary(trace)
            # plot_posterior(trace)
            # traceplot(trace)
            return pm.summary(trace)["mean"].values[0]

    def simple_model_error_dist(self, ymincentroid):
        """
        # Legacy code for Bayesian efforts - not used.
        pymc3 is not installed by default,
            if have plans to use it manual installation of pymc3
            package before using this method is needed
        """
        import pymc3 as pm
        # from pymc3 import traceplot, plot_posterior
        # import seaborn as sns
        # f, ax = pyplot.subplots(figsize=(6, 6))
        # sns.distplot(ymincentroid)
        # sns.kdeplot(ymincentroid, ax=ax, shade=True, color="g")
        # sns.rugplot(ymincentroid, color="black", ax=ax)
        # ax.set(xlabel= "Peak Minima Magnitude", ylabel= "Density")
        # pyplot.show()

        with pm.Model() as model:
            # mu = pm.Uniform('mu', lower=-1, upper=1)
            lower = ymincentroid.min()
            upper = ymincentroid.max()

            sd = pm.Uniform("sd", lower=lower, upper=upper)

            y = pm.HalfNormal("y", sd=sd, observed=ymincentroid)

            start = pm.find_MAP()
            step = pm.NUTS()  # Hamiltonian MCMC with No U-Turn Sampler
            trace = pm.sample(
                1000, step, start, random_seed=123, progressbar=True, tune=1000
            )

            return pm.summary(trace)["mean"].values[0]
