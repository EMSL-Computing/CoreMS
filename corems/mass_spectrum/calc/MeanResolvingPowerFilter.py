"""
Created on June 2nd 2023

@author: Will Kew

Module for mean resolving power filtration
Based upon the work in: 

Kanawati, B, Bader, TM, Wanczek, K-P, Li, Y, Schmitt-Kopplin, P. 
Fourier transform (FT)-artifacts and power-function resolution filter in Fourier transform mass spectrometry. 
Rapid Commun Mass Spectrom. 2017; 31: 1607- 1615. https://doi.org/10.1002/rcm.7940

Calculates a m/z normalised resolving power, fits a gaussian distribution to this, and then filters out peaks which are outside of the user defined number of standard deviations


"""

from lmfit.models import GaussianModel
from scipy import stats
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import copy

class MeanResolvingPowerFilter():
    '''
    This class covers HRMS filtering based on resolving power
    '''

    def __init__(self, mass_spectrum, ndeviations=3, plot=False, guess_pars = False):
        # we dont want the assignments made in this exploratory class to copy to the original object, so we make a copy of it.  
        # Possible future task - make mass spectrum base class copyable...
        #self.mass_spectrum
        self.mass_spectrum = mass_spectrum
        self.plot = plot
        self.ndeviations = ndeviations
        self.guess_pars = guess_pars

    def extract_peaks(self):
        '''
        Get the peak index, mz, and resolving power into a dataframe
        '''
        ids = []
        mzs = []
        rps = []
        for mspeak in self.mass_spectrum.mspeaks:
            ids.append(mspeak.index)
            mzs.append(mspeak.mz_exp)
            rps.append(mspeak.resolving_power)
        mzs = np.array(mzs)
        rps = np.array(rps)  

        tmpdf_ms = pd.DataFrame(index=ids,columns=['mz','rp','crp'])
        tmpdf_ms['mz'] = mzs
        tmpdf_ms['rp'] = rps
        return tmpdf_ms
    
    def normalise_rps(self,tmpdf_ms):
        '''
        Normalise the resolving powers to be indepent of m/z 
        '''
        if self.mass_spectrum.analyzer == 'ICR':
            tmpdf_ms['crp'] = tmpdf_ms['rp'] * np.sqrt(tmpdf_ms['mz']**2)
        else:
            print('Analyzer type not yet supported.')
        return tmpdf_ms

    def calculate_distribution(self,tmpdf_ms):
        '''
        Calculate the existing distribution of the resolving powers
        '''
        # Use Seaborn to create a KDE of the normalised resolving powers
        rps = sns.kdeplot(tmpdf_ms['crp']) 
        rps_data = rps.get_lines()[0].get_data()
        tmpdf = pd.Series(index=rps_data[0],data=rps_data[1])
        rps_apex_ppm = tmpdf.idxmax()
        rps_apex_val = tmpdf.max()
        plt.close(rps.figure)
        plt.close('all')

        # Use LMFIT to create a gaussian model of the distribution
        lmmodel = GaussianModel()
        lmpars = lmmodel.guess(rps_data[1], x=rps_data[0])
        if self.guess_pars:
            lmpars['sigma'].value = rps_data[0][-1]*0.01
            lmpars['center'].value = rps_apex_ppm
            lmpars['amplitude'].value = rps_apex_val
        lmout = lmmodel.fit(rps_data[1], lmpars, x=rps_data[0])

        if self.plot:
            fig,ax = plt.subplots(figsize=(8,4))
            lmout.plot_fit(ax=ax,data_kws ={'color':'tab:blue'},fit_kws ={'color':'tab:red'})
            ax.set_xlabel('Normalised Resolving Power')
            ax.set_ylabel('Density')
            plt.legend(facecolor='white', framealpha=0)

        mean_res = lmout.best_values['center']
        std_res = lmout.best_values['sigma']
        fwhm_res = std_res*np.sqrt(8*np.log(2))

        ndevs = self.ndeviations/2
        rps_thresh = [mean_res-(fwhm_res*ndevs),mean_res+(fwhm_res*ndevs)]
        return rps_thresh
    
    def create_index_list_to_remove(self,tmpdf_ms,rps_thresh):
        '''
        Subset the list of mspeaks to only the ones to keep, return an index list which can be passed back to the main 
        '''
        tmpdf_ms = tmpdf_ms[(tmpdf_ms['crp']<min(rps_thresh))|(tmpdf_ms['crp']>max(rps_thresh))]
        index_to_keep  = list(tmpdf_ms.index)
        return index_to_keep
    
    def main(self):
        tmpdf_ms = self.extract_peaks()
        tmpdf_ms  = self.normalise_rps(tmpdf_ms)
        rps_thresh = self.calculate_distribution(tmpdf_ms)
        index_to_remove = self.create_index_list_to_remove(tmpdf_ms,rps_thresh)
        return index_to_remove





