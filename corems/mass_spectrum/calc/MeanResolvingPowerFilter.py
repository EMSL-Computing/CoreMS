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

import warnings
from lmfit.models import GaussianModel
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class MeanResolvingPowerFilter:
    """Class for for mean resolving power filtration.

    This module implements a mean resolving power filter based on the work described [1]

    The MeanResolvingPowerFilter class provides methods to calculate the m/z normalized resolving power, fit a Gaussian distribution to it, and filter out peaks that are outside of the user-defined number of standard deviations.

    Attributes
    -------
    mass_spectrum (object): The mass spectrum object.
    ndeviations (int): The number of standard deviations used for filtering.
    plot (bool): Flag indicating whether to plot the results.
    guess_pars (bool): Flag indicating whether to guess the parameters for the Gaussian model.
    return_rps (bool): Flag indicating whether to return the calculated resolving powers, will return the NORMALISED resolving power, the mean, and the standard deviation.

    Methods
    ------
    * extract_peaks(): Extracts the peaks from the mass spectrum.
    * normalise_rps(tmpdf_ms): Normalizes the resolving powers to be independent of m/z.
    * calculate_distribution(tmpdf_ms): Calculates the distribution of the resolving powers.
    * create_index_list_to_remove(tmpdf_ms, rps_thresh): Creates an index list of peaks to remove based on the calculated thresholds.
    * main(): Executes the main filtering process and returns the index list of peaks to remove.

    References
    ----------
    1.  Kanawati, B, Bader, TM, Wanczek, K-P, Li, Y, Schmitt-Kopplin, P.
        Fourier transform (FT)-artifacts and power-function resolution filter in Fourier transform mass spectrometry.
        Rapid Commun Mass Spectrom. 2017; 31: 1607- 1615. https://doi.org/10.1002/rcm.7940
    """

    def __init__(
        self,
        mass_spectrum,
        ndeviations: float = 3,
        plot: bool = False,
        guess_pars: bool = False,
        return_rps: bool = False,
    ):
        # we dont want the assignments made in this exploratory class to copy to the original object, so we make a copy of it.
        # Possible future task - make mass spectrum base class copyable...
        # TODO see if there is redundancy in the AutoRecalibration function we can minimise here?
        # self.mass_spectrum
        self.mass_spectrum = mass_spectrum
        self.plot = plot
        self.ndeviations = ndeviations
        self.guess_pars = guess_pars
        self.return_rps = return_rps

    def extract_peaks(self):
        """Extracts the peaks from the mass spectrum.

        Returns
        ----------
        tmpdf_ms : Pandas DataFrame
            A DataFrame containing the extracted peaks.
        """
        ids = []
        mzs = []
        rps = []
        for mspeak in self.mass_spectrum.mspeaks:
            ids.append(mspeak.index)
            mzs.append(mspeak.mz_exp)
            rps.append(mspeak.resolving_power)
        mzs = np.array(mzs)
        rps = np.array(rps)

        tmpdf_ms = pd.DataFrame(index=ids, columns=["mz", "rp", "crp"])
        tmpdf_ms["mz"] = mzs
        tmpdf_ms["rp"] = rps
        return tmpdf_ms

    def normalise_rps(self, tmpdf_ms):
        """Normalizes the resolving powers to be independent of m/z.

        Parameters
        ------
        tmpdf_ms : Pandas DataFrame
            A DataFrame containing the extracted peaks.

        Returns
        --------
        tmpdf_ms : Pandas DataFrame
            A DataFrame with the resolving powers normalized.
        """

        if self.mass_spectrum.analyzer == "ICR":
            tmpdf_ms["crp"] = tmpdf_ms["rp"] * np.sqrt(tmpdf_ms["mz"] ** 2)
        else:
            warnings.warn(
                f"Analyzer type {self.mass_spectrum.analyzer} not yet supported.",
                UserWarning,
            )
        return tmpdf_ms

    def calculate_distribution(self, tmpdf_ms):
        """Calculates the distribution of the resolving powers.

        Parameters
        --------
        tmpdf_ms : Pandas DataFrame
            A DataFrame containing the extracted peaks with normalized resolving powers.

        Returns
        --------
        rps_thresh : list
            A list of the calculated thresholds for filtering.
        """

        # Use Seaborn to create a KDE of the normalised resolving powers
        rps = sns.kdeplot(tmpdf_ms["crp"])
        rps_data = rps.get_lines()[0].get_data()
        tmpdf = pd.Series(index=rps_data[0], data=rps_data[1])
        rps_apex_ppm = tmpdf.idxmax()
        rps_apex_val = tmpdf.max()
        plt.close(rps.figure)
        plt.close("all")

        # Use LMFIT to create a gaussian model of the distribution
        lmmodel = GaussianModel()
        lmpars = lmmodel.guess(rps_data[1], x=rps_data[0])
        if self.guess_pars:
            lmpars["sigma"].value = rps_data[0][-1] * 0.01
            lmpars["center"].value = rps_apex_ppm
            lmpars["amplitude"].value = rps_apex_val
        lmout = lmmodel.fit(rps_data[1], lmpars, x=rps_data[0])

        if self.plot:
            fig, ax = plt.subplots(figsize=(8, 4))
            lmout.plot_fit(
                ax=ax, data_kws={"color": "tab:blue"}, fit_kws={"color": "tab:red"}
            )
            ax.set_xlabel("Normalised Resolving Power")
            ax.set_ylabel("Density")
            plt.legend(facecolor="white", framealpha=0)

        mean_res = lmout.best_values["center"]
        std_res = lmout.best_values["sigma"]
        fwhm_res = std_res * np.sqrt(8 * np.log(2))

        ndevs = self.ndeviations / 2
        rps_thresh = [mean_res - (fwhm_res * ndevs), mean_res + (fwhm_res * ndevs)]
        if self.return_rps:
            return rps_thresh, mean_res, std_res
        return rps_thresh

    def create_index_list_to_remove(self, tmpdf_ms, rps_thresh: list):
        """Creates an index list of peaks to remove based on the calculated thresholds.

        Parameters
        ---------
        tmpdf_ms : Pandas DataFrame
            A DataFrame containing the extracted peaks with normalized resolving powers.
        rps_thresh : list
            A list of the calculated thresholds for filtering.

        Returns
        ----------
        index_to_keep :list
            A list of indices of peaks to keep.
        """
        # Subset the list of mspeaks to only the ones to keep, return an index list which can be passed back to the main

        tmpdf_ms = tmpdf_ms[
            (tmpdf_ms["crp"] < min(rps_thresh)) | (tmpdf_ms["crp"] > max(rps_thresh))
        ]
        index_to_keep = list(tmpdf_ms.index)
        return index_to_keep

    def main(self):
        """Executes the main filtering process and returns the index list of peaks to remove.

        Returns
        --------
        index_to_remove : list
            A list of indices of peaks to remove.
        """
        tmpdf_ms = self.extract_peaks()
        tmpdf_ms = self.normalise_rps(tmpdf_ms)
        if self.return_rps:
            rps_thresh, mean_res, std_res = self.calculate_distribution(tmpdf_ms)
        else:
            rps_thresh = self.calculate_distribution(tmpdf_ms)
        index_to_remove = self.create_index_list_to_remove(tmpdf_ms, rps_thresh)
        if self.return_rps:
            return index_to_remove, rps_thresh, mean_res, std_res
        else:
            return index_to_remove
