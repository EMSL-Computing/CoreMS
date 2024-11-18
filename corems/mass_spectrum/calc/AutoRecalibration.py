# -*- coding: utf-8 -*-
"""
Created on March 23 2023

@author: Will Kew

Modules for automatic mass internal recalibration
"""

from lmfit.models import GaussianModel
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
import copy


class HighResRecalibration:
    """
    This class is designed for high resolution (FTICR, Orbitrap) data of complex mixture, e.g. Organic matter

    The tool first does a broad mass range search for the most commonly expected ion type (i.e. CHO, deprotonated - for negative ESI)
    And then the assigned data mass error distribution is searched, with a gaussian fit to the most prominent range.
    This tool works when the data are of sufficient quality, and not outwith the typical expected range of the mass analyzer
    It presumes the mean error is out by 0-several ppm, but that the spread of error values is modest (<2ppm)

    Parameters
    ----------
    mass_spectrum : MassSpectrum
        CoreMS mass spectrum object
    plot : bool, optional
        Whether to plot the error distribution. The default is False.
    docker : bool, optional
        Whether to use the docker database. The default is True. If not, it uses a dynamically generated sqlite database.
    ppmFWHMprior : float, optional
        The FWHM of the prior distribution (ppm). The default is 3.
    ppmRangeprior : float, optional
        The range of the prior distribution (ppm). The default is 15.

    Methods
    --------
    * determine_error_boundaries(). Determine the error boundaries for recalibration space.

    Notes
    -----
    This initialisation function creates a copy of the MassSpectrum object to avoid over-writing assignments.
    Possible future task is to make the base class copyable.

    """

    def __init__(
        self,
        mass_spectrum,
        plot: bool = False,
        docker: bool = True,
        ppmFWHMprior: float = 3,
        ppmRangeprior: float = 15,
    ):
        self.mass_spectrum = copy.deepcopy(mass_spectrum)
        self.plot = plot
        self.docker = docker
        self.ppmFWHMprior = ppmFWHMprior
        self.ppmRangeprior = ppmRangeprior

    def set_uncal_settings(self):
        """Set uncalibrated formula search settings

        This function serves the uncalibrated data (hence broad error tolerance)
        It only allows CHO formula in deprotonated ion type- as most common for SRFA ESI negative mode

        This will not work for positive mode data, or for other ion types, or other expected elemental searches.

        """
        # TODO rework this.

        if self.docker:
            self.mass_spectrum.molecular_search_settings.url_database = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
        else:
            self.mass_spectrum.molecular_search_settings.url_database = None
        self.mass_spectrum.molecular_search_settings.error_method = None
        self.mass_spectrum.molecular_search_settings.score_method = "prob_score"

        self.mass_spectrum.molecular_search_settings.min_ppm_error = (
            -1 * self.ppmRangeprior / 2
        )  # -7.5
        self.mass_spectrum.molecular_search_settings.max_ppm_error = (
            self.ppmRangeprior / 2
        )  # 7.5

        self.mass_spectrum.molecular_search_settings.min_dbe = 0
        self.mass_spectrum.molecular_search_settings.max_dbe = 50

        self.mass_spectrum.molecular_search_settings.use_isotopologue_filter = False
        self.mass_spectrum.molecular_search_settings.min_abun_error = -30
        self.mass_spectrum.molecular_search_settings.max_abun_error = 70

        self.mass_spectrum.molecular_search_settings.use_min_peaks_filter = True
        self.mass_spectrum.molecular_search_settings.min_peaks_per_class = (
            10  # default is 15
        )

        self.mass_spectrum.molecular_search_settings.usedAtoms["C"] = (1, 90)
        self.mass_spectrum.molecular_search_settings.usedAtoms["H"] = (4, 200)
        self.mass_spectrum.molecular_search_settings.usedAtoms["O"] = (1, 23)
        self.mass_spectrum.molecular_search_settings.usedAtoms["N"] = (0, 0)
        self.mass_spectrum.molecular_search_settings.usedAtoms["S"] = (0, 0)
        self.mass_spectrum.molecular_search_settings.usedAtoms["P"] = (0, 0)

        self.mass_spectrum.molecular_search_settings.isProtonated = True
        self.mass_spectrum.molecular_search_settings.isRadical = False
        self.mass_spectrum.molecular_search_settings.isAdduct = False

    def positive_search_settings(self):
        """Set the positive mode elemental search settings"""
        self.mass_spectrum.molecular_search_settings.isProtonated = False
        self.mass_spectrum.molecular_search_settings.isAdduct = True
        self.mass_spectrum.molecular_search_settings.adduct_atoms_pos = ["Na"]

    @staticmethod
    def get_error_range(
        errors: list, ppmFWHMprior: float = 3, plot_logic: bool = False
    ):
        """Get the error range from the error distribution

        Using lmfit and seaborn kdeplot to extract the error range from the error distribution of assigned species.

        Parameters
        ----------
        errors : list
            list of the errors of the assigned species (ppm)
        ppmFWHMprior : float, optional
            The FWHM of the prior distribution (ppm). The default is 3.
        plot_logic : bool, optional
            Whether to plot the error distribution. The default is False.

        Returns
        -------
        mean_error : float
            mean mass error of the Gaussian distribution (ppm)
        fwhm_error : float
            full width half max of the gaussian error distribution (ppm)
        ppm_thresh : list
            recommended thresholds for the recalibration parameters (ppm)
            Consists of [mean_error-fwhm_error,mean_error+fwhm_error]

        """
        kde = sns.kdeplot(errors)

        kde_data = kde.get_lines()[0].get_data()

        tmpdf = pd.Series(index=kde_data[0], data=kde_data[1])
        kde_apex_ppm = tmpdf.idxmax()
        kde_apex_val = tmpdf.max()

        plt.close(kde.figure)
        plt.close("all")

        lmmodel = GaussianModel()
        lmpars = lmmodel.guess(kde_data[1], x=kde_data[0])
        lmpars["sigma"].value = 2.3548 / ppmFWHMprior
        lmpars["center"].value = kde_apex_ppm
        lmpars["amplitude"].value = kde_apex_val
        lmout = lmmodel.fit(kde_data[1], lmpars, x=kde_data[0])

        if plot_logic:
            fig, ax = plt.subplots(figsize=(8, 4))
            lmout.plot_fit(
                ax=ax, data_kws={"color": "tab:blue"}, fit_kws={"color": "tab:red"}
            )
            ax.set_xlabel("$m/z$ Error (ppm)")
            ax.set_ylabel("Density")
            plt.legend(facecolor="white", framealpha=0)

        mean_error = lmout.best_values["center"]
        std_error = lmout.best_values["sigma"]
        # FWHM from Sigma = approx. 2.355*sigma
        # fwhm_error = 2*np.sqrt(2*np.log(2))*std_error
        fwhm_error = std_error * np.sqrt(8 * np.log(2))

        ppm_thresh = [mean_error - fwhm_error, mean_error + fwhm_error]
        return mean_error, fwhm_error, ppm_thresh

    def determine_error_boundaries(self):
        """Determine the error boundaries for recalibration space

        This is the main function in this class
        Sets the Molecular Formulas search settings, performs the initial formula search
        Converts the data to a dataframe, and gets the error range
        Returns the error thresholds.

        Returns
        -------
        mean_error : float
            mean mass error of the Gaussian distribution (ppm)
        fwhm_error : float
            full width half max of the gaussian error distribution (ppm)
        ppm_thresh : list
            recommended thresholds for the recalibration parameters (ppm)
            Consists of [mean_error-fwhm_error,mean_error+fwhm_error]
        """

        # Set the search settings
        self.set_uncal_settings()

        # Set the positive mode settings
        # To do - have user defineable settings?
        if self.mass_spectrum.polarity == 1:
            self.positive_search_settings()

        # Search MFs
        SearchMolecularFormulas(
            self.mass_spectrum, first_hit=True
        ).run_worker_mass_spectrum()

        # Exporting to a DF is ~30x slower than just getting the errors, so this is fast.
        errors = []
        for mspeak in self.mass_spectrum.mspeaks:
            if len(mspeak.molecular_formulas) > 0:
                errors.append(mspeak.best_molecular_formula_candidate.mz_error)

        # If there are NO assignments, it'll fail on the next step. Need to check for that
        nassign = len(errors)
        # Here we say at least 5 features assigned are needed - it probably should be greater, but we are just trying to stop it breaking the code
        # We want to make sure the spectrum is capture in the database though - so we return the stats entries (0 assignments) and the number of assignments
        if nassign < 5:
            if self.mass_spectrum.parameters.mass_spectrum.verbose_processing:
                print("fewer than 5 peaks assigned, cannot determine error range")
            return np.nan, np.nan, [np.nan, np.nan]
        else:
            mean_error, fwhm_error, ppm_thresh = self.get_error_range(
                errors, self.ppmFWHMprior, self.plot
            )
            return mean_error, fwhm_error, ppm_thresh
