__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

from corems.chroma_peak.calc.ChromaPeakCalc import (
    GCPeakCalculation,
    LCMSMassFeatureCalculation,
)
from corems.mass_spectra.factory.chromat_data import EIC_Data
from corems.molecular_id.factory.EI_SQL import LowResCompoundRef


class ChromaPeakBase:
    """Base class for chromatographic peak (ChromaPeak) objects.

    Parameters
    -------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    start_index : int
        The start index of the peak.
    index : int
        The index of the peak.
    final_index : int
        The final index of the peak.

    Attributes
    --------
    start_scan : int
        The start scan of the peak.
    final_scan : int
        The final scan of the peak.
    apex_scan : int
        The apex scan of the peak.
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum : MassSpectrum
        The mass spectrum object.
    _area : float
        The area of the peak.

    Properties
    --------
    * retention_time : float.
        The retention time of the peak.
    * tic : float.
        The total ion current of the peak.
    * area : float.
        The area of the peak.
    * rt_list : list.
        The list of retention times within the peak.
    * tic_list : list.
        The list of total ion currents within the peak.

    Methods
    --------
    * None
    """

    def __init__(
        self, chromatogram_parent, mass_spectrum_obj, start_index, index, final_index
    ):
        self.start_scan = start_index
        self.final_scan = final_index
        self.apex_scan = int(index)
        self.chromatogram_parent = chromatogram_parent
        self.mass_spectrum = mass_spectrum_obj
        self._area = None

    @property
    def retention_time(self):
        """Retention Time"""
        return self.mass_spectrum.retention_time

    @property
    def tic(self):
        """Total Ion Current"""
        return self.mass_spectrum.tic

    @property
    def area(self):
        """Peak Area"""
        return self._area

    @property
    def rt_list(self):
        """Retention Time List"""
        return [
            self.chromatogram_parent.retention_time[i]
            for i in range(self.start_scan, self.final_scan + 1)
        ]

    @property
    def tic_list(self):
        """Total Ion Current List"""
        return [
            self.chromatogram_parent.tic[i]
            for i in range(self.start_scan, self.final_scan + 1)
        ]


class LCMSMassFeature(ChromaPeakBase, LCMSMassFeatureCalculation):
    """Class representing a mass feature in a liquid chromatography (LC) chromatogram.

    Parameters
    -------
    lcms_parent : LCMS
        The parent LCMSBase object.
    mz : float
        The observed mass to charge ratio of the feature.
    retention_time : float
        The retention time of the feature (in minutes), at the apex.
    intensity : float
        The intensity of the feature.
    apex_scan : int
        The scan number of the apex of the feature.
    persistence : float, optional
        The persistence of the feature. Default is None.

    Attributes
    --------
    _mz_exp : float
        The observed mass to charge ratio of the feature.
    _mz_cal : float
        The calibrated mass to charge ratio of the feature.
    _retention_time : float
        The retention time of the feature (in minutes), at the apex.
    _apex_scan : int
        The scan number of the apex of the feature.
    _intensity : float
        The intensity of the feature.
    _persistence : float
        The persistence of the feature.
    _eic_data : EIC_Data
        The EIC data object associated with the feature.
    _dispersity_index : float
        The dispersity index of the feature.
    _half_height_width : numpy.ndarray
        The half height width of the feature (in minutes, as an array of min and max values).
    _tailing_factor : float
        The tailing factor of the feature.
        > 1 indicates tailing, < 1 indicates fronting, = 1 indicates symmetrical peak.
    _ms_deconvoluted_idx : [int]
        The indexes of the mass_spectrum attribute in the deconvoluted mass spectrum.
    is_calibrated : bool
        If True, the feature has been calibrated. Default is False.
    monoisotopic_mf_id : int
        Mass feature id that is the monoisotopic version of self.
        If self.id, then self is the monoisotopic feature). Default is None.
    isotopologue_type : str
        The isotopic class of the feature, i.e. "13C1", "13C2", "13C1 37Cl1" etc.
        Default is None.
    ms2_scan_numbers : list
        List of scan numbers of the MS2 spectra associated with the feature.
        Default is an empty list.
    ms2_mass_spectra : dict
        Dictionary of MS2 spectra associated with the feature (key = scan number for DDA).
        Default is an empty dictionary.
    ms2_similarity_results : list
        List of MS2 similarity results associated with the mass feature.
        Default is an empty list.
    id : int
        The ID of the feature, also the key in the parent LCMS object's
        `mass_features` dictionary.
    mass_spectrum_deconvoluted_parent : bool
        If True, the mass feature corresponds to the most intense peak in the deconvoluted mass spectrum. Default is None.
    associated_mass_features_deconvoluted : list
        List of mass features associated with the deconvoluted mass spectrum. Default is an empty list.

    """

    def __init__(
        self,
        lcms_parent,
        mz: float,
        retention_time: float,
        intensity: float,
        apex_scan: int,
        persistence: float = None,
        id: int = None,
    ):
        super().__init__(
            chromatogram_parent=lcms_parent,
            mass_spectrum_obj=None,
            start_index=None,
            index=apex_scan,
            final_index=None,
        )
        # Core attributes, marked as private
        self._mz_exp: float = mz
        self._mz_cal: float = None
        self._retention_time: float = retention_time
        self._apex_scan: int = apex_scan
        self._intensity: float = intensity
        self._persistence: float = persistence
        self._eic_data: EIC_Data = None
        self._dispersity_index: float = None
        self._half_height_width: np.ndarray = None
        self._ms_deconvoluted_idx = None

        # Additional attributes
        self.monoisotopic_mf_id = None
        self.isotopologue_type = None
        self.ms2_scan_numbers = []
        self.ms2_mass_spectra = {}
        self.ms2_similarity_results = []
        self.mass_spectrum_deconvoluted_parent: bool = None
        self.associated_mass_features_deconvoluted = []

        if id:
            self.id = id
        else:
            # get the parent's mass feature keys and add 1 to the max value to get the new key
            self.id = (
                max(lcms_parent.mass_features.keys()) + 1
                if lcms_parent.mass_features.keys()
                else 0
            )

    def update_mz(self):
        """Update the mass to charge ratio from the mass spectrum object."""
        if self.mass_spectrum is None:
            raise ValueError(
                "The mass spectrum object is not set, cannot update the m/z from the MassSpectrum object"
            )
        if len(self.mass_spectrum.mz_exp) == 0:
            raise ValueError(
                "The mass spectrum object has no m/z values, cannot update the m/z from the MassSpectrum object until it is processed"
            )
        new_mz = self.ms1_peak.mz_exp

        # calculate the difference between the new and old m/z, only update if it is close
        mz_diff = new_mz - self.mz
        if abs(mz_diff) < 0.01:
            self._mz_exp = new_mz

    def plot(self, to_plot=["EIC", "MS1", "MS2"], return_fig=True):
        """Plot the mass feature.

        Parameters
        ----------
        to_plot : list, optional
            List of strings specifying what to plot, any iteration of
            "EIC", "MS2", and "MS1".
            Default is ["EIC", "MS1", "MS2"].
        return_fig : bool, optional
            If True, the figure is returned. Default is True.

        Returns
        -------
        matplotlib.figure.Figure or None
            The figure object if `return_fig` is True.
            Otherwise None and the figure is displayed.
        """

        # EIC plot preparation
        eic_buffer_time = self.chromatogram_parent.parameters.lc_ms.eic_buffer_time

        # Adjust to_plot list if there are not spectra added to the mass features
        if self.mass_spectrum is None:
            to_plot = [x for x in to_plot if x != "MS1"]
        if len(self.ms2_mass_spectra) == 0:
            to_plot = [x for x in to_plot if x != "MS2"]
        if self._eic_data is None:
            to_plot = [x for x in to_plot if x != "EIC"]
        if self._ms_deconvoluted_idx is not None:
            deconvoluted = True
        else:
            deconvoluted = False

        fig, axs = plt.subplots(
            len(to_plot), 1, figsize=(9, len(to_plot) * 4), squeeze=False
        )
        fig.suptitle(
            "Mass Feature "
            + str(self.id)
            + ": m/z = "
            + str(round(self.mz, ndigits=4))
            + "; time = "
            + str(round(self.retention_time, ndigits=1))
            + " minutes"
        )

        i = 0
        # EIC plot
        if "EIC" in to_plot:
            if self._eic_data is None:
                raise ValueError(
                    "EIC data is not available, cannot plot the mass feature's EIC"
                )
            axs[i][0].set_title("EIC", loc="left")
            axs[i][0].plot(self._eic_data.time, self._eic_data.eic)
            if self.start_scan is not None:
                axs[i][0].fill_between(
                    self.eic_rt_list, self.eic_list, color="b", alpha=0.2
                )
            else:
                if self.chromatogram_parent.parameters.lc_ms.verbose_processing:
                    print(
                        "No start and final scan numbers were provided for mass feature "
                        + str(self.id)
                    )
            axs[i][0].set_ylabel("Intensity")
            axs[i][0].set_xlabel("Time (minutes)")
            axs[i][0].set_ylim(0, self.eic_list.max() * 1.1)
            axs[i][0].set_xlim(
                self.retention_time - eic_buffer_time,
                self.retention_time + eic_buffer_time,
            )
            axs[i][0].axvline(
                x=self.retention_time, color="k", label="MS1 scan time (apex)"
            )
            if len(self.ms2_scan_numbers) > 0:
                axs[i][0].axvline(
                    x=self.chromatogram_parent.get_time_of_scan_id(
                        self.best_ms2.scan_number
                    ),
                    color="grey",
                    linestyle="--",
                    label="MS2 scan time",
                )
            axs[i][0].legend(loc="upper left")
            axs[i][0].yaxis.get_major_formatter().set_useOffset(False)
            i += 1

        # MS1 plot
        if "MS1" in to_plot:
            if deconvoluted:
                axs[i][0].set_title("MS1 (deconvoluted)", loc="left")
                axs[i][0].vlines(
                    self.mass_spectrum.mz_exp,
                    0,
                    self.mass_spectrum.abundance,
                    color="k",
                    alpha=0.2,
                    label="Raw MS1",
                )
                axs[i][0].vlines(
                    self.mass_spectrum_deconvoluted.mz_exp,
                    0,
                    self.mass_spectrum_deconvoluted.abundance,
                    color="k",
                    label="Deconvoluted MS1",
                )
                axs[i][0].set_xlim(
                    self.mass_spectrum_deconvoluted.mz_exp.min() * 0.8,
                    self.mass_spectrum_deconvoluted.mz_exp.max() * 1.1,
                )
                axs[i][0].set_ylim(
                    0, self.mass_spectrum_deconvoluted.abundance.max() * 1.1
                )
            else:
                axs[i][0].set_title("MS1 (raw)", loc="left")
                axs[i][0].vlines(
                    self.mass_spectrum.mz_exp,
                    0,
                    self.mass_spectrum.abundance,
                    color="k",
                    label="Raw MS1",
                )
                axs[i][0].set_xlim(
                    self.mass_spectrum.mz_exp.min() * 0.8,
                    self.mass_spectrum.mz_exp.max() * 1.1,
                )
                axs[i][0].set_ylim(bottom=0)

            if (self.ms1_peak.mz_exp - self.mz) < 0.01:
                axs[i][0].vlines(
                    self.ms1_peak.mz_exp,
                    0,
                    self.ms1_peak.abundance,
                    color="m",
                    label="Feature m/z",
                )

            else:
                if self.chromatogram_parent.parameters.lc_ms.verbose_processing:
                    print(
                        "The m/z of the mass feature "
                        + str(self.id)
                        + " is different from the m/z of MS1 peak, the MS1 peak will not be plotted"
                    )
            axs[i][0].legend(loc="upper left")
            axs[i][0].set_ylabel("Intensity")
            axs[i][0].set_xlabel("m/z")
            axs[i][0].yaxis.set_tick_params(labelleft=False)
            i += 1

        # MS2 plot
        if "MS2" in to_plot:
            axs[i][0].set_title("MS2", loc="left")
            axs[i][0].vlines(
                self.best_ms2.mz_exp, 0, self.best_ms2.abundance, color="k"
            )
            axs[i][0].set_ylabel("Intensity")
            axs[i][0].set_xlabel("m/z")
            axs[i][0].set_ylim(bottom=0)
            axs[i][0].yaxis.get_major_formatter().set_scientific(False)
            axs[i][0].yaxis.get_major_formatter().set_useOffset(False)
            axs[i][0].set_xlim(
                self.best_ms2.mz_exp.min() * 0.8, self.best_ms2.mz_exp.max() * 1.1
            )
            axs[i][0].yaxis.set_tick_params(labelleft=False)

        # Add space between subplots
        plt.tight_layout()

        if return_fig:
            # Close figure
            plt.close(fig)
            return fig

    @property
    def mz(self):
        """Mass to charge ratio of the mass feature"""
        # If the mass feature has been calibrated, return the calibrated m/z, otherwise return the measured m/z
        if self._mz_cal is not None:
            return self._mz_cal
        else:
            return self._mz_exp

    @property
    def mass_spectrum_deconvoluted(self):
        """Returns the deconvoluted mass spectrum object associated with the mass feature, if deconvolution has been performed."""
        if self._ms_deconvoluted_idx is not None:
            ms_deconvoluted = copy.deepcopy(self.mass_spectrum)
            ms_deconvoluted.set_indexes(self._ms_deconvoluted_idx)
            return ms_deconvoluted
        else:
            raise ValueError(
                "Deconvolution has not been performed for mass feature " + str(self.id)
            )

    @property
    def retention_time(self):
        """Retention time of the mass feature"""
        return self._retention_time

    @retention_time.setter
    def retention_time(self, value):
        """Set the retention time of the mass feature"""
        if not isinstance(value, float):
            raise ValueError("The retention time of the mass feature must be a float")
        self._retention_time = value

    @property
    def apex_scan(self):
        """Apex scan of the mass feature"""
        return self._apex_scan

    @apex_scan.setter
    def apex_scan(self, value):
        """Set the apex scan of the mass feature"""
        if not isinstance(value, int):
            raise ValueError("The apex scan of the mass feature must be an integer")
        self._apex_scan = value

    @property
    def intensity(self):
        """Intensity of the mass feature"""
        return self._intensity

    @intensity.setter
    def intensity(self, value):
        """Set the intensity of the mass feature"""
        if not isinstance(value, float):
            raise ValueError("The intensity of the mass feature must be a float")
        self._intensity = value

    @property
    def persistence(self):
        """Persistence of the mass feature"""
        return self._persistence

    @persistence.setter
    def persistence(self, value):
        """Set the persistence of the mass feature"""
        if not isinstance(value, float):
            raise ValueError("The persistence of the mass feature must be a float")
        self._persistence = value

    @property
    def eic_rt_list(self):
        """Retention time list between the beginning and end of the mass feature"""
        # Find index of the start and final scans in the EIC data
        start_index = self._eic_data.scans.tolist().index(self.start_scan)
        final_index = self._eic_data.scans.tolist().index(self.final_scan)

        # Get the retention time list
        rt_list = self._eic_data.time[start_index : final_index + 1]
        return rt_list

    @property
    def eic_list(self):
        """EIC List between the beginning and end of the mass feature"""
        # Find index of the start and final scans in the EIC data
        start_index = self._eic_data.scans.tolist().index(self.start_scan)
        final_index = self._eic_data.scans.tolist().index(self.final_scan)

        # Get the retention time list
        eic = self._eic_data.eic[start_index : final_index + 1]
        return eic

    @property
    def ms1_peak(self):
        """MS1 peak from associated mass spectrum that is closest to the mass feature's m/z"""
        # Find index array self.mass_spectrum.mz_exp that is closest to self.mz
        closest_mz = min(self.mass_spectrum.mz_exp, key=lambda x: abs(x - self.mz))
        closest_mz_index = self.mass_spectrum.mz_exp.tolist().index(closest_mz)

        return self.mass_spectrum._mspeaks[closest_mz_index]

    @property
    def tailing_factor(self):
        """Tailing factor of the mass feature"""
        return self._tailing_factor

    @tailing_factor.setter
    def tailing_factor(self, value):
        """Set the tailing factor of the mass feature"""
        if not isinstance(value, float):
            raise ValueError("The tailing factor of the mass feature must be a float")
        self._tailing_factor = value

    @property
    def dispersity_index(self):
        """Dispersity index of the mass feature"""
        return self._dispersity_index

    @dispersity_index.setter
    def dispersity_index(self, value):
        """Set the dispersity index of the mass feature"""
        if not isinstance(value, float):
            raise ValueError("The dispersity index of the mass feature must be a float")
        self._dispersity_index = value

    @property
    def half_height_width(self):
        """Half height width of the mass feature, average of min and max values, in minutes"""
        return np.mean(self._half_height_width)

    @property
    def best_ms2(self):
        """Points to the best representative MS2 mass spectrum

        Notes
        -----
        If there is only one MS2 mass spectrum, it will be returned
        If there are MS2 similarity results, this will return the MS2 mass spectrum with the highest entropy similarity score.
        If there are no MS2 similarity results, the best MS2 mass spectrum is determined by the closest scan time to the apex of the mass feature, with higher resolving power.  Checks for and disqualifies possible chimeric spectra.

        Returns
        -------
        MassSpectrum or None
            The best MS2 mass spectrum.
        """
        if len(self.ms2_similarity_results) > 0:
            # the scan number with the highest similarity score
            results_df = [x.to_dataframe() for x in self.ms2_similarity_results]
            results_df = pd.concat(results_df)
            results_df = results_df.sort_values(
                by="entropy_similarity", ascending=False
            )
            best_scan_number = results_df.iloc[0]["query_spectrum_id"]
            return self.ms2_mass_spectra[best_scan_number]

        ms2_scans = list(self.ms2_mass_spectra.keys())
        if len(ms2_scans) > 1:
            mz_diff_list = []  # List of mz difference between mz of mass feature and mass of nearest mz in each scan
            res_list = []  # List of maximum resolving power of peaks in each scan
            time_diff_list = []  # List of time difference between scan and apex scan in each scan
            for scan in ms2_scans:
                if len(self.ms2_mass_spectra[scan].mspeaks) > 0:
                    # Find mz closest to mass feature mz, return both the difference in mass and its resolution
                    closest_mz = min(
                        self.ms2_mass_spectra[scan].mz_exp,
                        key=lambda x: abs(x - self.mz),
                    )
                    if all(
                        np.isnan(self.ms2_mass_spectra[scan].resolving_power)
                    ):  # All NA for resolving power in peaks, not uncommon in CID spectra
                        res_list.append(2)  # Assumes very low resolving power
                    else:
                        res_list.append(
                            np.nanmax(self.ms2_mass_spectra[scan].resolving_power)
                        )
                    mz_diff_list.append(np.abs(closest_mz - self.mz))
                    time_diff_list.append(
                        np.abs(
                            self.chromatogram_parent.get_time_of_scan_id(scan)
                            - self.retention_time
                        )
                    )
                else:
                    res_list.append(np.nan)
                    mz_diff_list.append(np.nan)
                    time_diff_list.append(np.nan)
            # Convert diff_lists into logical scores (higher is better for each score)
            time_score = 1 - np.array(time_diff_list) / np.nanmax(
                np.array(time_diff_list)
            )
            res_score = np.array(res_list) / np.nanmax(np.array(res_list))
            # mz_score is 0 for possible chimerics, 1 for all others (already within mass tolerance before assigning)
            mz_score = np.zeros(len(ms2_scans))
            for i in np.arange(0, len(ms2_scans)):
                if mz_diff_list[i] < 0.8 and mz_diff_list[i] > 0.1:  # Possible chimeric
                    mz_score[i] = 0
                else:
                    mz_score[i] = 1
            # get the index of the best score and return the mass spectrum
            if len([np.nanargmax(time_score * res_score * mz_score)]) == 1:
                return self.ms2_mass_spectra[
                    ms2_scans[np.nanargmax(time_score * res_score * mz_score)]
                ]
            # remove the mz_score condition and try again
            elif len(np.argmax(time_score * res_score)) == 1:
                return self.ms2_mass_spectra[
                    ms2_scans[np.nanargmax(time_score * res_score)]
                ]
            else:
                raise ValueError(
                    "No best MS2 mass spectrum could be found for mass feature "
                    + str(self.id)
                )
        elif len(ms2_scans) == 1:  # if only one ms2 spectra, return it
            return self.ms2_mass_spectra[ms2_scans[0]]
        else:  # if no ms2 spectra, return None
            return None


class GCPeak(ChromaPeakBase, GCPeakCalculation):
    """Class representing a peak in a gas chromatography (GC) chromatogram.

    Parameters
    ----------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object associated with the peak.
    indexes : tuple
        The indexes of the peak in the chromatogram.

    Attributes
    ----------
    _compounds : list
        List of compounds associated with the peak.
    _ri : float or None
        Retention index of the peak.

    Methods
    -------
    * __len__(). Returns the number of compounds associated with the peak.
    * __getitem__(position).  Returns the compound at the specified position.
    * remove_compound(compounds_obj). Removes the specified compound from the peak.
    * clear_compounds(). Removes all compounds from the peak.
    * add_compound(compounds_dict, spectral_similarity_scores, ri_score=None, similarity_score=None). Adds a compound to the peak with the specified attributes.
    * ri().  Returns the retention index of the peak.
    * highest_ss_compound(). Returns the compound with the highest spectral similarity score.
    * highest_score_compound(). Returns the compound with the highest similarity score.
    * compound_names(). Returns a list of names of compounds associated with the peak.
    """

    def __init__(self, chromatogram_parent, mass_spectrum_obj, indexes):
        self._compounds = []
        self._ri = None
        super().__init__(chromatogram_parent, mass_spectrum_obj, *indexes)

    def __len__(self):
        return len(self._compounds)

    def __getitem__(self, position):
        return self._compounds[position]

    def remove_compound(self, compounds_obj):
        self._compounds.remove(compounds_obj)

    def clear_compounds(self):
        self._compounds = []

    def add_compound(
        self,
        compounds_dict,
        spectral_similarity_scores,
        ri_score=None,
        similarity_score=None,
    ):
        """Adds a compound to the peak with the specified attributes.

        Parameters
        ----------
        compounds_dict : dict
            Dictionary containing the compound information.
        spectral_similarity_scores : dict
            Dictionary containing the spectral similarity scores.
        ri_score : float or None, optional
            The retention index score of the compound. Default is None.
        similarity_score : float or None, optional
            The similarity score of the compound. Default is None.
        """
        compound_obj = LowResCompoundRef(compounds_dict)
        compound_obj.spectral_similarity_scores = spectral_similarity_scores
        compound_obj.spectral_similarity_score = spectral_similarity_scores.get(
            "cosine_correlation"
        )
        # TODO check is the above line correct?
        compound_obj.ri_score = ri_score
        compound_obj.similarity_score = similarity_score
        self._compounds.append(compound_obj)
        if similarity_score:
            self._compounds.sort(key=lambda c: c.similarity_score, reverse=True)
        else:
            self._compounds.sort(
                key=lambda c: c.spectral_similarity_score, reverse=True
            )

    @property
    def ri(self):
        """Returns the retention index of the peak.

        Returns
        -------
        float or None
            The retention index of the peak.
        """
        return self._ri

    @property
    def highest_ss_compound(self):
        """Returns the compound with the highest spectral similarity score.

        Returns
        -------
        LowResCompoundRef or None
            The compound with the highest spectral similarity score.
        """
        if self:
            return max(self, key=lambda c: c.spectral_similarity_score)
        else:
            return None

    @property
    def highest_score_compound(self):
        """Returns the compound with the highest similarity score.

        Returns
        -------
        LowResCompoundRef or None
            The compound with the highest similarity score.
        """
        if self:
            return max(self, key=lambda c: c.similarity_score)
        else:
            return None

    @property
    def compound_names(self):
        """Returns a list of names of compounds associated with the peak.

        Returns
        -------
        list
            List of names of compounds associated with the peak.
        """
        if self:
            return [c.name for c in self]
        else:
            return []


class GCPeakDeconvolved(GCPeak):
    """Represents a deconvolved peak in a chromatogram.

    Parameters
    ----------
    chromatogram_parent : Chromatogram
        The parent chromatogram object.
    mass_spectra : list
        List of mass spectra associated with the peak.
    apex_index : int
        Index of the apex mass spectrum in the `mass_spectra` list.
    rt_list : list
        List of retention times.
    tic_list : list
        List of total ion currents.
    """

    def __init__(
        self, chromatogram_parent, mass_spectra, apex_index, rt_list, tic_list
    ):
        self._ri = None
        self._rt_list = list(rt_list)
        self._tic_list = list(tic_list)
        self.mass_spectra = list(mass_spectra)
        super().__init__(
            chromatogram_parent,
            self.mass_spectra[apex_index],
            (0, apex_index, len(self.mass_spectra) - 1),
        )

    @property
    def rt_list(self):
        """Get the list of retention times.

        Returns
        -------
        list
            The list of retention times.
        """
        return self._rt_list

    @property
    def tic_list(self):
        """Get the list of total ion currents.

        Returns
        -------
        list
            The list of total ion currents.
        """
        return self._tic_list
