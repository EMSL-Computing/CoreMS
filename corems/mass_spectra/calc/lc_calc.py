import numpy as np
import pandas as pd
import warnings, scipy, multiprocessing
from ripser import ripser
from scipy import sparse
from scipy.spatial import KDTree
from sklearn.svm import SVR
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
from tqdm import tqdm

from corems.chroma_peak.factory.chroma_peak_classes import LCMSMassFeature
from corems.mass_spectra.calc import SignalProcessing as sp
from corems.mass_spectra.factory.chromat_data import EIC_Data
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile

warnings.filterwarnings("ignore", category=RuntimeWarning)

def find_closest(A, target):
    """Find the index of closest value in A to each value in target.

    Parameters
    ----------
    A : :obj:`~numpy.array`
        The array to search (blueprint). A must be sorted.
    target : :obj:`~numpy.array`
        The array of values to search for. target must be sorted.

    Returns
    -------
    :obj:`~numpy.array`
        The indices of the closest values in A to each value in target.
    """
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A) - 1)
    left = A[idx - 1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


class LCCalculations:
    """Methods for performing LC calculations on mass spectra data.

    Notes
    -----
    This class is intended to be used as a mixin for the LCMSBase class.

    Methods
    -------
    * get_max_eic(eic_data).
        Returns the maximum EIC value from the given EIC data. A static method.
    * smooth_tic(tic).
        Smooths the TIC data using the specified smoothing method and settings.
    * eic_centroid_detector(rt, eic, max_eic).
        Performs EIC centroid detection on the given EIC data.
    * find_nearest_scan(rt).
        Finds the nearest scan to the given retention time.
    * get_average_mass_spectrum(scan_list, apex_scan, spectrum_mode="profile", ms_level=1, auto_process=True, use_parser=False, perform_checks=True, polarity=None).
        Returns an averaged mass spectrum object.
    * find_mass_features(ms_level=1).
        Find regions of interest for a given MS level (default is MS1).
    * integrate_mass_features(drop_if_fail=False, ms_level=1).
        Integrate mass features of interest and extracts EICs.
    * find_c13_mass_features().
        Evaluate mass features and mark likely C13 isotopes.
    * deconvolute_ms1_mass_features().
        Deconvolute mass features' ms1 mass spectra.
    """

    @staticmethod
    def get_max_eic(eic_data: dict):
        """Returns the maximum EIC value from the given EIC data.

        Notes
        -----
        This is a static method.

        Parameters
        ----------
        eic_data : dict
            A dictionary containing EIC data.

        Returns
        -------
        float
            The maximum EIC value.
        """
        max_eic = 0
        for eic_data in eic_data.values():
            ind_max_eic = max(eic_data.get("EIC"))
            max_eic = ind_max_eic if ind_max_eic > max_eic else max_eic

        return max_eic

    def smooth_tic(self, tic):
        """Smooths the TIC or EIC data using the specified smoothing method and settings.

        Parameters
        ----------
        tic : numpy.ndarray
            The TIC (or EIC) data to be smoothed.

        Returns
        -------
        numpy.ndarray
            The smoothed TIC data.
        """
        implemented_smooth_method = self.parameters.lc_ms.implemented_smooth_method

        pol_order = self.parameters.lc_ms.savgol_pol_order

        window_len = self.parameters.lc_ms.smooth_window

        window = self.parameters.lc_ms.smooth_method

        return sp.smooth_signal(
            tic, window_len, window, pol_order, implemented_smooth_method
        )

    def eic_centroid_detector(self, rt, eic, max_eic, apex_indexes=[]):
        """Performs EIC centroid detection on the given EIC data.

        Parameters
        ----------
        rt : numpy.ndarray
            The retention time data.
        eic : numpy.ndarray
            The EIC data.
        max_eic : float
            The maximum EIC value.
        apex_indexes : list, optional
            The apexes of the EIC peaks. Defaults to [], which means that the apexes will be calculated by the function.

        Returns
        -------
        numpy.ndarray
            The indexes of left, apex, and right limits as a generator.
        """

        max_prominence = self.parameters.lc_ms.peak_max_prominence_percent

        max_height = self.parameters.lc_ms.peak_height_max_percent

        signal_threshold = self.parameters.lc_ms.eic_signal_threshold

        min_peak_datapoints = self.parameters.lc_ms.min_peak_datapoints

        peak_derivative_threshold = self.parameters.lc_ms.peak_derivative_threshold

        include_indexes = sp.peak_picking_first_derivative(
            domain=rt,
            signal=eic,
            max_height=max_height,
            max_prominence=max_prominence,
            max_signal=max_eic,
            min_peak_datapoints=min_peak_datapoints,
            peak_derivative_threshold=peak_derivative_threshold,
            signal_threshold=signal_threshold,
            correct_baseline=False,
            plot_res=False,
            apex_indexes=apex_indexes,
        )
        #include_indexes is a generator of tuples (left_index, apex_index, right_index)
        include_indexes = list(include_indexes)
        # Add check to make sure that there are at least 1/2 of min_peak_datapoints on either side of the apex
        indicies = [x for x in include_indexes]
        for idx in indicies:
            if (idx[1] - idx[0] < min_peak_datapoints / 2) or (
                idx[2] - idx[1] < min_peak_datapoints / 2
            ):
                include_indexes.remove(idx)
        return include_indexes

    def find_nearest_scan(self, rt):
        """Finds the nearest scan to the given retention time.

        Parameters
        ----------
        rt : float
            The retention time (in minutes) to find the nearest scan for.

        Returns
        -------
        int
            The scan number of the nearest scan.
        """
        array_rt = np.array(self.retention_time)

        scan_index = (np.abs(array_rt - rt)).argmin()

        real_scan = self.scans_number[scan_index]

        return real_scan

    def add_peak_metrics(self, remove_by_metrics=True, induced_features=False):
        """Add peak metrics to the mass features.

        This function calculates the peak metrics for each mass feature and adds them to the mass feature objects.

        Parameters
        ----------
        remove_by_metrics : bool, optional
            If True, remove mass features based on their peak metrics such as S/N, Gaussian similarity,
            dispersity index, and noise score. Default is True, which checks the setting in the processing parameters.
            If False, peak metrics are calculated but no mass features are removed, regardless of the setting in the processing parameters.
        induced_features : bool, optional
            Whether the mass features to be integrated were induced. Default is False.
        """
        # Check that at least some mass features have eic data
        if induced_features:
            mf_dict_values = self.induced_mass_features.values()
        else:
            mf_dict_values = self.mass_features.values()

        if not any([mf._eic_data is not None for mf in mf_dict_values]):
            raise ValueError(
                "No mass features have EIC data. Run integrate_mass_features first."
            )

        for mass_feature in mf_dict_values:
            # Check if the mass feature has been integrated
            if mass_feature._eic_data is not None and mass_feature.area is not None:
                # Calculate peak metrics
                mass_feature.calc_half_height_width()
                mass_feature.calc_tailing_factor()
                mass_feature.calc_dispersity_index()
                mass_feature.calc_gaussian_similarity()
                mass_feature.calc_noise_score()
        
        # Remove mass features by peak metrics if designated in parameters
        if self.parameters.lc_ms.remove_mass_features_by_peak_metrics and remove_by_metrics:
            self._remove_mass_features_by_peak_metrics(induced_features=induced_features)

    def get_average_mass_spectrum(
        self,
        scan_list,
        apex_scan,
        spectrum_mode="profile",
        ms_level=1,
        auto_process=True,
        use_parser=False,
        perform_checks=True,
        polarity=None,
        ms_params=None,
    ):
        """Returns an averaged mass spectrum object

        Parameters
        ----------
        scan_list : list
            List of scan numbers to average.
        apex_scan : int
            Number of the apex scan
        spectrum_mode : str, optional
            The spectrum mode to use. Defaults to "profile". Not that only "profile" mode is supported for averaging.
        ms_level : int, optional
            The MS level to use. Defaults to 1.
        auto_process : bool, optional
            If True, the averaged mass spectrum will be auto-processed. Defaults to True.
        use_parser : bool, optional
            If True, the mass spectra will be obtained from the parser. Defaults to False.
        perform_checks : bool, optional
            If True, the function will check if the data are within the ms_unprocessed dictionary and are the correct mode. Defaults to True. Only set to False if you are sure the data are profile, and (if not using the parser) are in the ms_unprocessed dictionary!  ms_unprocessed dictionary also must be indexed on scan
        polarity : int, optional
            The polarity of the mass spectra (1 or -1). If not set, the polarity will be determined from the dataset. Defaults to None. (fastest if set to -1 or 1)
        ms_params : MSParameters, optional
            The mass spectrum parameters to use. If not set (None), the globally set parameters will be used. Defaults to None.

        Returns
        -------
        MassSpectrumProfile
            The averaged mass spectrum object.

        Raises
        ------
        ValueError
            If the spectrum mode is not "profile".
            If the MS level is not found in the unprocessed mass spectra dictionary.
            If not all scan numbers are found in the unprocessed mass spectra dictionary.
        """
        if perform_checks:
            if spectrum_mode != "profile":
                raise ValueError("Averaging only supported for profile mode")

        if polarity is None:
            # set polarity to -1 if negative mode, 1 if positive mode (for mass spectrum creation)
            if self.polarity == "negative":
                polarity = -1
            elif self.polarity == "positive":
                polarity = 1
            else:
                raise ValueError(
                    "Polarity not set for dataset, must be a set containing either 'positive' or 'negative'"
                )

        # if not using_parser, check that scan numbers are in _ms_unprocessed
        if not use_parser:
            if perform_checks:
                # Set index to scan for faster lookup
                ms_df = (
                    self._ms_unprocessed[ms_level]
                    .copy()
                    .set_index("scan", drop=False)
                    .sort_index()
                )
                my_ms_df = ms_df.loc[scan_list]
                # Check that all scan numbers are in the ms_df
                if not all(np.isin(scan_list, ms_df.index)):
                    raise ValueError(
                        "Not all scan numbers found in the unprocessed mass spectra dictionary"
                    )
            else:
                my_ms_df = (
                    pd.DataFrame({"scan": scan_list})
                    .set_index("scan")
                    .join(self._ms_unprocessed[ms_level], how="left")
                )

        if use_parser:
            ms_list = [
                self.spectra_parser.get_mass_spectrum_from_scan(
                    x, spectrum_mode=spectrum_mode, auto_process=False
                )
                for x in scan_list
            ]
            ms_mz = [x._mz_exp for x in ms_list]
            ms_int = [x._abundance for x in ms_list]
            my_ms_df = []
            for i in np.arange(len(ms_mz)):
                my_ms_df.append(
                    pd.DataFrame(
                        {"mz": ms_mz[i], "intensity": ms_int[i], "scan": scan_list[i]}
                    )
                )
            my_ms_df = pd.concat(my_ms_df)

        if not self.check_if_grid(my_ms_df):
            my_ms_df = self.grid_data(my_ms_df)

        my_ms_ave = my_ms_df.groupby("mz")["intensity"].sum().reset_index()

        ms = ms_from_array_profile(
            my_ms_ave.mz,
            my_ms_ave.intensity,
            self.file_location,
            polarity=polarity,
            auto_process=False,
        )

        # Set the mass spectrum parameters, auto-process if auto_process is True, and add to the dataset
        if ms is not None:
            if ms_params is not None:
                ms.parameters = ms_params
            ms.scan_number = apex_scan
            if auto_process:
                ms.process_mass_spec()
        return ms

    def find_mass_features(self, ms_level=1, grid=True, assign_ms2_scans=False, ms2_scan_filter=None, 
                          targeted_search=False, target_search_dict=None, accumulate_features=False):
        """Find mass features within an LCMSBase object

        Note that this is a wrapper function that calls the find_mass_features_ph function, but can be extended to support other peak picking methods in the future.

        Parameters
        ----------
        ms_level : int, optional
            The MS level to use for peak picking Default is 1.
        grid : bool, optional
            If True, will regrid the data before running the persistent homology calculations (after checking if the data is gridded),
            used for persistent homology peak picking for profile data only. Default is True.
        assign_ms2_scans : bool, optional
            If True, assign MS2 scan numbers to mass features after peak picking.
            This populates the ms2_scan_numbers attribute on each mass feature, which enables
            choosing representative features based on MS2 availability. Default is False.
        ms2_scan_filter : str or None, optional
            Filter string for MS2 scans when assign_ms2_scans is True (e.g., 'hcd').
            If None, all MS2 scans are considered. Default is None.
        targeted_search : bool, optional
            If True, perform targeted mass feature search using the target_search_dict.
            This mode filters data to only m/z and RT windows of interest and bypasses
            intensity and persistence thresholds. Default is False.
        target_search_dict : dict or None, optional
            Dictionary containing target search parameters. Required if targeted_search is True.
            Must contain:
                - 'target_mz_list': list of target m/z values
                - 'target_rt_list': list of target retention times (in minutes)
                - 'mz_tolerance_ppm': m/z tolerance in ppm
                - 'rt_tolerance': retention time tolerance (in minutes)
            Optionally can contain:
                - 'type': type label for mass features (e.g., "internal standard")
                  If not provided, defaults to "targeted"
            Default is None.
        accumulate_features : bool, optional
            If True, new mass features will be added to existing features rather than replacing them.
            This allows multiple sequential calls to find_mass_features to build up a combined set.
            Default is False (replace existing features for backwards compatibility).

        Raises
        ------
        ValueError
            If no MS level data is found on the object.
            If persistent homology peak picking is attempted on non-profile mode data.
            If data is not gridded and grid is False.
            If peak picking method is not implemented.
            If targeted_search is True but target_search_dict is None or invalid.

        Returns
        -------
        None, but assigns the mass_features and eics attributes to the object.

        """
        # Validate targeted search parameters
        if targeted_search:
            if target_search_dict is None:
                raise ValueError("target_search_dict must be provided when targeted_search is True")
            required_keys = ['target_mz_list', 'target_rt_list', 'mz_tolerance_ppm', 'rt_tolerance']
            for key in required_keys:
                if key not in target_search_dict:
                    raise ValueError(f"target_search_dict must contain '{key}'")
            if len(target_search_dict['target_mz_list']) != len(target_search_dict['target_rt_list']):
                raise ValueError("target_mz_list and target_rt_list must have the same length")
        
        pp_method = self.parameters.lc_ms.peak_picking_method

        if pp_method == "persistent homology":
            msx_scan_df = self.scan_df[self.scan_df["ms_level"] == ms_level]
            if all(msx_scan_df["ms_format"] == "profile"):
                # Determine mass feature type
                if targeted_search:
                    mf_type = target_search_dict.get('type', 'targeted')
                else:
                    mf_type = 'untargeted'
                self.find_mass_features_ph(ms_level=ms_level, grid=grid, 
                                          targeted_search=targeted_search, 
                                          target_search_dict=target_search_dict,
                                          mf_type=mf_type,
                                          accumulate_features=accumulate_features)
            else:
                raise ValueError(
                    "MS{} scans are not profile mode, which is required for persistent homology peak picking.".format(
                        ms_level
                    )
                )
        elif pp_method == "centroided_persistent_homology":
            msx_scan_df = self.scan_df[self.scan_df["ms_level"] == ms_level]
            if all(msx_scan_df["ms_format"] == "centroid"):
                # Determine mass feature type
                if targeted_search:
                    mf_type = target_search_dict.get('type', 'targeted')
                else:
                    mf_type = 'untargeted'
                self.find_mass_features_ph_centroid(ms_level=ms_level, 
                                                    targeted_search=targeted_search, 
                                                    target_search_dict=target_search_dict,
                                                    mf_type=mf_type,
                                                    accumulate_features=accumulate_features)
            else:
                raise ValueError(
                    "MS{} scans are not centroid mode, which is required for persistent homology centroided peak picking.".format(
                        ms_level
                    )
                )
        else:
            raise ValueError("Peak picking method not implemented")
        
        # Cluster mass features to remove redundant features
        self.cluster_mass_features(drop_children=True)
        
        # Optionally assign MS2 scan numbers to mass features during peak picking
        # This helps with choosing representative features that have MS2 data
        if assign_ms2_scans:
            try:
                self._find_ms2_scans_for_mass_features(
                    mf_ids=None,  # Process all mass features
                    scan_filter=ms2_scan_filter
                )
            except ValueError:
                # No MS2 scans found - this is okay, just skip
                pass
        
        # Remove noisey mass features if designated in parameters
        if self.parameters.lc_ms.remove_redundant_mass_features and not targeted_search:
            self._remove_redundant_mass_features()

    def integrate_mass_features(
        self, drop_if_fail=True, drop_duplicates=True, ms_level=1, induced_features=False
    ):
        """Integrate mass features and extract EICs.

        Populates the _eics attribute on the LCMSBase object for each unique mz in the mass_features dataframe and adds data (start_scan, final_scan, area) to the mass_features attribute.

        Parameters
        ----------
        drop_if_fail : bool, optional
            Whether to drop mass features if the EIC limit calculations fail.
            Default is True.
        drop_duplicates : bool, optional
            Whether to mass features that appear to be duplicates
            (i.e., mz is similar to another mass feature and limits of the EIC are similar or encapsulating).
            Default is True.
        ms_level : int, optional
            The MS level to use. Default is 1.
        induced_features : bool, optional
            Whether the mass features to be intergrated were induced. Default is False.

        Raises
        ------
        ValueError
            If no mass features are found.
            If no MS level data is found for the given MS level (either in data or in the scan data)

        Returns
        -------
        None, but populates the eics attribute on the LCMSBase object and adds data (start_scan, final_scan, area) to the mass_features attribute.

        Notes
        -----
        drop_if_fail is useful for discarding mass features that do not have good shapes, usually due to a detection on a shoulder of a peak or a noisy region (especially if minimal smoothing is used during mass feature detection).
        """
        
        # Check if there is data
        if ms_level in self._ms_unprocessed.keys():
            raw_data = self._ms_unprocessed[ms_level].copy()
        else:
            raise ValueError("No MS level " + str(ms_level) + " data found")

        # Check if mass_spectrum exists on each mass feature
        if induced_features:
            mf_dict = self.induced_mass_features
            if len(mf_dict) == 0:
                raise ValueError(
                    "No induced mass features found, did you run fill_missing_cluster_features() first?"
                )

            ## remove not found induced mass features by mz <= 0 (-99 indicator)
            # also remove any where mz is nan
            mf_dict = {k:v for k, v in mf_dict.items() if v.mz > 0 and not np.isnan(v.mz)}

        else:
            mf_dict = self.mass_features
            if len(mf_dict) == 0:
                raise ValueError(
                    "No mass features found, did you run find_mass_features() first?"
                )

        # Subset scan data to only include correct ms_level
        scan_df_sub = self.scan_df[
            self.scan_df["ms_level"] == int(ms_level)
        ].reset_index(drop=True)
        if scan_df_sub.empty:
            raise ValueError("No MS level " + ms_level + " data found in scan data")
        scan_df_sub = scan_df_sub[["scan", "scan_time"]].copy()

        mzs_to_extract = np.unique([mf.mz for mf in mf_dict.values()])
        mzs_to_extract.sort()

        # Pre-sort raw_data by mz for faster filtering
        raw_data_sorted = raw_data.sort_values(["mz", "scan"]).reset_index(drop=True)
        raw_data_mz = raw_data_sorted["mz"].values

        # Get EICs for each unique mz in mass features list
        for mz in mzs_to_extract:
            mz_max = mz + self.parameters.lc_ms.eic_tolerance_ppm * mz / 1e6
            mz_min = mz - self.parameters.lc_ms.eic_tolerance_ppm * mz / 1e6

            # Use binary search for faster mz range filtering
            left_idx = np.searchsorted(raw_data_mz, mz_min, side="left")
            right_idx = np.searchsorted(raw_data_mz, mz_max, side="right")
            raw_data_sub = raw_data_sorted.iloc[left_idx:right_idx].copy()

            raw_data_sub = (
                raw_data_sub.groupby(["scan"])["intensity"].sum().reset_index()
            )
            raw_data_sub = scan_df_sub.merge(raw_data_sub, on="scan", how="left")
            raw_data_sub["intensity"] = raw_data_sub["intensity"].fillna(0)
            myEIC = EIC_Data(
                scans=raw_data_sub["scan"].values,
                time=raw_data_sub["scan_time"].values,
                eic=raw_data_sub["intensity"].values,
            )
            # Smooth EIC
            smoothed_eic = self.smooth_tic(myEIC.eic)
            smoothed_eic[smoothed_eic < 0] = 0
            myEIC.eic_smoothed = smoothed_eic
            self.eics[mz] = myEIC

        # Get limits of mass features using EIC centroid detector and integrate
        for idx, mass_feature in list(mf_dict.items()):
            mz = mass_feature.mz
            apex_scan = mass_feature.apex_scan

            # Pull EIC data and find apex scan index
            myEIC = self.eics[mz]
            mf_dict[idx]._eic_data = myEIC
            mf_dict[idx]._eic_mz = mz
            apex_index = np.searchsorted(myEIC.scans, apex_scan)

            # Find left and right limits of peak using EIC centroid detector, add to EICData
            centroid_eics = self.eic_centroid_detector(
                myEIC.time,
                myEIC.eic_smoothed,
                mass_feature.intensity * 1.1,
                apex_indexes=[int(apex_index)],
            )
            l_a_r_scan_idx = [i for i in centroid_eics]
            if len(l_a_r_scan_idx) > 0:
                # Calculate number of consecutive scans with intensity > 0 and check if it is above the minimum consecutive scans
                # Find the number of consecutive non-zero values in the EIC segment
                mask = myEIC.eic[l_a_r_scan_idx[0][0] : l_a_r_scan_idx[0][2] + 1] > 0
                # Find the longest run of consecutive True values
                if np.any(mask):
                    # Find indices where mask changes value
                    diff = np.diff(np.concatenate(([0], mask.astype(int), [0])))
                    starts = np.where(diff == 1)[0]
                    ends = np.where(diff == -1)[0]
                    consecutive_scans = (ends - starts).max()
                else:
                    consecutive_scans = 0
                if consecutive_scans < self.parameters.lc_ms.consecutive_scan_min:
                    mf_dict.pop(idx)
                    continue
                # Add start and final scan to mass_features and EICData
                left_scan, right_scan = (
                    myEIC.scans[l_a_r_scan_idx[0][0]],
                    myEIC.scans[l_a_r_scan_idx[0][2]],
                )
                mf_scan_apex = [(left_scan, int(apex_scan), right_scan)]
                myEIC.apexes = myEIC.apexes + mf_scan_apex
                mf_dict[idx].start_scan = left_scan
                mf_dict[idx].final_scan = right_scan

                # Find area under peak using limits from EIC centroid detector, add to mass_features and EICData
                area = np.trapz(
                    myEIC.eic_smoothed[l_a_r_scan_idx[0][0] : l_a_r_scan_idx[0][2] + 1],
                    myEIC.time[l_a_r_scan_idx[0][0] : l_a_r_scan_idx[0][2] + 1],
                )
                myEIC.areas = myEIC.areas + [area]
                self.eics[mz] = myEIC
                mf_dict[idx]._area = area
            else:
                if drop_if_fail is True:
                    mf_dict.pop(idx)

        if drop_duplicates:
            # Prepare mass feature dataframe
            if induced_features:
                mf_df = self.mass_features_to_df(induced_features = True).copy()
                mf_df = mf_df[mf_df.start_scan.notna()]        
            else:
                mf_df = self.mass_features_to_df(induced_features = False).copy()

            # For each mass feature, find all mass features within the clustering tolerance ppm and drop if their start and end times are within another mass feature
            # Keep the first mass feature (highest persistence)
            for idx, mass_feature in mf_df.iterrows():
                mz = mass_feature.mz
                apex_scan = mass_feature.apex_scan

                mf_df["mz_diff_ppm"] = np.abs(mf_df["mz"] - mz) / mz * 10**6
                mf_df_sub = mf_df[
                    mf_df["mz_diff_ppm"]
                    < self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel
                    * 10**6
                ].copy()

                # For all mass features within the clustering tolerance, check if the start and end times are within the start and end times of the mass feature
                for idx2, mass_feature2 in mf_df_sub.iterrows():
                    if idx2 != idx:
                        if (
                            mass_feature2.start_scan >= mass_feature.start_scan
                            and mass_feature2.final_scan <= mass_feature.final_scan
                        ):
                            if idx2 in self.mass_features.keys():
                                self.mass_features.pop(idx2)
        
        # Filter MS2 scans to only include those within integration bounds
        # This ensures MS2 scans outside start_scan to final_scan are removed
        if induced_features:
            self._filter_ms2_scans_by_integration_bounds(mf_dict=self.induced_mass_features)
        else:
            self._filter_ms2_scans_by_integration_bounds(mf_dict=self.mass_features)

    def find_c13_mass_features(self):
        """Mark likely C13 isotopes and connect to monoisoitopic mass features.

        Returns
        -------
        None, but populates the monoisotopic_mf_id and isotopologue_type attributes to the indivual LCMSMassFeatures within the mass_features attribute of the LCMSBase object.

        Raises
        ------
        ValueError
            If no mass features are found.
        """
        verbose = self.parameters.lc_ms.verbose_processing
        if verbose:
            print("evaluating mass features for C13 isotopes")
        if self.mass_features is None:
            raise ValueError("No mass features found, run find_mass_features() first")

        # Data prep fo sparse distance matrix
        dims = ["mz", "scan_time"]
        mf_df = self.mass_features_to_df().copy()
        # Drop mass features that have no area (these are likely to be noise)
        mf_df = mf_df[mf_df["area"].notnull()]
        mf_df["mf_id"] = mf_df.index.values
        dims = ["mz", "scan_time"]

        # Sort my ascending mz so we always get the monoisotopic mass first, regardless of the order/intensity of the mass features
        mf_df = mf_df.sort_values(by=["mz"]).reset_index(drop=True).copy()

        mz_diff = 1.003355  # C13-C12 mass difference
        tol = [
            mf_df["mz"].median()
            * self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
            self.parameters.lc_ms.mass_feature_cluster_rt_tolerance * 0.5,
        ]  # mz, in relative; scan_time in minutes

        # Compute inter-feature distances
        distances = None
        for i in range(len(dims)):
            # Construct k-d tree
            values = mf_df[dims[i]].values
            tree = KDTree(values.reshape(-1, 1))

            max_tol = tol[i]
            if dims[i] == "mz":
                # Maximum absolute tolerance
                max_tol = mz_diff + tol[i]

            # Compute sparse distance matrix
            # the larger the max_tol, the slower this operation is
            sdm = tree.sparse_distance_matrix(tree, max_tol, output_type="coo_matrix")

            # Only consider forward case, exclude diagonal
            sdm = sparse.triu(sdm, k=1)

            if dims[i] == "mz":
                min_tol = mz_diff - tol[i]
                # Get only the ones that are above the min tol
                idx = sdm.data > min_tol

                # Reconstruct sparse distance matrix
                sdm = sparse.coo_matrix(
                    (sdm.data[idx], (sdm.row[idx], sdm.col[idx])),
                    shape=(len(values), len(values)),
                )

            # Cast as binary matrix
            sdm.data = np.ones_like(sdm.data)

            # Stack distances
            if distances is None:
                distances = sdm
            else:
                distances = distances.multiply(sdm)

        # Extract indices of within-tolerance points
        distances = distances.tocoo()
        pairs = np.stack((distances.row, distances.col), axis=1)  # C12 to C13 pairs

        # Turn pairs (which are index of mf_df) into mf_id and then into two dataframes to join to mf_df
        pairs_mf = pairs.copy()
        pairs_mf[:, 0] = mf_df.iloc[pairs[:, 0]].mf_id.values
        pairs_mf[:, 1] = mf_df.iloc[pairs[:, 1]].mf_id.values

        # Connect monoisotopic masses with isotopologes within mass_features
        monos = np.setdiff1d(np.unique(pairs_mf[:, 0]), np.unique(pairs_mf[:, 1]))
        for mono in monos:
            self.mass_features[mono].monoisotopic_mf_id = mono
        pairs_iso_df = pd.DataFrame(pairs_mf, columns=["parent", "child"])
        while not pairs_iso_df.empty:
            pairs_iso_df = pairs_iso_df.set_index("parent", drop=False)
            m1_isos = pairs_iso_df.loc[monos, "child"].unique()
            for iso in m1_isos:
                # Set monoisotopic_mf_id and isotopologue_type for isotopologues
                parent = pairs_mf[pairs_mf[:, 1] == iso, 0]
                if len(parent) > 1:
                    # Choose the parent that is closest in time to the isotopologue
                    parent_time = [self.mass_features[p].retention_time for p in parent]
                    time_diff = [
                        np.abs(self.mass_features[iso].retention_time - x)
                        for x in parent_time
                    ]
                    parent = parent[np.argmin(time_diff)]
                else:
                    parent = parent[0]
                self.mass_features[iso].monoisotopic_mf_id = self.mass_features[
                    parent
                ].monoisotopic_mf_id
                if self.mass_features[iso].monoisotopic_mf_id is not None:
                    mass_diff = (
                        self.mass_features[iso].mz
                        - self.mass_features[
                            self.mass_features[iso].monoisotopic_mf_id
                        ].mz
                    )
                    self.mass_features[iso].isotopologue_type = "13C" + str(
                        int(round(mass_diff, 0))
                    )

            # Drop the mono and iso from the pairs_iso_df
            pairs_iso_df = pairs_iso_df.drop(
                index=monos, errors="ignore"
            )  # Drop pairs where the parent is a child that is a child of a root
            pairs_iso_df = pairs_iso_df.set_index("child", drop=False)
            pairs_iso_df = pairs_iso_df.drop(index=m1_isos, errors="ignore")

            if not pairs_iso_df.empty:
                # Get new monos, recognizing that these are just 13C isotopologues that are connected to other 13C isotopologues to repeat the process
                monos = np.setdiff1d(
                    np.unique(pairs_iso_df.parent), np.unique(pairs_iso_df.child)
                )
        if verbose:
            # Report fraction of compounds annotated with isotopes
            mf_df["c13_flag"] = np.where(
                np.logical_or(
                    np.isin(mf_df["mf_id"], pairs_mf[:, 0]),
                    np.isin(mf_df["mf_id"], pairs_mf[:, 1]),
                ),
                1,
                0,
            )
            print(
                str(round(len(mf_df[mf_df["c13_flag"] == 1]) / len(mf_df), ndigits=3))
                + " of mass features have or are C13 isotopes"
            )

    def deconvolute_ms1_mass_features(self):
        """Deconvolute MS1 mass features

        Deconvolute mass features ms1 spectrum based on the correlation of all masses within a spectrum over the EIC of the mass features

        Parameters
        ----------
        None

        Returns
        -------
        None, but assigns the _ms_deconvoluted_idx, mass_spectrum_deconvoluted_parent,
        and associated_mass_features_deconvoluted attributes to the mass features in the
        mass_features attribute of the LCMSBase object.

        Raises
        ------
        ValueError
            If no mass features are found, must run find_mass_features() first.
            If no EICs are found, did you run integrate_mass_features() first?

        """
        # Checks for set mass_features and eics
        if self.mass_features is None:
            raise ValueError(
                "No mass features found, did you run find_mass_features() first?"
            )

        if self.eics == {}:
            raise ValueError(
                "No EICs found, did you run integrate_mass_features() first?"
            )

        if 1 not in self._ms_unprocessed.keys():
            raise ValueError("No unprocessed MS1 spectra found.")

        # Prep ms1 data
        ms1_data = self._ms_unprocessed[1].copy()
        ms1_data = ms1_data.set_index("scan")

        # Prep mass feature summary
        mass_feature_df = self.mass_features_to_df()

        # Loop through each mass feature
        for mf_id, mass_feature in self.mass_features.items():
            # Check that the mass_feature.mz attribute == the mz of the mass feature in the mass_feature_df
            if mass_feature.mz != mass_feature.ms1_peak.mz_exp:
                continue

            # Get the left and right limits of the EIC of the mass feature
            l_scan, _, r_scan = mass_feature._eic_data.apexes[0]

            # Pull from the _ms1_unprocessed data the scan range of interest and sort by mz
            ms1_data_sub = ms1_data.loc[l_scan:r_scan].copy()
            ms1_data_sub = ms1_data_sub.sort_values(by=["mz"]).reset_index(drop=False)

            # Get the centroided masses of the mass feature
            mf_mspeak_mzs = mass_feature.mass_spectrum.mz_exp

            # Find the closest mz in the ms1 data to the centroided masses of the mass feature
            ms1_data_sub["mass_feature_mz"] = mf_mspeak_mzs[
                find_closest(mf_mspeak_mzs, ms1_data_sub.mz.values)
            ]

            # Drop rows with mz_diff > 0.01 between the mass feature mz and the ms1 data mz
            ms1_data_sub["mz_diff_rel"] = (
                np.abs(ms1_data_sub["mass_feature_mz"] - ms1_data_sub["mz"])
                / ms1_data_sub["mz"]
            )
            ms1_data_sub = ms1_data_sub[
                ms1_data_sub["mz_diff_rel"]
                < self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel
            ].reset_index(drop=True)

            # Group by mass_feature_mz and scan and sum intensity
            ms1_data_sub_group = (
                ms1_data_sub.groupby(["mass_feature_mz", "scan"])["intensity"]
                .sum()
                .reset_index()
            )

            # Calculate the correlation of the intensities of the mass feature and the ms1 data (set to 0 if no intensity)
            corr = (
                ms1_data_sub_group.pivot(
                    index="scan", columns="mass_feature_mz", values="intensity"
                )
                .fillna(0)
                .corr()
            )

            # Subset the correlation matrix to only include the masses of the mass feature and those with a correlation > 0.8
            decon_corr_min = self.parameters.lc_ms.ms1_deconvolution_corr_min

            # Try catch for KeyError in case the mass feature mz is not in the correlation matrix
            try:
                corr_subset = corr.loc[mass_feature.mz,]
            except KeyError:
                # If the mass feature mz is not in the correlation matrix, skip to the next mass feature
                continue

            corr_subset = corr_subset[corr_subset > decon_corr_min]

            # Get the masses from the mass spectrum that are the result of the deconvolution
            mzs_decon = corr_subset.index.values

            # Get the indices of the mzs_decon in mass_feature.mass_spectrum.mz_exp and assign to the mass feature
            mzs_decon_idx = [
                id
                for id, mz in enumerate(mass_feature.mass_spectrum.mz_exp)
                if mz in mzs_decon
            ]
            mass_feature._ms_deconvoluted_idx = mzs_decon_idx

            # Check if the mass feature's ms1 peak is the largest in the deconvoluted mass spectrum
            if (
                mass_feature.ms1_peak.abundance
                == mass_feature.mass_spectrum.abundance[mzs_decon_idx].max()
            ):
                mass_feature.mass_spectrum_deconvoluted_parent = True
            else:
                mass_feature.mass_spectrum_deconvoluted_parent = False

            # Check for other mass features that are in the deconvoluted mass spectrum and add the deconvoluted mass spectrum to the mass feature
            # Subset mass_feature_df to only include mass features that are within the clustering tolerance
            mass_feature_df_sub = mass_feature_df[
                abs(mass_feature.retention_time - mass_feature_df["scan_time"])
                < self.parameters.lc_ms.mass_feature_cluster_rt_tolerance
            ].copy()
            # Calculate the mz difference in ppm between the mass feature and the peaks in the deconvoluted mass spectrum
            mass_feature_df_sub["mz_diff_ppm"] = [
                np.abs(mzs_decon - mz).min() / mz * 10**6
                for mz in mass_feature_df_sub["mz"]
            ]
            # Subset mass_feature_df to only include mass features that are within 1 ppm of the deconvoluted masses
            mfs_associated_decon = mass_feature_df_sub[
                mass_feature_df_sub["mz_diff_ppm"]
                < self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel * 10**6
            ].index.values

            mass_feature.associated_mass_features_deconvoluted = mfs_associated_decon

    def _remove_redundant_mass_features(
        self,
        ) -> None:
        """
        Identify and remove redundant mass features that are likely contaminants based on their m/z values and scan frequency. 
        Especially useful for HILIC data where signals do not return to baseline between peaks or for data with significant background noise.
        
        Contaminants are characterized by:
        1. Similar m/z values (within ppm_tolerance)
        2. High frequency across scan numbers (ubiquitous presence)
    
        Notes
        -----
        Depends on self.mass_features being populated, uses the parameters in self.parameters.lc_ms for tolerances (mass_feature_cluster_mz_tolerance_rel)
        """
        ppm_tolerance = self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel*1e6
        min_scan_frequency = self.parameters.lc_ms.redundant_scan_frequency_min
        n_retain = self.parameters.lc_ms.redundant_feature_retain_n

        df = self.mass_features_to_df()

        if df.empty:
            return pd.DataFrame()
        # df index should be mf_id
        if 'mf_id' not in df.columns:
            if 'mf_id' in df.index.names:
                df = df.reset_index()
            else:
                raise ValueError("DataFrame must contain 'mf_id' column or index.")
        
        # Sort by m/z for efficient grouping
        df_sorted = df.sort_values('mz').reset_index(drop=True)
        
        # Calculate total number of unique scans for frequency calculation
        # Calculate total possible scans (check the cluster rt tolerance and the min rt and max rt of the data)
        total_time = self.scan_df['scan_time'].max() - self.scan_df['scan_time'].min()
        cluster_rt_tolerance = self.parameters.lc_ms.mass_feature_cluster_rt_tolerance
        # If the feature was detected in every possible scan (and then rolled up), it would be in this many scans
        total_scans =  int(total_time / cluster_rt_tolerance) + 1

        # Group similar m/z values using ppm tolerance
        mz_groups = []
        current_group = []
        
        for i, row in df_sorted.iterrows():
            current_mz = row['mz']
            
            if not current_group:
                # Start first group
                current_group = [i]
            else:
                # Check if current m/z is within tolerance of group representative
                group_representative_mz = df_sorted.iloc[current_group[0]]['mz']
                ppm_diff = abs(current_mz - group_representative_mz) / group_representative_mz * 1e6
                
                if ppm_diff <= ppm_tolerance:
                    # Add to current group
                    current_group.append(i)
                else:
                    # Start new group, but first process current group
                    if len(current_group) > 0:
                        mz_groups.append(current_group)
                    current_group = [i]
        
        # Don't forget the last group
        if current_group:
            mz_groups.append(current_group)
        
        # Analyze each m/z group for contaminant characteristics
        
        for group_indices in mz_groups:
            group_data = df_sorted.iloc[group_indices]
            
            # Calculate group statistics
            unique_scans = group_data['apex_scan'].nunique()
            scan_frequency = unique_scans / total_scans
            
            # Check if this group meets contaminant criteria
            if scan_frequency >= min_scan_frequency:
                group_data = group_data.sort_values('intensity', ascending=False)
                non_representative_mf_id = group_data.iloc[n_retain:]['mf_id'].tolist()  # These will be removed

                self.mass_features = {
                    k: v for k, v in self.mass_features.items() if k not in non_representative_mf_id
                }

    def _remove_mass_features_by_peak_metrics(self, induced_features=False) -> None:
        """Remove mass features based on peak metrics defined in mass_feature_attribute_filter_dict.
        
        This method filters mass features based on various peak shape metrics and quality indicators
        such as noise scores, Gaussian similarity, tailing factors, dispersity index, etc.
        
        The filtering criteria are defined in the mass_feature_attribute_filter_dict parameter,
        which should contain attribute names as keys and filter specifications as values.
        
        Filter specification format:
        {attribute_name: {'value': threshold, 'operator': comparison}}
        
        Available operators:
        - '>' or 'greater': Keep features where attribute > threshold
        - '<' or 'less': Keep features where attribute < threshold  
        - '>=' or 'greater_equal': Keep features where attribute >= threshold
        - '<=' or 'less_equal': Keep features where attribute <= threshold
        
        Examples:
        - {'noise_score_max': {'value': 0.5, 'operator': '>='}} - Keep features with noise_score_max >= 0.5
        - {'dispersity_index': {'value': 0.1, 'operator': '<'}} - Keep features with dispersity_index < 0.1
        - {'gaussian_similarity': {'value': 0.7, 'operator': '>='}} - Keep features with gaussian_similarity >= 0.7
        
        Parameters
        ----------
        induced_features : bool, optional
            If True, filter induced_mass_features instead of regular mass_features. Default is False.
        
        Returns
        -------
        None
            Modifies self.mass_features or self.induced_mass_features in place by removing filtered features.
            
        Raises
        ------
        ValueError
            If no mass features are found, if an invalid attribute is specified, or if filter specification is malformed.
        """
        # Select the appropriate mass features dictionary
        if induced_features:
            mf_dict = self.induced_mass_features
            mf_type = "induced mass features"
        else:
            mf_dict = self.mass_features
            mf_type = "mass features"
            
        if mf_dict is None or len(mf_dict) == 0:
            raise ValueError(f"No {mf_type} found, run {'gap filling' if induced_features else 'find_mass_features()'} first")
            
        filter_dict = self.parameters.lc_ms.mass_feature_attribute_filter_dict
        
        if not filter_dict:
            # No filtering criteria specified, return early
            return
            
        verbose = self.parameters.lc_ms.verbose_processing
        initial_count = len(mf_dict)
        
        if verbose:
            print(f"Filtering {mf_type} using peak metrics. Initial count: {initial_count}")
            
        # List to collect IDs of mass features to remove
        features_to_remove = []
        
        for mf_id, mass_feature in mf_dict.items():
            should_remove = False
            
            for attribute_name, filter_spec in filter_dict.items():
                # Validate filter specification structure
                if not isinstance(filter_spec, dict):
                    raise ValueError(f"Filter specification for '{attribute_name}' must be a dictionary with 'value' and 'operator' keys")
                
                if 'value' not in filter_spec or 'operator' not in filter_spec:
                    raise ValueError(f"Filter specification for '{attribute_name}' must contain both 'value' and 'operator' keys")
                
                threshold_value = filter_spec['value']
                operator = filter_spec['operator'].lower().strip()
                
                # Validate operator
                valid_operators = {'>', '<', '>=', '<=', 'greater', 'less', 'greater_equal', 'less_equal'}
                if operator not in valid_operators:
                    raise ValueError(f"Invalid operator '{operator}' for attribute '{attribute_name}'. Valid operators: {valid_operators}")
                
                # Normalize operator names
                operator_map = {
                    'greater': '>',
                    'less': '<', 
                    'greater_equal': '>=',
                    'less_equal': '<='
                }
                operator = operator_map.get(operator, operator)
                
                # Get the attribute value from the mass feature
                try:
                    if hasattr(mass_feature, attribute_name):
                        attribute_value = getattr(mass_feature, attribute_name)
                    else:
                        raise ValueError(f"Mass feature does not have attribute '{attribute_name}'")
                        
                    # Handle None values or attributes that haven't been calculated
                    if attribute_value is None:
                        if verbose:
                            print(f"Warning: Mass feature {mf_id} has None value for '{attribute_name}'. Removing feature.")
                        should_remove = True
                        break
                        
                    # Handle numpy arrays (like half_height_width which returns mean)
                    if hasattr(attribute_value, '__len__') and not isinstance(attribute_value, str):
                        # For arrays, we use the mean or appropriate summary statistic
                        if attribute_name == 'half_height_width':
                            # half_height_width property already returns the mean
                            pass
                        else:
                            attribute_value = float(np.mean(attribute_value))
                    
                    # Handle NaN values
                    if np.isnan(float(attribute_value)):
                        if verbose:
                            print(f"Warning: Mass feature {mf_id} has NaN value for '{attribute_name}'. Removing feature.")
                        should_remove = True
                        break
                    
                    # Apply the threshold comparison based on operator
                    attribute_value = float(attribute_value)
                    threshold_value = float(threshold_value)
                    
                    if operator == '>' and not (attribute_value > threshold_value):
                        should_remove = True
                        break
                    elif operator == '<' and not (attribute_value < threshold_value):
                        should_remove = True
                        break
                    elif operator == '>=' and not (attribute_value >= threshold_value):
                        should_remove = True
                        break  
                    elif operator == '<=' and not (attribute_value <= threshold_value):
                        should_remove = True
                        break
                        
                except (AttributeError, ValueError, TypeError) as e:
                    if verbose:
                        print(f"Error evaluating filter '{attribute_name}' for mass feature {mf_id}: {e}")
                    should_remove = True
                    break
            
            if should_remove:
                features_to_remove.append(mf_id)
        
        # Remove filtered mass features
        for mf_id in features_to_remove:
            del mf_dict[mf_id]
        
        if verbose and len(features_to_remove) > 0:
            print(f"Removed {len(features_to_remove)} {mf_type} based on peak metrics. Remaining: {len(mf_dict)}")
        
        # Update the appropriate dictionary
        if induced_features:
            self.induced_mass_features = mf_dict
        else:
            self.mass_features = mf_dict
        
        # Clean up unassociated EICs and ms1 data (only for regular features)
        self._remove_unassociated_eics()
        self._remove_unassociated_ms1_spectra()
            
    def _remove_unassociated_eics(self) -> None:
        """Remove EICs that are not associated with any mass features.

        This method cleans up the eics attribute by removing any EICs that do not correspond to
        any mass features currently stored in the mass_features attribute. This is useful for
        freeing up memory and ensuring that only relevant EICs are retained.

        Returns
        -------
        None
            Modifies self.eics in place by removing unassociated EICs.
        """
        if self.mass_features is None or len(self.mass_features) == 0:
            self.eics = {}
            return

        # Get the set of m/z values associated with current mass features
        associated_mzs = {mf.mz for mf in self.mass_features.values()}

        # Remove EICs that are not associated with any mass features
        self.eics = {mz: eic for mz, eic in self.eics.items() if mz in associated_mzs}
    
    def _remove_unassociated_ms1_spectra(self) -> None:
        """Remove MS1 spectra that are not associated with any mass features.
        This method cleans up the _ms_unprocessed attribute by removing any MS1 spectra that do not correspond to
        any mass features currently stored in the mass_features attribute. This is useful for freeing up memory
        and ensuring that only relevant MS1 spectra are retained.

        Returns
        -------
        None
        """
        if self.mass_features is None or len(self.mass_features) == 0:
            self._ms_unprocessed = {}
            return

        # Get the set of m/z values associated with current mass features
        associated_ms1_scans = {mf.apex_scan for mf in self.mass_features.values()}
        associated_ms1_scans = [int(scan) for scan in associated_ms1_scans]
        
        # Get keys within the _ms attribute (these are individual MassSpectrum objects)
        current_stored_spectra = list(set(self._ms.keys()))
        if len(current_stored_spectra) == 0:
            return
        current_stored_spectra = [int(scan) for scan in current_stored_spectra]

        # Filter the current_stored_spectra to only ms1 scans
        current_stored_spectra_ms1 = [ scan for scan in current_stored_spectra if scan in self.ms1_scans ]

        # Remove MS1 spectra that are not associated with any mass features
        scans_to_drop = [scan for scan in current_stored_spectra_ms1 if scan not in associated_ms1_scans]
        for scan in scans_to_drop:
            if scan in self._ms:
                del self._ms[scan]

class PHCalculations:
    """Methods for performing calculations related to 2D peak picking via persistent homology on LCMS data.

    Notes
    -----
    This class is intended to be used as a mixin for the LCMSBase class.

    Methods
    -------
    * sparse_mean_filter(idx, V, radius=[0, 1, 1]).
        Sparse implementation of a mean filter.
    * embed_unique_indices(a).
        Creates an array of indices, sorted by unique element.
    * sparse_upper_star(idx, V).
        Sparse implementation of an upper star filtration.
    * check_if_grid(data).
        Check if the data is gridded in mz space.
    * grid_data(data).
        Grid the data in the mz dimension.
    * find_mass_features_ph(ms_level=1, grid=True).
        Find mass features within an LCMSBase object using persistent homology.
    * cluster_mass_features(drop_children=True).
        Cluster regions of interest.
    """

    @staticmethod
    def sparse_mean_filter(idx, V, radius=[0, 1, 1]):
        """Sparse implementation of a mean filter.

        Parameters
        ----------
        idx : :obj:`~numpy.array`
            Edge indices for each dimension (MxN).
        V : :obj:`~numpy.array`
            Array of intensity data (Mx1).
        radius : float or list
            Radius of the sparse filter in each dimension. Values less than
            zero indicate no connectivity in that dimension.

        Returns
        -------
        :obj:`~numpy.array`
            Filtered intensities (Mx1).

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package: https://github.com/pnnl/deimos.
        This is a static method.
        """

        # Copy indices
        idx = idx.copy().astype(V.dtype)

        # Scale
        for i, r in enumerate(radius):
            # Increase inter-index distance
            if r < 1:
                idx[:, i] *= 2

            # Do nothing
            elif r == 1:
                pass

            # Decrease inter-index distance
            else:
                idx[:, i] /= r

        # Connectivity matrix
        cmat = KDTree(idx)
        cmat = cmat.sparse_distance_matrix(cmat, 1, p=np.inf, output_type="coo_matrix")
        cmat.setdiag(1)

        # Pair indices
        I, J = cmat.nonzero()

        # Delete cmat
        cmat_shape = cmat.shape
        del cmat

        # Sum over columns
        V_sum = sparse.bsr_matrix(
            (V[J], (I, I)), shape=cmat_shape, dtype=V.dtype
        ).diagonal(0)

        # Count over columns
        V_count = sparse.bsr_matrix(
            (np.ones_like(J), (I, I)), shape=cmat_shape, dtype=V.dtype
        ).diagonal(0)

        return V_sum / V_count

    @staticmethod
    def embed_unique_indices(a):
        """Creates an array of indices, sorted by unique element.

        Parameters
        ----------
        a : :obj:`~numpy.array`
            Array of unique elements (Mx1).

        Returns
        -------
        :obj:`~numpy.array`
            Array of indices (Mx1).

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package: https://github.com/pnnl/deimos
        This is a static method.
        """

        def count_tens(n):
            # Count tens
            ntens = (n - 1) // 10

            while True:
                ntens_test = (ntens + n - 1) // 10

                if ntens_test == ntens:
                    return ntens
                else:
                    ntens = ntens_test

        def arange_exclude_10s(n):
            # How many 10s will there be?
            ntens = count_tens(n)

            # Base array
            arr = np.arange(0, n + ntens)

            # Exclude 10s
            arr = arr[(arr == 0) | (arr % 10 != 0)][:n]

            return arr

        # Creates an array of indices, sorted by unique element
        idx_sort = np.argsort(a)
        idx_unsort = np.argsort(idx_sort)

        # Sorts records array so all unique elements are together
        sorted_a = a[idx_sort]

        # Returns the unique values, the index of the first occurrence,
        # and the count for each element
        vals, idx_start, count = np.unique(
            sorted_a, return_index=True, return_counts=True
        )

        # Splits the indices into separate arrays
        splits = np.split(idx_sort, idx_start[1:])

        # Creates unique indices for each split
        idx_unq = np.concatenate([arange_exclude_10s(len(x)) for x in splits])

        # Reorders according to input array
        idx_unq = idx_unq[idx_unsort]

        # Magnitude of each index
        exp = np.log10(
            idx_unq, where=idx_unq > 0, out=np.zeros_like(idx_unq, dtype=np.float64)
        )
        idx_unq_mag = np.power(10, np.floor(exp) + 1)

        # Result
        return a + idx_unq / idx_unq_mag

    @staticmethod
    def roll_up_dataframe(
        df: pd.DataFrame,
        sort_by: str,
        tol: list,
        relative: list,
        dims: list,
        memory_opt_threshold: int = 10000,
    ):
        """Subset data by rolling up into apex in appropriate dimensions.

        Parameters
        ----------
        data : pd.DataFrame
            The input data containing "dims" columns and the "sort_by" column.
        sort_by : str
            The column to sort the data by, this will determine which mass features get rolled up into a parent mass feature
            (i.e., the mass feature with the highest value in the sort_by column).
        dims : list
            A list of dimension names (column names in the data DataFrame) to roll up the mass features by.
        tol : list
            A list of tolerances for each dimension. The length of the list must match the number of dimensions.
            The tolerances can be relative (as a fraction of the maximum value in that dimension) or absolute (in the units of that dimension).
            If relative is True, the tolerance will be multiplied by the maximum value in that dimension.
            If relative is False, the tolerance will be used as is.
        relative : list
            A list of booleans indicating whether the tolerance for each dimension is relative (True) or absolute (False).
        memory_opt_threshold : int, optional
            Minimum number of rows to trigger memory-optimized processing. Default is 10000.

        Returns
        -------
        pd.DataFrame
            A DataFrame with only the rolled up mass features, with the original index and columns.


        Raises
        ------
        ValueError
            If the input data is not a pandas DataFrame.
            If the input data does not have columns for each of the dimensions in "dims".
            If the length of "dims", "tol", and "relative" do not match.
        """
        og_columns = df.columns.copy()

        # Unindex the data, but keep the original index
        if df.index.name is not None:
            og_index = df.index.name
        else:
            og_index = "index"
        df = df.reset_index(drop=False)

        # Sort data by sort_by column, and reindex
        df = df.sort_values(by=sort_by, ascending=False).reset_index(drop=True)

        # Check that data is a DataFrame and has columns for each of the dims
        if not isinstance(df, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame")
        for dim in dims:
            if dim not in df.columns:
                raise ValueError(f"Data must have a column for {dim}")
        if len(dims) != len(tol) or len(dims) != len(relative):
            raise ValueError(
                "Dimensions, tolerances, and relative flags must be the same length"
            )

        # Pre-compute all values arrays
        all_values = [df[dim].values for dim in dims]

        # Choose processing method based on dataframe size
        if len(df) >= memory_opt_threshold:
            # Memory-optimized approach for large dataframes
            distances = PHCalculations._compute_distances_memory_optimized(
                all_values, tol, relative
            )
        else:
            # Faster approach for smaller dataframes
            distances = PHCalculations._compute_distances_original(
                all_values, tol, relative
            )

        # Process pairs with original logic but memory optimizations
        distances = distances.tocoo()
        pairs = np.stack((distances.row, distances.col), axis=1)
        pairs_df = pd.DataFrame(pairs, columns=["parent", "child"]).set_index("parent")
        del distances, pairs  # Free memory immediately

        to_drop = []
        while not pairs_df.empty:
            # Find root_parents and their children (original logic preserved)
            root_parents = np.setdiff1d(
                np.unique(pairs_df.index.values), np.unique(pairs_df.child.values)
            )
            children_of_roots = pairs_df.loc[root_parents, "child"].unique()
            to_drop.extend(children_of_roots)  # Use extend instead of append

            # Remove root_children as possible parents from pairs_df for next iteration
            pairs_df = pairs_df.drop(index=children_of_roots, errors="ignore")
            pairs_df = pairs_df.reset_index().set_index("child")
            # Remove root_children as possible children from pairs_df for next iteration
            pairs_df = pairs_df.drop(index=children_of_roots)

            # Prepare for next iteration
            pairs_df = pairs_df.reset_index().set_index("parent")

        # Convert to numpy array for efficient dropping
        to_drop = np.array(to_drop)

        # Drop mass features that are not cluster parents
        df_sub = df.drop(index=to_drop)

        # Set index back to og_index and only keep original columns
        df_sub = df_sub.set_index(og_index).sort_index()[og_columns]

        return df_sub

    @staticmethod
    def _compute_distances_original(all_values, tol, relative):
        """Original distance computation method for smaller datasets.

        This method computes the pairwise distances between features in the dataset
        using a straightforward approach. It is suitable for smaller datasets where
        memory usage is not a primary concern.

        Parameters
        ----------
        all_values : list of :obj:`~numpy.array`
            List of arrays containing the values for each dimension.
        tol : list of float
            List of tolerances for each dimension.
        relative : list of bool
            List of booleans indicating whether the tolerance for each dimension is relative (True) or absolute (False).

        Returns
        -------
        :obj:`~scipy.sparse.coo_matrix`
            Sparse matrix indicating pairwise distances within tolerances.
        """
        # Compute inter-feature distances with memory optimization
        distances = None
        for i in range(len(all_values)):
            values = all_values[i]
            # Use single precision if possible to reduce memory
            tree = KDTree(values.reshape(-1, 1).astype(np.float32))

            max_tol = tol[i]
            if relative[i] is True:
                max_tol = tol[i] * values.max()

            # Compute sparse distance matrix with smaller chunks if memory is an issue
            sdm = tree.sparse_distance_matrix(tree, max_tol, output_type="coo_matrix")

            # Only consider forward case, exclude diagonal
            sdm = sparse.triu(sdm, k=1)

            # Process relative distances more efficiently
            if relative[i] is True:
                # Vectorized computation without creating intermediate arrays
                row_values = values[sdm.row]
                valid_idx = sdm.data <= tol[i] * row_values

                # Reconstruct sparse matrix more efficiently
                sdm = sparse.coo_matrix(
                    (
                        np.ones(valid_idx.sum(), dtype=np.uint8),
                        (sdm.row[valid_idx], sdm.col[valid_idx]),
                    ),
                    shape=(len(values), len(values)),
                )
            else:
                # Cast as binary matrix with smaller data type
                sdm.data = np.ones(len(sdm.data), dtype=np.uint8)

            # Stack distances with memory-efficient multiplication
            if distances is None:
                distances = sdm
            else:
                # Use in-place operations where possible
                distances = distances.multiply(sdm)
                del sdm  # Free memory immediately

        return distances

    @staticmethod
    def _compute_distances_memory_optimized(all_values, tol, relative):
        """Memory-optimized distance computation for large datasets.

        This method computes the pairwise distances between features in the dataset
        using a more memory-efficient approach. It is suitable for larger datasets
        where memory usage is a primary concern.

        Parameters
        ----------
        all_values : list of :obj:`~numpy.array`
            List of arrays containing the values for each dimension.
        tol : list of float
            List of tolerances for each dimension.
        relative : list of bool
            List of booleans indicating whether the tolerance for each dimension is relative (True) or absolute (False).

        Returns
        -------
        :obj:`~scipy.sparse.coo_matrix`
            Sparse matrix indicating pairwise distances within tolerances.
        """
        # Compute distance matrix for first dimension (full matrix as before)
        values_0 = all_values[0].astype(np.float32)
        tree_0 = KDTree(values_0.reshape(-1, 1))

        max_tol_0 = tol[0]
        if relative[0] is True:
            max_tol_0 = tol[0] * values_0.max()

        # Compute sparse distance matrix for first dimension
        distances = tree_0.sparse_distance_matrix(
            tree_0, max_tol_0, output_type="coo_matrix"
        )
        distances = sparse.triu(distances, k=1)

        # Process relative distances for first dimension
        if relative[0] is True:
            row_values = values_0[distances.row]
            valid_idx = distances.data <= tol[0] * row_values
            distances = sparse.coo_matrix(
                (
                    np.ones(valid_idx.sum(), dtype=np.uint8),
                    (distances.row[valid_idx], distances.col[valid_idx]),
                ),
                shape=(len(values_0), len(values_0)),
            )
        else:
            distances.data = np.ones(len(distances.data), dtype=np.uint8)

        # For remaining dimensions, work only on chunks defined by first dimension pairs
        if len(all_values) > 1:
            distances_coo = distances.tocoo()
            valid_pairs = []

            # Process each pair from first dimension
            for idx in range(len(distances_coo.data)):
                i, j = distances_coo.row[idx], distances_coo.col[idx]
                is_valid_pair = True

                # Check remaining dimensions for this specific pair
                for dim_idx in range(1, len(all_values)):
                    values = all_values[dim_idx]
                    val_i, val_j = values[i], values[j]

                    max_tol = tol[dim_idx]
                    if relative[dim_idx] is True:
                        max_tol = tol[dim_idx] * values.max()

                    distance_ij = abs(val_i - val_j)

                    # Check if this pair satisfies the tolerance for this dimension
                    if relative[dim_idx] is True:
                        if distance_ij > tol[dim_idx] * val_i:
                            is_valid_pair = False
                            break
                    else:
                        if distance_ij > max_tol:
                            is_valid_pair = False
                            break

                if is_valid_pair:
                    valid_pairs.append((i, j))

            # Rebuild distances matrix with only valid pairs
            if valid_pairs:
                valid_pairs = np.array(valid_pairs)
                distances = sparse.coo_matrix(
                    (
                        np.ones(len(valid_pairs), dtype=np.uint8),
                        (valid_pairs[:, 0], valid_pairs[:, 1]),
                    ),
                    shape=(len(values_0), len(values_0)),
                )
            else:
                # No valid pairs found
                distances = sparse.coo_matrix(
                    (len(values_0), len(values_0)), dtype=np.uint8
                )

        return distances

    def sparse_upper_star(self, idx, V):
        """Sparse implementation of an upper star filtration.

        Parameters
        ----------
        idx : :obj:`~numpy.array`
            Edge indices for each dimension (MxN).
        V : :obj:`~numpy.array`
            Array of intensity data (Mx1).
        Returns
        -------
        idx : :obj:`~numpy.array`
            Index of filtered points (Mx1).
        persistence : :obj:`~numpy.array`
            Persistence of each filtered point (Mx1).

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package: https://github.com/pnnl/deimos
        """

        # Invert
        V = -1 * V.copy().astype(int)

        # Embed indices
        V = self.embed_unique_indices(V)

        # Connectivity matrix
        cmat = KDTree(idx)
        cmat = cmat.sparse_distance_matrix(cmat, 1, p=np.inf, output_type="coo_matrix")
        cmat.setdiag(1)
        cmat = sparse.triu(cmat)

        # Pairwise minimums
        I, J = cmat.nonzero()
        d = np.maximum(V[I], V[J])

        # Delete connectiity matrix
        cmat_shape = cmat.shape
        del cmat

        # Sparse distance matrix
        sdm = sparse.coo_matrix((d, (I, J)), shape=cmat_shape)

        # Delete pairwise mins
        del d, I, J

        # Persistence homology
        ph = ripser(sdm, distance_matrix=True, maxdim=0)["dgms"][0]

        # Bound death values
        ph[ph[:, 1] == np.inf, 1] = np.max(V)

        # Construct tree to query against
        tree = KDTree(V.reshape((-1, 1)))

        # Get the indexes of the first nearest neighbor by birth
        _, nn = tree.query(ph[:, 0].reshape((-1, 1)), k=1, workers=-1)

        return nn, -(ph[:, 0] // 1 - ph[:, 1] // 1)

    def check_if_grid(self, data):
        """Check if the data are gridded in mz space.

        Parameters
        ----------
        data : DataFrame
            DataFrame containing the mass spectrometry data.  Needs to have mz and scan columns.

        Returns
        -------
        bool
            True if the data is gridded in the mz direction, False otherwise.

        Notes
        -----
        This function is used within the grid_data function and the find_mass_features function and is not intended to be called directly.
        """
        # Calculate the difference between consecutive mz values in a single scan
        dat_check = data.copy().reset_index(drop=True)
        dat_check["mz_diff"] = np.abs(dat_check["mz"].diff())
        mz_diff_min = (
            dat_check.groupby("scan")["mz_diff"].min().min()
        )  # within each scan, what is the smallest mz difference between consecutive mz values

        # Find the mininum mz difference between mz values in the data; regardless of scan
        dat_check_mz = dat_check[["mz"]].drop_duplicates().copy()
        dat_check_mz = dat_check_mz.sort_values(by=["mz"]).reset_index(drop=True)
        dat_check_mz["mz_diff"] = np.abs(dat_check_mz["mz"].diff())

        # Get minimum mz_diff between mz values in the data
        mz_diff_min_raw = dat_check_mz["mz_diff"].min()

        # If the minimum mz difference between mz values in the data is less than the minimum mz difference between mz values within a single scan, then the data is not gridded
        if mz_diff_min_raw < mz_diff_min:
            return False
        else:
            return True

    def grid_data(self, data, attempts=5):
        """Grid the data in the mz dimension.

        Data must be gridded prior to persistent homology calculations and computing average mass spectrum

        Parameters
        ----------
        data : DataFrame
            The input data containing mz, scan, scan_time, and intensity columns.
        attempts : int, optional
            The number of attempts to grid the data. Default is 5.

        Returns
        -------
        DataFrame
            The gridded data with mz, scan, scan_time, and intensity columns.

        Raises
        ------
        ValueError
            If gridding fails after the specified number of attempts.
        """
        attempt_i = 0
        while attempt_i < attempts:
            attempt_i += 1
            data = self._grid_data(data)

            if self.check_if_grid(data):
                return data

        if not self.check_if_grid(data):
            raise ValueError(
                "Gridding failed after "
                + str(attempt_i)
                + " attempts. Please check the data."
            )
        else:
            return data

    def _grid_data(self, data):
        """Internal method to grid the data in the mz dimension.

        Notes
        -----
        This method is called by the grid_data method and should not be called directly.
        It will attempt to grid the data in the mz dimension by creating a grid of mz values based on the minimum mz difference within each scan,
        but it does not check if the data is already gridded or if the gridding is successful.

        Parameters
        ----------
        data : pd.DataFrame or pl.DataFrame
            The input data to grid.

        Returns
        -------
        pd.DataFrame or pl.DataFrame
            The data after attempting to grid it in the mz dimension.
        """
        # Calculate the difference between consecutive mz values in a single scan for grid spacing
        data_w = data.copy().reset_index(drop=True)
        data_w["mz_diff"] = np.abs(data_w["mz"].diff())
        mz_diff_min = data_w.groupby("scan")["mz_diff"].min().min() * 0.99999

        # Need high intensity mz values first so they are parents in the output pairs stack
        dat_mz = data_w[["mz", "intensity"]].sort_values(
            by=["intensity"], ascending=False
        )
        dat_mz = dat_mz[["mz"]].drop_duplicates().reset_index(drop=True).copy()

        # Construct KD tree
        tree = KDTree(dat_mz.mz.values.reshape(-1, 1))
        sdm = tree.sparse_distance_matrix(tree, mz_diff_min, output_type="coo_matrix")
        sdm = sparse.triu(sdm, k=1)
        sdm.data = np.ones_like(sdm.data)
        distances = sdm.tocoo()
        pairs = np.stack((distances.row, distances.col), axis=1)

        # Cull pairs to just get root
        to_drop = []
        while len(pairs) > 0:
            root_parents = np.setdiff1d(np.unique(pairs[:, 0]), np.unique(pairs[:, 1]))
            id_root_parents = np.isin(pairs[:, 0], root_parents)
            children_of_roots = np.unique(pairs[id_root_parents, 1])
            to_drop = np.append(to_drop, children_of_roots)

            # Set up pairs array for next iteration by removing pairs with children or parents already dropped
            pairs = pairs[~np.isin(pairs[:, 1], to_drop), :]
            pairs = pairs[~np.isin(pairs[:, 0], to_drop), :]
        dat_mz = dat_mz.reset_index(drop=True).drop(index=np.array(to_drop))
        mz_dat_np = (
            dat_mz[["mz"]]
            .sort_values(by=["mz"])
            .reset_index(drop=True)
            .values.flatten()
        )

        # Sort data by mz and recast mz to nearest value in mz_dat_np
        data_w = data_w.sort_values(by=["mz"]).reset_index(drop=True).copy()
        data_w["mz_new"] = mz_dat_np[find_closest(mz_dat_np, data_w["mz"].values)]
        data_w["mz_diff"] = np.abs(data_w["mz"] - data_w["mz_new"])

        # Rename mz_new as mz; drop mz_diff; groupby scan and mz and sum intensity
        new_data_w = data_w.rename(columns={"mz": "mz_orig", "mz_new": "mz"}).copy()
        new_data_w = (
            new_data_w.drop(columns=["mz_diff", "mz_orig"])
            .groupby(["scan", "mz"])["intensity"]
            .sum()
            .reset_index()
        )
        new_data_w = (
            new_data_w.sort_values(by=["scan", "mz"], ascending=[True, True])
            .reset_index(drop=True)
            .copy()
        )

        return new_data_w

    def _filter_data_by_targets(self, data, target_search_dict):
        """Filter MS data to only include m/z and RT windows around target values.
        
        Parameters
        ----------
        data : pd.DataFrame
            MS data with 'mz' and 'scan_time' columns
        target_search_dict : dict
            Dictionary with target_mz_list, target_rt_list, mz_tolerance_ppm, rt_tolerance
            
        Returns
        -------
        pd.DataFrame
            Filtered data containing only points within target windows
        """
        target_mz_list = target_search_dict['target_mz_list']
        target_rt_list = target_search_dict['target_rt_list']
        mz_tolerance_ppm = target_search_dict['mz_tolerance_ppm']
        rt_tolerance = target_search_dict['rt_tolerance']
        
        # Create a mask for data points that fall within any target window
        mask = np.zeros(len(data), dtype=bool)
        
        for target_mz, target_rt in zip(target_mz_list, target_rt_list):
            # Calculate m/z window
            mz_tol = target_mz * mz_tolerance_ppm / 1e6
            mz_min = target_mz - mz_tol
            mz_max = target_mz + mz_tol
            
            # Calculate RT window
            rt_min = target_rt - rt_tolerance
            rt_max = target_rt + rt_tolerance
            
            # Create mask for this target
            target_mask = (
                (data['mz'] >= mz_min) & (data['mz'] <= mz_max) &
                (data['scan_time'] >= rt_min) & (data['scan_time'] <= rt_max)
            )
            
            # Combine with overall mask
            mask |= target_mask
        
        return data[mask].reset_index(drop=True)
    
    def find_mass_features_ph(self, ms_level=1, grid=True, targeted_search=False, target_search_dict=None, mf_type="untargeted", accumulate_features=False):
        """Find mass features within an LCMSBase object using persistent homology.

        Assigns the mass_features attribute to the object (a dictionary of LCMSMassFeature objects, keyed by mass feature id)

        Parameters
        ----------
        ms_level : int, optional
            The MS level to use. Default is 1.
        grid : bool, optional
            If True, will regrid the data before running the persistent homology calculations (after checking if the data is gridded). Default is True.
        targeted_search : bool, optional
            If True, perform targeted search mode. Default is False.
        target_search_dict : dict or None, optional
            Dictionary with target parameters for targeted search. Default is None.
        mf_type : str, optional
            Type label for the mass features. Default is "untargeted".
        accumulate_features : bool, optional
            If True, add to existing features rather than replacing them. Default is False.

        Raises
        ------
        ValueError
            If no MS level data is found on the object.
            If data is not gridded and grid is False.

        Returns
        -------
        None, but assigns the mass_features attribute to the object.

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package: https://github.com/pnnl/deimos
        """
        # Check that ms_level is a key in self._ms_uprocessed
        if ms_level not in self._ms_unprocessed.keys():
            raise ValueError(
                "No MS level "
                + str(ms_level)
                + " data found, did you instantiate with parser specific to MS level?"
            )

        # Get ms data
        data = self._ms_unprocessed[ms_level].copy()

        # Drop rows with missing intensity values and reset index
        data = data.dropna(subset=["intensity"]).reset_index(drop=True)
        
        # Add scan_time for filtering if in targeted mode
        if targeted_search:
            data = data.merge(self.scan_df[["scan", "scan_time"]], on="scan", how="left")

        # Threshold data (bypass thresholds in targeted mode)
        dims = ["mz", "scan_time"]
        if targeted_search:
            # In targeted mode, bypass intensity and persistence thresholds
            threshold = 0
            persistence_threshold = 0
            # Filter data to only target windows
            data_thres = self._filter_data_by_targets(data, target_search_dict)
            if len(data_thres) == 0:
                if self.parameters.lc_ms.verbose_processing:
                    print("No data found in target windows")
                self.mass_features = {}
                return
        else:
            threshold = self.parameters.lc_ms.ph_inten_min_rel * data.intensity.max()
            persistence_threshold = (
                self.parameters.lc_ms.ph_persis_min_rel * data.intensity.max()
            )
            data_thres = data[data["intensity"] > threshold].reset_index(drop=True).copy()

        # Check if gridded, if not, grid
        gridded_mz = self.check_if_grid(data_thres)
        if gridded_mz is False:
            if grid is False:
                raise ValueError(
                    "Data are not gridded in mz dimension, try reprocessing with a different params or grid data before running this function"
                )
            else:
                data_thres = self.grid_data(data_thres)

        # Add scan_time (skip if already present from targeted mode)
        if 'scan_time' not in data_thres.columns:
            data_thres = data_thres.merge(self.scan_df[["scan", "scan_time"]], on="scan")
        # Process in chunks if required
        if len(data_thres) > 10000:
            return self._find_mass_features_ph_partition(
                data_thres, dims, persistence_threshold, mf_type, accumulate_features
            )
        else:
            # Process all at once
            return self._find_mass_features_ph_single(
                data_thres, dims, persistence_threshold, mf_type, accumulate_features
            )
            return self._find_mass_features_ph_single(
                data_thres, dims, persistence_threshold, mf_type
            )

    def _find_mass_features_ph_single(self, data_thres, dims, persistence_threshold, mf_type="untargeted", accumulate_features=False):
        """Process all data at once (original logic)."""
        # Build factors
        factors = {
            dim: pd.factorize(data_thres[dim], sort=True)[1].astype(np.float32)
            for dim in dims
        }

        # Build indexes
        index = {
            dim: np.searchsorted(factors[dim], data_thres[dim]).astype(np.float32)
            for dim in factors
        }

        # Smooth and process
        mass_features_df = self._process_partition_ph(
            data_thres, index, dims, persistence_threshold
        )

        # Roll up within chunk to remove duplicates
        mass_features_df = self.roll_up_dataframe(
            df=mass_features_df,
            sort_by="persistence",
            dims=["mz", "scan_time"],
            tol=[
                self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
                self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
            ],
            relative=[True, False],
        )
        mass_features_df = mass_features_df.reset_index(drop=True)

        # Populate mass_features attribute
        self._populate_mass_features(mass_features_df, mf_type, accumulate_features)

    def _find_mass_features_ph_partition(self, data_thres, dims, persistence_threshold, mf_type="untargeted", accumulate_features=False):
        """Partition the persistent homology mass feature detection for large datasets.

        This method splits the input data into overlapping scan partitions, processes each partition to detect mass features
        using persistent homology, rolls up duplicates within and across partitions, and populates the mass_features attribute.

        Parameters
        ----------
        data_thres : pd.DataFrame
            The thresholded input data containing mass spectrometry information.
        dims : list
            List of dimension names (e.g., ["mz", "scan_time"]) used for feature detection.
        persistence_threshold : float
            Minimum persistence value required for a detected mass feature to be retained.
        mf_type : str, optional
            Type label for the mass features. Default is "untargeted".
        accumulate_features : bool, optional
            If True, add to existing features rather than replacing them. Default is False.

        Returns
        -------
        None
            Populates the mass_features attribute of the object with detected mass features.
        """
        all_mass_features = []

        # Split scans into partitions
        unique_scans = sorted(data_thres["scan"].unique())
        unique_scans_n = len(unique_scans)

        # Calculate partition size in scans based on goal
        partition_size_goal = 5000
        scans_per_partition = max(
            1, partition_size_goal // (len(data_thres) // unique_scans_n)
        )
        if scans_per_partition == 0:
            scans_per_partition = 1

        # Make partitions based on scans, with overlapping in partitioned scans
        scan_overlap = 4
        partition_scans = []
        for i in range(0, unique_scans_n, scans_per_partition):
            start_idx = max(0, i - scan_overlap)
            end_idx = min(
                unique_scans_n - 1, i + scans_per_partition - 1 + scan_overlap
            )
            scans_group = [int(s) for s in unique_scans[start_idx : end_idx + 1]]
            partition_scans.append(scans_group)

        # Set index to scan for faster filtering
        data_thres = data_thres.set_index("scan")
        for scans in partition_scans:
            # Determine start and end scan for partition, with 5 scans overlap
            partition_data = data_thres.loc[scans].reset_index(drop=False).copy()

            if len(partition_data) == 0:
                continue

            # Build factors for this partition
            factors = {
                dim: pd.factorize(partition_data[dim], sort=True)[1].astype(np.float32)
                for dim in dims
            }

            # Build indexes
            index = {
                dim: np.searchsorted(factors[dim], partition_data[dim]).astype(
                    np.float32
                )
                for dim in factors
            }

            # Process partition
            partition_features = self._process_partition_ph(
                partition_data, index, dims, persistence_threshold
            )

            if len(partition_features) == 0:
                continue

            # Roll up within partition
            partition_features = self.roll_up_dataframe(
                df=partition_features,
                sort_by="persistence",
                dims=["mz", "scan_time"],
                tol=[
                    self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
                    self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
                ],
                relative=[True, False],
            )
            partition_features = partition_features.reset_index(drop=True)

            if len(partition_features) > 0:
                all_mass_features.append(partition_features)

        # Combine results from all partitions
        if all_mass_features:
            combined_features = pd.concat(all_mass_features, ignore_index=True)

            # Sort by persistence
            combined_features = combined_features.sort_values(
                by="persistence", ascending=False
            ).reset_index(drop=True)

            # Remove duplicates from overlapping regions
            combined_features = self.roll_up_dataframe(
                df=combined_features,
                sort_by="persistence",
                dims=["mz", "scan_time"],
                tol=[
                    self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
                    self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
                ],
                relative=[True, False],
            )

            # resort by persistence and reset index
            combined_features = combined_features.reset_index(drop=True)

            # Populate mass_features attribute
            self._populate_mass_features(combined_features, mf_type, accumulate_features)
        else:
            self.mass_features = {}

    def _process_partition_ph(self, partition_data, index, dims, persistence_threshold):
        """Process a single partition with persistent homology."""
        # Smooth data
        iterations = self.parameters.lc_ms.ph_smooth_it
        smooth_radius = [
            self.parameters.lc_ms.ph_smooth_radius_mz,
            self.parameters.lc_ms.ph_smooth_radius_scan,
        ]

        index_array = np.vstack([index[dim] for dim in dims]).T
        V = partition_data["intensity"].values
        resid = np.inf

        for i in range(iterations):
            # Previous iteration
            V_prev = V.copy()
            resid_prev = resid
            V = self.sparse_mean_filter(index_array, V, radius=smooth_radius)

            # Calculate residual with previous iteration
            resid = np.sqrt(np.mean(np.square(V - V_prev)))

            # Evaluate convergence
            if i > 0:
                # Percent change in residual
                test = np.abs(resid - resid_prev) / resid_prev

                # Exit criteria
                if test <= 0:
                    break

        # Overwrite values
        partition_data = partition_data.copy()
        partition_data["intensity"] = V

        # Use persistent homology to find regions of interest
        pidx, pers = self.sparse_upper_star(index_array, V)
        pidx = pidx[pers > 1]
        pers = pers[pers > 1]

        if len(pidx) == 0:
            return pd.DataFrame()

        # Get peaks
        peaks = partition_data.iloc[pidx, :].reset_index(drop=True)

        # Add persistence column
        peaks["persistence"] = pers
        mass_features = peaks.sort_values(
            by="persistence", ascending=False
        ).reset_index(drop=True)

        # Filter by persistence threshold
        mass_features = mass_features.loc[
            mass_features["persistence"] > persistence_threshold, :
        ].reset_index(drop=True)

        return mass_features

    def _populate_mass_features(self, mass_features_df, mf_type="untargeted", accumulate_features=False):
        """Populate the mass_features attribute from a DataFrame.

        Parameters
        ----------
        mass_features_df : pd.DataFrame
            DataFrame containing mass feature information.
            Note that the order of this DataFrame will determine the order of mass features in the mass_features attribute.
        mf_type : str, optional
            Type label for the mass features. Default is "untargeted".
        accumulate_features : bool, optional
            If True, new features will be added to existing features rather than replacing them.
            Mass feature IDs will be offset to avoid conflicts. Default is False.

        Returns
        -------
        None, but assigns or updates the mass_features attribute to the object.
        """
        # Rename scan column to apex_scan
        mass_features_df = mass_features_df.rename(
            columns={"scan": "apex_scan", "scan_time": "retention_time"}
        )

        # Initialize or preserve existing mass_features attribute
        if accumulate_features and self.mass_features is not None and len(self.mass_features) > 0:
            # Find the maximum existing ID to offset new IDs and avoid conflicts
            id_offset = max(self.mass_features.keys()) + 1
            initial_count = len(self.mass_features)
        else:
            # Replace mode (default/backwards compatible)
            self.mass_features = {}
            id_offset = 0
            initial_count = 0
        
        # Add new mass features
        for idx, row in enumerate(mass_features_df.itertuples()):
            row_dict = mass_features_df.iloc[row.Index].to_dict()
            lcms_feature = LCMSMassFeature(self, **row_dict)
            lcms_feature.type = mf_type
            # Use sequential ID starting from id_offset to avoid conflicts with existing features
            new_id = idx + id_offset
            lcms_feature._id = new_id  # Update the internal ID
            self.mass_features[new_id] = lcms_feature

        if self.parameters.lc_ms.verbose_processing:
            if accumulate_features and initial_count > 0:
                print(f"Found {len(mass_features_df)} new mass features (total: {len(self.mass_features)})")
            else:
                print("Found " + str(len(mass_features_df)) + " initial mass features")

    def find_mass_features_ph_centroid(self, ms_level=1, targeted_search=False, target_search_dict=None, mf_type="untargeted", accumulate_features=False):
        """Find mass features within an LCMSBase object using persistent homology-type approach but with centroided data.

        Parameters
        ----------
        ms_level : int, optional
            The MS level to use. Default is 1.
        targeted_search : bool, optional
            If True, perform targeted search mode. Default is False.
        target_search_dict : dict or None, optional
            Dictionary with target parameters for targeted search. Default is None.
        mf_type : str, optional
            Type label for the mass features. Default is "untargeted".
        accumulate_features : bool, optional
            If True, add to existing features rather than replacing them. Default is False.

        Raises
        ------
        ValueError
            If no MS level data is found on the object.

        Returns
        -------
        None, but assigns the mass_features attribute to the object.
        """
        # Check that ms_level is a key in self._ms_uprocessed
        if ms_level not in self._ms_unprocessed.keys():
            raise ValueError(
                "No MS level "
                + str(ms_level)
                + " data found, did you instantiate with parser specific to MS level?"
            )

        # Work with reference instead of copy
        data = self._ms_unprocessed[ms_level]

        # Merge with scan data first (needed for filtering in targeted mode)
        scan_subset = self.scan_df[["scan", "scan_time"]]
        data_with_time = data.merge(scan_subset, on="scan", how="inner")
        
        # Calculate threshold and filter (bypass in targeted mode)
        if targeted_search:
            # In targeted mode, bypass intensity threshold
            threshold = 0
            valid_mask = data_with_time["intensity"].notna()
            required_cols = ["mz", "intensity", "scan", "scan_time"]
            data_thres = data_with_time.loc[valid_mask, required_cols].copy()
            
            # Filter to target windows
            data_thres = self._filter_data_by_targets(data_thres, target_search_dict)
            
            if len(data_thres) == 0:
                if self.parameters.lc_ms.verbose_processing:
                    print("No data found in target windows")
                self.mass_features = {}
                return
        else:
            # Normal mode with threshold
            max_intensity = data_with_time["intensity"].max()
            threshold = self.parameters.lc_ms.ph_inten_min_rel * max_intensity
            valid_mask = data_with_time["intensity"].notna() & (data_with_time["intensity"] > threshold)
            required_cols = ["mz", "intensity", "scan", "scan_time"]
            data_thres = data_with_time.loc[valid_mask, required_cols].copy()
        
        data_thres["persistence"] = data_thres["intensity"]
        mf_df = data_thres
        del data_thres, scan_subset, data_with_time

        # Order by scan_time and then mz to ensure features near in rt are processed together
        # It's ok that different scans are in different partitions; we will roll up later
        mf_df = mf_df.sort_values(
            by=["scan_time", "mz"], ascending=[True, True]
        ).reset_index(drop=True)
        partition_size = 10000
        partitions = [
            mf_df.iloc[i : i + partition_size].reset_index(drop=True)
            for i in range(0, len(mf_df), partition_size)
        ]
        del mf_df

        # Run roll_up_dataframe on each partition
        rolled_partitions = []
        for part in partitions:
            rolled = self.roll_up_dataframe(
                df=part,
                sort_by="persistence",
                dims=["mz", "scan_time"],
                tol=[
                    self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
                    self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
                ],
                relative=[True, False],
            )
            rolled_partitions.append(rolled)
        del partitions

        # Run roll_up_dataframe on the rolled_up partitions to merge features near partition boundaries

        # Combine results and run a final roll-up to merge features near partition boundaries
        mf_df_final = pd.concat(rolled_partitions, ignore_index=True)
        del rolled_partitions

        # Reorder by persistence before final roll-up
        mf_df_final = mf_df_final.sort_values(
            by="persistence", ascending=False
        ).reset_index(drop=True)

        mf_df_final = self.roll_up_dataframe(
            df=mf_df_final,
            sort_by="persistence",
            dims=["mz", "scan_time"],
            tol=[
                self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
                self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
            ],
            relative=[True, False],
        )
        # reset index
        mf_df_final = mf_df_final.reset_index(drop=True)

        # Combine rename and sort operations
        mass_features = (
            mf_df_final.rename(
                columns={"scan": "apex_scan", "scan_time": "retention_time"}
            )
            .sort_values(by="persistence", ascending=False)
            .reset_index(drop=True)
        )
        del mf_df_final  # Free memory

        # Order by persistence and reset index
        mass_features = mass_features.sort_values(
            by="persistence", ascending=False
        ).reset_index(drop=True)

        self.mass_features = {}
        for idx, row in mass_features.iterrows():
            row_dict = row.to_dict()
            lcms_feature = LCMSMassFeature(self, **row_dict)
            lcms_feature.type = mf_type
            self.mass_features[lcms_feature.id] = lcms_feature

        if self.parameters.lc_ms.verbose_processing:
            print("Found " + str(len(mass_features)) + " initial mass features")
    
    def cluster_mass_features(self, drop_children=True, sort_by="persistence"):
        """Cluster mass features

        Based on their proximity in the mz and scan_time dimensions, priorizies the mass features with the highest persistence.

        Parameters
        ----------
        drop_children : bool, optional
            Whether to drop the mass features that are not cluster parents. Default is True.
        sort_by : str, optional
            The column to sort the mass features by, this will determine which mass features get rolled up into a parent mass feature. Default is "persistence".

        Raises
        ------
        ValueError
            If no mass features are found.
            If too many mass features are found.

        Returns
        -------
        None if drop_children is True, otherwise returns a list of mass feature ids that are not cluster parents.
        """
        if self.mass_features is None:
            raise ValueError("No mass features found, run find_mass_features() first")
        if len(self.mass_features) > 400000:
            raise ValueError(
                "Too many mass features of interest found, run find_mass_features() with a higher intensity threshold"
            )
        dims = ["mz", "scan_time"]
        mf_df_og = self.mass_features_to_df()
        mf_df = mf_df_og.copy()

        tol = [
            self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
            self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
        ]  # mz, in relative; scan_time in minutes
        relative = [True, False]

        # Roll up mass features based on their proximity in the declared dimensions
        mf_df_new = self.roll_up_dataframe(
            df=mf_df, sort_by=sort_by, dims=dims, tol=tol, relative=relative
        )

        mf_df["cluster_parent"] = np.where(
            np.isin(mf_df.index, mf_df_new.index), True, False
        )

        # get mass feature ids of features that are not cluster parents
        cluster_daughters = mf_df[~mf_df["cluster_parent"]].index.values
        if drop_children is True:
            # Drop mass features that are not cluster parents from self
            self.mass_features = {
                k: v
                for k, v in self.mass_features.items()
                if k not in cluster_daughters
            }
        else:
            return cluster_daughters


class LCMSCollectionCalculations:
    """Methods for performing calculations related to LCMSCollection objects.

    Notes
    -----
    This class is intended as a mixin for the LCMSCollection class.
    """

    @staticmethod
    def _plot_multiple_eics(ax, cluster_mfs, induced_cluster_mfs, rep_sample_id, rep_mf_id,
                           median_rt, eic_buffer_time, plot_smoothed=False, 
                           plot_datapoints=False, label_samples=False, lcms_collection=None):
        """Internal method to plot multiple EICs from different samples on a given axis.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to plot on.
        cluster_mfs : pd.DataFrame
            DataFrame containing cluster mass features (non-induced).
        induced_cluster_mfs : pd.DataFrame or None
            DataFrame containing induced (gap-filled) mass features.
        rep_sample_id : int
            Sample ID of the representative mass feature.
        rep_mf_id : int
            Mass feature ID of the representative mass feature.
        median_rt : float
            Median retention time for the cluster.
        eic_buffer_time : float
            Time buffer around the peak (minutes).
        plot_smoothed : bool, optional
            If True, plot smoothed EICs. Default is False.
        plot_datapoints : bool, optional
            If True, plot EIC datapoints. Default is False.
        label_samples : bool, optional
            If True, label each sample individually. Default is False.
        lcms_collection : LCMSCollection, optional
            The parent collection object for accessing samples. Required.
        """
        ax.set_title("EICs from all samples", loc="left")
        
        # Track if we've added labels for legend (to avoid duplicates)
        rep_labeled = False
        regular_labeled = False
        induced_labeled = False
        
        # Plot regular (non-induced) mass features
        for _, row in cluster_mfs.iterrows():
            sample_id = int(row['sample_id'])
            mf_id = row['mf_id']
            sample = lcms_collection[sample_id]
            sample_name = row['sample_name']
            
            # Get EIC using eic_mz column from dataframe
            eic_mz = row.get('_eic_mz')
            if eic_mz is not None and not pd.isna(eic_mz) and hasattr(sample, 'eics') and sample.eics:
                eic_data = sample.eics.get(eic_mz)
            else:
                eic_data = None
            
            if eic_data:
                # Determine line style and width
                if sample_id == rep_sample_id and mf_id == rep_mf_id:
                    # Representative feature - bold line
                    linewidth = 2.5
                    alpha = 1.0
                    color = 'tab:blue'
                    if label_samples:
                        label = f"{sample_name} (representative)"
                    else:
                        label = "Representative" if not rep_labeled else None
                        rep_labeled = True
                else:
                    # Other features - thinner line
                    linewidth = 1.0
                    alpha = 0.5
                    color = 'tab:blue'
                    if label_samples:
                        label = sample_name
                    else:
                        label = "Regular features" if not regular_labeled else None
                        regular_labeled = True
                
                ax.plot(
                    eic_data.time,
                    eic_data.eic,
                    c=color,
                    linewidth=linewidth,
                    alpha=alpha,
                    linestyle='-',
                    label=label
                )
                
                if plot_datapoints:
                    ax.scatter(
                        eic_data.time,
                        eic_data.eic,
                        c=color,
                        alpha=alpha,
                        s=10
                    )
                
                if plot_smoothed and hasattr(eic_data, 'eic_smoothed'):
                    ax.plot(
                        eic_data.time,
                        eic_data.eic_smoothed,
                        c=color,
                        linestyle='--',
                        alpha=alpha * 0.8,
                        linewidth=linewidth * 0.8
                    )
        
        # Plot induced (gap-filled) mass features if available
        if induced_cluster_mfs is not None and not induced_cluster_mfs.empty:
            for _, row in induced_cluster_mfs.iterrows():
                sample_id = int(row['sample_id'])
                mf_id = row['mf_id']
                sample = lcms_collection[sample_id]
                sample_name = row['sample_name']
                
                # Get EIC using eic_mz column from dataframe
                eic_mz = row.get('_eic_mz')
                if eic_mz is not None and not pd.isna(eic_mz) and hasattr(sample, 'eics') and sample.eics:
                    eic_data = sample.eics.get(eic_mz)
                else:
                    eic_data = None
                
                if eic_data:
                    # Induced features - even thinner line
                    linewidth = 0.5
                    alpha = 0.4
                    color = 'tab:orange'
                    
                    if label_samples:
                        label = f"{sample_name} (induced)"
                    else:
                        label = "Gap-filled features" if not induced_labeled else None
                        induced_labeled = True
                    
                    ax.plot(
                        eic_data.time,
                        eic_data.eic,
                        c=color,
                        linewidth=linewidth,
                        alpha=alpha,
                        linestyle='-',
                        label=label
                    )
                    
                    if plot_datapoints:
                        ax.scatter(
                            eic_data.time,
                            eic_data.eic,
                            c=color,
                            alpha=alpha,
                            s=5
                        )
                    
                    if plot_smoothed and hasattr(eic_data, 'eic_smoothed'):
                        ax.plot(
                            eic_data.time,
                            eic_data.eic_smoothed,
                            c=color,
                            linestyle='--',
                            alpha=alpha * 0.8,
                            linewidth=linewidth * 0.8
                        )
        
        # Add vertical line at median RT
        ax.axvline(
            x=median_rt,
            color='k',
            linestyle='--',
            alpha=0.7,
            label='Median RT'
        )
        
        ax.set_ylabel("Intensity")
        ax.set_xlabel("Time (minutes)")
        ax.set_xlim(
            median_rt - eic_buffer_time,
            median_rt + eic_buffer_time,
        )
        ax.legend(loc='upper left', fontsize=8)
        ax.yaxis.get_major_formatter().set_useOffset(False)

    def clean_sparse_matrix(self, sparse_matrix):
        """Clean a sparse matrix by removing duplicates and sorting.

        Parameters
        ----------
        sparse_matrix : :obj:`~numpy.array`
            A sparse matrix to clean.

        Returns
        -------
        :obj:`~numpy.array`
            A cleaned sparse matrix.
        """
        for match in sparse_matrix:
            match.sort()
        sparse_matrix.sort()
        dereplicated_sparse_matrix = np.unique(sparse_matrix, axis=0)
        return dereplicated_sparse_matrix

    def match_mfs(self, mf_c, mf_i):
        """Match mass features between two LCMS objects.

        Parameters
        ----------
        mf_c : :obj:`~pandas.DataFrame`
            The mass features to match against.
        mf_i : :obj:`~pandas.DataFrame`
            The mass features to match.

        Returns
        -------
        :obj:`~pandas.DataFrame`
            The matched mass features from mf_c.
        :obj:`~pandas.DataFrame`
            The matched mass features from mf_i.

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package:
        https://github.com/pnnl/deimos
        """
        if mf_c is None or mf_i is None or len(mf_c.index) < 1 or len(mf_i.index) < 1:
            return None, None

        # Prepare dataframes
        mf_c = mf_c.copy()
        mf_c["id_i"] = 0
        mf_i = mf_i.copy()
        mf_i["id_i"] = 1

        # Set dimensions for matching
        dims = ["mz", "scan_time"]
        relative = [True, False]
        mz_tol = self.parameters.lcms_collection.alignment_mz_tol_ppm * 1e-6
        rt_tol = self.parameters.lcms_collection.alignment_rt_tol
        tol = [mz_tol, rt_tol]

        # Compute inter-feature distances
        idx = []
        for i, f in enumerate(dims):
            # vectors
            v1 = mf_c[f].values.reshape(-1, 1)
            v2 = mf_i[f].values.reshape(-1, 1)

            # Distances
            d = scipy.spatial.distance.cdist(v1, v2)

            if relative[i] is True:
                # Divisor
                basis = np.repeat(v1, v2.shape[0], axis=1)
                fix = np.repeat(v2, v1.shape[0], axis=1).T
                basis = np.where(basis == 0, fix, basis)

                # Divide
                d = np.divide(d, basis, out=np.zeros_like(basis), where=basis != 0)

            # Check tol
            idx.append(d <= tol[i])

        # Stack truth arrays
        idx = np.prod(np.dstack(idx), axis=-1, dtype=bool)

        # Compute normalized 3d distance
        v1 = mf_c[dims].values / tol
        v2 = mf_i[dims].values / tol
        dist3d = scipy.spatial.distance.cdist(v1, v2, "cityblock")
        
        # Separate features within tolerance from those outside
        # Features outside tolerance should be inf, features within tolerance keep their distance
        # Use idx mask: True for within tolerance, False for outside
        dist3d_within_tol = np.where(idx, dist3d, np.inf)

        # Normalize to 0-1 (only affects within-tolerance distances)
        mx = np.max(dist3d_within_tol[idx]) if np.sum(idx) > 0 else 0
        if mx > 0:
            # Lower distance is better - normalize only the within-tolerance values
            dist3d_within_tol = np.where(idx, dist3d_within_tol / mx, np.inf)
        else:
            # All matches are perfect (distance=0), assign tiny value to within-tolerance pairs
            dist3d_within_tol = np.where(idx, 1e-10, np.inf)
        
        # Use the masked distance matrix
        dist3d = dist3d_within_tol

        # Min over dims
        mincols = np.min(dist3d, axis=0, keepdims=True)

        # Zero out mincols over dims
        dist3d[dist3d != mincols] = np.inf

        # Min over clusters
        minrows = np.min(dist3d, axis=1, keepdims=True)

        # Where max and nonzero
        ii, jj = np.where((dist3d == minrows) & (dist3d < np.inf))

        # Reorder
        mf_c = mf_c.iloc[ii]
        mf_i = mf_i.iloc[jj]

        if len(mf_c.index) < 1 or len(mf_i.index) < 1:
            return None, None

        return mf_c, mf_i

    def fit_rts(self, a, b, align="scan_time", **kwargs):
        """
        Fit a support vector regressor to matched features.

        Parameters
        ----------
        a : :obj:`~pandas.DataFrame`
            First set of input feature coordinates and intensities; the center object and the object to align to.
        b : :obj:`~pandas.DataFrame`
            Second set of input feature coordinates and intensities; the object to align to the center object.
        align : str
            Dimension to align.
        kwargs
            Keyword arguments for support vector regressor
            (:class:`sklearn.svm.SVR`).

        Returns
        -------
        :obj:`~function`
            An interpolation function where one can input a retention time and get the predicted retention time.

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package:
        https://github.com/pnnl/deimos

        """

        # Uniqueify
        x = a[align].values
        y = b[align].values
        arr = np.vstack((x, y)).T
        arr = np.unique(arr, axis=0)

        # Safety check: ensure we have data to work with
        if len(arr) == 0:
            warnings.warn("No data points available for retention time fitting. Returning identity function.")
            return lambda x: x

        # Check kwargs
        if "kernel" in kwargs:
            kernel = kwargs.get("kernel")
        else:
            kernel = "linear"

        # Construct interpolation axis
        newx = np.linspace(arr[:, 0].min(), arr[:, 0].max(), 1000)

        # Linear kernel
        if kernel == "linear":
            reg = scipy.stats.linregress(x, y)
            newy = reg.slope * newx + reg.intercept

        # Other kernels
        else:
            # Fit
            svr = SVR(**kwargs)
            svr.fit(arr[:, 1].reshape(-1, 1), arr[:, 0])

            # Predict
            newy = svr.predict(newx.reshape(-1, 1))

        # Pad x and y_pred with zeros to force interpolation to start at 0
        newx = np.concatenate(([0], newx))
        newy = np.concatenate(([0], newy))

        # Pad x and y_pred with max time to force interpolation to end at max time to force interpolation to match at end max time
        max_time = self[0].scan_df["scan_time"].max()
        newx = np.concatenate((newx, [max_time]))
        newy = np.concatenate((newy, [max_time]))

        # Return an interpolation function for the x and y_pred
        def interp(x):
            pred_y = np.interp(x, newx, newy)
            return pred_y

        return interp

    def get_anchor_mass_features(self, mf_df):
        """
        Get the anchor mass features from a DataFrame of mass features.

        Parameters
        ----------
        mf_df : :obj:`~pandas.DataFrame`
            The mass features to filter to just the anchor mass features.

        Returns
        -------
        :obj:`~pandas.DataFrame`
            The anchor mass features dataframe.
        """
        mf_df = mf_df.copy()

        if (
            "deconvoluted_mass_spectra"
            in self.parameters.lcms_collection.mass_feature_anchor_technique
        ):
            # Drop features that are not mass_spectrum_deconvoluted_parent or are NA as mass_spectrum_deconvoluted_parent
            mf_df = mf_df.dropna(subset=["mass_spectrum_deconvoluted_parent"])
            mf_df = mf_df[mf_df["mass_spectrum_deconvoluted_parent"]]

        if (
            "absolute_intensity"
            in self.parameters.lcms_collection.mass_feature_anchor_technique
        ):
            # Drop features that have an intensity lower than the threshold
            threshold = self.parameters.lcms_collection.mass_feature_anchor_absolute_intensity_threshold
            mf_df = mf_df[mf_df["intensity"] > threshold]

        if (
            "relative_intensity"
            in self.parameters.lcms_collection.mass_feature_anchor_technique
        ):
            # Drop features in the lower fraction of intensities
            threshold_quantile = self.parameters.lcms_collection.mass_feature_anchor_relative_intensity_threshold
            intensity_threshold = mf_df["intensity"].quantile(threshold_quantile)
            mf_df = mf_df[mf_df["intensity"] >= intensity_threshold]

        return mf_df

    def attempt_alignment(self, matches_c, matches_i):
        """
        Check if alignment is needed for the LCMS objects in the collection.
        """

        # Hold out a subset of matches_c and matches_i for spline fitting
        matches_c.reset_index(drop=False, inplace=True)
        matches_i.reset_index(drop=False, inplace=True)

        # Check if there are enough matches to attempt alignment
        minimum_matches = self.parameters.lcms_collection.alignment_minimum_matches
        if len(matches_c) < minimum_matches:
            # Return False (no alignment) and identity function (returns original time) 
            # which isn't used but is a placeholder to avoid errors in downstream code since 
            # the function expects a callable to be returned
            return False, lambda x: x

        # Rearrange matches_c and matches_i to be in the order of the scan_time of matches_c
        matches_c = matches_c.sort_values(by="scan_time")
        matches_i = matches_i.iloc[matches_c.index.values]

        hold_out_fraction = self.parameters.lcms_collection.alignment_hold_out_fraction
        # starting with an array of length len(matches_c), select equally spaced indices to hold out
        idx_holdout = matches_c.index.values[
            np.arange(0, len(matches_c), int(1 / hold_out_fraction))
        ]

        matches_c_holdout = matches_c.loc[idx_holdout].copy()
        matches_i_holdout = matches_i.loc[idx_holdout].copy()

        # Remove the holdout matches from the matches_c and matches_i DataFrames and reset the index
        matches_c = matches_c.drop(index=idx_holdout).set_index("sample_name")
        matches_i = matches_i.drop(index=idx_holdout).set_index("sample_name")

        # Reset the scan_time to the original scan_time
        matches_i = matches_i.copy()
        matches_i["scan_time"] = matches_i["scan_time_og"]

        # Fit the retention times of the LCMS object to the center LCMS object using the matched mass features
        spl = self.fit_rts(matches_c, matches_i, kernel="rbf", C=1000)

        # Check if the spline fitting improved the alignment for the holdout matches
        matches_i_holdout["scan_time_fit"] = spl(matches_i_holdout["scan_time"])
        og_diff = np.abs(
            matches_i_holdout["scan_time"] - matches_c_holdout["scan_time"]
        )
        fit_diff = np.abs(
            matches_i_holdout["scan_time_fit"] - matches_c_holdout["scan_time"]
        )

        if (
            "fraction_improved"
            in self.parameters.lcms_collection.alignment_acceptance_technique
        ):
            fraction_improved = np.sum(fit_diff < og_diff) / len(og_diff)
            use_spline_alignment = (
                fraction_improved
                > self.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold
            )
        if (
            "mean_squared_error_improved"
            in self.parameters.lcms_collection.alignment_acceptance_technique
        ):
            mse_og = np.mean(og_diff**2)
            mse = np.mean(fit_diff**2)
            use_spline_alignment = mse < mse_og
            # Convert to boolean
            use_spline_alignment = bool(use_spline_alignment)

        return use_spline_alignment, spl

    def align_lcms_objects(self, overwrite=False):
        """
        Align LCMS objects in the collection.

        Aligns the LCMS objects in the collection by aligning the retention times of the mass features in the LCMS objects.
        First, the mass features in the center LCMS object are matched to the mass features in the other LCMS objects,
        starting with the LCMS object immediately following the center LCMS object. The retention times of the LCMS objects
        are then fit to the center LCMS object using the matched mass features.

        Returns
        -------
        None, but aligns the LCMS objects in the collection and sets the scan_time_aligned column in the scan_df attribute of each LCMS object.

        Notes
        -----
        This function has been adapted from the original implementation in the Deimos package:
        https://github.com/pnnl/deimos
        """
       
        # Prepare the center LCMS object
        center_obj_ids = self.manifest_dataframe[
            self.manifest_dataframe["center"]
        ].collection_id.values

        full_mf_df = self.mass_features_dataframe
        # re-index to sample_name for faster lookups
        full_mf_df = full_mf_df.reset_index().set_index("sample_name")

        if "scan_time_aligned" in full_mf_df.columns and not overwrite:
            raise ValueError("Mass features have already been aligned")

        anchor_mf_dfs = []
        for center_obj_id in center_obj_ids:
            # Get the anchor mass features from the center LCMS object
            mf_df_c = full_mf_df.loc[self.samples[center_obj_id]]
            mf_df_c = self.get_anchor_mass_features(mf_df_c)
            anchor_mf_dfs.append(mf_df_c)

            # Set scan_time_aligned to scan_time for the center LCMS object
            center_scan_df = self[center_obj_id].scan_df.copy()
            center_scan_df["scan_time_aligned"] = center_scan_df["scan_time"]
            self[center_obj_id].scan_df = center_scan_df
            
            # Store alignment data for center object (identity mapping)
            center_sample_name = self.samples[center_obj_id]

            index_steps = (1, -1)
            # Run this twice, once going forward (+1 indexing) and once going backward (-1 indexing)
            for index_step in index_steps:
                # Initialize spline for propagation to samples without features
                spl = None
                use_spline_alignment = False
                
                # Loop through the other LCMS objects in the collection (going forward)
                i = center_obj_id + index_step
                if i < len(self) and i >= 0:
                    # Check if this sample has any features in the dataframe
                    sample_name = self.samples[i]
                    if sample_name not in full_mf_df.index.get_level_values('sample_name'):
                        # For samples with no mass features, use the same alignment as the previous sample
                        if use_spline_alignment and spl is not None:
                            # Use the spline from the adjacent sample
                            self[i]._scan_info["scan_time_aligned"] = {k: spl(v) for k, v in self[i]._scan_info["scan_time"].items()}
                        else:
                            # No spline available, use original times
                            self[i]._scan_info["scan_time_aligned"] = self[i]._scan_info["scan_time"].copy()
                        self.rt_alignment_attempted = True
                        
                        # Move to next sample
                        i += index_step
                        while i < len(self) and i >= 0:
                            sample_name = self.samples[i]
                            if sample_name not in full_mf_df.index.get_level_values('sample_name'):
                                # Apply same alignment to this empty sample
                                if use_spline_alignment and spl is not None:
                                    self[i]._scan_info["scan_time_aligned"] = {k: spl(v) for k, v in self[i]._scan_info["scan_time"].items()}
                                else:
                                    self[i]._scan_info["scan_time_aligned"] = self[i]._scan_info["scan_time"].copy()
                                i += index_step
                            else:
                                # Found a sample with features, exit inner loop to process it
                                break
                        
                        # If we've processed all remaining empty samples, continue to next index_step
                        if i >= len(self) or i < 0:
                            continue
                    
                    # Grab the first LCMS object after the center object
                    mf_df_i = full_mf_df.loc[sample_name].copy()
                    mf_df_i["scan_time_og"] = mf_df_i["scan_time"]

                    while mf_df_i is not None:
                        mf_df_i = self.get_anchor_mass_features(mf_df_i)

                        # Match the mass features in the LCMS object to the anchor mass features in the center LCMS object.
                        matches_c, matches_i = self.match_mfs(mf_df_c, mf_df_i)

                        if matches_c is not None:
                            use_spline_alignment, spl = self.attempt_alignment(
                                matches_c, matches_i
                            )

                            # Record if we used alignment for this sample
                            sample_name = self.samples[i]
                            self._manifest_dict[sample_name]["use_rt_alignment"] = (
                                use_spline_alignment
                            )

                            if use_spline_alignment:
                                # Set new retention times on scan_df for lc_obj using the spline fitting
                                matches_i["scan_time_fit"] = spl(matches_i["scan_time"])

                                # Add "scan_time_aligned" to LCMSObject's _scan_info dict
                                self[i]._scan_info["scan_time_aligned"] = {k: spl(v) for k, v in self[i]._scan_info["scan_time"].items()}

                                # Retrieve the new aligned times for all scans in the LCMS object
                                new_times = [x for k, x in sorted(self[i]._scan_info["scan_time_aligned"].items())]
                                
                                # Switch the rt_aligned flag to True and attempted to True
                                self.rt_aligned = True
                                self.rt_alignment_attempted = True
                            else:
                                # Set aligned retention times on scan_df for lc_obj using the original retention times
                                self[i]._scan_info["scan_time_aligned"] = self[i]._scan_info["scan_time"].copy()
                                # Switch the rt_attempted flag to True
                                self.rt_aligned = False
                                self.rt_alignment_attempted = True

                            i += index_step
                            if i >= len(self) or i < 0:
                                mf_df_i = None
                            else:
                                # Check if this sample has any features in the dataframe
                                sample_name_next = self.samples[i]
                                if sample_name_next not in full_mf_df.index.get_level_values('sample_name'):
                                    # For samples with no mass features, apply the same alignment transformation
                                    if use_spline_alignment and spl is not None:
                                        self[i]._scan_info["scan_time_aligned"] = {k: spl(v) for k, v in self[i]._scan_info["scan_time"].items()}
                                    else:
                                        self[i]._scan_info["scan_time_aligned"] = self[i]._scan_info["scan_time"].copy()
                                    self.rt_alignment_attempted = True
                                    
                                    # Continue to the next sample
                                    i += index_step
                                    while i < len(self) and i >= 0:
                                        sample_name = self.samples[i]
                                        if sample_name not in full_mf_df.index.get_level_values('sample_name'):
                                            # Apply same alignment to consecutive empty samples
                                            if use_spline_alignment and spl is not None:
                                                self[i]._scan_info["scan_time_aligned"] = {k: spl(v) for k, v in self[i]._scan_info["scan_time"].items()}
                                            else:
                                                self[i]._scan_info["scan_time_aligned"] = self[i]._scan_info["scan_time"].copy()
                                            i += index_step
                                        else:
                                            # Found a sample with features
                                            break
                                    
                                    # Check if we're done
                                    if i >= len(self) or i < 0:
                                        mf_df_i = None
                                    else:
                                        # Grab the next LCMS object with features
                                        mf_df_i = full_mf_df.loc[self.samples[i]].copy()
                                        mf_df_i["scan_time_og"] = mf_df_i["scan_time"]
                                        mf_df_i = mf_df_i.reset_index(drop=False)
                                        if use_spline_alignment:
                                            # Set scan_time to previous sample's predicted scan_time to find closer matches
                                            mf_df_i["scan_time"] = spl(mf_df_i["scan_time"])
                                else:
                                    # Grab the next LCMS object and use the previous spline fitting to get a better starting point
                                    mf_df_i = full_mf_df.loc[sample_name_next].copy()
                                    mf_df_i["scan_time_og"] = mf_df_i["scan_time"]
                                    mf_df_i = mf_df_i.reset_index(drop=False)
                                    if use_spline_alignment:
                                        # Set scan_time to previous sample's predicted scan_time to find closer matches
                                        mf_df_i["scan_time"] = spl(mf_df_i["scan_time"])
                        else:
                            raise ValueError(
                                f"No matches found between the center object and {self.samples[i]}"
                            )

        # Now align each batch using the center objects as anchors with the other batches
        mf_df_c = anchor_mf_dfs[0]
        for i in center_obj_ids[1:]:
            mf_df_i = full_mf_df.loc[self.samples[i]].copy()
            mf_df_i["scan_time_og"] = mf_df_i["scan_time"]
            mf_df_i = self.get_anchor_mass_features(mf_df_i)

            matches_c, matches_i = self.match_mfs(mf_df_c, mf_df_i)
            if matches_c is not None:
                use_spline_alignment, spl = self.attempt_alignment(matches_c, matches_i)

                # Record if we used alignment for this sample
                sample_name = self.samples[i]
                self._manifest_dict[sample_name]["use_rt_alignment"] = (
                    use_spline_alignment
                )

                if use_spline_alignment:
                    # Set new retention times on all this object's
                    new_times = spl(self[i].scan_df["scan_time"])
                    new_scan_info = self[i].scan_df.copy()
                    new_scan_info["scan_time_aligned"] = new_times
                    self[i].scan_df = new_scan_info
                    

                    # Get the batch that this object belongs to
                    batch = self.manifest[self.samples[i]]["batch"]

                    for j in range(len(self)):
                        if self.manifest[self.samples[j]]["batch"] == batch:
                            if j != i:
                                sample_name_j = self.samples[j]
                                self._manifest_dict[sample_name_j]["use_rt_alignment"] = (
                                    use_spline_alignment
                                )
                                new_scan_info = self[j].scan_df.copy()
                                aligned_times = spl(self[j].scan_df["scan_time_aligned"])
                                new_scan_info["scan_time_aligned"] = aligned_times
                                self[j].scan_df = new_scan_info
                                
        # Set final mass_features_dataframe with the aligned scan_time
        center_sample_name = self.samples[center_obj_ids[0]]
        self._manifest_dict[center_sample_name]["use_rt_alignment"] = False
        new_scan_info = self[center_obj_ids[0]].scan_df.copy()
        new_scan_info["scan_time_aligned"] = new_scan_info["scan_time"]

    def add_consensus_mass_features(self):
        """
        Create consensus mass features by clustering aligned features across samples.
        
        This method clusters mass features from all samples in the collection based on
        their m/z and aligned retention time proximity. Features that cluster together
        across samples are assigned a common cluster ID, creating consensus features
        that represent the same compound detected across multiple samples.
        
        The clustering process:
        1. Partitions features by m/z to avoid large sparse matrices and enable parallelization
        2. Clusters features within each partition using hierarchical clustering
        3. Merges partition-boundary clusters that represent the same feature
        4. Filters out clusters not present in minimum fraction of samples
        
        Must be run after align_lcms_objects(). Results are stored in the 
        mass_features_dataframe with a 'cluster' column added.
        
        Parameters
        ----------
        None
            Uses parameters from self.parameters.lcms_collection:
            - consensus_mz_tol_ppm: m/z tolerance for clustering (ppm)
            - consensus_rt_tol: retention time tolerance for clustering (minutes)
            - consensus_partition_size: target partition size for managing memory and parallelization
            - consensus_min_sample_fraction: minimum fraction of samples a cluster
              must appear in to be retained (0-1)
            - cores: number of CPU cores to use for parallel partition processing
            
        Returns
        -------
        None
            Updates self.mass_features_dataframe in place by adding 'cluster' column
            and filtering to retain only clusters meeting minimum sample presence.
            
        Raises
        ------
        ValueError
            If mass features have not been aligned (run align_lcms_objects() first).
            
        Notes
        -----
        - Partitioning prevents memory issues with large sparse distance matrices
        - Each partition is processed in parallel (up to cores limit)
        - Clusters not meeting consensus_min_sample_fraction are automatically removed
        - Access cluster_summary_dataframe property for summary statistics
        - Use fill_missing_cluster_features() for gap-filling after clustering
        
        See Also
        --------
        align_lcms_objects : Aligns retention times before consensus clustering
        cluster_summary_dataframe : Property that generates summary statistics for clusters
        fill_missing_cluster_features : Gap-fill missing features in clusters
        """
        # Get the combined mass features from all LCMS objects, keep the original index as a separate column
        combined_mfs = self.mass_features_dataframe.copy()
        combined_mfs["coll_mf_id"] = combined_mfs.index

        # Check if the mass features have been aligned
        if "scan_time_aligned" not in combined_mfs.columns:
            raise ValueError(
                "Mass features have not been aligned, run align_lcms_objects() first"
            )

        # Partition the mass features by mz so we can parallelize the matching before clustering
        from corems.chroma_peak.calc import subset as corems_subset

        # get max mz from combined_mfs and calculate tolerance from ppm
        mz_tol = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
        n_partition_size = self.parameters.lcms_collection.consensus_partition_size
        lazy_partitions = corems_subset.multi_sample_partition(
            combined_mfs,
            split_on="mz",
            size=n_partition_size,
            tol=mz_tol,
            relative=True,
        )

        # If any of lazy_partitions._counts is 2xn_partition_size, issue a warning
        if np.array(lazy_partitions._counts).max() > 2 * n_partition_size:
            warnings.warn(
                "Some partitions are larger than 2x the goal partition size. Consider increasing the partition or decreasing the mz_tol."
            )

        # Cluster the mass features within each partition
        if self.parameters.lcms_collection.cores > lazy_partitions.n_partitions:
            cores_to_use = lazy_partitions.n_partitions
        else:
            cores_to_use = self.parameters.lcms_collection.cores
        # mfs_with_clusters = lazy_partitions.map(self.cluster_mass_features, processes=cores_to_use)
        mfs_with_clusters = lazy_partitions.map(
            self.cluster_mass_features_agg_cluster, processes=cores_to_use
        )

        # Clean up cluster id names after partitioning
        new_cluster_ids = (
            mfs_with_clusters[["cluster", "partition_idx"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        new_cluster_ids["cluster_unqiue"] = new_cluster_ids.index
        mfs_with_clusters = mfs_with_clusters.merge(
            new_cluster_ids, on=["cluster", "partition_idx"]
        )
        mfs_with_clusters["cluster"] = mfs_with_clusters["cluster_unqiue"]
        mfs_with_clusters = mfs_with_clusters.drop(columns=["cluster_unqiue"])

        # Embed a new cluster id into the mass features dataframe and set as index
        mfs_with_clusters["idx"] = mfs_with_clusters.index

        try:
            # Check if any clusters can be merged into a single cluster
            eval_dict = self.evaluate_clusters_for_repeats(mfs_with_clusters)

            # Merge clusters identified in eval_dict
            while len(eval_dict["merge_these_clusters"]) > 0:
                list_of_clusters_to_merge = [
                    [x[0], x[1]] for x in eval_dict["merge_these_clusters"]
                ]
                # Convert to a dataframe with columns "new_cluster" and "cluster"
                df = pd.DataFrame(
                    np.array(list_of_clusters_to_merge), columns=["new_cluster", "cluster"]
                )
                # Drop duplicates of "child" clusters
                df = df.drop_duplicates("cluster", keep="first")
                df = df.drop_duplicates("new_cluster", keep="first")
                mfs_with_clusters = mfs_with_clusters.merge(df, on="cluster", how="left")
                mfs_with_clusters["cluster"] = mfs_with_clusters["new_cluster"].fillna(
                    mfs_with_clusters["cluster"]
                )
                mfs_with_clusters = mfs_with_clusters.drop(columns=["new_cluster"])

                # Re-evaluate clusters for repeats
                eval_dict = self.evaluate_clusters_for_repeats(mfs_with_clusters)
                self.mass_features_dataframe = mfs_with_clusters

        except:
            mfs_with_clusters.set_index('coll_mf_id', inplace = True)
            self.mass_features_dataframe = mfs_with_clusters
            
        # Filter out clusters that don't meet minimum sample fraction
        self._filter_clusters_by_sample_presence()
            
        # TODO KRH: Deal with isomers better? Pool them together and then split them out using samples with 2 as the template?
    
    def _filter_clusters_by_sample_presence(self):
        """
        Filter out clusters that don't meet the minimum sample fraction threshold.
        
        Removes clusters (and their associated mass features) from the mass_features_dataframe
        if they don't appear in at least consensus_min_sample_fraction of samples.
        
        This is called automatically at the end of add_consensus_mass_features().
        
        Returns
        -------
        None
            Updates self.mass_features_dataframe in place by removing clusters that don't
            meet the minimum sample presence threshold.
        """
        if self.mass_features_dataframe is None or len(self.mass_features_dataframe) == 0:
            return
        
        min_sample_fraction = self.parameters.lcms_collection.consensus_min_sample_fraction
        
        # Validate parameter
        if not 0 <= min_sample_fraction <= 1:
            raise ValueError("consensus_min_sample_fraction must be between 0 and 1")
        
        # Calculate minimum number of samples required
        total_samples = len(self.samples)
        min_samples_required = min_sample_fraction * total_samples
        
        # Count unique samples per cluster
        cluster_sample_counts = (
            self.mass_features_dataframe.groupby('cluster')['sample_id']
            .nunique()
            .reset_index(name='sample_count')
        )
        
        # Identify clusters to keep
        clusters_to_keep = cluster_sample_counts[
            cluster_sample_counts['sample_count'] > min_samples_required
        ]['cluster'].values
        
        # Filter mass features dataframe
        self.mass_features_dataframe = self.mass_features_dataframe[
            self.mass_features_dataframe['cluster'].isin(clusters_to_keep)
        ]
        
    def summarize_clusters(self):
        """
        Generate summary statistics for consensus mass feature clusters.
        
        Computes aggregate statistics (median, mean, std, min, max) for each cluster
        across all samples. Combines both regular mass features and induced mass features
        (from gap-filling) when available to provide complete cluster statistics.
        
        Must be run after add_consensus_mass_features() which creates the cluster assignments.
        Results are stored in cluster_summary_dataframe property and used by plotting methods.
        
        Parameters
        ----------
        None
            Operates on self.mass_features_dataframe and self.induced_mass_features_dataframe.
            Both must contain 'cluster' column.
            
        Returns
        -------
        :obj:`~pandas.DataFrame` or None
            DataFrame with one row per cluster containing summary statistics:
            - cluster: cluster ID
            - mz_{median,mean,std,max,min}: m/z statistics
            - scan_time_aligned_{median,mean,std,max,min}: aligned RT statistics
            - half_height_width_{median,mean,std,max,min}: peak width statistics
            - tailing_factor_{median,mean,std,max,min}: peak shape statistics
            - dispersity_index_{median,mean,std,max,min}: peak quality statistics
            - sample_id_nunique: number of unique samples containing the cluster
            - intensity_{max,median,mean,std,min}: intensity statistics
            - persistence_{max,median,mean,std,min}: persistence statistics
            
            Returns None if mass_features_dataframe is empty.
            
        Notes
        -----
        - Summary DataFrame is automatically stored in cluster_summary_dataframe property
        - Includes both regular and induced (gap-filled) mass features when available
        - Used by plotting methods: plot_consensus_mz_features, plot_mz_features_per_cluster
        - Sample count (sample_id_nunique) indicates cluster prevalence across samples
        - Filters applied by consensus_min_sample_fraction affect which clusters appear
        
        See Also
        --------
        add_consensus_mass_features : Creates clusters before summarization
        fill_missing_cluster_features : Creates induced mass features via gap-filling
        plot_consensus_mz_features : Visualizes cluster summaries
        plot_mz_features_per_cluster : Shows cluster size distribution
        """
        # First check if there are minimum columns in the features dataframe
        if len(self.mass_features_dataframe.columns) < 1:
            return None

        # Combine regular and induced mass features
        mf_df = self.mass_features_dataframe.copy()
        mf_df = mf_df.reset_index(drop=False)
        
        # Check if induced mass features are available and combine them
        if self.induced_mass_features_dataframe is not None and len(self.induced_mass_features_dataframe) > 0:
            imf_df = self.induced_mass_features_dataframe.copy()
            imf_df = imf_df.reset_index(drop=False)
            # Cluster column extracted from mf_id in _prepare_lcms_mass_features_for_combination
            # Combine regular and induced features
            mf_df = pd.concat([mf_df, imf_df], axis=0)
            mf_df = mf_df.reset_index(drop=True)
        
        # Filter out any rows with NaN cluster values before converting to int
        if 'cluster' in mf_df.columns:
            mf_df = mf_df.dropna(subset=['cluster'])
            mf_df['cluster'] = mf_df['cluster'].astype(int)

        # Build aggregation dictionary based on available columns
        agg_dict = {
            "mz": ["median", "mean", "std", "max", "min"],
            "scan_time_aligned": ["median", "mean", "std", "max", "min"],
            "sample_id": ["nunique"],
            "intensity": ["max", "median", "mean", "std", "min"],
        }
        
        # Add optional columns if they exist
        optional_columns = {
            "half_height_width": ["median", "mean", "std", "max", "min"],
            "tailing_factor": ["median", "mean", "std", "max", "min"],
            "dispersity_index": ["median", "mean", "std", "max", "min"],
            "persistence": ["max", "median", "mean", "std", "min"],
        }
        
        for col, funcs in optional_columns.items():
            if col in mf_df.columns:
                agg_dict[col] = funcs

        summary_df = (
            mf_df.groupby("cluster")
            .agg(agg_dict)
            .reset_index()
        )

        # Fix the column names
        summary_df.columns = [
            "_".join(col).strip()
            for col in summary_df.columns.values
            if col != "cluster"
        ]
        summary_df = summary_df.rename(columns={"cluster_": "cluster"})
        # Set cluster as the index for easy lookup
        summary_df = summary_df.set_index('cluster')
        return summary_df

    def plot_mz_features_per_cluster(self, return_fig = False):
        """
        Plot the number of mass features in a cluster against how many clusters
        contain that number of mass features

        Parameters
        -----------
        return_fig : boolean
            Indicates whether to plot composite feature map (False) or return figure object (True). Defaults to False.

        Returns
        --------
        matplotlib.pyplot.Figure
            A figure displaying the frequency with which clusters contain the given number of m/z features

        Raises
        ------
        Warning
            If consensus features haven't been added to the object yet
        """

        if not hasattr(self, 'cluster_summary_dataframe'):
            raise ValueError(
                'cluster_summary_dataframe is not set, must run add_consensus_mass_features() first'
            )
        else:
            sum_data = self.cluster_summary_dataframe
            fig, ax = plt.subplots()
            sum_data.sample_id_nunique.value_counts().sort_index().plot(ax = ax, kind = 'bar')
            plt.xlabel('Number of mass features in a cluster')
            plt.ylabel('Number of clusters with this many mass features')
            if return_fig:
                plt.close(fig)
                return fig
            else:
                plt.show()
        
    def plot_mz_features_across_samples(self, alpha = 0.75, s = 0.005, return_fig = False):
        """
        Generate Scan Time vs m/z plot of all the mass features across all 
        samples in collection where intensity of color on the plot indicates
        density of mass features, NOT INTENSITY

        Parameters
        -----------
        alpha :  float
            Desired transparency for plotted m/z features.  Defaults to 0.75.
        s : float
            Desired size of plotted m/z features. Defaults to 0.005.
        return_fig : boolean
            Indicates whether to plot composite feature map (False) or return figure object (True). Defaults to False.

        Returns
        --------
        matplotlib.pyplot.Figure
            A figure displaying a scan time vs m/z scatterplot of all the m/z features identified in the collection.
            Parameters alpha (transparency) and s (marker size) allow the user to emphasize the density of features.
            Intensity of features is not represented.
        """
        df = self.mass_features_dataframe.copy()
        fig = plt.figure()
        plt.scatter(
            df.scan_time_aligned,
            df.mz,
            c = 'tab:gray',
            alpha = alpha,
            s = s
        )

        plt.xlabel('Scan time')
        plt.ylabel('m/z')
        plt.ylim(0, np.ceil(np.max(df.mz)))
        plt.xlim(0, np.ceil(np.max(df.scan_time)))
        plt.title('All mass features, all samples')
        
        if return_fig:
            plt.close(fig)
            return fig
        else:
            plt.show()

    def plot_consensus_mz_features(self, xb = 'xb', xt = 'xt', yb = 'yb', yt = 'yt', show_all = True, return_fig = False):
        """
        Generate Scan Time vs m/z plot of the consensus features scaled by size
        with option ('show_all') of leaving the individual m/z features in the figure.

        Parameters
        -----------
        xb :  float
            Desired starting scan time value for the x-axis. Defaults to 0.
        xt : float
            Desired ending scan time for the x-axis. Defaults to the maximum scan time value in the provided data.
        yb :  float
            Desired starting m/z value for the y-axis. Defaults to 0.
        yt : float
            Desired ending m/z for the y-axis. Defaults to the maximum m/z value in the provided data.
        show_all : boolean
            Indicates whether to display all identified m/z features (True) or just the consensus features (False). Defaults to True.
        return_fig : boolean
            Indicates whether to plot composite feature map (False) or return figure object (True). Defaults to False.

        Returns
        --------
        matplotlib.pyplot.Figure
            A scalable figure that overlays the consensus features over all the m/z features identified in the collection.
            Consensus features are scaled by how many m/z features are represented in the consensus. Figure can be scaled by
            inputting desired boundaries on the scan time (xb, xt) and m/z values (yb, yt).
        """
        df = self.cluster_summary_dataframe.copy()
        mfdf = self.mass_features_dataframe.copy()

        fig = plt.figure()
        if show_all:
            plt.scatter(
                mfdf.scan_time_aligned,
                mfdf.mz,
                c = 'tab:gray',
                s = 1
            )

        m = plt.scatter(
            df.scan_time_aligned_median,
            df.mz_median, 
            c = 'tab:orange',
            alpha = 0.7, 
            s = (df.sample_id_nunique**2)/5
        )

        plt.xlabel('Scan time')
        plt.ylabel('m/z')
        
        if xt == 'xt':
            xt = np.ceil(np.max(mfdf.mz))
        if yt == 'yt':
            yt = np.ceil(np.max(mfdf.scan_time))
        if xb == 'xb':
            xb = 0
        if yb == 'yb':
            yb = 0
        plt.ylim(xb, xt)
        plt.xlim(yb, yt)

        kw = dict(
            prop = 'sizes',
            num = max(1, int(len(df.sample_id_nunique.unique())/3)),
            color = 'tab:orange',
            alpha = 0.7,
            func = lambda s: np.sqrt(s*5)
        )

        plt.legend(
            *m.legend_elements(**kw), 
            title = 'Features\nper cluster',
            bbox_to_anchor = (1.01, 0.4, 0.225, 0.5)
        )
        plt.tight_layout()
        plt.title('Consensus Features')

        if return_fig:
            plt.close(fig)
            return fig
        else:
            plt.show()
    
    def plot_cluster(
        self,
        cluster_id,
        to_plot=["EIC", "MS1", "MS2"],
        return_fig=False,
        plot_smoothed_eic=False,
        plot_eic_datapoints=False,
        eic_buffer_time=None,
        label_samples=False,
        molecular_metadata=None,
        spectral_library=None,
    ):
        """
        Plot a consensus mass feature cluster across all samples.
        
        Similar to LCMSMassFeature.plot() but shows EICs from all samples in the cluster,
        highlighting the representative mass feature.
        
        Parameters
        ----------
        cluster_id : int
            The cluster ID to plot
        to_plot : list, optional
            List of strings specifying what to plot: "EIC", "MS1", "MS2", "MS2_mirror".
            Default is ["EIC", "MS1", "MS2"].
        return_fig : bool, optional
            If True, returns the figure object. Default is False.
        plot_smoothed_eic : bool, optional
            If True, plots smoothed EICs. Default is False.
        plot_eic_datapoints : bool, optional
            If True, plots EIC data points. Default is False.
        eic_buffer_time : float, optional
            Time buffer around the peak for EIC plotting (minutes).
            If None, uses parameter setting. Default is None.
        label_samples : bool, optional
            If True, labels each sample in the legend. Default is False.
        molecular_metadata : dict, optional
            Dictionary mapping molecular IDs to MetaboliteMetadata objects.
            Required for MS2_mirror plots. Default is None.
        spectral_library : FlashEntropySearch, optional
            FlashEntropy spectral library containing MS2 spectra.
            Required for MS2_mirror plots to retrieve library spectra. Default is None.
            
        Returns
        -------
        matplotlib.figure.Figure or None
            The figure object if return_fig=True, otherwise None
            
        Raises
        ------
        ValueError
            If cluster_id is not found or if required data is not loaded
        """
        import matplotlib.pyplot as plt
        
        # Get cluster summary for median values
        if cluster_id not in self.cluster_summary_dataframe.index:
            raise ValueError(
                f"Cluster {cluster_id} not found in cluster_summary_dataframe. "
                f"Run add_consensus_mass_features() first."
            )
        
        cluster_summary = self.cluster_summary_dataframe.loc[cluster_id]
        
        # Get representative mass feature info
        rep_info = self.get_most_representative_sample_for_cluster(cluster_id)
        rep_sample_id = rep_info['sample_id']
        rep_mf_id = rep_info['mf_id']
        rep_sample = self[rep_sample_id]
        
        # Check if representative mass feature is loaded
        if rep_mf_id not in rep_sample.mass_features:
            raise ValueError(
                f"Representative mass feature {rep_mf_id} not loaded in sample {rep_sample.sample_name}. "
                f"Run reload_representative_mass_features() or process_consensus_features() first."
            )
        
        rep_mf = rep_sample.mass_features[rep_mf_id]
        
        # Get eic buffer time
        if eic_buffer_time is None:
            eic_buffer_time = self[0].parameters.lc_ms.eic_buffer_time
        
        # Adjust to_plot based on available data
        if rep_mf.mass_spectrum is None:
            to_plot = [x for x in to_plot if x != "MS1"]
        if len(rep_mf.ms2_mass_spectra) == 0:
            to_plot = [x for x in to_plot if x not in ["MS2", "MS2_mirror"]]
        
        # Check if EICs are available
        cluster_mfs = self.mass_features_dataframe[
            self.mass_features_dataframe['cluster'] == cluster_id
        ]
        
        has_eics = False
        # Check regular features
        for _, row in cluster_mfs.iterrows():
            sample_id = int(row['sample_id'])
            sample = self[sample_id]
            if hasattr(sample, 'eics') and sample.eics:
                has_eics = True
                break
        
        # Also check induced features if available
        induced_cluster_mfs = None
        if not has_eics and self.induced_mass_features_dataframe is not None:
            induced_cluster_mfs = self.induced_mass_features_dataframe[
                self.induced_mass_features_dataframe['cluster'] == cluster_id
            ]
            for _, row in induced_cluster_mfs.iterrows():
                sample_id = int(row['sample_id'])
                sample = self[sample_id]
                if hasattr(sample, 'eics') and sample.eics:
                    has_eics = True
                    break
        
        if not has_eics:
            to_plot = [x for x in to_plot if x != "EIC"]
            if len(to_plot) == 0:
                raise ValueError(
                    f"No plottable data available for cluster {cluster_id}. "
                    f"Run process_consensus_features(gather_eics=True, add_ms1=True, add_ms2=True) first."
                )
        
        # Get induced features if not already retrieved
        if induced_cluster_mfs is None and self.induced_mass_features_dataframe is not None:
            induced_cluster_mfs = self.induced_mass_features_dataframe[
                self.induced_mass_features_dataframe['cluster'] == cluster_id
            ]
        
        # Check if MS1 is deconvoluted
        deconvoluted = rep_mf._ms_deconvoluted_idx is not None
        
        # Create figure
        fig, axs = plt.subplots(
            len(to_plot), 1, figsize=(10, len(to_plot) * 4), squeeze=False
        )
        
        fig.suptitle(
            f"Consensus Cluster {cluster_id}: "
            f"m/z = {cluster_summary['mz_median']:.4f} "
            f"(±{cluster_summary['mz_std']:.4f}); "
            f"RT = {cluster_summary['scan_time_aligned_median']:.2f} min "
            f"(±{cluster_summary['scan_time_aligned_std']:.2f}); "
            f"{int(cluster_summary['sample_id_nunique'])} samples"
        )
        
        i = 0
        
        # EIC plot - show all samples using helper method
        if "EIC" in to_plot:
            self._plot_multiple_eics(
                axs[i][0],
                cluster_mfs,
                induced_cluster_mfs,
                rep_sample_id,
                rep_mf_id,
                cluster_summary['scan_time_aligned_median'],
                eic_buffer_time,
                plot_smoothed=plot_smoothed_eic,
                plot_datapoints=plot_eic_datapoints,
                label_samples=label_samples,
                lcms_collection=self
            )
            i += 1
        
        # MS1 plot - from representative using helper method
        if "MS1" in to_plot:
            rep_mf._plot_ms1_spectrum(
                axs[i][0], 
                deconvoluted=deconvoluted, 
                sample_name=rep_sample.sample_name
            )
            i += 1
        
        # MS2 plot - from representative using helper method
        if "MS2" in to_plot:
            rep_mf._plot_ms2_spectrum(axs[i][0], sample_name=rep_sample.sample_name)
            i += 1
        
        # MS2 mirror plot - from representative using helper method
        if "MS2_mirror" in to_plot:
            rep_mf._plot_ms2_mirror(axs[i][0], molecular_metadata=molecular_metadata, spectral_library=spectral_library)
            i += 1
        
        plt.tight_layout()
        
        if return_fig:
            plt.close(fig)
            return fig
        else:
            plt.show()
            return None
    
    def get_representative_mass_features_for_all_clusters(self, representative_metric=None):
        """
        Get the most representative mass feature for all clusters in bulk.
        
        This is much more efficient than calling get_most_representative_sample_for_cluster
        in a loop, as it processes all clusters in a single pass over the dataframe.
        
        Parameters
        ----------
        representative_metric : str, optional
            The metric to use to determine the most representative sample.
            If None, uses the value from self.parameters.lcms_collection.consensus_representative_metric.
            Options:
            - 'intensity': Selects the mass feature with the highest intensity
            - 'intensity_prefer_ms2': Selects the highest intensity feature that has MS2 scans,
              or the highest intensity overall if none have MS2
            Default is None (uses parameter setting).
            
        Returns
        -------
        :obj:`~pandas.DataFrame`
            DataFrame with one row per cluster containing:
            - cluster: cluster ID
            - sample_id: The sample ID of the most representative sample
            - mf_id: The mass feature ID in the sample
            - coll_mf_id: The collection-level mass feature ID (index)
            - has_ms2: Whether this mass feature has MS2 scan numbers
            - intensity: The intensity value of the representative mass feature
        """
        # Use default from parameters if not specified
        if representative_metric is None:
            representative_metric = self.parameters.lcms_collection.consensus_representative_metric
        
        mf_df = self.mass_features_dataframe.copy()
        # Reset index to make coll_mf_id a column we can work with
        mf_df = mf_df.reset_index(drop=False)
        
        # Handle special metric 'intensity_prefer_ms2'
        if representative_metric == 'intensity_prefer_ms2':
            if 'intensity' not in mf_df.columns:
                raise ValueError(
                    f"'intensity' column not found in mass_features_dataframe. "
                    f"Available columns: {mf_df.columns.tolist()}"
                )
            
            # Add has_ms2 flag if ms2_scan_numbers column exists
            if 'ms2_scan_numbers' in mf_df.columns:
                def has_ms2_scans(val):
                    if val is None:
                        return False
                    try:
                        return len(val) > 0
                    except (TypeError, ValueError):
                        return False
                
                mf_df['has_ms2'] = mf_df['ms2_scan_numbers'].apply(has_ms2_scans)
                
                # Sort by has_ms2 (descending) then intensity (descending)
                # This ensures features with MS2 are preferred when intensities are equal
                mf_df = mf_df.sort_values(['has_ms2', 'intensity'], ascending=[False, False])
            else:
                mf_df['has_ms2'] = False
                mf_df = mf_df.sort_values('intensity', ascending=False)
            
            # Group by cluster and take the first (highest intensity, preferring MS2)
            representatives = mf_df.groupby('cluster').first().reset_index()
            
        else:
            # Standard metric - check if it exists
            if representative_metric not in mf_df.columns:
                raise ValueError(
                    f"Metric '{representative_metric}' not found. Available columns: {mf_df.columns.tolist()}"
                )
            
            # Add has_ms2 flag for consistency
            if 'ms2_scan_numbers' in mf_df.columns:
                def has_ms2_scans(val):
                    if val is None:
                        return False
                    try:
                        return len(val) > 0
                    except (TypeError, ValueError):
                        return False
                mf_df['has_ms2'] = mf_df['ms2_scan_numbers'].apply(has_ms2_scans)
            else:
                mf_df['has_ms2'] = False
            
            # Get the index of max value for each cluster
            idx = mf_df.groupby('cluster')[representative_metric].idxmax()
            representatives = mf_df.loc[idx].copy()
        
        # Select only the columns we need
        result_cols = ['cluster', 'sample_id', 'mf_id', 'coll_mf_id', 'has_ms2', 'intensity']
        representatives = representatives[result_cols]
        
        return representatives
    
    def get_sample_mf_map_for_representatives(self, representative_metric=None, include_cluster_id=True):
        """
        Build a mapping of sample_id -> list of representative mass feature IDs to load.
        
        This is a DRY helper method used by both process_consensus_features() and
        ReadSavedLCMSCollection to determine which mass features should be loaded
        for each sample when loading representatives.
        
        Parameters
        ----------
        representative_metric : str, optional
            The metric to use to determine the most representative sample.
            If None, uses the value from self.parameters.lcms_collection.consensus_representative_metric.
            Default is None.
        include_cluster_id : bool, optional
            If True, returns tuples of (mf_id, cluster_id). If False, returns just mf_id.
            Default is True.
        
        Returns
        -------
        dict
            Dictionary mapping sample_id (int) to list of mass feature identifiers.
            If include_cluster_id=True: list of tuples (mf_id, cluster_id)
            If include_cluster_id=False: list of mf_id integers
        
        Examples
        --------
        >>> # Get map with cluster IDs for loading
        >>> sample_mf_map = collection.get_sample_mf_map_for_representatives()
        >>> # sample_mf_map = {0: [(123, 0), (456, 1)], 1: [(789, 2)], ...}
        >>> 
        >>> # Get map without cluster IDs for pipeline
        >>> sample_mf_map = collection.get_sample_mf_map_for_representatives(include_cluster_id=False)
        >>> # sample_mf_map = {0: [123, 456], 1: [789], ...}
        """
        # Get all representative mass features in bulk (much faster than looping)
        representatives = self.get_representative_mass_features_for_all_clusters(
            representative_metric=representative_metric
        )
        
        # Build sample_mf_map
        sample_mf_map = {}
        for _, row in representatives.iterrows():
            sample_id = row['sample_id']
            mf_id = row['mf_id']
            cluster_id = row['cluster']
            
            if sample_id not in sample_mf_map:
                sample_mf_map[sample_id] = []
            
            if include_cluster_id:
                sample_mf_map[sample_id].append((mf_id, cluster_id))
            else:
                sample_mf_map[sample_id].append(mf_id)
        
        return sample_mf_map
    
    def get_most_representative_sample_for_cluster(self, cluster_id, representative_metric=None):
        """
        Get the most representative sample for a given cluster based on a metric.
        
        Parameters
        ----------
        cluster_id : int
            The cluster ID to find the representative sample for.
        representative_metric : str, optional
            The metric to use to determine the most representative sample.
            If None, uses the value from self.parameters.lcms_collection.consensus_representative_metric.
            Options:
            - 'intensity': Selects the mass feature with the highest intensity
            - 'intensity_prefer_ms2': Selects the highest intensity feature that has MS2 scans,
              or the highest intensity overall if none have MS2
            Default is None (uses parameter setting).
            
        Returns
        -------
        dict
            Dictionary containing:
            - 'sample_id': The sample ID of the most representative sample
            - 'sample_name': The sample name of the most representative sample
            - 'mf_id': The mass feature ID in the sample
            - 'coll_mf_id': The collection-level mass feature ID (index)
            - 'has_ms2': Whether this mass feature has MS2 scan numbers
            - 'intensity': The intensity value of the representative mass feature
        
        Raises
        ------
        ValueError
            If cluster_id is not found or if representative_metric is not a valid column.
        """
        # Use the bulk method to get all representatives, then filter to this cluster
        # This follows DRY principle and ensures consistency
        all_representatives = self.get_representative_mass_features_for_all_clusters(
            representative_metric=representative_metric
        )
        
        # Filter to the requested cluster
        cluster_rep = all_representatives[all_representatives['cluster'] == cluster_id]
        
        if len(cluster_rep) == 0:
            # Try to provide helpful error message
            available_clusters = self.mass_features_dataframe['cluster'].unique()
            raise ValueError(
                f"Cluster {cluster_id} not found in mass_features_dataframe. "
                f"Available clusters: {sorted(available_clusters[:10].tolist())}... "
                f"(showing first 10 of {len(available_clusters)} total clusters)"
            )
        
        # Get the representative row (should only be one)
        rep_row = cluster_rep.iloc[0]
        
        # Get sample name from sample_id (convert to int for list indexing)
        sample_id = int(rep_row['sample_id'])
        sample_name = self.samples[sample_id]
        
        return {
            'sample_id': sample_id,
            'sample_name': sample_name,
            'mf_id': rep_row['mf_id'],
            'coll_mf_id': rep_row['coll_mf_id'],
            'has_ms2': rep_row['has_ms2'],
            'intensity': rep_row['intensity']
        }
    
    def reload_representative_mass_features(self, add_ms2=False, auto_process_ms2=True, ms2_spectrum_mode=None, ms2_scan_filter=None):
        """
        Reload mass features for all representative samples in the cluster summary.
        
        This method is useful when the collection was loaded with load_light=True,
        which stores mass features only in the collection dataframe. This reloads
        the specific mass features that are representatives for each cluster,
        allowing them to be accessed as LCMSMassFeature objects.
        
        Parameters
        ----------
        add_ms2 : bool, optional
            If True, also loads and associates MS2 spectra with mass features. Default is False.
        auto_process_ms2 : bool, optional
            If True and add_ms2=True, auto-processes MS2 spectra. Default is True.
        ms2_spectrum_mode : str or None, optional
            Spectrum mode for MS2 spectra. If None, determines from parser. Default is None.
        ms2_scan_filter : str or None, optional
            Filter string for MS2 scans (e.g., 'hcd'). Default is None.
        
        Returns
        -------
        dict
            Dictionary mapping sample_id to list of reloaded mf_ids.
            
        Raises
        ------
        ValueError
            If cluster_summary_dataframe is not set (run add_consensus_mass_features first).
            
        Notes
        -----
        - Only reloads mass features that are cluster representatives
        - Uses get_most_representative_sample_for_cluster() to determine which to reload
        - More memory-efficient than reloading all mass features
        - Parallelized based on lcms_collection.cores parameter
        - MS2 association uses same logic as add_associated_ms2_dda()
        
        See Also
        --------
        _reload_sample_mass_features : Low-level method to reload specific mass features
        get_most_representative_sample_for_cluster : Gets representative sample for cluster
        """
        # Validate prerequisites
        if not hasattr(self, 'cluster_summary_dataframe') or self.cluster_summary_dataframe is None:
            raise ValueError(
                "cluster_summary_dataframe not found. Must run add_consensus_mass_features() first."
            )
        
        # Get all representative mass features in bulk (much faster than looping)
        representatives = self.get_representative_mass_features_for_all_clusters()
        
        # Build a dictionary of sample_id -> list of mf_ids that are representatives
        sample_mf_map = {}
        for _, row in representatives.iterrows():
            sample_id = row['sample_id']
            mf_id = row['mf_id']
            
            if sample_id not in sample_mf_map:
                sample_mf_map[sample_id] = []
            sample_mf_map[sample_id].append(mf_id)
        
        # Reload mass features for each sample (parallelized)
        if self.parameters.lcms_collection.cores == 1:
            # Serial processing
            from tqdm import tqdm
            for sample_id in tqdm(sample_mf_map.keys(), desc="Reloading representative mass features", unit="sample"):
                mf_ids = sample_mf_map[sample_id]
                self._reload_sample_mass_features(sample_id, mf_ids_to_load=mf_ids, add_ms2=add_ms2, 
                                                  auto_process_ms2=auto_process_ms2, ms2_spectrum_mode=ms2_spectrum_mode,
                                                  ms2_scan_filter=ms2_scan_filter)
        else:
            # Parallel processing
            import multiprocessing
            from tqdm import tqdm
            
            if self.parameters.lcms_collection.cores > len(sample_mf_map):
                ncores = len(sample_mf_map)
            else:
                ncores = self.parameters.lcms_collection.cores
            
            pool = multiprocessing.Pool(ncores)
            
            # Build arguments list for starmap
            args_list = [
                (sample_id, sample_mf_map[sample_id], add_ms2, auto_process_ms2, 
                 ms2_spectrum_mode, ms2_scan_filter, False)
                for sample_id in sample_mf_map.keys()
            ]
            
            # Execute in parallel
            mp_result = pool.starmap(self._reload_sample_mass_features, args_list)
            pool.close()
            pool.join()
            
            # Collect results back into samples
            for i, sample_id in enumerate(tqdm(sample_mf_map.keys(), desc="Collecting reloaded mass features", unit="sample")):
                self[sample_id].mass_features = mp_result[i]
        
        return sample_mf_map
    
    def _associate_ms2_with_mass_features(self, sample, local_mf_ids, auto_process=True, 
                                          spectrum_mode=None, scan_filter=None):
        """
        Associate MS2 spectra with specific mass features in a sample.
        
        Uses the LCMSBase helper method to find and load MS2 scans for the specified mass features.
        
        Parameters
        ----------
        sample : LCMSBase
            The sample object containing mass features and scan data.
        local_mf_ids : list of int
            List of local (sample-level) mass feature IDs to find MS2 for.
        auto_process : bool, optional
            If True, auto-processes the MS2 spectra. Default is True.
        spectrum_mode : str or None, optional
            Spectrum mode for MS2 spectra. If None, determines from parser. Default is None.
        scan_filter : str or None, optional
            Filter string for MS2 scans (e.g., 'hcd'). Default is None.
            
        Returns
        -------
        dict
            Dictionary of scan_number -> MassSpectrum objects for the loaded MS2 spectra.
        """
        # Check if we have scan data
        if not hasattr(sample, 'scan_df') or sample.scan_df is None:
            return {}
        
        # Separate mass features into those that need scan finding vs those that already have scans
        mfs_needing_scan_finding = []
        unique_dda_scans = set()
        
        for mf_id in local_mf_ids:
            if mf_id not in sample.mass_features:
                continue
            mf = sample.mass_features[mf_id]
            # If this mass feature already has MS2 scans, add them to our set
            if mf.ms2_scan_numbers is not None and len(mf.ms2_scan_numbers) > 0:
                # Convert to integers in case they come from HDF5 as numpy types
                unique_dda_scans.update([int(scan) for scan in mf.ms2_scan_numbers])
            else:
                # Otherwise, we need to find scans for this mass feature
                mfs_needing_scan_finding.append(mf_id)
        
        # Only run the scan finding for mass features that need it
        if mfs_needing_scan_finding:
            found_scans = sample._find_ms2_scans_for_mass_features(
                mf_ids=mfs_needing_scan_finding,
                scan_filter=scan_filter
            )
            unique_dda_scans.update(found_scans)

        if len(unique_dda_scans) == 0:
            return {}
        
        # Get ms2 parameters from sample
        #TODO KRH: deal with different ms2 scan types here (CID vs HCD), may need to add scan translator to the initializeion
        ms_params = sample.parameters.mass_spectrum['ms2']

        # Load MS2 spectra (convert set to list)
        sample.add_mass_spectra(
            scan_list=list(unique_dda_scans),
            auto_process=auto_process,
            spectrum_mode=spectrum_mode,
            ms_level=2,
            use_parser=True,
            ms_params=ms_params,
        )
        
        # Associate MS2 spectra with mass features
        for mf_id in local_mf_ids:
            if mf_id not in sample.mass_features:
                continue
            if sample.mass_features[mf_id].ms2_scan_numbers is not None and len(sample.mass_features[mf_id].ms2_scan_numbers) > 0:
                for dda_scan in sample.mass_features[mf_id].ms2_scan_numbers:
                    if dda_scan in sample._ms:
                        sample.mass_features[mf_id].ms2_mass_spectra[dda_scan] = sample._ms[dda_scan]
        
        # Return only the MS2 spectra we loaded (for parallel processing)
        return {scan: sample._ms[scan] for scan in unique_dda_scans if scan in sample._ms}
    
    def _reload_sample_mass_features(self, sample_id, mf_ids_to_load=None, add_ms2=False, 
                                     auto_process_ms2=True, ms2_spectrum_mode=None, ms2_scan_filter=None,
                                     inplace=True):
        """
        Reload specific mass features for a sample from HDF5.
        
        This is useful when the collection was loaded with load_light=True,
        which stores mass features only in the collection dataframe and not
        as LCMSMassFeature objects in individual samples.
        
        Parameters
        ----------
        sample_id : int
            The sample ID to reload mass features for.
        mf_ids_to_load : list of str, optional
            List of collection-level mf_ids (format: '{sample_id}_{local_mf_id}') to load.
            If None, loads all mass features for the sample.
        add_ms2 : bool, optional
            If True, also loads and associates MS2 spectra. Default is False.
        auto_process_ms2 : bool, optional
            If True, auto-processes MS2 spectra. Default is True.
        ms2_spectrum_mode : str or None, optional
            Spectrum mode for MS2 spectra. Default is None.
        ms2_scan_filter : str or None, optional
            Filter string for MS2 scans. Default is None.
        inplace : bool, optional
            If True, updates the sample's mass_features in place. If False, returns the
            mass_features dictionary (for multiprocessing). Default is True.
            
        Returns
        -------
        dict or None
            If inplace=False, returns dictionary of mass features.
            Otherwise returns None and updates object in place.
        """
        sample = self[sample_id]
        sample_name = self.samples[sample_id]
        
        # Check if we have a collection parser that can reload
        if not hasattr(self, 'collection_parser') or self.collection_parser is None:
            print("Warning: Cannot reload mass features - no collection_parser available")
            if not inplace:
                return {}
            return
        
        # Get the HDF5 file for this sample
        hdf5_file = self.collection_parser.folder_location / f"{sample_name}.corems/{sample_name}.hdf5"
        
        if not hdf5_file.exists():
            print(f"Warning: HDF5 file not found for sample {sample_name}: {hdf5_file}")
            if not inplace:
                return {}
            return
        
        # Import here to avoid circular imports
        from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
        
        # If specific mf_ids requested, use them directly
        local_mf_ids_to_load = None
        if mf_ids_to_load is not None:
            # mf_ids_to_load is already a list of sample-level mf_ids (integers)
            # No parsing needed - they come from the mf_id column in the dataframe
            local_mf_ids_to_load = set(mf_ids_to_load)
        
        # Reload mass features from HDF5
        with ReadCoreMSHDFMassSpectra(hdf5_file) as parser:
            # Load mass features - if specific IDs requested, only load those
            parser.import_mass_features(sample, mf_ids=local_mf_ids_to_load)
        
        # If add_ms2, associate MS2 spectra with the loaded mass features
        if add_ms2 and local_mf_ids_to_load is not None:
            self._associate_ms2_with_mass_features(
                sample, 
                list(local_mf_ids_to_load),
                auto_process=auto_process_ms2,
                spectrum_mode=ms2_spectrum_mode,
                scan_filter=ms2_scan_filter
            )
        
        # Return mass features if not inplace (for multiprocessing)
        if not inplace:
            return sample.mass_features
        
    def add_sparse_distance_matrix(self, features):
        if features is None:
            return None
        else:
            features = features.copy()

        # Parameters for calculating distance between features
        dims = ["mz", "scan_time_aligned"]
        relative = [True, False]
        mz_tol_relative = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
        tol = [mz_tol_relative, self.parameters.lcms_collection.consensus_rt_tol]
        dist_weight = [1, 1]

        # Check that the dimensions and tolerances are the same length
        if (
            len(dims) != len(tol)
            or len(dims) != len(relative)
            or len(dims) != len(dist_weight)
        ):
            raise ValueError(
                "The dimensions, tolerances, relative, dist_weight, and na_allow lists must be the same length"
            )

        # Make connectivity matrix for masking within sample mass features
        ## Masking matrix cmat will mark all features from the same sample as 0
        ## To mask, a matrix can be multiplied by cmat and features from same
        ## samples are multiplied by 0, while features from different samples 
        ## are multiplied by 1

        if "sample_id" not in features.columns:
            cmat = None
        else:
            vals = features["sample_id"].values.reshape(-1, 1)
            cmat = scipy.spatial.distance.cdist(vals, vals)
            # Convert to binary (0 if same sample, 1 if different)
            cmat = np.where(cmat == 0, 0, 1)
            # Convert to coordinate matrix for sparse operations later
            cmat = sparse.coo_matrix(cmat)

        # Compute inter-feature distances using sparse matrix approach
        distances = None # clear the distances object before starting
        for i in range(len(dims)): # iterate through all dimensions to be considered
            # Construct k-d tree
            values = features[dims[i]].values

            tree = KDTree(values.reshape(-1, 1))

            max_tol = tol[i]
            if relative[i] is True:
                # Maximum absolute tolerance
                max_tol = tol[i] * values.max()

            # Compute sparse distance matrix
            # the larger the max_tol, the slower this operation is
            sdm = tree.sparse_distance_matrix(tree, max_tol, output_type="coo_matrix")

            # Only consider forward case, exclude diagonal
            sdm = sparse.triu(sdm, k=1)

            # Filter relative distances
            if relative[i] is True:
                # Compute relative distances
                rel_dists = sdm.data / values[sdm.row]

                # Indices of relative distances less than tolerance
                idx = rel_dists <= tol[i]

                # Reconstruct sparse distance matrix
                sdm = sparse.coo_matrix(
                    (rel_dists[idx], (sdm.row[idx], sdm.col[idx])),
                    shape=(len(values), len(values)),
                )

            # Scaled distances wrt the maximum tolerance for the dimension
            sdm.data = sdm.data / tol[i]

            # Stack distances for dimensions where na_allow is False
            if distances is None:
                sdm.data = sdm.data * dist_weight[i]
                # Replace zeros with epsilon to handle perfect matches
                sdm.data[sdm.data == 0] = 1e-10
                distances = sdm
            else:
                # Prepare sdm to match shape of existing distances
                distances_truth = distances.copy()
                # make new sparse matrix with same positions as previous 
                # distance matrix but all ones for values
                distances_truth.data = np.ones_like(distances_truth.data)
                
                # Replace zeros with epsilon BEFORE multiply to prevent sparse matrix from dropping them
                sdm.data[sdm.data == 0] = 1e-10
                
                # multiply the new sparse matrix (sdm) by this mask to remove 
                # data that doesn't exist in original sparse matrix
                sdm = distances_truth.multiply(sdm)
                
                sdm.data = sdm.data * dist_weight[i]
                # Replace zeros with epsilon to handle perfect matches
                sdm.data[sdm.data == 0] = 1e-10

                # use same process as before to remove data from previous
                # distances matrix that isn't in new distances matrix
                sdm_truth = sdm.copy()
                sdm_truth.data = np.ones_like(sdm_truth.data)

                # remove the distances that are not sdm
                distances = distances.multiply(sdm_truth)

                # Sum the new distances
                distances = distances + sdm

        # Multiply by connectivity matrix for more masking
        distances = distances.multiply(cmat)

        # Set attribute holding distance matrix
        self._sparse_distance_matrix = distances

    def evaluate_clusters_for_repeats(self, features):
        raise NotImplementedError('evaluate_clusters_for_repeats not implemented yet')
        summary_df = self.cluster_summary_dataframe.copy()

        # Arrange by decreasing median intensity
        summary_df = summary_df.sort_values(
            by="intensity_median", ascending=False
        ).reset_index(drop=True)

        # Find clusters that are within the mz_tol and rt_tol of each other (on the medians)
        # Create a distance matrix
        # Define how to calculate the distance between features
        dims = ["mz_median", "scan_time_aligned_median"]
        relative = [True, False]
        mz_tol_relative = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
        tol = [mz_tol_relative, self.parameters.lcms_collection.consensus_rt_tol]

        # Compute inter-feature distances
        distances = None
        for i in range(len(dims)):
            # Construct k-d tree
            values = summary_df[dims[i]].values
            tree = KDTree(values.reshape(-1, 1))

            max_tol = tol[i]
            if relative[i] is True:
                # Maximum absolute tolerance
                max_tol = tol[i] * values.max()

            # Compute sparse distance matrix
            # the larger the max_tol, the slower this operation is
            sdm = tree.sparse_distance_matrix(tree, max_tol, output_type="coo_matrix")

            # Only consider forward case, exclude diagonal
            sdm = sparse.triu(sdm, k=1)

            # Filter relative distances
            if relative[i] is True:
                # Compute relative distances
                rel_dists = sdm.data / values[sdm.row]  # or col?

                # Indices of relative distances less than tolerance
                idx = rel_dists <= tol[i]

                # Reconstruct sparse distance matrix
                sdm = sparse.coo_matrix(
                    (rel_dists[idx], (sdm.row[idx], sdm.col[idx])),
                    shape=(len(values), len(values)),
                )

            # Cast as binary matrix
            sdm.data = np.ones_like(sdm.data)

            # Stack distances
            if distances is None:
                distances = sdm
            else:
                distances = distances.multiply(sdm)

        # Roll up features
        # Extract indices of within-tolerance points
        distances = distances.tocoo()
        pairs = np.stack(
            (distances.row, distances.col), axis=1
        )  # These are the index values of the clusters, not the cluster ids
        # Conver to cluster ids
        pairs_df = pd.DataFrame(pairs, columns=["parent", "child"])
        pairs_df["parent"] = summary_df.loc[pairs[:, 0]]["cluster"].values
        pairs_df["child"] = summary_df.loc[pairs[:, 1]]["cluster"].values
        pairs_df = pairs_df.set_index("parent")

        merge_these_clusters = []
        possible_overlaps = []
        root_parents = np.setdiff1d(
            np.unique(pairs_df.index.values), np.unique(pairs_df.child.values)
        )
        for parent in root_parents:
            parent_features = features[features["cluster"] == parent]
            children = pairs_df.loc[[parent], "child"].tolist()
            for child in children:
                overlap = self.check_merge(parent_features, child, features)
                if len(overlap) == 0:
                    merge_these_clusters.append((parent, child, len(overlap)))
                else:
                    possible_overlaps.append((parent, child, len(overlap)))

        result_dict = {}
        result_dict["merge_these_clusters"] = merge_these_clusters
        result_dict["possible_overlaps"] = possible_overlaps

        return result_dict

    def check_merge(self, parent_features, child, features):
        # Grab the features of the parent and children
        child_features = features[features["cluster"] == child]

        # Check if there is an overlap between mf_coll_id in the parent and child clusters
        overlap = np.intersect1d(
            parent_features["sample_id"].values, child_features["sample_id"].values
        )

        return overlap

    def cluster_mass_features_agg_cluster(self, features):
        if features is None:
            return None

        features = features.copy()

        self.add_sparse_distance_matrix(features)

        distances = self._sparse_distance_matrix

        # Convert to full matrix
        distances = distances.todense()
        
        # Cast all 0s to 1s for a distance matrix
        distances[distances == 0] = 1
        distances = np.asarray(distances)

        # Perform clustering
        try:
            clustering = AgglomerativeClustering(
                n_clusters=None,
                linkage="complete",
                # using complete linkage will prevent one sample from being assigned to multiple clusters
                metric="precomputed",
                distance_threshold=1,
            ).fit(distances)
            features["cluster"] = clustering.labels_

        # All data points are singleton clusters
        except:
            features["cluster"] = np.arange(len(features.index))

        return features

    def cluster_inspection_plot(self, clu, return_fig = False):        
        """
        Generate Scan Time vs m/z plot for a narrow range around a given 
        cluster. This tool is meant to support the user in fine tuning the
        tolerances used for the clustering algorithm. The user-provided cluster
        ID is highlighted in larger, magenta marker and the ten largest of the
        remaining clusters are idenfitied with different colors while the
        smallest clusters are light gray.

        Parameters
        -----------
        clu :  integer
            A cluster ID that exists in self.mass_features_dataframe
        return_fig : boolean
            Indicates whether to plot cluster inspection figure (False) or 
            return figure object (True). Defaults to False.

        Returns
        --------
        matplotlib.pyplot.Figure
            A figure displaying a scan time vs m/z scatterplot of small region
            around a given cluster with the ten largest clusters in the region
            distinctly identified

        Raises
        ------
        Warning
            If cluster data haven't been added to the object yet
        """

        if 'cluster' not in self.mass_features_dataframe.columns:
            raise ValueError(
            'Cluster information is not yet added to mass_features_dataframe, must run add_consensus_mass_features() first'
            )
        
        else:
            mztol = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
            rttol = self.parameters.lcms_collection.consensus_rt_tol
            clu_features = self.mass_features_dataframe.copy()

            inclu = clu_features[clu_features.cluster == clu]
            exclu = clu_features[clu_features.cluster != clu]

            dt_ymin = np.floor(min(inclu.mz)) - 1
            dt_ymax = np.ceil(max(inclu.mz)) + 1
            dt_xmin = np.floor(min(inclu.scan_time_aligned)) - 1
            dt_xmax = np.ceil(max(inclu.scan_time_aligned)) + 1

            exclu = exclu[
                (
                    exclu.mz.between(dt_ymin, dt_ymax, inclusive = 'both')
                ) & (
                    exclu.scan_time_aligned.between(dt_xmin, dt_xmax, inclusive = 'both')
                )
            ]

            bigclulist = list(exclu.cluster.value_counts()[:10].index)
            bigclu = exclu[exclu.cluster.isin(bigclulist)]
            smclu = exclu[~exclu.cluster.isin(bigclulist)]

            colors = np.arange(0, 10)
            colordict = dict(zip(bigclulist, colors))
            bigclu['color'] = bigclu.cluster.apply(lambda x: colordict[x])

            fig = plt.figure(figsize = (7.5, 5))

            plt.scatter(
                inclu.scan_time_aligned,
                inclu.mz,
                c = 'm',
                s = 3,
                label = 'Cluster ' + str(clu)
            )

            plt.scatter(
                bigclu.scan_time_aligned,
                bigclu.mz,
                c = bigclu.color,
                cmap = 'tab10',
                s = 1.5
            )

            plt.scatter(
                smclu.scan_time_aligned,
                smclu.mz,
                c = 'silver',
                s = 2,
                label = 'Small clusters'
            )

            plt.ylim(dt_ymin, dt_ymax)
            plt.xlim(dt_xmin, dt_xmax)
            plt.legend(ncol = 2, bbox_to_anchor = (0.8, -0.1))
            plt.xlabel('Scan time')
            plt.ylabel('m/z')
            title_str = 'Cluster ' + str(clu)
            title_str += ': representing ' + str(len(inclu.sample_id.unique())) 
            title_str += ' of ' + str(len(clu_features.sample_id.unique())) 
            title_str += ' samples\n'
            title_str += 'M/Z tolerance: ' + str(mztol) + '\n'
            title_str += 'Scan Time tolerance: ' + str(rttol)
            plt.title(title_str, fontsize = 10)

            if return_fig:
                plt.close(fig)
                return fig
            else:
                plt.show()

    def plot_cluster_outlier_frequency(self, dim_list = ['mz', 'scan_time_aligned'], clu_size_thresh = 0.5, return_fig = False):
        """
        Generate histogram showing the frequency of outlier occurrences by
        clustering dimension across all clusters

        Parameters
        -----------
        dim_list :  list
            List of strings describing dimensions that can be used in 
            clustering. Available list items:
                - 'mz'
                - 'scan_time_aligned'
                - 'half_height_width'
                - 'tailing_factor'
                - 'dispersity_index'
                - 'intensity'
                - 'persistence'
        clu_size_thresh : float
            Value between 0 and 1 that indicates what percentage of samples 
            need to be present in a cluster before it's evaluated for outliers.
            Defaults to 0.5.
        return_fig : boolean
            Indicates whether to plot cluster inspection figure (False) or 
            return figure object (True). Defaults to False.

        Returns
        --------
        matplotlib.pyplot.Figure
            A figure displaying the frequency of outlier occurrences across all
            clusters in the provided measurement dimensions

        Raises
        ------
        Warning
            If cluster data haven't been added to the object yet
        """

        if not hasattr(self, 'cluster_summary_dataframe'):
            raise ValueError(
                'cluster_summary_dataframe is not yet added, must run add_consensus_mass_features() first'
            )

        mfdf = self.mass_features_dataframe.copy()
        summarydf = self.cluster_summary_dataframe

        numsamples = len(self)
        sumdf = summarydf[summarydf.sample_id_nunique > numsamples * clu_size_thresh].reset_index(drop = True).copy()

        ## find the ranges for non-outlier values and add them to sumdf
        mergelist = ['cluster']
        for dim in dim_list:
            maxtag = dim + '_outmax'
            mintag = dim + '_outmin'
            mergelist.append(maxtag)
            mergelist.append(mintag)
            # Calculate outlier thresholds using vectorized operations
            sumdf[mintag] = sumdf[dim + '_mean'] - 3*sumdf[dim + '_std']
            sumdf[maxtag] = sumdf[dim + '_mean'] + 3*sumdf[dim + '_std']
            ## If NaN shows up anywhere in dim_min, dim_max calculations, value is set to NaN and it's 
            ## not flagged. This happens when there's not enough values to compute median/std for that 
            ## dimension therefore can't have outliers

        ## add ranges to mfdf and identify mass features that fall outside the ranges
        # Merge without dropping NaN - we'll handle it per-dimension
        outdf = pd.merge(mfdf, sumdf[mergelist], on = 'cluster')

        outtags = ['cluster']
        for dim in dim_list:
            dimtag = dim + '_outlier'
            outtags.append(dimtag)
            maxtag = dim + '_outmax'
            mintag = dim + '_outmin'
            # Only flag as outlier if thresholds are valid (not NaN)
            outdf[dimtag] = np.where(
                (outdf[maxtag].notna() & outdf[mintag].notna()) &
                (((outdf[dim] > outdf[maxtag])) | ((outdf[dim] < outdf[mintag]))), 
                True, 
                False
            )

        ## identify number of outliers in each cluster
        outliers = outdf[outtags]
        outliers = outliers.groupby(['cluster']).sum()

        ## plot number of clusters that contain any outliers
        fig = plt.figure()
        plt.bar(dim_list, outliers.sum().values, width = 0.5)
        plt.xticks(rotation = 90)
        plt.title('Frequency of outliers across all clusters by category')
        
        if return_fig:
            plt.close(fig)
            return fig
        else:
            plt.show()
            
    def _search_for_targeted_mass_features_in_sample(self, obj_idx, missingdf, cluster_dict, expand_on_miss=False, inplace=True):
        """
        Helper method to search for missing mass features in a single sample.
        
        Internal method called by fill_missing_cluster_features() to perform
        gap-filling for one sample in the collection.
        
        Parameters
        ----------
        obj_idx : int
            Index of the sample being processed
        missingdf : pd.DataFrame
            DataFrame containing cluster information with columns:
            'cluster', 'sample_id_nunique', 'mz_min', 'mz_max', 
            'scan_time_aligned_min', 'scan_time_aligned_max', 'mz_min_allowed', 
            'mz_max_allowed', 'scan_time_aligned_min_allowed', 
            'scan_time_aligned_max_allowed', 'missing_samples'
        cluster_dict : dict
            Pre-computed cluster feature dictionary to avoid recomputation
        expand_on_miss : bool
            If True, expands search window when no peak found initially
        inplace : bool
            If True, assigns induced_mass_features in place. If False, returns the
            induced features dictionary (for multiprocessing)
            
        Returns
        -------
        dict or None
            If inplace=False, returns dictionary of induced mass features.
            Otherwise returns None and updates object in place.
        """
        ## Use the pre-computed cluster dictionary passed as parameter
        
        ## to get clusters missing data based on sample index:
        sampledf = missingdf[
            missingdf.missing_samples.apply(lambda x: obj_idx in x)
        ].reset_index(drop = True).copy()

        # Skip if no missing features for this sample
        if len(sampledf) == 0:
            if not inplace:
                return {}
            return

        self.load_raw_data(obj_idx, 1)
               
        ## this is the line that bugs due to _ms_unprocessed not having key 1
        ms1df = self[obj_idx]._ms_unprocessed[1].copy()
        scan_df = self[obj_idx].scan_df[['scan', 'scan_time_aligned']]
        ms1df = pd.merge(ms1df, scan_df, on = 'scan')

        # Pre-extract all values from sampledf to avoid repeated .iloc calls
        clusters = sampledf.cluster.values
        mz_mins = sampledf.mz_min.values
        mz_maxs = sampledf.mz_max.values
        st_mins = sampledf.scan_time_aligned_min.values
        st_maxs = sampledf.scan_time_aligned_max.values
        
        if expand_on_miss:
            mz_mins_allowed = sampledf.mz_min_allowed.values
            mz_maxs_allowed = sampledf.mz_max_allowed.values
            st_mins_allowed = sampledf.sta_min_allowed.values
            st_maxs_allowed = sampledf.sta_max_allowed.values

        # Pre-filter ms1df to reduce search space
        mz_global_min = mz_mins.min()
        mz_global_max = mz_maxs.max()
        st_global_min = st_mins.min()
        st_global_max = st_maxs.max()
        
        if expand_on_miss:
            mz_global_min = min(mz_global_min, mz_mins_allowed.min())
            mz_global_max = max(mz_global_max, mz_maxs_allowed.max())
            st_global_min = min(st_global_min, st_mins_allowed.min())
            st_global_max = max(st_global_max, st_maxs_allowed.max())
        
        ms1df_filtered = ms1df[
            (ms1df.mz >= mz_global_min) & 
            (ms1df.mz <= mz_global_max) &
            (ms1df.scan_time_aligned >= st_global_min) &
            (ms1df.scan_time_aligned <= st_global_max)
        ].copy()

        # Generate set_ids for all features
        set_ids = [f'c{clusters[i]}_{i}_i' for i in range(len(sampledf))]
        
        # Use batch method to process all features at once
        if expand_on_miss:
            # First try with normal bounds
            peaks_dict = self[obj_idx].search_for_targeted_mass_features_batch(
                ms1df_filtered,
                mz_mins,
                mz_maxs,
                st_mins,
                st_maxs,
                set_ids,
                obj_idx=obj_idx,
                st_aligned=True
            )
            
            # Retry failed features with expanded bounds
            failed_indices = [i for i, sid in enumerate(set_ids) if peaks_dict[sid].apex_scan == -99]
            if failed_indices:
                failed_ids = [set_ids[i] for i in failed_indices]
                retry_peaks = self[obj_idx].search_for_targeted_mass_features_batch(
                    ms1df_filtered,
                    mz_mins_allowed[failed_indices],
                    mz_maxs_allowed[failed_indices],
                    st_mins_allowed[failed_indices],
                    st_maxs_allowed[failed_indices],
                    failed_ids,
                    obj_idx=obj_idx,
                    st_aligned=True
                )
                peaks_dict.update(retry_peaks)
        else:
            peaks_dict = self[obj_idx].search_for_targeted_mass_features_batch(
                ms1df_filtered,
                mz_mins,
                mz_maxs,
                st_mins,
                st_maxs,
                set_ids,
                obj_idx=obj_idx,
                st_aligned=True
            )
        
        # Assign peaks to induced_mass_features and cluster_dict
        for i in range(len(sampledf)):
            peak = peaks_dict[set_ids[i]]
            self[obj_idx].induced_mass_features[peak.id] = peak
            cluster_dict[clusters[i]] += [set_ids[i]]

        # TODO KRH: Let's try to avoid these steps unless asked for by parameters to pick up speed
        if False:
            self[obj_idx].add_associated_ms1(induced_features = True)
            # need to set drop_if_fail to false for induced features as they will fail
            self[obj_idx].add_peak_metrics(induced_features = True)
            
        self[obj_idx].integrate_mass_features(drop_if_fail = False, induced_features = True)

        if not inplace:
            return self[obj_idx].induced_mass_features
    
    def fill_missing_cluster_features(self):
        """
        Gap-filling for consensus mass features across collection samples.
        
        For clusters present in multiple samples but missing from others, searches
        raw MS1 data to find peaks in expected m/z and retention time windows. This
        creates "induced" mass features for peaks that exist in the data but weren't
        detected in the initial peak detection.
        
        Must be run after add_consensus_mass_features(). Results are accessible via
        induced_mass_features_dataframe property and included in collection_pivot_table
        and collection_consensus_report outputs.

        Parameters
        ----------
        None
            Uses parameters from self.parameters.lcms_collection:
            - consensus_min_sample_fraction: Minimum fraction of samples (0-1) that must contain
              a cluster before gap-filling is attempted
            - gap_fill_expand_on_miss: If True, expands search window when no peak is found
            
        Returns
        -------
        None
            Updates induced_mass_features attribute for each LCMSBase object and
            combines them into induced_mass_features_dataframe.
            
        Raises
        ------
        ValueError
            If cluster_summary_dataframe is not set (must run add_consensus_mass_features first).
            
        Notes
        -----
        - Loads raw MS1 data for each sample, which may be memory intensive
        - Induced features are integrated and metrics calculated automatically
        - Processing can be parallelized using parameters.lcms_collection.cores
        
        See Also
        --------
        add_consensus_mass_features : Creates consensus features before gap-filling
        collection_pivot_table : Includes both regular and induced features
        collection_consensus_report : Reports on complete feature matrix
        """
        
        # Validate prerequisites
        if not hasattr(self, 'cluster_summary_dataframe') or self.cluster_summary_dataframe is None:
            raise ValueError(
                "cluster_summary_dataframe not found. Must run add_consensus_mass_features() first."
            )
        
        # Get parameters from settings
        min_cluster_presence = self.parameters.lcms_collection.consensus_min_sample_fraction
        expand_on_miss = self.parameters.lcms_collection.gap_fill_expand_on_miss
        
        # Validate parameters
        if not 0 <= min_cluster_presence <= 1:
            raise ValueError("consensus_min_sample_fraction must be between 0 and 1")
        
        summarydf = self.cluster_summary_dataframe
        mfdf = self.mass_features_dataframe
        
        sample_ct = len(self.samples)
        
        # Identify clusters present in sufficient samples but not all samples
        missingdf = summarydf[[
            'cluster', 
            'sample_id_nunique', 
            'mz_min', 
            'mz_max', 
            'scan_time_aligned_min', 
            'scan_time_aligned_max'
        ]]
        missingdf = missingdf[missingdf.sample_id_nunique > min_cluster_presence * sample_ct]
        missingdf = missingdf[missingdf.sample_id_nunique != sample_ct]
        
        # Check if there are any clusters to gap-fill
        if len(missingdf) == 0:
            return

        # Find which samples are missing for each cluster
        # Use range(sample_ct) to include all samples, even those with no mass features
        all_sample_ids = list(range(sample_ct))
        missing_samples_list = []
        for c in missingdf.cluster.to_numpy():
            cludf = mfdf[mfdf.cluster == c]
            missing = [x for x in all_sample_ids if x not in cludf.sample_id.unique()]
            missing_samples_list.append(missing)
        missingdf['missing_samples'] = missing_samples_list
        
        # Calculate expanded search windows for expand_on_miss option
        mz_clu_tol = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
        rt_clu_tol = self.parameters.lcms_collection.consensus_rt_tol
        missingdf['mz_max_allowed'] = missingdf.mz_max + mz_clu_tol * missingdf.mz_max
        missingdf['mz_min_allowed'] = missingdf.mz_min - mz_clu_tol * missingdf.mz_min
        missingdf['sta_max_allowed'] = missingdf.scan_time_aligned_max + rt_clu_tol * missingdf.scan_time_aligned_max
        missingdf['sta_min_allowed'] = missingdf.scan_time_aligned_min - rt_clu_tol * missingdf.scan_time_aligned_min

        # Compute cluster dictionary once to avoid recomputing for each sample
        cluster_dict = self.cluster_feature_dictionary
        
        # Process each sample to search for missing features
        if self.parameters.lcms_collection.cores == 1:
            for i in tqdm(range(sample_ct), desc="Gap-filling samples", unit="sample"):
                self._search_for_targeted_mass_features_in_sample(i, missingdf, cluster_dict, expand_on_miss)

        if self.parameters.lcms_collection.cores > 1:
            if self.parameters.lcms_collection.cores > len(self):
                ncores = len(self)
            else:
                ncores = self.parameters.lcms_collection.cores
            pool = multiprocessing.Pool(ncores)
            mp_result = pool.starmap(
                self._search_for_targeted_mass_features_in_sample, 
                [(x, missingdf, cluster_dict, expand_on_miss, False) for x in range(sample_ct)]
            )

            for i in tqdm(range(sample_ct), desc="Collecting gap-filled features", unit="sample"):
                self[i].induced_mass_features = mp_result[i]
                
        self._combine_mass_features(induced_features = True)
        
        # Mark that gap-filling has been performed
        self.missing_mass_features_searched = True
        
        for sample_name in self.samples:
            self._lcms[sample_name].mass_features = {}
    
    def process_samples_pipeline(self, operations, description=None, keep_raw_data=False, show_progress=True):
        """
        Execute a pipeline of operations on all samples in parallel.
        
        This method provides a flexible framework for performing multiple
        sample-level operations in a single parallelized pass, which is more
        efficient than calling separate methods sequentially.
        
        Parameters
        ----------
        operations : list of SampleOperation
            List of operations to perform on each sample, in order.
            Each operation should be an instance of a class derived from
            SampleOperation (see lc_calc_operations module).
        description : str or None, optional
            Progress bar description. If None, automatically generates description
            from operation descriptions (e.g., "gap-filling, reloading features").
            Default is None.
        keep_raw_data : bool, optional
            If True, keeps raw MS data loaded in memory after pipeline completes.
            If False, cleans up raw data to free memory. Default is False.
        show_progress : bool, optional
            If True, displays progress bars during processing. If False, runs silently.
            Default is True.
            
        Returns
        -------
        dict
            Dictionary with results from pipeline execution, keyed by operation name.
            Structure: {operation_name: {sample_id: result, ...}, ...}
            
        Raises
        ------
        ValueError
            If operations list is empty or contains invalid operations.
            
        Notes
        -----
        - Operations are executed sequentially within each sample
        - Samples are processed in parallel based on parameters.lcms_collection.cores
        - Each operation can have conditional execution via can_execute()
        - Results are collected back via collect_results() method of each operation
        - Failed operations for a sample are logged but don't halt processing
        - Raw MS data loaded by operations is automatically cleaned up unless keep_raw_data=True
        
        Examples
        --------
        >>> from corems.mass_spectra.calc.lc_calc_operations import (
        ...     GapFillOperation, ReloadFeaturesOperation
        ... )
        >>> ops = [
        ...     GapFillOperation('gap_fill', expand_on_miss=True),
        ...     ReloadFeaturesOperation('reload', add_ms2=True)
        ... ]
        >>> results = lcms_collection.process_samples_pipeline(ops)
        
        See Also
        --------
        lc_calc_operations : Module containing built-in operation classes
        fill_and_process_features : Convenience method combining common operations
        """
        from corems.mass_spectra.calc.lc_calc_operations import SampleOperation
        
        # Validate operations
        if not operations or len(operations) == 0:
            raise ValueError("operations list cannot be empty")
        
        for op in operations:
            if not isinstance(op, SampleOperation):
                raise ValueError(f"All operations must be SampleOperation instances, got {type(op)}")
        
        # Generate description from operations if not provided
        if description is None:
            operation_descriptions = [op.description for op in operations]
            description = ", ".join(operation_descriptions).capitalize()
        
        # Prepare runtime parameters for each operation
        # This is where we gather collection-level data that operations need
        runtime_params = self._prepare_pipeline_runtime_params(operations)
        runtime_params['keep_raw_data'] = keep_raw_data
        
        # Execute pipeline
        sample_ct = len(self.samples)
        
        if self.parameters.lcms_collection.cores == 1:
            # Serial processing
            results_by_operation = {op.name: {} for op in operations}
            
            if show_progress:
                from tqdm import tqdm
                # Print description on its own line before progress bar
                print(f"\n{description.capitalize()}:")
                iterator = tqdm(range(sample_ct), unit="sample", ncols=80)
            else:
                iterator = range(sample_ct)
            
            for sample_id in iterator:
                sample_results = self._execute_sample_pipeline(
                    sample_id, operations, runtime_params, inplace=True
                )
                # Collect results (collect_results already called in _execute_sample_pipeline when inplace=True)
                # Skip 'sample_id' key which is added for tracking
                for op_name, result in sample_results.items():
                    if op_name != 'sample_id':
                        results_by_operation[op_name][sample_id] = result
        else:
            # Parallel processing
            import multiprocessing
            
            if self.parameters.lcms_collection.cores > sample_ct:
                ncores = sample_ct
            else:
                ncores = self.parameters.lcms_collection.cores
            
            pool = multiprocessing.Pool(ncores)
            
            # Build arguments for each sample
            args_list = [
                (sample_id, operations, runtime_params, False)
                for sample_id in range(sample_ct)
            ]
            
            # Execute in parallel with progress tracking
            results_by_operation = {op.name: {} for op in operations}
            
            if show_progress:
                from tqdm import tqdm
                import time
                
                # Use starmap_async for parallel execution with progress tracking
                async_result = pool.starmap_async(self._execute_sample_pipeline, args_list)
                
                # Poll for completion and update progress bar
                print(description)
                pbar = tqdm(
                    total=sample_ct, 
                    desc="",
                    unit="sample", 
                    position=0,
                    leave=True,
                    dynamic_ncols=True
                )
                prev_completed = 0
                while not async_result.ready():
                    # Get number of completed tasks by checking remaining
                    completed = sample_ct - async_result._number_left
                    if completed > prev_completed:
                        pbar.update(completed - prev_completed)
                        prev_completed = completed
                    time.sleep(0.5)  # Poll every 500ms to avoid spam
                
                # Final update to 100%
                if prev_completed < sample_ct:
                    pbar.update(sample_ct - prev_completed)
                pbar.close()
                
                # Get all results
                mp_results = async_result.get()
            else:
                # Execute without progress
                mp_results = pool.starmap(self._execute_sample_pipeline, args_list)
            
            pool.close()
            pool.join()
            
            # Collect results back into collection
            for result in mp_results:
                sample_id = result.get('sample_id')
                for op in operations:
                    op_result = result.get(op.name)
                    if op_result is not None:
                        op.collect_results(sample_id, op_result, self)
                        results_by_operation[op.name][sample_id] = op_result
        
        return results_by_operation
    
    def _prepare_pipeline_runtime_params(self, operations):
        """
        Prepare runtime parameters needed by operations in the pipeline.
        
        This method gathers collection-level data that operations need,
        such as cluster information for gap-filling or mf_ids for reloading.
        
        Parameters
        ----------
        operations : list of SampleOperation
            List of operations that will be executed
            
        Returns
        -------
        dict
            Dictionary of runtime parameters for operations
        """
        from corems.mass_spectra.calc.lc_calc_operations import (
            GapFillOperation, ReloadFeaturesOperation, MS2SpectralSearchOperation,
            LoadEICsOperation
        )
        
        runtime_params = {}
        
        # Check if any operation needs gap-fill parameters
        needs_gap_fill = any(isinstance(op, GapFillOperation) for op in operations)
        if needs_gap_fill:
            # Prepare gap-fill parameters (same as fill_missing_cluster_features)
            min_cluster_presence = self.parameters.lcms_collection.consensus_min_sample_fraction
            expand_on_miss = self.parameters.lcms_collection.gap_fill_expand_on_miss
            
            summarydf = self.cluster_summary_dataframe
            mfdf = self.mass_features_dataframe
            sample_ct = len(self.samples)
            
            # Identify clusters needing gap-filling
            # Note: cluster_summary_dataframe has 'cluster' as index, need to reset it
            missingdf = summarydf.reset_index()[[
                'cluster', 
                'sample_id_nunique', 
                'mz_min', 
                'mz_max', 
                'scan_time_aligned_min', 
                'scan_time_aligned_max'
            ]].copy()
            missingdf = missingdf[missingdf.sample_id_nunique > min_cluster_presence * sample_ct]
            missingdf = missingdf[missingdf.sample_id_nunique != sample_ct]
            
            if len(missingdf) > 0:
                # Find which samples are missing for each cluster
                # Use range(sample_ct) to include all samples, even those with no mass features
                all_sample_ids = list(range(sample_ct))
                missing_samples_list = []
                for c in missingdf.cluster.to_numpy():
                    cludf = mfdf[mfdf.cluster == c]
                    missing = [x for x in all_sample_ids if x not in cludf.sample_id.unique()]
                    missing_samples_list.append(missing)
                missingdf['missing_samples'] = missing_samples_list
                
                # Calculate expanded search windows
                mz_clu_tol = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
                rt_clu_tol = self.parameters.lcms_collection.consensus_rt_tol
                missingdf['mz_max_allowed'] = missingdf.mz_max + mz_clu_tol * missingdf.mz_max
                missingdf['mz_min_allowed'] = missingdf.mz_min - mz_clu_tol * missingdf.mz_min
                missingdf['sta_max_allowed'] = missingdf.scan_time_aligned_max + rt_clu_tol * missingdf.scan_time_aligned_max
                missingdf['sta_min_allowed'] = missingdf.scan_time_aligned_min - rt_clu_tol * missingdf.scan_time_aligned_min
                
                runtime_params['missingdf'] = missingdf
                runtime_params['cluster_dict'] = self.cluster_feature_dictionary
                runtime_params['expand_on_miss'] = expand_on_miss
        
        # Check if any operation needs reload parameters
        needs_reload = any(isinstance(op, ReloadFeaturesOperation) for op in operations)
        if needs_reload:
            # Use DRY helper method to build sample_mf_map
            sample_mf_map = self.get_sample_mf_map_for_representatives(include_cluster_id=False)
            runtime_params['sample_mf_map'] = sample_mf_map
        
        # Check if any operation needs MS2 spectral search parameters
        needs_ms2_search = any(isinstance(op, MS2SpectralSearchOperation) for op in operations)
        if needs_ms2_search:
            # Pass through pre-prepared spectral library
            if hasattr(self, '_spectral_lib') and self._spectral_lib is not None:
                runtime_params['fe_lib'] = self._spectral_lib
            if hasattr(self, '_spectral_search_molecular_metadata'):
                runtime_params['molecular_metadata'] = self._spectral_search_molecular_metadata
        
        # Check if any operation needs EIC loading parameters
        needs_eic_loading = any(isinstance(op, LoadEICsOperation) for op in operations)
        if needs_eic_loading:
            # Build cluster_mz_dict: map of sample_id -> list of m/z values in clusters
            mfdf = self.mass_features_dataframe
            cluster_mz_dict = {}
            
            # Get all mass features that belong to clusters (cluster is not NaN)
            clustered_mf = mfdf[mfdf['cluster'].notna()]
            
            # Group by sample_id and collect all m/z values associated with eics
            for sample_id in clustered_mf['sample_id'].unique():
                sample_df = clustered_mf[clustered_mf['sample_id'] == sample_id]
                sample = self[sample_id]  # Get the LCMS object for this sample
                
                # Extract _eic_mz from actual mass feature objects, not from dataframe
                eic_mz_list = []
                for mf_id in sample_df['mf_id'].values:
                    if mf_id in sample.mass_features:
                        mf = sample.mass_features[mf_id]
                        if hasattr(mf, '_eic_mz') and mf._eic_mz is not None:
                            eic_mz_list.append(mf._eic_mz)
                
                # Use the collected m/z values, or fallback to empty list if none found
                cluster_mz_dict[sample_id] = list(set(eic_mz_list)) if eic_mz_list else []
            
            runtime_params['cluster_mz_dict'] = cluster_mz_dict
        
        return runtime_params
    
    def _execute_sample_pipeline(self, sample_id, operations, runtime_params, inplace=True):
        """
        Execute a pipeline of operations on a single sample.
        
        This is the worker function called (potentially in parallel) for each sample.
        
        Parameters
        ----------
        sample_id : int
            Sample ID to process
        operations : list of SampleOperation
            Operations to execute in order
        runtime_params : dict
            Runtime parameters prepared by _prepare_pipeline_runtime_params
        inplace : bool, optional
            If True, updates sample in place. If False, returns results for
            multiprocessing. Default is True.
            
        Returns
        -------
        dict
            Dictionary with results from each operation, keyed by operation name.
            If inplace=True, returns results that need to be collected.
            If inplace=False, returns all results for multiprocessing collection.
        """
        results = {}
        
        # Check if any operations need raw MS data
        needs_raw_data = {}  # {ms_level: True/False}
        for op in operations:
            needs_raw, ms_level = op.needs_raw_ms_data()
            if needs_raw and ms_level:
                needs_raw_data[ms_level] = True
        
        # Load raw data once if any operations need it
        # Note: For gap-filling, it loads data internally, so we just track it here
        for ms_level in needs_raw_data.keys():
            # Gap-filling loads its own data, but we want to keep track that it's loaded
            # Other operations can then use the loaded data
            pass
        
        for op in operations:
            # Check if operation can execute on this sample
            sample = self[sample_id]
            if not op.can_execute(sample, self):
                # Skip this operation for this sample if prerequisites aren't met
                # This allows processing to continue for samples that don't have
                # all required data (e.g., MS2 spectra)
                results[op.name] = None
                continue
            
            # Prepare operation-specific runtime params
            op_runtime_params = {}
            
            # Add gap-fill params if this is a gap-fill operation
            from corems.mass_spectra.calc.lc_calc_operations import (
                GapFillOperation, ReloadFeaturesOperation, MS2SpectralSearchOperation, LoadEICsOperation
            )
            
            if isinstance(op, GapFillOperation):
                if 'missingdf' in runtime_params:
                    op_runtime_params['missingdf'] = runtime_params['missingdf']
                    op_runtime_params['cluster_dict'] = runtime_params['cluster_dict']
                    op_runtime_params['expand_on_miss'] = runtime_params['expand_on_miss']
            
            elif isinstance(op, ReloadFeaturesOperation):
                if 'sample_mf_map' in runtime_params:
                    sample_mf_map = runtime_params['sample_mf_map']
                    # Always pass mf_ids_to_load to ensure we only load what's needed
                    # If sample not in map, it has no representatives - pass empty list
                    op_runtime_params['mf_ids_to_load'] = sample_mf_map.get(sample_id, [])
            
            elif isinstance(op, MS2SpectralSearchOperation):
                # Add MS2 spectral search parameters
                if 'fe_lib' in runtime_params:
                    op_runtime_params['fe_lib'] = runtime_params['fe_lib']
                if 'molecular_metadata' in runtime_params:
                    op_runtime_params['molecular_metadata'] = runtime_params['molecular_metadata']
            
            elif isinstance(op, LoadEICsOperation):
                # Add EIC loading parameters
                if 'cluster_mz_dict' in runtime_params:
                    op_runtime_params['cluster_mz_dict'] = runtime_params['cluster_mz_dict']
            
            # Execute the operation
            result = op.execute(sample_id, self, **op_runtime_params)
            results[op.name] = result
            
            # If inplace, collect immediately
            if inplace and result is not None:
                op.collect_results(sample_id, result, self)
        
        # Clean up raw data if requested
        keep_raw_data = runtime_params.get('keep_raw_data', False)
        if not keep_raw_data:
            for ms_level in needs_raw_data.keys():
                if ms_level in self[sample_id]._ms_unprocessed:
                    del self[sample_id]._ms_unprocessed[ms_level]
        
        # Include sample_id in results for tracking (especially important for imap_unordered)
        results['sample_id'] = sample_id
        return results
    
    def process_consensus_features(self, load_representatives=True, perform_gap_filling=True,
                                   add_ms1=False, add_ms2=False,
                                   ms2_scan_filter=None, molecular_formula_search=False,
                                   ms2_spectral_search=False, spectral_lib=None,
                                   molecular_metadata=None,
                                   gather_eics=False,
                                   keep_raw_data=False,
                                   show_progress=True):
        """
        Process consensus mass features across the collection in a single parallelized pass.
        
        This method provides a convenient interface to the sample processing pipeline,
        allowing multiple operations (gap-filling, feature reloading, MS1/MS2 association,
        molecular formula search, and MS2 spectral search) to be performed efficiently in 
        a single pass through all samples.
        
        Parameters
        ----------
        load_representatives : bool, optional
            If True, loads representative mass features from HDF5. Default is True.
        perform_gap_filling : bool, optional
            If True, performs gap-filling for missing cluster features. Default is True.
            This operation loads raw MS1 data which can be reused by subsequent operations.
        add_ms1 : bool, optional
            If True and load_representatives=True, associates MS1 spectra with
            loaded features. Automatically uses raw data from gap-filling if available,
            otherwise uses parser. Spectrum mode is auto-detected. Default is False.
        add_ms2 : bool, optional
            If True and load_representatives=True, associates MS2 spectra with
            loaded features and automatically processes them. Spectrum mode is auto-detected. Default is False.
        ms2_scan_filter : str or None, optional
            Filter string for MS2 scans (e.g., 'hcd'). Default is None.
        molecular_formula_search : bool, optional
            If True, performs molecular formula search on mass features using
            associated MS1 spectra. Requires add_ms1=True or that MS1 spectra
            are already associated. Uses parameters from 
            parameters.mass_spectrum["ms1"].molecular_search. Default is False.
        ms2_spectral_search : bool, optional
            If True, performs MS2 spectral library search using FlashEntropy.
            Requires add_ms2=True and spectral_lib to be provided. Default is False.
        spectral_lib : FlashEntropy library, optional
            Pre-prepared FlashEntropy spectral library for MS2 search.
            Create using MSPInterface.get_metabolomics_spectra_library().
            Required if ms2_spectral_search=True. Default is None.
        molecular_metadata : pd.DataFrame, optional
            Molecular metadata corresponding to spectral_lib.
            Returned from MSPInterface.get_metabolomics_spectra_library().
            Stored as self.spectral_search_molecular_metadata for later export.
            Default is None.
        gather_eics : bool, optional
            If True, loads extracted ion chromatograms (EICs) from HDF5 for all
            mass features with assigned cluster_index (including gap-filled features).
            Enables access to EICs via get_eics_for_cluster(cluster_id) method.
            Requires that EICs were previously exported with export_eics=True.
            Default is False.
        keep_raw_data : bool, optional
            If True, keeps raw MS data loaded in memory after pipeline completes.
            If False, cleans up raw data to free memory. Default is False.
        show_progress : bool, optional
            If True, displays progress bars during processing. If False, runs silently.
            Default is True.
            
        Returns
        -------
        dict
            Dictionary with pipeline results. Keys include:
            - 'gap_fill': dict mapping sample_id to induced mass features (if gap-filling)
            - 'reload': dict mapping sample_id to reloaded mass features (if reloading)
            - 'mf_search': dict mapping sample_id to number of features searched (if molecular formula search)
            - 'ms2_search': dict mapping sample_id to number of spectra searched (if MS2 spectral search)
            
        Raises
        ------
        ValueError
            If neither operation is enabled, or if required parameters are missing.
            
        Notes
        -----
        - Must run add_consensus_mass_features() before calling this method
        - Processes samples in parallel based on parameters.lcms_collection.cores
        - Raw MS1 data loaded by gap-filling is automatically reused by MS1 association
        - MS2 spectral search requires add_ms2=True and msp_file_path
        - FlashEntropy library is created once and reused across all samples
        - More efficient than calling individual methods separately
        - After gap-filling, sets missing_mass_features_searched = True
        - Mass features remain loaded in memory for downstream processing
        - For more advanced workflows, use process_samples_pipeline() directly
        
        Examples
        --------
        >>> # Prepare spectral library for MS2 search
        >>> from corems.molecular_id.search.database_interfaces import MSPInterface
        >>> my_msp = MSPInterface(file_path='path/to/library.msp')
        >>> spectral_lib, molecular_metadata = my_msp.get_metabolomics_spectra_library(
        ...     polarity='negative',
        ...     format='flashentropy',
        ...     normalize=True,
        ...     fe_kwargs={
        ...         'normalize_intensity': True,
        ...         'min_ms2_difference_in_da': 0.02,
        ...         'max_ms2_tolerance_in_da': 0.01,
        ...         'max_indexed_mz': 3000,
        ...         'precursor_ions_removal_da': None,
        ...         'noise_threshold': 0,
        ...     }
        ... )
        >>> 
        >>> # Gap-fill, reload with MS1/MS2, perform molecular formula and spectral search
        >>> results = lcms_collection.process_consensus_features(
        ...     load_representatives=True,
        ...     perform_gap_filling=True,
        ...     add_ms1=True,
        ...     add_ms2=True,
        ...     molecular_formula_search=True,
        ...     ms2_spectral_search=True,
        ...     spectral_lib=spectral_lib,
        ...     molecular_metadata=molecular_metadata
        ... )
        
        See Also
        --------
        process_samples_pipeline : Generic pipeline executor for custom workflows
        fill_missing_cluster_features : Original gap-filling method
        reload_representative_mass_features : Original reload method
        """
        from corems.mass_spectra.calc.lc_calc_operations import (
            GapFillOperation, ReloadFeaturesOperation, MolecularFormulaSearchOperation,
            MS2SpectralSearchOperation, LoadEICsOperation
        )
        
        # Validate that at least one meaningful operation is enabled
        has_operations = (
            perform_gap_filling or 
            load_representatives or 
            molecular_formula_search or 
            ms2_spectral_search or 
            gather_eics or
            add_ms1 or
            add_ms2
        )
        
        if not has_operations:
            raise ValueError(
                "At least one operation must be enabled: perform_gap_filling, load_representatives, "
                "molecular_formula_search, ms2_spectral_search, gather_eics, add_ms1, or add_ms2"
            )
        
        # Validate prerequisites for gap-filling
        if perform_gap_filling:
            if not hasattr(self, 'cluster_summary_dataframe') or self.cluster_summary_dataframe is None:
                raise ValueError(
                    "Cannot perform gap-filling: cluster_summary_dataframe not set. "
                    "You must run add_consensus_mass_features() before calling process_consensus_features()."
                )
        
        # Validate prerequisites for MS2 spectral search
        if ms2_spectral_search:
            if spectral_lib is None:
                raise ValueError(
                    "MS2 spectral search requires spectral_lib to be provided. "
                    "Create it using MSPInterface.get_metabolomics_spectra_library() before calling this method."
                )
            # Check if mass features will be loaded OR are already loaded
            # (The operation's can_execute will check if MS2 spectra are actually present)
            if not load_representatives and not perform_gap_filling:
                # Check if at least one sample has mass features loaded
                # This allows MS2 search on already-loaded features
                has_loaded_features = any(
                    len(self[i].mass_features) > 0 if hasattr(self[i], 'mass_features') and self[i].mass_features is not None else False
                    for i in range(len(self.samples))
                )
                if not has_loaded_features:
                    raise ValueError(
                        "MS2 spectral search requires mass features to be loaded. "
                        "Either set load_representatives=True or perform_gap_filling=True to load them, "
                        "or load them in a previous call to process_consensus_features() before calling "
                        "with ms2_spectral_search=True."
                    )
        
        # Build pipeline
        operations = []
        
        if perform_gap_filling:
            expand_on_miss = self.parameters.lcms_collection.gap_fill_expand_on_miss
            operations.append(GapFillOperation('gap_fill', expand_on_miss=expand_on_miss))
        
        if load_representatives:
            operations.append(ReloadFeaturesOperation(
                'reload',
                add_ms1=add_ms1,
                add_ms2=add_ms2,
                auto_process_ms2=add_ms2,  # Auto-process MS2 if add_ms2 is enabled
                ms2_scan_filter=ms2_scan_filter
            ))
        
        if molecular_formula_search:
            operations.append(MolecularFormulaSearchOperation('mf_search'))
        
        if ms2_spectral_search:
            operations.append(MS2SpectralSearchOperation(
                'ms2_search',
                ms2_scan_filter=ms2_scan_filter
            ))
            # Store spectral library and metadata for runtime preparation
            self._spectral_lib = spectral_lib
            self._spectral_search_molecular_metadata = molecular_metadata
        
        if gather_eics:
            operations.append(LoadEICsOperation('load_eics'))
        
        # Execute pipeline (description auto-generated from operations)
        results = self.process_samples_pipeline(
            operations,
            keep_raw_data=keep_raw_data,
            show_progress=show_progress
        )
        
        # Store molecular metadata if spectral search was performed
        if ms2_spectral_search and hasattr(self, '_spectral_search_molecular_metadata'):
            # This allows users to access the metadata for reporting
            self.spectral_search_molecular_metadata = self._spectral_search_molecular_metadata
        # Post-processing
        if perform_gap_filling:
            # Combine induced mass features into dataframe
            self._combine_mass_features(induced_features=True)
            # Mark that gap-filling has been performed
            self.missing_mass_features_searched = True

            # Add ._eic_mz to induced_mass_features_dataframe if it exists
            if self.induced_mass_features_dataframe is not None and len(self.induced_mass_features_dataframe) > 0:
                eics_mz = []
                for i, row in self.induced_mass_features_dataframe.iterrows():
                    sample_id = row['sample_id']
                    sample = self[sample_id]
                    if row['mf_id'] in sample.induced_mass_features.keys():
                        eic_mz = sample.induced_mass_features[row['mf_id']]._eic_mz
                        eics_mz.append(eic_mz)
                    else:
                        eics_mz.append(None)
                self.induced_mass_features_dataframe['_eic_mz'] = eics_mz

            # Clear mass features from samples to free memory
            for sample_name in self.samples:
                self._lcms[sample_name].induced_mass_features = {}
        
        # Associate EICs with mass features if they were loaded
        # This must happen after all operations complete to work on the actual sample objects
        if gather_eics:
            print("\nAssociating EICs with mass features:")
            from tqdm import tqdm
            
            for sample_id in tqdm(range(len(self.samples)), unit="sample", ncols=80):
                sample = self[sample_id]
                if sample.eics:  # Only if EICs were loaded
                    # Associate EICs with regular mass features
                    sample.associate_eics_with_mass_features(induced=False)
                    # Associate EICs with induced mass features
                    sample.associate_eics_with_mass_features(induced=True)
                
        return results