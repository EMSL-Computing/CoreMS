import numpy as np
import pandas as pd
import warnings
from ripser import ripser
import scipy
from scipy import sparse
from scipy.spatial import KDTree
from sklearn.svm import SVR
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
from ast import literal_eval

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

    def add_peak_metrics(self, induced_features = False):
        """Add peak metrics to the mass features.

        This function calculates the peak metrics for each mass feature and adds them to the mass feature objects.
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

    def find_mass_features(self, ms_level=1, grid=True):
        """Find mass features within an LCMSBase object

        Note that this is a wrapper function that calls the find_mass_features_ph function, but can be extended to support other peak picking methods in the future.

        Parameters
        ----------
        ms_level : int, optional
            The MS level to use for peak picking Default is 1.
        grid : bool, optional
            If True, will regrid the data before running the persistent homology calculations (after checking if the data is gridded),
            used for persistent homology peak picking for profile data only. Default is True.

        Raises
        ------
        ValueError
            If no MS level data is found on the object.
            If persistent homology peak picking is attempted on non-profile mode data.
            If data is not gridded and grid is False.
            If peak picking method is not implemented.

        Returns
        -------
        None, but assigns the mass_features and eics attributes to the object.

        """
        pp_method = self.parameters.lc_ms.peak_picking_method

        if pp_method == "persistent homology":
            msx_scan_df = self.scan_df[self.scan_df["ms_level"] == ms_level]
            if all(msx_scan_df["ms_format"] == "profile"):
                self.find_mass_features_ph(ms_level=ms_level, grid=grid)
            else:
                raise ValueError(
                    "MS{} scans are not profile mode, which is required for persistent homology peak picking.".format(
                        ms_level
                    )
                )
        elif pp_method == "centroided_persistent_homology":
            msx_scan_df = self.scan_df[self.scan_df["ms_level"] == ms_level]
            if all(msx_scan_df["ms_format"] == "centroid"):
                self.find_mass_features_ph_centroid(ms_level=ms_level)
            else:
                raise ValueError(
                    "MS{} scans are not centroid mode, which is required for persistent homology centroided peak picking.".format(
                        ms_level
                    )
                )
        else:
            raise ValueError("Peak picking method not implemented")

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
                    "No mass features found, did you run search_for_targeted_mass_feature() first?"
                )
            if not any(
                [mf.mass_spectrum is not None for mf in mf_dict.values()]
            ):
                raise ValueError(
                    "Mass spectrum must be associated with induced mass features, did you run add_associated_ms1() first?"
                )
            ## remove not found induced mass features
            mf_dict = {k:v for k, v in mf_dict.items() if v.mass_spectrum is not None}

        else:
            mf_dict = self.mass_features
            if len(mf_dict) == 0:
                raise ValueError(
                    "No mass features found, did you run find_mass_features() first?"
                )
            if not all(
                [mf.mass_spectrum is not None for mf in mf_dict.values()]
            ):
                raise ValueError(
                    "Mass spectrum must be associated with each mass feature, did you run add_associated_ms1() first?"
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
            apex_index = np.where(myEIC.scans == apex_scan)[0][0]

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

    def find_mass_features_ph(self, ms_level=1, grid=True):
        """Find mass features within an LCMSBase object using persistent homology.

        Assigns the mass_features attribute to the object (a dictionary of LCMSMassFeature objects, keyed by mass feature id)

        Parameters
        ----------
        ms_level : int, optional
            The MS level to use. Default is 1.
        grid : bool, optional
            If True, will regrid the data before running the persistent homology calculations (after checking if the data is gridded). Default is True.

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

        # Threshold data
        dims = ["mz", "scan_time"]
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

        # Add scan_time
        data_thres = data_thres.merge(self.scan_df[["scan", "scan_time"]], on="scan")
        # Process in chunks if required
        if len(data_thres) > 10000:
            return self._find_mass_features_ph_partition(
                data_thres, dims, persistence_threshold
            )
        else:
            # Process all at once
            return self._find_mass_features_ph_single(
                data_thres, dims, persistence_threshold
            )

    def _find_mass_features_ph_single(self, data_thres, dims, persistence_threshold):
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
        self._populate_mass_features(mass_features_df)

    def _find_mass_features_ph_partition(self, data_thres, dims, persistence_threshold):
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
            self._populate_mass_features(combined_features)
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

    def _populate_mass_features(self, mass_features_df):
        """Populate the mass_features attribute from a DataFrame.

        Parameters
        ----------
        mass_features_df : pd.DataFrame
            DataFrame containing mass feature information.
            Note that the order of this DataFrame will determine the order of mass features in the mass_features attribute.

        Returns
        -------
        None, but assigns the mass_features attribute to the object.
        """
        # Rename scan column to apex_scan
        mass_features_df = mass_features_df.rename(
            columns={"scan": "apex_scan", "scan_time": "retention_time"}
        )

        # Populate mass_features attribute
        self.mass_features = {}
        for row in mass_features_df.itertuples():
            row_dict = mass_features_df.iloc[row.Index].to_dict()
            lcms_feature = LCMSMassFeature(self, **row_dict)
            self.mass_features[lcms_feature.id] = lcms_feature

        if self.parameters.lc_ms.verbose_processing:
            print("Found " + str(len(mass_features_df)) + " initial mass features")

    def find_mass_features_ph_centroid(self, ms_level=1):
        """Find mass features within an LCMSBase object using persistent homology-type approach but with centroided data.

        Parameters
        ----------
        ms_level : int, optional
            The MS level to use. Default is 1.

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

        # Calculate threshold first to avoid multiple operations
        max_intensity = data["intensity"].max()
        threshold = self.parameters.lc_ms.ph_inten_min_rel * max_intensity

        # Create single filter condition and apply to required columns only
        valid_mask = data["intensity"].notna() & (data["intensity"] > threshold)
        required_cols = ["mz", "intensity", "scan"]
        data_thres = data.loc[valid_mask, required_cols].copy()
        data_thres["persistence"] = data_thres["intensity"]

        # Merge with required scan data
        scan_subset = self.scan_df[["scan", "scan_time"]]
        mf_df = data_thres.merge(scan_subset, on="scan", how="inner")
        del data_thres, scan_subset

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
        if mf_c is None or mf_i is None:
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
        dist3d = np.multiply(dist3d, idx)

        # Normalize to 0-1
        mx = dist3d.max()
        if mx > 0:
            # Lower distance is better
            dist3d = dist3d / dist3d.max()

        # Turn zeros to inf (no match)
        dist3d[dist3d == 0] = np.inf

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
            # Drop features that have an absolute_intensity lower than the threshold
            threshold = self.parameters.lcms_collection.mass_feature_anchor_aboslute_intensity_threshold
            mf_df = mf_df[mf_df["absolute_intensity"] > threshold]

        return mf_df

    def attempt_alignment(self, matches_c, matches_i):
        """
        Check if alignment is needed for the LCMS objects in the collection.
        """

        # Hold out a subset of matches_c and matches_i for spline fitting
        matches_c.reset_index(drop=False, inplace=True)
        matches_i.reset_index(drop=False, inplace=True)

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
            in self.parameters.lcms_collection.alignment_acceptance_techinque
        ):
            fraction_improved = np.sum(fit_diff < og_diff) / len(og_diff)
            use_spline_alignment = (
                fraction_improved
                > self.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold
            )
        if (
            "mean_squared_error_improved"
            in self.parameters.lcms_collection.alignment_acceptance_techinque
        ):
            mse_og = np.mean(og_diff**2)
            mse = np.mean(fit_diff**2)
            use_spline_alignment = mse < mse_og

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

            index_steps = (1, -1)
            # Run this twice, once going forward (+1 indexing) and once going backward (-1 indexing)
            for index_step in index_steps:
                # Loop through the other LCMS objects in the collection (going forward)
                i = center_obj_id + index_step
                if i < len(self) and i >= 0:
                    # Grab the first LCMS object after the center object
                    mf_df_i = full_mf_df.loc[self.samples[i]].copy()
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
                                new_times = spl(self[i].scan_df["scan_time"])
                                new_scan_info = self[i].scan_df.copy()
                                new_scan_info["scan_time_aligned"] = new_times
                                self[i].scan_df = new_scan_info
                            else:
                                # Set aligned retention times on scan_df for lc_obj using the original retention times
                                new_scan_info = self[i].scan_df.copy()
                                new_scan_info["scan_time_aligned"] = new_scan_info[
                                    "scan_time"
                                ]
                                self[i].scan_df = new_scan_info

                            i += index_step
                            if i >= len(self) or i < 0:
                                mf_df_i = None
                            else:
                                # Grab the next LCMS object and use the previous spline fitting to get a better starting point
                                mf_df_i = full_mf_df.loc[self.samples[i]].copy()
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
                                sample_name = self.samples[j]
                                self._manifest_dict[sample_name]["use_rt_alignment"] = (
                                    use_spline_alignment
                                )
                                new_scan_info = self[j].scan_df.copy()
                                new_scan_info["scan_time_aligned"] = spl(
                                    self[j].scan_df["scan_time_aligned"]
                                )
                                self[j].scan_df = new_scan_info

        # Set final mass_features_dataframe with the aligned scan_time
        center_sample_name = self.samples[center_obj_ids[0]]
        self._manifest_dict[center_sample_name]["use_rt_alignment"] = False
        new_scan_info = self[center_obj_ids[0]].scan_df.copy()
        new_scan_info["scan_time_aligned"] = new_scan_info["scan_time"]

    def add_consensus_mass_features(self):
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
            
        # TODO KRH: Deal with isomers better? Pool them together and then split them out using samples with 2 as the template?
        
    def summarize_clusters(self):
        """
        Summarize the clusters of mass features by median attributes
        """
        # First check if there are minimum columns in the features dataframe
        if len(self.mass_features_dataframe.columns) < 1:
            return None

        summary_df = (
            self.mass_features_dataframe.groupby("cluster")
            .agg(
                {
                    "mz": ["median", "mean", "std", "max", "min"],
                    "scan_time_aligned": ["median", "mean", "std", "max", "min"],
                    "half_height_width": ["median", "mean", "std", "max", "min"],
                    "tailing_factor": ["median", "mean", "std", "max", "min"],
                    "dispersity_index": ["median", "mean", "std", "max", "min"],
                    "sample_id": ["nunique"],
                    "intensity": ["max", "median", "mean", "std", "max", "min"],
                    "persistence": ["max", "median", "mean", "std", "max", "min"],
                }
            )
            .reset_index()
        )

        # Fix the column names
        summary_df.columns = [
            "_".join(col).strip()
            for col in summary_df.columns.values
            if col != "cluster"
        ]
        summary_df = summary_df.rename(columns={"cluster_": "cluster"})
        summary_df = summary_df.reset_index(drop=True)
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

        fig = plt.figure()
        if show_all:
            plt.scatter(
                df.scan_time_aligned,
                df.mz,
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
            xt = np.ceil(np.max(df.mz))
        if yt == 'yt':
            yt = np.ceil(np.max(df.scan_time))
        if xb == 'xb':
            xb = 0
        if yb == 'yb':
            yb = 0
        plt.ylim(xb, xt)
        plt.xlim(yb, yt)

        kw = dict(
            prop = 'sizes',
            num = len(df.sample_id_nunique.unique())/3,
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
                distances = sdm
            else:
                # Prepare sdm to match shape of existing distances
                distances_truth = distances.copy()
                # make new sparse matrix with same positions as previous 
                # distance matrix but all ones for values
                distances_truth.data = np.ones_like(distances_truth.data)
                # multiply the new sparse matrix (sdm) by this mask to remove 
                # data that doesn't exist in original sparse matrix
                sdm = distances_truth.multiply(sdm)
                sdm.data = sdm.data * dist_weight[i]

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

        if not hasattr(self, 'mass_features_dataframe.cluster'):
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
        sumdf = self.cluster_summary_dataframe.copy()

        numsamples = mfdf.sample_id.max() + 1
        sumdf = sumdf[sumdf.sample_id_nunique > numsamples * clu_size_thresh]

        ## find the ranges for non-outlier values and add them to sumdf
        mergelist = ['cluster']
        for dim in dim_list:
            maxtag = dim + '_outmax'
            mintag = dim + '_outmin'
            mergelist.append(maxtag)
            mergelist.append(mintag)
            sumdf[maxtag] = None
            sumdf[mintag] = None
            for i in range(len(sumdf)):
                sumdf.loc[i, mintag] = sumdf[dim + '_mean'].iloc[i] - 3*sumdf[dim + '_std'].iloc[i]
                sumdf.loc[i, maxtag] = sumdf[dim + '_mean'].iloc[i] + 3*sumdf[dim + '_std'].iloc[i]
                ## If NaN shows up anywhere in dim_min, dim_max calculations, value is set to NaN and it's 
                ## not flagged. This happens when there's not enough values to compute median/std for that 
                ## dimension therefore can't have outliers

        ## add ranges to mfdf and identify mass features that fall outside the ranges
        outdf = pd.merge(mfdf, sumdf[mergelist], on = 'cluster')
        outtags = ['cluster']
        for dim in dim_list:
            dimtag = dim + '_outlier'
            outtags.append(dimtag)
            outdf[dimtag] = np.where(
                ((outdf[dim] > outdf[dim + '_outmax'])) | ((outdf[dim] < outdf[dim + '_outmin'])), 
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
            
    def search_for_missing_mass_features_in_one_sample(self, samplename, threshold = 0.5, tol_flag = 0):
        '''
        Work in progress temporary code
        threshold default to 0.5 --> 
            only consider clusters that contain at least 50% of the sample
        tol_flag default to 0 -->
            don't check for possible mass features on the edges of the cluster
            for sampleindex == 7, tol_flag == 1 picks up 6 more mass features
        '''
        summarydf = self.cluster_summary_dataframe
        mfdf = self.mass_features_dataframe
        sampleindex = self.samples.index(samplename)
        self.load_raw_data(sampleindex, 1)
        
        sample_ct = len(self.samples)
        missingdf = summarydf[[
            'cluster', 
            'sample_id_nunique', 
            'mz_min', 
            'mz_max', 
            'scan_time_aligned_min', 
            'scan_time_aligned_max'
        ]]
        missingdf = missingdf[missingdf.sample_id_nunique > threshold*(sample_ct)]
        missingdf = missingdf[missingdf.sample_id_nunique != sample_ct]

        missingdf['missing_samples'] = None
        for c in missingdf.cluster.to_numpy():
            cludf = mfdf[mfdf.cluster == c]
            missingdf.loc[c, 'missing_samples'] = str(
                [x for x in mfdf.sample_name.unique() if x not in cludf.sample_name.unique()]
            )
        missingdf['missing_samples'] = missingdf.missing_samples.apply(literal_eval)

        ## to get clusters missing data based on sample name:
        sampledf = missingdf[
            missingdf.missing_samples.apply(lambda x: samplename in x)
        ].reset_index(drop = True).copy()
        
        if tol_flag == 1:
            mz_clu_tol = self.parameters.lcms_collection.consensus_mz_tol_ppm * 1e-6
            rt_clu_tol = self.parameters.lcms_collection.consensus_rt_tol
            sampledf['mz_max_allowed'] = sampledf.mz_max + mz_clu_tol*sampledf.mz_max
            sampledf['mz_min_allowed'] = sampledf.mz_min - mz_clu_tol*sampledf.mz_min
            sampledf['sta_max_allowed'] = sampledf.scan_time_aligned_max + rt_clu_tol*sampledf.scan_time_aligned_max
            sampledf['sta_min_allowed'] = sampledf.scan_time_aligned_min - rt_clu_tol*sampledf.scan_time_aligned_min

        ms1df = self[sampleindex]._ms_unprocessed[1].copy()
        scan_df = self[sampleindex].scan_df[['scan', 'scan_time_aligned']]
        ms1df = pd.merge(ms1df, scan_df, on = 'scan')

        for i in range(len(sampledf)):
            mz_min = sampledf.mz_min.iloc[i]
            mz_max = sampledf.mz_max.iloc[i]
            st_min = sampledf.scan_time_aligned_min.iloc[i]
            st_max = sampledf.scan_time_aligned_max.iloc[i]
            found_feature = self[sampleindex].search_for_targeted_mass_feature(            
                ms1df, 
                mz_min, 
                mz_max, 
                st_min, 
                st_max,
                set_id = 'c' + str(sampledf.cluster.iloc[i]) + '_' + str(i) + '_i',
                obj_idx = sampleindex,
                st_aligned = True
            )

            if tol_flag == 1 and found_feature.apex_scan == -99:
                mz_min = sampledf.mz_min_allowed.iloc[i] 
                mz_max = sampledf.mz_max_allowed.iloc[i]
                st_min = sampledf.sta_min_allowed.iloc[i] 
                st_max = sampledf.sta_max_allowed.iloc[i] 

                found_feature = self[sampleindex].search_for_targeted_mass_feature(            
                    ms1df, 
                    mz_min, 
                    mz_max, 
                    st_min, 
                    st_max,
                    set_id = 'c' + str(sampledf.cluster.iloc[i]) + '_' + str(i) + '_i',
                    obj_idx = sampleindex,
                    st_aligned = True
                )
        
            self[sampleindex].induced_mass_features[found_feature.id] = found_feature

        self[sampleindex].add_associated_ms1(induced_features = True)
        # need to set drop_if_fail to false for induced features as they will fail
        self[sampleindex].integrate_mass_features(drop_if_fail = False, induced_features = True)
        self[sampleindex].add_peak_metrics(induced_features = True)
