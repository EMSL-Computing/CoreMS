import numpy as np
import pandas as pd
from ripser import ripser
from scipy import sparse
from scipy.spatial import KDTree

from corems.chroma_peak.factory.chroma_peak_classes import LCMSMassFeature
from corems.mass_spectra.calc import SignalProcessing as sp
from corems.mass_spectra.factory.chromat_data import EIC_Data
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile


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

    def add_peak_metrics(self):
        """Add peak metrics to the mass features.

        This function calculates the peak metrics for each mass feature and adds them to the mass feature objects.
        """
        # Check that at least some mass features have eic data
        if not any([mf._eic_data is not None for mf in self.mass_features.values()]):
            raise ValueError(
                "No mass features have EIC data. Run integrate_mass_features first."
            )

        for mass_feature in self.mass_features.values():
            # Check if the mass feature has been integrated
            if mass_feature._eic_data is not None:
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
            If True, will regrid the data before running the persistent homology calculations (after checking if the data is gridded, used for persistent homology peak picking. Default is True.

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
                self.cluster_mass_features(drop_children=True, sort_by="persistence")
            else:
                raise ValueError(
                    "MS{} scans are not profile mode, which is required for persistent homology peak picking.".format(
                        ms_level
                    )
                )
        else:
            raise ValueError("Peak picking method not implemented")

    def integrate_mass_features(
        self, drop_if_fail=True, drop_duplicates=True, ms_level=1
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
        if self.mass_features is not None:
            mf_df = self.mass_features_to_df().copy()
        else:
            raise ValueError(
                "No mass features found, did you run find_mass_features() first?"
            )
        # Check if mass_spectrum exists on each mass feature
        if not all(
            [mf.mass_spectrum is not None for mf in self.mass_features.values()]
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

        mzs_to_extract = np.unique(mf_df["mz"].values)
        mzs_to_extract.sort()

        # Get EICs for each unique mz in mass features list
        for mz in mzs_to_extract:
            mz_max = mz + self.parameters.lc_ms.eic_tolerance_ppm * mz / 1e6
            mz_min = mz - self.parameters.lc_ms.eic_tolerance_ppm * mz / 1e6
            raw_data_sub = raw_data[
                (raw_data["mz"] >= mz_min) & (raw_data["mz"] <= mz_max)
            ].reset_index(drop=True)
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
        mf_df["area"] = np.nan
        for idx, mass_feature in mf_df.iterrows():
            mz = mass_feature.mz
            apex_scan = mass_feature.apex_scan

            # Pull EIC data and find apex scan index
            myEIC = self.eics[mz]
            self.mass_features[idx]._eic_data = myEIC
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
                # Add start and final scan to mass_features and EICData
                left_scan, right_scan = (
                    myEIC.scans[l_a_r_scan_idx[0][0]],
                    myEIC.scans[l_a_r_scan_idx[0][2]],
                )
                mf_scan_apex = [(left_scan, int(apex_scan), right_scan)]
                myEIC.apexes = myEIC.apexes + mf_scan_apex
                self.mass_features[idx].start_scan = left_scan
                self.mass_features[idx].final_scan = right_scan

                # Find area under peak using limits from EIC centroid detector, add to mass_features and EICData
                area = np.trapz(
                    myEIC.eic_smoothed[l_a_r_scan_idx[0][0] : l_a_r_scan_idx[0][2] + 1],
                    myEIC.time[l_a_r_scan_idx[0][0] : l_a_r_scan_idx[0][2] + 1],
                )
                mf_df.at[idx, "area"] = area
                myEIC.areas = myEIC.areas + [area]
                self.eics[mz] = myEIC
                self.mass_features[idx]._area = area
            else:
                if drop_if_fail is True:
                    self.mass_features.pop(idx)

        if drop_duplicates:
            # Prepare mass feature dataframe
            mf_df = self.mass_features_to_df().copy()

            # For each mass feature, find all mass features within the clustering tolerance ppm and drop if their start and end times are within another mass feature
            # Kepp the first mass fea
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
            corr_subset = corr.loc[mass_feature.mz,]
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

    def grid_data(self, data):
        """Grid the data in the mz dimension.

        Data must be gridded prior to persistent homology calculations.

        Parameters
        ----------
        data : DataFrame
            The input data containing mz, scan, scan_time, and intensity columns.

        Returns
        -------
        DataFrame
            The gridded data with mz, scan, scan_time, and intensity columns.

        Raises
        ------
        ValueError
            If gridding fails.
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

        # Check if grid worked and return
        if self.check_if_grid(new_data_w):
            return new_data_w
        else:
            raise ValueError("Gridding failed")

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

        # Add build factors and add scan_time
        data_thres = data_thres.merge(self.scan_df[["scan", "scan_time"]], on="scan")
        factors = {
            dim: pd.factorize(data_thres[dim], sort=True)[1].astype(np.float32)
            for dim in dims
        }  # this is return a float64 index

        # Build indexes
        index = {
            dim: np.searchsorted(factors[dim], data_thres[dim]).astype(np.float32)
            for dim in factors
        }

        # Smooth data
        iterations = self.parameters.lc_ms.ph_smooth_it
        smooth_radius = [
            self.parameters.lc_ms.ph_smooth_radius_mz,
            self.parameters.lc_ms.ph_smooth_radius_scan,
        ]  # mz, scan_time smoothing radius (in steps)

        index = np.vstack([index[dim] for dim in dims]).T
        V = data_thres["intensity"].values
        resid = np.inf
        for i in range(iterations):
            # Previous iteration
            V_prev = V.copy()
            resid_prev = resid
            V = self.sparse_mean_filter(index, V, radius=smooth_radius)

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
        data_thres["intensity"] = V

        # Use persistent homology to find regions of interest
        pidx, pers = self.sparse_upper_star(index, V)
        pidx = pidx[pers > 1]
        pers = pers[pers > 1]

        # Get peaks
        peaks = data_thres.iloc[pidx, :].reset_index(drop=True)

        # Add persistence column
        peaks["persistence"] = pers
        mass_features = peaks.sort_values(
            by="persistence", ascending=False
        ).reset_index(drop=True)

        # Filter by persistence threshold
        persistence_threshold = (
            self.parameters.lc_ms.ph_persis_min_rel * data.intensity.max()
        )
        mass_features = mass_features.loc[
            mass_features["persistence"] > persistence_threshold, :
        ].reset_index(drop=True)

        # Rename scan column to apex_scan
        mass_features = mass_features.rename(
            columns={"scan": "apex_scan", "scan_time": "retention_time"}
        )

        # Populate mass_features attribute
        self.mass_features = {}
        for row in mass_features.itertuples():
            row_dict = mass_features.iloc[row.Index].to_dict()
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
        verbose = self.parameters.lc_ms.verbose_processing

        if self.mass_features is None:
            raise ValueError("No mass features found, run find_mass_features() first")
        if len(self.mass_features) > 400000:
            raise ValueError(
                "Too many mass featuers of interest found, run find_mass_features() with a higher intensity threshold"
            )
        dims = ["mz", "scan_time"]
        mf_df_og = self.mass_features_to_df()
        mf_df = mf_df_og.copy()

        # Sort mass features by sort_by column, make mf_id its own column for easier bookkeeping
        mf_df = mf_df.sort_values(by=sort_by, ascending=False).reset_index(drop=False)

        tol = [
            self.parameters.lc_ms.mass_feature_cluster_mz_tolerance_rel,
            self.parameters.lc_ms.mass_feature_cluster_rt_tolerance,
        ]  # mz, in relative; scan_time in minutes
        relative = [True, False]

        # Compute inter-feature distances
        distances = None
        for i in range(len(dims)):
            # Construct k-d tree
            values = mf_df[dims[i]].values
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

        # Extract indices of within-tolerance points
        distances = distances.tocoo()
        pairs = np.stack((distances.row, distances.col), axis=1)
        pairs_df = pd.DataFrame(pairs, columns=["parent", "child"])
        pairs_df = pairs_df.set_index("parent")

        to_drop = []
        while not pairs_df.empty:
            # Find root_parents and their children
            root_parents = np.setdiff1d(
                np.unique(pairs_df.index.values), np.unique(pairs_df.child.values)
            )
            children_of_roots = pairs_df.loc[root_parents, "child"].unique()
            to_drop = np.append(to_drop, children_of_roots)

            # Remove root_children as possible parents from pairs_df for next iteration
            pairs_df = pairs_df.drop(index=children_of_roots, errors="ignore")
            pairs_df = pairs_df.reset_index().set_index("child")
            # Remove root_children as possible children from pairs_df for next iteration
            pairs_df = pairs_df.drop(index=children_of_roots)

            # Prepare for next iteration
            pairs_df = pairs_df.reset_index().set_index("parent")

        # Drop mass features that are not cluster parents
        mf_df = mf_df.drop(index=np.array(to_drop))

        # Set index back to mf_id
        mf_df = mf_df.set_index("mf_id")
        if verbose:
            print(str(len(mf_df)) + " mass features remaining")

        mf_df_new = mf_df_og.copy()
        mf_df_new["cluster_parent"] = np.where(
            np.isin(mf_df_new.index, mf_df.index), True, False
        )

        # get mass feature ids of features that are not cluster parents
        cluster_daughters = mf_df_new[mf_df_new["cluster_parent"] == False].index.values
        if drop_children is True:
            # Drop mass features that are not cluster parents from self
            self.mass_features = {
                k: v
                for k, v in self.mass_features.items()
                if k not in cluster_daughters
            }
        else:
            return cluster_daughters
