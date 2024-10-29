from pathlib import Path

import numpy as np
import pandas as pd
import warnings

from corems.encapsulation.factory.parameters import LCMSParameters
from corems.mass_spectra.calc.lc_calc import LCCalculations, PHCalculations
from corems.molecular_id.search.lcms_spectral_search import LCMSSpectralSearch
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile


class MassSpectraBase:
    """Base class for mass spectra objects.

    Parameters
    -----------
    file_location : str or Path
        The location of the file containing the mass spectra data.
    analyzer : str, optional
        The type of analyzer used to generate the mass spectra data. Defaults to 'Unknown'.
    instrument_label : str, optional
        The type of instrument used to generate the mass spectra data. Defaults to 'Unknown'.
    sample_name : str, optional
        The name of the sample; defaults to the file name if not provided to the parser. Defaults to None.
    spectra_parser : object, optional
        The spectra parser object used to create the mass spectra object. Defaults to None.

    Attributes
    -----------
    spectra_parser_class : class
        The class of the spectra parser used to create the mass spectra object.
    file_location : str or Path
        The location of the file containing the mass spectra data.
    sample_name : str
        The name of the sample; defaults to the file name if not provided to the parser.
    analyzer : str
        The type of analyzer used to generate the mass spectra data. Derived from the spectra parser.
    instrument_label : str
        The type of instrument used to generate the mass spectra data. Derived from the spectra parser.
    _scan_info : dict
        A dictionary containing the scan data with columns for scan number, scan time, ms level, precursor m/z,
        scan text, and scan window (lower and upper).
        Associated with the property scan_df, which returns a pandas DataFrame or can set this attribute from a pandas DataFrame.
    _ms : dict
        A dictionary containing mass spectra for the dataset, keys of dictionary are scan numbers. Initialized as an empty dictionary.
    _ms_unprocessed: dictionary of pandas.DataFrames or None
        A dictionary of unprocssed mass spectra data, as an (optional) intermediate data product for peak picking.
        Key is ms_level, and value is dataframe with columns for scan number, m/z, and intensity. Default is None.

    Methods
    --------
    * add_mass_spectra(scan_list, spectrum_mode: str = 'profile', use_parser = True, auto_process=True).
        Add mass spectra (or singlel mass spectrum) to _ms slot, from a list of scans
    * get_time_of_scan_id(scan).
        Returns the scan time for the specified scan number.
    """

    def __init__(
        self,
        file_location,
        analyzer="Unknown",
        instrument_label="Unknown",
        sample_name=None,
        spectra_parser=None,
    ):
        if isinstance(file_location, str):
            file_location = Path(file_location)
        else:
            file_location = file_location
        if not file_location.exists():
            raise FileExistsError("File does not exist: " + str(file_location))

        if sample_name:
            self.sample_name = sample_name
        else:
            self.sample_name = file_location.stem

        self.file_location = file_location
        self.analyzer = analyzer
        self.instrument_label = instrument_label

        # Add the spectra parser class to the object if it is not None
        if spectra_parser is not None:
            self.spectra_parser_class = spectra_parser.__class__
            self.spectra_parser = spectra_parser
            # Check that spectra_pasrser.sample_name is same as sample_name etc, raise warning if not
            if (
                self.sample_name is not None
                and self.sample_name != self.spectra_parser.sample_name
            ):
                warnings.warn(
                    "sample_name provided to MassSpectraBase object does not match sample_name provided to spectra parser object",
                    UserWarning,
                )
            if self.analyzer != self.spectra_parser.analyzer:
                warnings.warn(
                    "analyzer provided to MassSpectraBase object does not match analyzer provided to spectra parser object",
                    UserWarning,
                )
            if self.instrument_label != self.spectra_parser.instrument_label:
                warnings.warn(
                    "instrument provided to MassSpectraBase object does not match instrument provided to spectra parser object",
                    UserWarning,
                )
            if self.file_location != self.spectra_parser.file_location:
                warnings.warn(
                    "file_location provided to MassSpectraBase object does not match file_location provided to spectra parser object",
                    UserWarning,
                )

        # Instantiate empty dictionaries for scan information and mass spectra
        self._scan_info = {}
        self._ms = {}
        self._ms_unprocessed = {}

    def add_mass_spectrum(self, mass_spec):
        """Adds a mass spectrum to the dataset.

        Parameters
        -----------
        mass_spec : MassSpectrum
            The corems MassSpectrum object to be added to the dataset.

        Notes
        -----
        This is a helper function for the add_mass_spectra() method, and is not intended to be called directly.
        """
        # check if mass_spec has a scan_number attribute
        if not hasattr(mass_spec, "scan_number"):
            raise ValueError(
                "Mass spectrum must have a scan_number attribute to be added to the dataset correctly"
            )
        self._ms[mass_spec.scan_number] = mass_spec

    def add_mass_spectra(
        self,
        scan_list,
        spectrum_mode=None,
        ms_level=1,
        use_parser=True,
        auto_process=True,
        ms_params=None,
    ):
        """Add mass spectra to _ms dictionary, from a list of scans or single scan

        Notes
        -----
        The mass spectra will inherit the mass_spectrum, ms_peak, and molecular_search parameters from the LCMSBase object.


        Parameters
        -----------
        scan_list : list of ints
            List of scans to use to populate _ms slot
        spectrum_mode : str or None
            The spectrum mode to use for the mass spectra.
            If None, method will use the spectrum mode from the spectra parser to ascertain the spectrum mode (this allows for mixed types).
            Defaults to None.
        ms_level : int, optional
            The MS level to use for the mass spectra.
            This is used to pass the molecular_search parameters from the LCMS object to the individual MassSpectrum objects.
            Defaults to 1.
        using_parser : bool
            Whether to use the mass spectra parser to get the mass spectra.  Defaults to True.
        auto_process : bool
            Whether to auto-process the mass spectra.  Defaults to True.
        ms_params : MSParameters or None
            The mass spectrum parameters to use for the mass spectra.  If None, uses the globally set MSParameters.

        Raises
        ------
        TypeError
            If scan_list is not a list of ints
        ValueError
            If polarity is not 'positive' or 'negative'
            If ms_level is not 1 or 2
        """

        # check if scan_list is a list or a single int; if single int, convert to list
        if isinstance(scan_list, int):
            scan_list = [scan_list]
        if not isinstance(scan_list, list):
            raise TypeError("scan_list must be a list of integers")
        for scan in scan_list:
            if not isinstance(scan, int):
                raise TypeError("scan_list must be a list of integers")

        # set polarity to -1 if negative mode, 1 if positive mode (for mass spectrum creation)
        if self.polarity == "negative":
            polarity = -1
        elif self.polarity == "positive":
            polarity = 1
        else:
            raise ValueError(
                "Polarity not set for dataset, must be a either 'positive' or 'negative'"
            )

        # is not using_parser, check that ms1 and ms2 are not None
        if not use_parser:
            if ms_level not in self._ms_unprocessed.keys():
                raise ValueError(
                    "ms_level {} not found in _ms_unprocessed dictionary".format(
                        ms_level
                    )
                )

        scan_list = list(set(scan_list))
        scan_list.sort()
        if not use_parser:
            if self._ms_unprocessed[ms_level] is None:
                raise ValueError(
                    "No unprocessed data found for ms_level {}".format(ms_level)
                )
            if (
                len(
                    np.setdiff1d(
                        scan_list, self._ms_unprocessed[ms_level].scan.tolist()
                    )
                )
                > 0
            ):
                raise ValueError(
                    "Not all scans in scan_list are present in the unprocessed data"
                )
            # Prepare the ms_df for parsing
            ms_df = self._ms_unprocessed[ms_level].copy().set_index("scan", drop=False)

        for scan in scan_list:
            ms = None
            if spectrum_mode is None:
                # get spectrum mode from _scan_info
                spectrum_mode_scan = self.scan_df.loc[scan, "ms_format"]
            else:
                spectrum_mode_scan = spectrum_mode
            # Instantiate the mass spectrum object using the parser or the unprocessed data
            if not use_parser:
                my_ms_df = ms_df.loc[scan]
                if spectrum_mode_scan == "profile":
                    # Check this - it might be better to use the MassSpectrumProfile class to instantiate the mass spectrum
                    ms = ms_from_array_profile(
                        my_ms_df.mz,
                        my_ms_df.intensity,
                        self.file_location,
                        polarity=polarity,
                        auto_process=False,
                    )
                else:
                    raise ValueError(
                        "Only profile mode is supported for unprocessed data"
                    )
            if use_parser:
                ms = self.spectra_parser.get_mass_spectrum_from_scan(
                    scan_number=scan,
                    spectrum_mode=spectrum_mode_scan,
                    auto_process=False,
                )

            # Set the mass spectrum parameters, auto-process if auto_process is True, and add to the dataset
            if ms is not None:
                if ms_params is not None:
                    ms.parameters = ms_params
                ms.scan_number = scan
                if auto_process:
                    ms.process_mass_spec()
                self.add_mass_spectrum(ms)

    def get_time_of_scan_id(self, scan):
        """Returns the scan time for the specified scan number.

        Parameters
        -----------
        scan : int
            The scan number of the desired scan time.

        Returns
        --------
        float
            The scan time for the specified scan number (in minutes).

        Raises
        ------
        ValueError
            If no scan time is found for the specified scan number.
        """
        # Check if _retenion_time_list is empty and raise error if so
        if len(self._retention_time_list) == 0:
            raise ValueError("No retention times found in dataset")
        rt = self._retention_time_list[self._scans_number_list.index(scan)]
        return rt

    @property
    def scan_df(self):
        """
        pandas.DataFrame : A pandas DataFrame containing the scan info data with columns for scan number, scan time, ms level, precursor m/z, scan text, and scan window (lower and upper).
        """
        scan_df = pd.DataFrame.from_dict(self._scan_info)
        return scan_df

    @scan_df.setter
    def scan_df(self, df):
        """
        Sets the scan data for the dataset.

        Parameters
        -----------
        df : pandas.DataFrame
            A pandas DataFrame containing the scan data with columns for scan number, scan time, ms level,
            precursor m/z, scan text, and scan window (lower and upper).
        """
        self._scan_info = df.to_dict()

    def __getitem__(self, scan_number):
        return self._ms.get(scan_number)


class LCMSBase(MassSpectraBase, LCCalculations, PHCalculations, LCMSSpectralSearch):
    """A class representing a liquid chromatography-mass spectrometry (LC-MS) data object.

    This class is not intended to be instantiated directly, but rather to be instantiated by an appropriate mass spectra parser using the get_lcms_obj() method.

    Parameters
    -----------
    file_location : str or Path
        The location of the file containing the mass spectra data.
    analyzer : str, optional
        The type of analyzer used to generate the mass spectra data. Defaults to 'Unknown'.
    instrument_label : str, optional
        The type of instrument used to generate the mass spectra data. Defaults to 'Unknown'.
    sample_name : str, optional
        The name of the sample; defaults to the file name if not provided to the parser. Defaults to None.
    spectra_parser : object, optional
        The spectra parser object used to create the mass spectra object. Defaults to None.

    Attributes
    -----------
    polarity : str
        The polarity of the ionization mode used for the dataset.
    _parameters : LCMSParameters
        The parameters used for all methods called on the LCMSBase object. Set upon instantiation from LCMSParameters.
    _retention_time_list : numpy.ndarray
        An array of retention times for the dataset.
    _scans_number_list : list
        A list of scan numbers for the dataset.
    _tic_list : numpy.ndarray
        An array of total ion current (TIC) values for the dataset.
    eics : dict
        A dictionary containing extracted ion chromatograms (EICs) for the dataset.
        Key is the mz of the EIC. Initialized as an empty dictionary.
    mass_features : dictionary of LCMSMassFeature objects
        A dictionary containing mass features for the dataset.
        Key is mass feature ID. Initialized as an empty dictionary.
    spectral_search_results : dictionary of MS2SearchResults objects
        A dictionary containing spectral search results for the dataset.
        Key is scan number : precursor mz. Initialized as an empty dictionary.

    Methods
    --------
    * get_parameters_json().
        Returns the parameters used for the LC-MS analysis in JSON format.
    * add_associated_ms2_dda(add_to_lcmsobj=True, auto_process=True, use_parser=True)
        Adds which MS2 scans are associated with each mass feature to the
        mass_features dictionary and optionally adds the MS2 spectra to the _ms dictionary.
    * add_associated_ms1(add_to_lcmsobj=True, auto_process=True, use_parser=True)
        Adds the MS1 spectra associated with each mass feature to the
        mass_features dictionary and adds the MS1 spectra to the _ms dictionary.
    * mass_features_to_df()
        Returns a pandas dataframe summarizing the mass features in the dataset.
    * set_tic_list_from_data(overwrite=False)
        Sets the TIC list from the mass spectrum objects within the _ms dictionary.
    * set_retention_time_from_data(overwrite=False)
        Sets the retention time list from the data in the _ms dictionary.
    * set_scans_number_from_data(overwrite=False)
        Sets the scan number list from the data in the _ms dictionary.
    """

    def __init__(
        self,
        file_location,
        analyzer="Unknown",
        instrument_label="Unknown",
        sample_name=None,
        spectra_parser=None,
    ):
        super().__init__(
            file_location, analyzer, instrument_label, sample_name, spectra_parser
        )
        self.polarity = ""
        self._parameters = LCMSParameters()
        self._retention_time_list = []
        self._scans_number_list = []
        self._tic_list = []
        self.eics = {}
        self.mass_features = {}
        self.spectral_search_results = {}

    def get_parameters_json(self):
        """Returns the parameters stored for the LC-MS object in JSON format.

        Returns
        --------
        str
            The parameters used for the LC-MS analysis in JSON format.
        """
        return self.parameters.to_json()

    def remove_unprocessed_data(self, ms_level=None):
        """Removes the unprocessed data from the LCMSBase object.

        Parameters
        -----------
        ms_level : int, optional
            The MS level to remove the unprocessed data for. If None, removes unprocessed data for all MS levels.

        Raises
        ------
        ValueError
            If ms_level is not 1 or 2.

        Notes
        -----
        This method is useful for freeing up memory after the data has been processed.
        """
        if ms_level is None:
            for ms_level in self._ms_unprocessed.keys():
                self._ms_unprocessed[ms_level] = None
        if ms_level not in [1, 2]:
            raise ValueError("ms_level must be 1 or 2")
        self._ms_unprocessed[ms_level] = None

    def add_associated_ms2_dda(
        self,
        auto_process=True,
        use_parser=True,
        spectrum_mode=None,
        ms_params_key="ms2",
        scan_filter=None,
    ):
        """Add MS2 spectra associated with mass features to the dataset.

        Populates the mass_features ms2_scan_numbers attribute (on mass_features dictionary on LCMSObject)

        Parameters
        -----------
        auto_process : bool, optional
            If True, auto-processes the MS2 spectra before adding it to the object's _ms dictionary. Default is True.
        use_parser : bool, optional
            If True, envoke the spectra parser to get the MS2 spectra. Default is True.
        spectrum_mode : str or None, optional
            The spectrum mode to use for the mass spectra.  If None, method will use the spectrum mode
            from the spectra parser to ascertain the spectrum mode (this allows for mixed types).
            Defaults to None. (faster if defined, otherwise will check each scan)
        ms_params_key : string, optional
            The key of the mass spectrum parameters to use for the mass spectra, accessed from the LCMSObject.parameters.mass_spectrum attribute.
            Defaults to 'ms2'.
        scan_filter : str
            A string to filter the scans to add to the _ms dictionary.  If None, all scans are added.  Defaults to None.
            "hcd" will pull out only HCD scans.

        Raises
        ------
        ValueError
            If mass_features is not set, must run find_mass_features() first.
            If no MS2 scans are found in the dataset.
            If no precursor m/z values are found in MS2 scans, not a DDA dataset.
        """
        # Check if mass_features is set, raise error if not
        if self.mass_features is None:
            raise ValueError(
                "mass_features not set, must run find_mass_features() first"
            )

        # reconfigure ms_params to get the correct mass spectrum parameters from the key
        ms_params = self.parameters.mass_spectrum[ms_params_key]

        mf_df = self.mass_features_to_df().copy()
        # Find ms2 scans that have a precursor m/z value
        ms2_scans = self.scan_df[self.scan_df.ms_level == 2]
        ms2_scans = ms2_scans[~ms2_scans.precursor_mz.isna()]
        # drop ms2 scans that have no tic
        ms2_scans = ms2_scans[ms2_scans.tic > 0]
        if ms2_scans is None:
            raise ValueError("No DDA scans found in dataset")

        if scan_filter is not None:
            ms2_scans = ms2_scans[ms2_scans.scan_text.str.contains(scan_filter)]
        # set tolerance in rt space (in minutes) and mz space (in daltons)
        time_tol = self.parameters.lc_ms.ms2_dda_rt_tolerance
        mz_tol = self.parameters.lc_ms.ms2_dda_mz_tolerance

        # for each mass feature, find the ms2 scans that are within the roi scan time and mz range
        dda_scans = []
        for i, row in mf_df.iterrows():
            ms2_scans_filtered = ms2_scans[
                ms2_scans.scan_time.between(
                    row.scan_time - time_tol, row.scan_time + time_tol
                )
            ]
            ms2_scans_filtered = ms2_scans_filtered[
                ms2_scans_filtered.precursor_mz.between(
                    row.mz - mz_tol, row.mz + mz_tol
                )
            ]
            dda_scans = dda_scans + ms2_scans_filtered.scan.tolist()
            self.mass_features[i].ms2_scan_numbers = (
                ms2_scans_filtered.scan.tolist()
                + self.mass_features[i].ms2_scan_numbers
            )
        # add to _ms attribute
        self.add_mass_spectra(
            scan_list=list(set(dda_scans)),
            auto_process=auto_process,
            spectrum_mode=spectrum_mode,
            use_parser=use_parser,
            ms_params=ms_params,
        )
        # associate appropriate _ms attribute to appropriate mass feature's ms2_mass_spectra attribute
        for mf_id in self.mass_features:
            if self.mass_features[mf_id].ms2_scan_numbers is not None:
                for dda_scan in self.mass_features[mf_id].ms2_scan_numbers:
                    if dda_scan in self._ms.keys():
                        self.mass_features[mf_id].ms2_mass_spectra[dda_scan] = self._ms[
                            dda_scan
                        ]

    def add_associated_ms1(
        self, auto_process=True, use_parser=True, spectrum_mode=None
    ):
        """Add MS1 spectra associated with mass features to the dataset.

        Parameters
        -----------
        auto_process : bool, optional
            If True, auto-processes the MS1 spectra before adding it to the object's _ms dictionary. Default is True.
        use_parser : bool, optional
            If True, envoke the spectra parser to get the MS1 spectra. Default is True.
        spectrum_mode : str or None, optional
            The spectrum mode to use for the mass spectra.  If None, method will use the spectrum mode
            from the spectra parser to ascertain the spectrum mode (this allows for mixed types).
            Defaults to None. (faster if defined, otherwise will check each scan)

        Raises
        ------
        ValueError
            If mass_features is not set, must run find_mass_features() first.
            If apex scans are not profile mode, all apex scans must be profile mode for averaging.
            If number of scans to average is not  1 or an integer with an integer median (i.e. 3, 5, 7, 9).
            If deconvolute is True and no EICs are found, did you run integrate_mass_features() first?
        """
        # Check if mass_features is set, raise error if not
        if self.mass_features is None:
            raise ValueError(
                "mass_features not set, must run find_mass_features() first"
            )
        scans_to_average = self.parameters.lc_ms.ms1_scans_to_average

        if scans_to_average == 1:
            # Add to LCMSobj
            self.add_mass_spectra(
                scan_list=[
                    int(x) for x in self.mass_features_to_df().apex_scan.tolist()
                ],
                auto_process=auto_process,
                use_parser=use_parser,
                spectrum_mode=spectrum_mode,
                ms_params=self.parameters.mass_spectrum["ms1"],
            )

        elif (
            (scans_to_average - 1) % 2
        ) == 0:  # scans_to_average = 3, 5, 7 etc, mirror l/r around apex
            apex_scans = list(set(self.mass_features_to_df().apex_scan.tolist()))
            # Check if all apex scans are profile mode, raise error if not
            if not all(self.scan_df.loc[apex_scans, "ms_format"] == "profile"):
                raise ValueError("All apex scans must be profile mode for averaging")

            # First get sets of scans to average
            def get_scans_from_apex(ms1_scans, apex_scan, scans_to_average):
                ms1_idx_start = ms1_scans.index(apex_scan) - int(
                    (scans_to_average - 1) / 2
                )
                if ms1_idx_start < 0:
                    ms1_idx_start = 0
                ms1_idx_end = (
                    ms1_scans.index(apex_scan) + int((scans_to_average - 1) / 2) + 1
                )
                if ms1_idx_end > (len(ms1_scans) - 1):
                    ms1_idx_end = len(ms1_scans) - 1
                scan_list = ms1_scans[ms1_idx_start:ms1_idx_end]
                return scan_list

            ms1_scans = self.ms1_scans
            scans_lists = [
                get_scans_from_apex(ms1_scans, apex_scan, scans_to_average)
                for apex_scan in apex_scans
            ]

            # set polarity to -1 if negative mode, 1 if positive mode (for mass spectrum creation)
            if self.polarity == "negative":
                polarity = -1
            elif self.polarity == "positive":
                polarity = 1

            if not use_parser:
                # Perform checks and prepare _ms_unprocessed dictionary if use_parser is False (saves time to do this once)
                ms1_unprocessed = self._ms_unprocessed[1].copy()
                # Set the index on _ms_unprocessed[1] to scan number
                ms1_unprocessed = ms1_unprocessed.set_index("scan", drop=False)
                self._ms_unprocessed[1] = ms1_unprocessed

                # Check that all the scans in scan_lists are indexs in self._ms_unprocessed[1]
                scans_lists_flat = list(
                    set([scan for sublist in scans_lists for scan in sublist])
                )
                if (
                    len(
                        np.setdiff1d(
                            np.sort(scans_lists_flat),
                            np.sort(ms1_unprocessed.index.values),
                        )
                    )
                    > 0
                ):
                    raise ValueError(
                        "Not all scans to average are present in the unprocessed data"
                    )

            for scan_list_average, apex_scan in zip(scans_lists, apex_scans):
                # Get unprocessed mass spectrum from scans
                ms = self.get_average_mass_spectrum(
                    scan_list=scan_list_average,
                    apex_scan=apex_scan,
                    spectrum_mode="profile",
                    ms_level=1,
                    auto_process=auto_process,
                    use_parser=use_parser,
                    perform_checks=False,
                    polarity=polarity,
                    ms_params=self.parameters.mass_spectrum["ms1"],
                )
                # Add mass spectrum to LCMS object and associated with mass feature
                self.add_mass_spectrum(ms)

            if not use_parser:
                # Reset the index on _ms_unprocessed[1] to not be scan number
                ms1_unprocessed = ms1_unprocessed.reset_index(drop=True)
                self._ms_unprocessed[1] = ms1_unprocessed
        else:
            raise ValueError(
                "Number of scans to average must be 1 or an integer with an integer median (i.e. 3, 5, 7, 9)"
            )

        # Associate the ms1 spectra with the mass features
        for mf_id in self.mass_features:
            self.mass_features[mf_id].mass_spectrum = self._ms[
                self.mass_features[mf_id].apex_scan
            ]
            self.mass_features[mf_id].update_mz()

        # Re-process clustering if persistent homology is selected to remove duplicate mass features after adding and processing MS1 spectra
        if self.parameters.lc_ms.peak_picking_method == "persistent homology":
            self.cluster_mass_features(drop_children=True, sort_by="persistence")

    def mass_features_to_df(self):
        """Returns a pandas dataframe summarizing the mass features.

        The dataframe contains the following columns: mf_id, mz, apex_scan, scan_time, intensity,
        persistence, area, monoisotopic_mf_id, and isotopologue_type.  The index is set to mf_id (mass feature ID).


        Returns
        --------
        pandas.DataFrame
            A pandas dataframe of mass features with the following columns:
            mf_id, mz, apex_scan, scan_time, intensity, persistence, area.
        """

        def mass_spectrum_to_string(
            mass_spec, normalize=True, min_normalized_abun=0.01
        ):
            """Converts a mass spectrum to a string of m/z:abundance pairs.

            Parameters
            -----------
            mass_spec : MassSpectrum
                A MassSpectrum object to be converted to a string.
            normalize : bool, optional
                If True, normalizes the abundance values to a maximum of 1. Defaults to True.
            min_normalized_abun : float, optional
                The minimum normalized abundance value to include in the string, only used if normalize is True. Defaults to 0.01.

            Returns
            --------
            str
                A string of m/z:abundance pairs from the mass spectrum, separated by a semicolon.
            """
            mz_np = mass_spec.to_dataframe()["m/z"].values
            abun_np = mass_spec.to_dataframe()["Peak Height"].values
            if normalize:
                abun_np = abun_np / abun_np.max()
            mz_abun = np.column_stack((mz_np, abun_np))
            if normalize:
                mz_abun = mz_abun[mz_abun[:, 1] > min_normalized_abun]
            mz_abun_str = [
                str(round(mz, ndigits=4)) + ":" + str(round(abun, ndigits=2))
                for mz, abun in mz_abun
            ]
            return "; ".join(mz_abun_str)

        cols_in_df = [
            "id",
            "_apex_scan",
            "start_scan",
            "final_scan",
            "_retention_time",
            "_intensity",
            "_persistence",
            "_area",
            "_dispersity_index",
            "_tailing_factor",
            "monoisotopic_mf_id",
            "isotopologue_type",
            "mass_spectrum_deconvoluted_parent",
        ]
        df_mf_list = []
        for mf_id in self.mass_features.keys():
            # Find cols_in_df that are in single_mf
            df_keys = list(
                set(cols_in_df).intersection(self.mass_features[mf_id].__dir__())
            )
            dict_mf = {}
            for key in df_keys:
                dict_mf[key] = getattr(self.mass_features[mf_id], key)
            if len(self.mass_features[mf_id].ms2_scan_numbers) > 0:
                # Add MS2 spectra info
                best_ms2_spectrum = self.mass_features[mf_id].best_ms2
                dict_mf["ms2_spectrum"] = mass_spectrum_to_string(best_ms2_spectrum)
            if len(self.mass_features[mf_id].associated_mass_features_deconvoluted) > 0:
                dict_mf["associated_mass_features"] = ", ".join(
                    map(
                        str,
                        self.mass_features[mf_id].associated_mass_features_deconvoluted,
                    )
                )
            if self.mass_features[mf_id]._half_height_width is not None:
                dict_mf["half_height_width"] = self.mass_features[
                    mf_id
                ].half_height_width
            # Check if EIC for mass feature is set
            df_mf_single = pd.DataFrame(dict_mf, index=[mf_id])
            df_mf_single["mz"] = self.mass_features[mf_id].mz
            df_mf_list.append(df_mf_single)
        df_mf = pd.concat(df_mf_list)

        # rename _area to area and id to mf_id
        df_mf = df_mf.rename(
            columns={
                "_area": "area",
                "id": "mf_id",
                "_apex_scan": "apex_scan",
                "_retention_time": "scan_time",
                "_intensity": "intensity",
                "_persistence": "persistence",
                "_dispersity_index": "dispersity_index",
                "_tailing_factor": "tailing_factor",
            }
        )

        # reorder columns
        col_order = [
            "mf_id",
            "scan_time",
            "mz",
            "apex_scan",
            "start_scan",
            "final_scan",
            "intensity",
            "persistence",
            "area",
            "half_height_width",
            "tailing_factor",
            "dispersity_index",
            "monoisotopic_mf_id",
            "isotopologue_type",
            "mass_spectrum_deconvoluted_parent",
            "associated_mass_features",
            "ms2_spectrum",
        ]
        # drop columns that are not in col_order
        cols_to_order = [col for col in col_order if col in df_mf.columns]
        df_mf = df_mf[cols_to_order]

        # reset index to mf_id
        df_mf = df_mf.set_index("mf_id")
        df_mf.index.name = "mf_id"

        return df_mf

    def mass_features_ms1_annot_to_df(self):
        """Returns a pandas dataframe summarizing the MS1 annotations for the mass features in the dataset.

        Returns
        --------
        pandas.DataFrame
            A pandas dataframe of MS1 annotations for the mass features in the dataset.
            The index is set to mf_id (mass feature ID)

        Raises
        ------
        Warning
            If no MS1 annotations were found for the mass features in the dataset.
        """
        annot_df_list_ms1 = []
        for mf_id in self.mass_features.keys():
            if self.mass_features[mf_id].mass_spectrum is None:
                pass
            else:
                # Add ms1 annotations to ms1 annotation list
                if (
                    np.abs(
                        (
                            self.mass_features[mf_id].ms1_peak.mz_exp
                            - self.mass_features[mf_id].mz
                        )
                    )
                    < 0.01
                ):
                    # Get the molecular formula from the mass spectrum
                    annot_df = self.mass_features[mf_id].mass_spectrum.to_dataframe()
                    # Subset to pull out only the peak associated with the mass feature
                    annot_df = annot_df[
                        annot_df["Index"] == self.mass_features[mf_id].ms1_peak.index
                    ].copy()

                    # Remove the index column and add column for mf_id
                    annot_df = annot_df.drop(columns=["Index"])
                    annot_df["mf_id"] = mf_id
                    annot_df_list_ms1.append(annot_df)

        if len(annot_df_list_ms1) > 0:
            annot_ms1_df_full = pd.concat(annot_df_list_ms1)
            annot_ms1_df_full = annot_ms1_df_full.set_index("mf_id")
            annot_ms1_df_full.index.name = "mf_id"

        else:
            annot_ms1_df_full = None
            # Warn that no ms1 annotations were found
            warnings.warn(
                "No MS1 annotations found for mass features in dataset, were MS1 spectra added and processed within the dataset?",
                UserWarning,
            )

        return annot_ms1_df_full

    def mass_features_ms2_annot_to_df(self, molecular_metadata=None):
        """Returns a pandas dataframe summarizing the MS2 annotations for the mass features in the dataset.

        Parameters
        -----------
        molecular_metadata :  dict of MolecularMetadata objects
            A dictionary of MolecularMetadata objects, keyed by metabref_mol_id.  Defaults to None.

        Returns
        --------
        pandas.DataFrame
            A pandas dataframe of MS2 annotations for the mass features in the dataset,
            and optionally molecular metadata. The index is set to mf_id (mass feature ID)

        Raises
        ------
        Warning
            If no MS2 annotations were found for the mass features in the dataset.
        """
        annot_df_list_ms2 = []
        for mf_id in self.mass_features.keys():
            if len(self.mass_features[mf_id].ms2_similarity_results) > 0:
                # Add ms2 annotations to ms2 annotation list
                for result in self.mass_features[mf_id].ms2_similarity_results:
                    annot_df_ms2 = result.to_dataframe()
                    annot_df_ms2["mf_id"] = mf_id
                    annot_df_list_ms2.append(annot_df_ms2)

        if len(annot_df_list_ms2) > 0:
            annot_ms2_df_full = pd.concat(annot_df_list_ms2)
            if molecular_metadata is not None:
                molecular_metadata_df = pd.concat(
                    [
                        pd.DataFrame.from_dict(v.__dict__, orient="index").transpose()
                        for k, v in molecular_metadata.items()
                    ],
                    ignore_index=True,
                )
                molecular_metadata_df = molecular_metadata_df.rename(
                    columns={"id": "ref_mol_id"}
                )
                annot_ms2_df_full = annot_ms2_df_full.merge(
                    molecular_metadata_df, on="ref_mol_id", how="left"
                )
            annot_ms2_df_full = annot_ms2_df_full.drop_duplicates(
                subset=["mf_id", "query_spectrum_id", "ref_ms_id"]
            ).copy()
            annot_ms2_df_full = annot_ms2_df_full.set_index("mf_id")
            annot_ms2_df_full.index.name = "mf_id"
        else:
            annot_ms2_df_full = None
            # Warn that no ms2 annotations were found
            warnings.warn(
                "No MS2 annotations found for mass features in dataset, were MS2 spectra added and searched against a database?",
                UserWarning,
            )

        return annot_ms2_df_full

    def __len__(self):
        """
        Returns the number of mass spectra in the dataset.

        Returns
        --------
        int
            The number of mass spectra in the dataset.
        """
        return len(self._ms)

    def __getitem__(self, scan_number):
        """
        Returns the mass spectrum corresponding to the specified scan number.

        Parameters
        -----------
        scan_number : int
            The scan number of the desired mass spectrum.

        Returns
        --------
        MassSpectrum
            The mass spectrum corresponding to the specified scan number.
        """
        return self._ms.get(scan_number)

    def __iter__(self):
        """Returns an iterator over the mass spectra in the dataset.

        Returns
        --------
        iterator
            An iterator over the mass spectra in the dataset.
        """
        return iter(self._ms.values())

    def set_tic_list_from_data(self, overwrite=False):
        """Sets the TIC list from the mass spectrum objects within the _ms dictionary.

        Parameters
        -----------
        overwrite : bool, optional
            If True, overwrites the TIC list if it is already set. Defaults to False.

        Notes
        -----
        If the _ms dictionary is incomplete, sets the TIC list to an empty list.

        Raises
        ------
        ValueError
            If no mass spectra are found in the dataset.
            If the TIC list is already set and overwrite is False.
        """
        # Check if _ms is empty and raise error if so
        if len(self._ms) == 0:
            raise ValueError("No mass spectra found in dataset")

        # Check if tic_list is already set and raise error if so
        if len(self.tic) > 0 and not overwrite:
            raise ValueError("TIC list already set, use overwrite=True to overwrite")

        self.tic = [self._ms.get(i).tic for i in self.scans_number]

    def set_retention_time_from_data(self, overwrite=False):
        """Sets the retention time list from the data in the _ms dictionary.

        Parameters
        -----------
        overwrite : bool, optional
            If True, overwrites the retention time list if it is already set. Defaults to False.

        Notes
        -----
        If the _ms dictionary is empty or incomplete, sets the retention time list to an empty list.

        Raises
        ------
        ValueError
            If no mass spectra are found in the dataset.
            If the retention time list is already set and overwrite is False.
        """
        # Check if _ms is empty and raise error if so
        if len(self._ms) == 0:
            raise ValueError("No mass spectra found in dataset")

        # Check if retention_time_list is already set and raise error if so
        if len(self.retention_time) > 0 and not overwrite:
            raise ValueError(
                "Retention time list already set, use overwrite=True to overwrite"
            )

        retention_time_list = []
        for key_ms in sorted(self._ms.keys()):
            retention_time_list.append(self._ms.get(key_ms).retention_time)
        self.retention_time = retention_time_list

    def set_scans_number_from_data(self, overwrite=False):
        """Sets the scan number list from the data in the _ms dictionary.

        Notes
        -----
        If the _ms dictionary is empty or incomplete, sets the scan number list to an empty list.

        Raises
        ------
        ValueError
            If no mass spectra are found in the dataset.
            If the scan number list is already set and overwrite is False.
        """
        # Check if _ms is empty and raise error if so
        if len(self._ms) == 0:
            raise ValueError("No mass spectra found in dataset")

        # Check if scans_number_list is already set and raise error if so
        if len(self.scans_number) > 0 and not overwrite:
            raise ValueError(
                "Scan number list already set, use overwrite=True to overwrite"
            )

        self.scans_number = sorted(self._ms.keys())

    @property
    def ms1_scans(self):
        """
        list : A list of MS1 scan numbers for the dataset.
        """
        return self.scan_df[self.scan_df.ms_level == 1].index.tolist()

    @property
    def parameters(self):
        """
        LCMSParameters : The parameters used for the LC-MS analysis.
        """
        return self._parameters

    @parameters.setter
    def parameters(self, paramsinstance):
        """
        Sets the parameters used for the LC-MS analysis.

        Parameters
        -----------
        paramsinstance : LCMSParameters
            The parameters used for the LC-MS analysis.
        """
        self._parameters = paramsinstance

    @property
    def scans_number(self):
        """
        list : A list of scan numbers for the dataset.
        """
        return self._scans_number_list

    @scans_number.setter
    def scans_number(self, scan_numbers_list):
        """
        Sets the scan numbers for the dataset.

        Parameters
        -----------
        scan_numbers_list : list
            A list of scan numbers for the dataset.
        """
        self._scans_number_list = scan_numbers_list

    @property
    def retention_time(self):
        """
        numpy.ndarray : An array of retention times for the dataset.
        """
        return self._retention_time_list

    @retention_time.setter
    def retention_time(self, rt_list):
        """
        Sets the retention times for the dataset.

        Parameters
        -----------
        rt_list : list
            A list of retention times for the dataset.
        """
        self._retention_time_list = np.array(rt_list)

    @property
    def tic(self):
        """
        numpy.ndarray : An array of TIC values for the dataset.
        """
        return self._tic_list

    @tic.setter
    def tic(self, tic_list):
        """
        Sets the TIC values for the dataset.

        Parameters
        -----------
        tic_list : list
            A list of TIC values for the dataset.
        """
        self._tic_list = np.array(tic_list)
