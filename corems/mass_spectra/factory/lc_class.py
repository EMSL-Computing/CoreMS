from pathlib import Path

import numpy as np
import pandas as pd
import warnings
import multiprocessing

import matplotlib.pyplot as plt

from corems.encapsulation.factory.parameters import LCMSParameters, LCMSCollectionParameters
from corems.mass_spectra.calc.lc_calc import LCCalculations, PHCalculations, LCMSCollectionCalculations
from corems.molecular_id.search.lcms_spectral_search import LCMSSpectralSearch
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile, ms_from_array_centroid
from corems.mass_spectra.calc.lc_calc import find_closest
from corems.chroma_peak.factory.chroma_peak_classes import LCMSMassFeature


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
        self._raw_file_location = None

        # Add the spectra parser class to the object if it is not None
        if spectra_parser is not None:
            self.spectra_parser_class = spectra_parser.__class__
            if self.spectra_parser_class.__name__ == "ReadCoreMSHDFMassSpectra":
                self.raw_file_location = spectra_parser.get_raw_file_location()

            # Check that spectra_parser.sample_name is same as sample_name etc, raise warning if not
            if (
                self.sample_name is not None
                and self.sample_name != self.spectra_parser.sample_name
                and self.spectra_parser_class.__name__ != "ReadCoreMSHDFMassSpectra"
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

    @property
    def spectra_parser(self):
        """Returns an instance of the spectra parser class."""
        # Check if a file exists at the raw_file_location
        if not Path(self.raw_file_location).exists():
            raise FileNotFoundError(
                f"Raw file not found at location: {self.raw_file_location}, update raw_file_location property to point to correct location."
            )        
        return self.spectra_parser_class(self.raw_file_location)

    @property
    def raw_file_location(self):
        """Returns the file_location unless the _raw_file_location is not None."""
        return self._raw_file_location if self._raw_file_location is not None else self.file_location   
    
    @raw_file_location.setter
    def raw_file_location(self, value):
        self._raw_file_location = value

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

        if use_parser:
            # Use batch function to get all mass spectra at once
            if spectrum_mode is None:
                # get spectrum mode from _scan_info for each scan
                spectrum_modes = [self.scan_df.loc[scan, "ms_format"] for scan in scan_list]
                spectrum_mode_batch = spectrum_modes[0] if len(set(spectrum_modes)) == 1 else None
            else:
                spectrum_mode_batch = spectrum_mode
            
            ms_list = self.spectra_parser.get_mass_spectra_from_scan_list(
                scan_list=scan_list,
                spectrum_mode=spectrum_mode_batch,
                auto_process=False,
            )
            
            # Process each mass spectrum
            for i, scan in enumerate(scan_list):
                ms = ms_list[i] if i < len(ms_list) else None
                if ms is not None:
                    if ms_params is not None:
                        ms.parameters = ms_params
                    ms.scan_number = scan
                    if auto_process:
                        ms.process_mass_spec()
                    self.add_mass_spectrum(ms)
        else:
            # Original non-parser logic remains unchanged
            for scan in scan_list:
                ms = None
                if spectrum_mode is None:
                    # get spectrum mode from _scan_info
                    spectrum_mode_scan = self.scan_df.loc[scan, "ms_format"]
                else:
                    spectrum_mode_scan = spectrum_mode
                
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
                   ms = ms_from_array_centroid(
                        mz = my_ms_df.mz,
                        abundance = my_ms_df.intensity,
                        rp = [np.nan] * len(my_ms_df.mz),
                        s2n = [np.nan] * len(my_ms_df.mz),
                        dataname = self.file_location,
                        polarity=polarity,
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
        
    @property
    def ms(self):
        """
        dictionary : contains the key associated with mass spectra and values are the associated MassSpecProfiles
        """
        return self._ms

    
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
    induced_mass_features: dictionary of LCMSMassFeature objects
        A dictionary containing mass features from a collection that don't
        satisfy criteria for initial mass features. Key is mass feature ID.
        Initialized as an empty dictionary.
    missing_mass_features: pandas.DataFrame
        A table of clusters in a given sample for which a mass feature was
        sought and not found
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
    * plot_composite_mz_features(binsize = 1e-4, ph_int_min_thresh = 0.001, mf_plot = True, ms2_plot = True, return_fig = False)
        Generates plot of M/Z features comparing scan time vs M/Z value
    * search_for_targeted_mass_feature(ms1df: pd.DataFrame, sample: pd.Series, tol_flag = 0)
        Searches for mass features in specific M/Z and scan time windows that
        were missed by the persistent homology search
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
        self.induced_mass_features = {}
        self.spectral_search_results = {}

    def get_eic_mz_for_mass_feature(self, mf_mz, tolerance=0.0001):
        """Get the EIC dictionary key (m/z) that best matches a mass feature's m/z.
        
        Finds the closest EIC m/z key within the specified tolerance.
        
        Parameters
        ----------
        mf_mz : float
            The m/z value of the mass feature to match.
        tolerance : float, optional
            Maximum m/z difference for matching. Default is 0.0001 Da.
            
        Returns
        -------
        float or None
            The EIC dictionary key (m/z) of the closest matching EIC,
            or None if no EIC is within tolerance.
        """
        if not hasattr(self, 'eics') or not self.eics:
            return None
        
        best_eic_mz = None
        best_diff = tolerance
        for eic_mz in self.eics.keys():
            diff = abs(mf_mz - eic_mz)
            if diff < best_diff:
                best_diff = diff
                best_eic_mz = eic_mz
        return best_eic_mz
    
    def associate_eics_with_mass_features(self, tolerance=0.0001, induced=False):
        """Associate EICs with mass features using tolerance-based m/z matching.
        
        Associates EIC_Data objects from self.eics with mass features by finding
        the closest EIC within the specified m/z tolerance. This is more robust
        than exact matching which can fail due to floating point precision issues.
        
        Parameters
        ----------
        tolerance : float, optional
            Maximum m/z difference for matching EICs to mass features. Default is 0.0001 Da.
        induced : bool, optional
            If True, associates EICs with induced_mass_features instead of mass_features.
            Default is False.
            
        Notes
        -----
        For each mass feature, this method finds the EIC with the closest m/z value
        within the tolerance window and assigns it to the mass feature's _eic_data attribute.
        If multiple EICs are within tolerance, the one with the smallest m/z difference is chosen.
        """
        # Select which mass features dictionary to use
        mf_dict = self.induced_mass_features if induced else self.mass_features
        
        # Use the _eic_mz attribute on each mass_feature to find the closest matching EIC
        for idx in mf_dict.keys():
            mf_mz = mf_dict[idx]._eic_mz
            # Find closest EIC within tolerance
            best_match = None
            best_diff = tolerance
            for eic_mz, eic_data in self.eics.items():
                diff = abs(mf_mz - eic_mz)
                if diff < best_diff:
                    best_diff = diff
                    best_match = eic_data
            if best_match is not None:
                mf_dict[idx]._eic_data = best_match

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

    def _find_ms2_scans_for_mass_features(self, mf_ids=None, scan_filter=None):
        """Find MS2 scans associated with mass features.
        
        This helper method finds MS2 scans that match mass features based on RT and m/z tolerances.
        It updates the ms2_scan_numbers attribute on each mass feature.
        
        Parameters
        ----------
        mf_ids : list of int, optional
            List of mass feature IDs to find MS2 for. If None, finds for all mass features.
        scan_filter : str, optional
            Filter string for MS2 scans (e.g., 'hcd'). Default is None.
            
        Returns
        -------
        list
            List of unique MS2 scan numbers found across all mass features.
            
        Raises
        ------
        ValueError
            If no MS2 scans are found in the dataset.
        """
        # Get mass features to process
        if mf_ids is None:
            mf_ids = list(self.mass_features.keys())
        
        # Get mass features dataframe
        mf_df = self.mass_features_to_df()
        mf_df = mf_df.loc[mf_ids].copy()
        
        # Find ms2 scans that have a precursor m/z value
        ms2_scans = self.scan_df[self.scan_df.ms_level == 2]
        ms2_scans = ms2_scans[~ms2_scans.precursor_mz.isna()]
        ms2_scans = ms2_scans[ms2_scans.tic > 0]
        
        if len(ms2_scans) == 0:
            raise ValueError("No DDA scans found in dataset")
        
        if scan_filter is not None:
            ms2_scans = ms2_scans[ms2_scans.scan_text.str.contains(scan_filter)]
        
        # Get tolerances from parameters
        time_tol = self.parameters.lc_ms.ms2_dda_rt_tolerance
        mz_tol = self.parameters.lc_ms.ms2_dda_mz_tolerance
        
        # For each mass feature, find the ms2 scans that are within the roi scan time and mz range
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
            scan_list = ms2_scans_filtered.scan.tolist()
            if scan_list:
                self.mass_features[i].ms2_scan_numbers = (
                    scan_list + list(self.mass_features[i].ms2_scan_numbers)
                )
                dda_scans.extend(scan_list)
        
        return list(set(dda_scans))
    
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

        # Find MS2 scans for all mass features
        dda_scans = self._find_ms2_scans_for_mass_features(scan_filter=scan_filter)
        
        # Load MS2 spectra
        self.add_mass_spectra(
            scan_list=dda_scans,
            auto_process=auto_process,
            spectrum_mode=spectrum_mode,
            use_parser=use_parser,
            ms_params=ms_params,
        )
        
        # Associate appropriate _ms attribute to appropriate mass feature's ms2_mass_spectra attribute
        for mf_id in self.mass_features:
            if self.mass_features[mf_id].ms2_scan_numbers is not None:
                for dda_scan in self.mass_features[mf_id].ms2_scan_numbers:
                    if dda_scan in self._ms.keys():
                        self.mass_features[mf_id].ms2_mass_spectra[dda_scan] = self._ms[
                            dda_scan
                        ]

    def add_associated_ms1(
        self, auto_process=True, use_parser=True, spectrum_mode=None, induced_features=False
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
        induced_features : bool, optional
            If True, add associated MS1 of the induced mass features instead of the primary mass features

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
            
        if induced_features:
            mf_dict = self.induced_mass_features
        else:
            mf_dict = self.mass_features

        scans_to_average = self.parameters.lc_ms.ms1_scans_to_average
        
        ## sketchy work around for induced mass features
        scan_list = [
            int(mf_dict[x].apex_scan) for x in mf_dict if int(mf_dict[x].apex_scan) != -99
        ]

        if scans_to_average == 1:
            # Add to LCMSobj
            self.add_mass_spectra(
                scan_list = scan_list,
                auto_process=auto_process,
                use_parser=use_parser,
                spectrum_mode=spectrum_mode,
                ms_params=self.parameters.mass_spectrum["ms1"],
            )

        elif (
            (scans_to_average - 1) % 2
        ) == 0:  # scans_to_average = 3, 5, 7 etc, mirror l/r around apex
            apex_scans = list(set(scan_list))
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
        for k in mf_dict.keys():
            ## another induced feature work around
            if mf_dict[k].apex_scan != -99:
                mf_dict[k].mass_spectrum = self._ms[
                    mf_dict[k].apex_scan
                ]
                mf_dict[k].update_mz()

    def mass_features_to_df(self, induced_features=False, drop_na_cols=False, include_cols=None):
        """Returns a pandas dataframe summarizing the mass features.

        The dataframe contains the following columns: mf_id, mz, apex_scan, scan_time, intensity,
        persistence, area, monoisotopic_mf_id, and isotopologue_type.  The index is set to mf_id (mass feature ID).
        Parameters
        -----------
        induced_features : bool, optional
            If True, calls the induced_mass_features dictionary. Defaults to False.
        drop_na_cols : bool, optional
            If True, drops columns that are entirely NA. Defaults to False.
        include_cols : list of str, optional
            If provided, only includes the specified columns in the output (in addition to 'mf_id' which is always included as the index).
            If None, includes all available columns. Defaults to None.

        Raises
        --------
        ValueError
            If the sample provided doesn't contain the mass feature data.        

        Returns
        --------
        pandas.DataFrame
            A pandas dataframe of mass features with the following columns:
            mf_id, mz, apex_scan, scan_time, intensity, persistence, area.
        """
        import pandas as pd

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

        if induced_features:
            mf_dict = self.induced_mass_features
        else:
            mf_dict = self.mass_features
        
        if len(mf_dict) == 0:
            # Return an empty dataframe with the expected structure
            # This allows collection processing to continue even if some samples have no features
            return pd.DataFrame()
            
        cols_in_df = [
            "id",
            "apex_scan",
            "start_scan",
            "final_scan",
            "retention_time",
            "intensity",
            "persistence",
            "area",
            "dispersity_index",
            "normalized_dispersity_index",
            "tailing_factor",
            "gaussian_similarity",
            "noise_score",
            "noise_score_min",
            "noise_score_max",
            "monoisotopic_mf_id",
            "isotopologue_type",
            "mass_spectrum_deconvoluted_parent",
            "ms2_scan_numbers",
            "type"
        ]

        df_mf_list = []
        for mf_id in mf_dict.keys():
            # Find cols_in_df that are in single_mf
            df_keys = list(
                set(cols_in_df).intersection(mf_dict[mf_id].__dir__())
            )
            dict_mf = {}
            # Get the values for each key in df_keys from the mass feature object
            for key in df_keys:
                value = getattr(mf_dict[mf_id], key)
                # Wrap list/array values in a list so pandas treats them as single cell values
                if key == 'ms2_scan_numbers' and isinstance(value, (list, np.ndarray)):
                    dict_mf[key] = [value]
                else:
                    dict_mf[key] = value
            if len(mf_dict[mf_id].ms2_scan_numbers) > 0:
                # Add MS2 spectra info
                best_ms2_spectrum = mf_dict[mf_id].best_ms2
                if best_ms2_spectrum is not None:
                    dict_mf["ms2_spectrum"] = mass_spectrum_to_string(best_ms2_spectrum)
            if len(mf_dict[mf_id].associated_mass_features_deconvoluted) > 0:
                dict_mf["associated_mass_features"] = ", ".join(
                    map(
                        str,
                        mf_dict[mf_id].associated_mass_features_deconvoluted,
                    )
                )
            if mf_dict[mf_id]._half_height_width is not None:
                dict_mf["half_height_width"] = mf_dict[
                    mf_id
                ].half_height_width
            # Check if EIC for mass feature is set
            df_mf_single = pd.DataFrame(dict_mf, index=[mf_id])
            df_mf_single["mz"] = mf_dict[mf_id].mz
            df_mf_list.append(df_mf_single)
        df_mf = pd.concat(df_mf_list)

        # rename _area to area and id to mf_id
        df_mf = df_mf.rename(
            columns={
                "id": "mf_id",
                "retention_time": "scan_time",            
            }
        )

        # reorder columns
        col_order = [
            "mf_id",
            "type",
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
            "normalized_dispersity_index",
            "gaussian_similarity",
            "noise_score",
            "noise_score_min",
            "noise_score_max",
            "monoisotopic_mf_id",
            "isotopologue_type",
            "mass_spectrum_deconvoluted_parent",
            "associated_mass_features",
            "ms2_scan_numbers",
            "ms2_spectrum",
        ]
        # drop columns that are not in col_order
        cols_to_order = [col for col in col_order if col in df_mf.columns]
        df_mf = df_mf[cols_to_order]

        # reset index to mf_id
        df_mf = df_mf.set_index("mf_id")
        df_mf.index.name = "mf_id"
        
        if 'half_height_width' in df_mf.columns:
            df_mf["half_height_width"] = df_mf["half_height_width"].astype('float64')
        if 'tailing_factor' in df_mf.columns:
            df_mf["tailing_factor"] = df_mf["tailing_factor"].astype('float64')
        if 'dispersity_index' in df_mf.columns:
            df_mf["dispersity_index"] = df_mf["dispersity_index"].astype('float64')
        if 'normalized_dispersity_index' in df_mf.columns:
            df_mf["normalized_dispersity_index"] = df_mf["normalized_dispersity_index"].astype('float64')
        
        # Filter columns if include_cols is specified
        if include_cols is not None:
            # Ensure include_cols is a list
            if not isinstance(include_cols, list):
                raise ValueError("include_cols must be a list of column names")
            # Keep only requested columns that exist in the dataframe
            available_cols = [col for col in include_cols if col in df_mf.columns]
            df_mf = df_mf[available_cols]
        
        # Drop columns that are entirely NA if requested
        if drop_na_cols:
            df_mf = df_mf.dropna(axis=1, how='all')
        
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

                    # If there are more than 1 row, remove any rows without a molecular formula
                    if len(annot_df) > 1:
                        annot_df = annot_df[~annot_df["Molecular Formula"].isna()]

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

    def plot_composite_mz_features(self, binsize = 1e-4, ph_int_min_thresh = 0.001, mf_plot = True, ms2_plot = True, return_fig = False):
        """Returns a figure displaying 
            (1) thresholded, unprocessed data
            (2) the m/z features
            (3) which m/z features are associated with MS2 spectra

        Parameters
        -----------
        binsize :  float
            Desired binsize for the m/z axis of the composite feature map.  Defaults to 1e-4.
        mf_plot : boolean
            Indicates whether to plot the m/z features. Defaults to True.
        ms2_plot : boolean
            Indicates whether to identify m/z features with associated MS2 spectra. Defaults to True.
        return_fig : boolean
            Indicates whether to plot composite feature map (False) or return figure object (True). Defaults to False.

        Returns
        --------
        matplotlib.pyplot.Figure
            A figure with the thresholded, unprocessed data on an axis of m/z value with respect to 
            scan time. Unprocessed data is displayed in gray scale with darker colors indicating 
            higher intensities. If m/z features are plotted, they are displayed in cyan. If m/z
            features with associated with MS2 spectra are plotted, they are displayed in red.

        Raises
        ------
        Warning
            If m/z features are set to be plot but aren't in the dataset.
            If m/z features with associated MS2 data are set to be plot but no MS2 annotations 
            were found for the m/z features in the dataset.
        """
        if mf_plot:
            # Check if mass_features is set, raise error if not
            if self.mass_features is None:
                raise ValueError(
                    "mass_features not set, must run find_mass_features() first"
                )
            ## call mass feature data
            mf_df = self.mass_features_to_df()

        if ms2_plot:
            if not mf_plot:
                # Check if mass_features is set, raise error if not
                if self.mass_features is None:
                    raise ValueError(
                        "mass_features not set, must run find_mass_features() first"
                    )

            ## call m/z feature data
            mf_df = self.mass_features_to_df()

            # Check if ms2_spectrum is set, raise error if not
            if 'ms2_spectrum' not in mf_df.columns:
                raise ValueError(                
                    "ms2_spectrum not set, must run add_associated_ms2_dda() first"            
                )
    
        ## threshold and grid unprocessed data
        df = self._ms_unprocessed[1].copy()
        df = df.dropna(subset=['intensity']).reset_index(drop = True)
        threshold = ph_int_min_thresh * df.intensity.max()
        df_thres = df[df["intensity"] > threshold].reset_index(drop = True).copy()
        df = self.grid_data(df_thres)
    
        ## format unprocessed data for plotting
        df = df.merge(self.scan_df[['scan', 'scan_time']], on = 'scan')
        mz_grid = np.arange(0, np.max(df.mz), binsize)
        mz_data = np.array(df.mz)
        df['mz_bin'] = find_closest(mz_grid, mz_data)
        df['ab_bin'] = df.groupby(['mz_bin', 'scan_time']).intensity.transform(sum)
        unproc_df = df[['scan_time', 'mz_bin', 'ab_bin']].drop_duplicates(ignore_index = True)

        ## generate figure
        fig = plt.figure()
        plt.scatter(
            unproc_df.scan_time,
            unproc_df.mz_bin*binsize,
            c = unproc_df.ab_bin/np.max(unproc_df.ab_bin),
            alpha = unproc_df.ab_bin/np.max(unproc_df.ab_bin), 
            cmap = 'Greys_r',
            s = 1
        )

        if mf_plot:
            if ms2_plot:
                plt.scatter(
                    mf_df[mf_df.ms2_spectrum.isna()].scan_time,
                    mf_df[mf_df.ms2_spectrum.isna()].mz,
                    c = 'c',
                    s = 4,
                    label = 'M/Z features without MS2'
                )
            else:
                plt.scatter(
                    mf_df.scan_time,
                    mf_df.mz,
                    c = 'c',
                    s = 4,
                    label = 'M/Z features'
                )

        if ms2_plot: 
            plt.scatter(
                mf_df[~mf_df.ms2_spectrum.isna()].scan_time,
                mf_df[~mf_df.ms2_spectrum.isna()].mz,
                c = 'r',
                s = 2,
                label = 'M/Z features with MS2'
            )

        if mf_plot == True or ms2_plot == True:
            plt.legend(loc = 'lower center', bbox_to_anchor = (0.5, -0.25), ncol = 2)
        plt.xlabel('Scan time')
        plt.ylabel('m/z')
        plt.ylim(0, np.ceil(np.max(df.mz)))
        plt.xlim(0, np.ceil(np.max(df.scan_time)))
        plt.title('Composite Feature Map')

        if return_fig:
            plt.close(fig)
            return fig

        else:
            plt.show()
            
    def search_for_targeted_mass_features_batch(
            self,
            ms1df,
            mz_mins,
            mz_maxs,
            st_mins,
            st_maxs,
            set_ids,
            obj_idx=0,
            st_aligned=False
            ):
        """
        Returns multiple LCMSMassFeatures from a specific sample within specific mass and time ranges.
        Vectorized batch version of search_for_targeted_mass_feature for improved performance.

        Parameters
        -----------
        ms1df : pd.DataFrame
            Dataframe containing all the possible MS1 values to consider, collected by calling _ms_unprocessed[1] on the sample.
        mz_mins : np.ndarray
            Array of lower bounds of m/z values to use to find peaks.
        mz_maxs : np.ndarray
            Array of upper bounds of m/z values to use to find peaks.
        st_mins : np.ndarray
            Array of lower bounds of scan times to use to find peaks.
        st_maxs : np.ndarray
            Array of upper bounds of scan times to use to find peaks.
        set_ids : np.ndarray or list
            Array of strings used as IDs in LCMSMassFeatures.
        obj_idx : int
            Identifies index of sample in a collection. Defaults to 0.
        st_aligned : bool
            Whether to use scan_time_aligned or scan_time. Defaults to False.

        Returns
        --------
        dict
            Dictionary mapping set_id to LCMSMassFeature objects.

        Raises
        ------
        ValueError
            If appropriate scan time data is not contained in ms1df or if array lengths don't match.
        """
        # Validate inputs
        n_features = len(mz_mins)
        if not all(len(arr) == n_features for arr in [mz_maxs, st_mins, st_maxs, set_ids]):
            raise ValueError("All input arrays must have the same length")

        # Validate scan time column
        time_col = 'scan_time_aligned' if st_aligned else 'scan_time'
        if time_col not in ms1df.columns:
            raise ValueError(f"{time_col} not contained in ms1df")

        # Pre-extract columns for faster access
        mz_vals = ms1df.mz.values
        st_vals = ms1df[time_col].values
        scan_vals = ms1df.scan.values
        intensity_vals = ms1df.intensity.values

        # Process all features
        results = {}
        for i in range(n_features):
            # Vectorized filtering
            mask = (
                (mz_vals >= mz_mins[i]) & (mz_vals <= mz_maxs[i]) &
                (st_vals >= st_mins[i]) & (st_vals <= st_maxs[i])
            )
            
            if not mask.any():
                row_dict = {
                    'apex_scan': -99,
                    'mz': np.nan,
                    'intensity': np.nan,
                    'retention_time': np.nan,
                    'persistence': np.nan,
                    'id': set_ids[i]
                }
            else:
                # Find max intensity within filtered region
                filtered_intensities = intensity_vals[mask]
                max_idx = np.argmax(filtered_intensities)
                
                # Get indices of filtered data
                filtered_indices = np.where(mask)[0]
                peak_idx = filtered_indices[max_idx]
                
                row_dict = {
                    'apex_scan': scan_vals[peak_idx],
                    'mz': mz_vals[peak_idx],
                    'intensity': intensity_vals[peak_idx],
                    'retention_time': st_vals[peak_idx],
                    'persistence': np.nan,
                    'id': set_ids[i]
                }

            results[set_ids[i]] = LCMSMassFeature(self, **row_dict)

        return results

    def search_for_targeted_mass_feature(
            self,
            ms1df, 
            mz_min,
            mz_max, 
            st_min, 
            st_max,
            set_id,
            obj_idx = 0,
            st_aligned = False
            ):
        """
        Returns an LCMSMassFeature from a specific sample within a specific mass and time range. Returns an empty
        LCMSMassFeature if no satisfactory peak is found in the given window.

        Parameters
        -----------
        ms1df :  Pandas DataFrame
            Dataframe containing all the possible MS1 values to consider, collected by calling _ms_unprocessed[1] on the sample.
        mz_min : float
            Identifies lower bound of the weights to use to find a peak.
        mz_max : float
            Identifies upper bound of the weights to use to find a peak.
        st_min : float
            Identifies lower bound of the scan times to use to find a peak.
        st_max : float
            Identifies upper bound of the scan times to use to find a peak.
        set_id : str
            Indicates string used as ID in LCMSMassFeature.
        obj_idx : int
            Identifies index of sample in a collection that LCMSMassFeature should be assigned to. Defaults to 0 and is not used
            if data provided is an LCMSBase instead of an LCMSCollection.
        st_aligned : boolean
            Indicates whether to call scan time from scan_time or from scan_time_aligned if using a collection. Defaults to False.

        Returns
        --------
        LCMSMassFeature
            Object from ChromaPeak that contains data on selected MS1 peak. If no peak is found, will contain missing 
            information and list the apex scan value as -99.

        Raises
        ------
        Warning
            If appropriate scan time data is not contained in ms1df.
        """
        # Convert single feature to arrays and call batch method
        results = self.search_for_targeted_mass_features_batch(
            ms1df,
            np.array([mz_min]),
            np.array([mz_max]),
            np.array([st_min]),
            np.array([st_max]),
            [set_id],
            obj_idx=obj_idx,
            st_aligned=st_aligned
        )
        return results[set_id]


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

class LCMSCollection(LCMSCollectionCalculations):
    """A class representing a collection of liquid chromatography-mass spectrometry (LC-MS) runs.
    These runs can be from the same or different samples, but must be from the same instrument and have the same parameters 
    for the initial processing steps.  The LCMS objects are stored in an ordered dictionary with the sample name as the key.

    Parameters
    -----------

    Attributes
    -----------

    Methods
    --------

    Notes
    ------
    This class is not intended to be instantiated directly, but rather instantiated using a parser object and then interacted with.
    """

    def __init__(
            self,
            collection_location,
            manifest,
            collection_parser=None
    ):
        self.collection_location = collection_location
        self._manifest_dict = manifest
        self.collection_parser = collection_parser
        self.raw_files_relocated = False

        # These attributes are generally set by the parser during instantiation of this class
        self._lcms = {}
        self._combined_mass_features = None
        self._combined_induced_mass_features = None
        self.consensus_mass_features = {}
        self._parameters = LCMSCollectionParameters()
        self.isotopes_dropped = False
        self._mass_features_locked = False  # Prevents rebuilding mass_features_dataframe from samples

        # These attributes are set during processing
        self.rt_aligned = False
        self.rt_alignment_attempted = False
        self.missing_mass_features_searched = False

    def _reorder_lcms_objects(self):
        """
        Reorders the LCMS objects in the collection based on the order in the manifest.
        """
        ordered_samples = self.samples
        self._lcms = {k: self._lcms[k] for k in ordered_samples}

    def __getitem__(self, index):
        samp_name = self.samples[index]
        self._lcms[samp_name]
        return self._lcms[samp_name]
    
    def __len__(self):
        return len(self.samples)
    
    def _prepare_lcms_mass_features_for_combination(self, lcms_obj, induced_features = False):
        """
        Prepares the mass features in the LCMS objects in the collection for combination.
        """        
        if induced_features == True:
            mf_df = lcms_obj.mass_features_to_df(induced_features = True)
        # Check if lcms_obj has attribute light_mf_df
        elif hasattr(lcms_obj, "light_mf_df"):
            mf_df = lcms_obj.light_mf_df
        else:
            mf_df = lcms_obj.mass_features_to_df()
        
        # If dataframe is empty, add minimal required columns and return
        if len(mf_df) == 0:
            import pandas as pd
            mf_df["sample_name"] = []
            mf_df["sample_id"] = []
            mf_df["coll_mf_id"] = []
            mf_df["mf_id"] = []
            mf_df["_eic_mz"] = []  # Include _eic_mz for consistency with non-empty dataframes
            if induced_features:
                mf_df["cluster"] = []
            return mf_df
        
        # Remove index
        mf_df = mf_df.reset_index(drop=False)
        # Add sample name and sample id to the dataframe
        mf_df["sample_name"] = lcms_obj.sample_name
        mf_df["sample_id"] = self.manifest[lcms_obj.sample_name]["collection_id"]
        mf_df["coll_mf_id"] = mf_df["sample_id"].astype(str) + "_" + mf_df["mf_id"].astype(str)

        # For induced features, extract cluster from mf_id (format: c{cluster}_{index}_i)
        # and add as a column since cluster_index attribute may not be set on the object
        if induced_features:
            def extract_cluster(mf_id):
                # mf_id format: c{cluster}_{index}_i
                # Example: c123_5_i -> cluster 123
                if isinstance(mf_id, str) and mf_id.startswith('c') and '_i' in mf_id:
                    parts = mf_id.split('_')
                    if len(parts) >= 2:
                        cluster_str = parts[0][1:]  # Remove 'c' prefix
                        try:
                            return int(cluster_str)
                        except ValueError:
                            return None
                return None
            
            mf_df['cluster'] = mf_df['mf_id'].apply(extract_cluster)

        # Check if scan_df has scan_time_aligned and add to mf_df if so
        if "scan_time_aligned" in lcms_obj.scan_df.columns:
            scan_df = lcms_obj.scan_df[["scan", "scan_time_aligned"]].copy()
            scan_df = scan_df.rename(columns={"scan": "apex_scan"})
            mf_df = mf_df.merge(scan_df, on="apex_scan")
        
        return mf_df
       
    def _combine_mass_features(self, induced_features = False):
        """
        Concatenates the mass features from all the LCMS objects in the collection.

        Returns
        --------
        None, sets the _combined_mass_features or _combined_induced_mass_feature attribute.
        
        Notes
        -----
        If _mass_features_locked is True (e.g., when only representative features are loaded),
        this method will skip rebuilding the regular mass features dataframe to preserve
        the full collection-level dataframe. Induced features are always rebuilt since they
        are created during processing.
        """
        
        # Skip rebuilding regular mass features if locked (preserves full dataframe)
        if not induced_features and self._mass_features_locked:
            return

        ## TODO: See why this function runs slower on multiprocessing,
        ## especially for induced features
        ## has only been considered so far on ~20 samples
#        if self.parameters.lcms_collection.cores == 1:
#            # Prepare mass features for combination sequentially
#            mf_df_list = []
#            for lcms_obj in self:
#                mf_df = self._prepare_lcms_mass_features_for_combination(lcms_obj, induced_features)
#                mf_df_list.append(mf_df)

#        if self.parameters.lcms_collection.cores > 1:
#            # Parallelize the mass feature preparation
#            if self.parameters.lcms_collection.cores > len(self):
#                ncores = len(self)
#            else:
#                ncores = self.parameters.lcms_collection.cores
#            pool = multiprocessing.Pool(ncores)
#            mf_df_list = pool.starmap(
#                self._prepare_lcms_mass_features_for_combination, 
#                [(lcms_obj, induced_features) for lcms_obj in self]
#            )

        # Prepare mass features for combination sequentially
        mf_df_list = []
        for lcms_obj in self:
            # Skip samples with no induced mass features if processing induced features
            if induced_features:
                if not hasattr(lcms_obj, 'induced_mass_features') or len(lcms_obj.induced_mass_features) == 0:
                    continue
            mf_df = self._prepare_lcms_mass_features_for_combination(lcms_obj, induced_features)
            mf_df_list.append(mf_df)

        # If no mass features were collected (e.g., no induced features exist), return early
        if len(mf_df_list) == 0:
            # Add a warning here, not sure how one might reach this state, clearly saying if they are induced features or not
            warnings.warn("No mass features found to combine in the collection.", UserWarning)
            if induced_features:
                self._combined_induced_mass_features = None
            else:
                self._combined_mass_features = None
            return

        combined_mass_features = pd.concat(mf_df_list)
        # Move coll_mf_id, sample_name, sample_id, and mf_id to front
        cols = combined_mass_features.columns.tolist()
        top_cols = ["coll_mf_id", "sample_name", "sample_id", "mf_id", "mz", "scan_time_aligned", "cluster"]
        cols = [x for x in top_cols + [col for col in cols if col not in top_cols] if x in cols]
        combined_mass_features = combined_mass_features[cols]
        # Make coll_mf_id the index
        combined_mass_features = combined_mass_features.set_index("coll_mf_id")
        if induced_features == True:
            self._combined_induced_mass_features = combined_mass_features
        else:
            self._combined_mass_features = combined_mass_features

    def _check_mass_features_df(self, induced_features = False):
        """Checks if the mass features dataframe has expected columns.  If not, adds them.
        
        Returns
        --------
        pandas.DataFrame
            A pandas dataframe of mass features in the collection.

        Notes
        ------
        If scan_time_aligned is not in the _combined_mass_features or 
        _combined_induced_mass_features, tries to add it.

        """
        
        if induced_features:
            cmf_df = self._combined_induced_mass_features
        else:
            cmf_df = self._combined_mass_features
        # Check if parameters are set to drop isotopologues and drop if so
        if self.parameters.lcms_collection.drop_isotopologues:
            if not self.isotopes_dropped:
                self._drop_isotopologues()
        # Check if scan_time_aligned is in combined_mass_features, try to add if not
        if cmf_df is not None and "scan_time_aligned" not in cmf_df.columns:
            cmb_mf = cmf_df.copy()
            cmb_mf = cmb_mf.reset_index(drop=False)
            lcms_aligned = [True for x in self if "scan_time_aligned" in x.scan_df.columns]
            if len(lcms_aligned) == len(self):
                # Add scan_time_aligned to combined_mass_features dataframe
                scan_time_aligned_list = []
                for lcms_obj in self:
                    scan_time_df_i = lcms_obj.scan_df[["scan", "scan_time_aligned"]]
                    scan_time_df_i["sample_name"] = lcms_obj.sample_name
                    scan_time_aligned_list.append(scan_time_df_i)
                scan_time_aligned_df = pd.concat(scan_time_aligned_list)
                # Rename scan to apex_scan
                scan_time_aligned_df = scan_time_aligned_df.rename(columns={"scan": "apex_scan"})
                cmb_mf_merged = cmb_mf.merge(scan_time_aligned_df, on=["apex_scan", "sample_name"])
                cmb_mf_merged = cmb_mf_merged.set_index("coll_mf_id")
                # Merge scan_time_aligned_df with combined_mass_features on apex_scan and sample_name
                if induced_features:
                    self._combined_induced_mass_features = cmb_mf_merged
                else:
                    self._combined_mass_features = cmb_mf_merged
    
    def plot_tics(self, ms_level=1, type = "raw", plot_legend=False):
        """Plots the TICs for all the LCMS objects in the collection.
        
        Parameters
        -----------
        ms_level : int, optional
            The MS level to plot the TICs for. Defaults to 1.
        type : str, optional
            The type of TIC to plot, either "raw" or "corrected" or "both". Defaults to "raw".
        plot_legend : bool, optional
            If True, plots a legend on the TIC plot that labels each sample. Defaults to False.
        """
        to_plot = []
        if type == "both":
            to_plot = ["raw", "corrected"]
        else:
            to_plot = [type]

        fig, axs = plt.subplots(
            len(to_plot), 1, figsize=(10, 5 * len(to_plot)), sharex=True, squeeze=False
        )
        
        for i, plot_type in enumerate(to_plot):
            ax = axs[i, 0]
            colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(self))))
            for lcms_obj in self:
                c = next(colors)
                # check if lcms_obj is the center of the collection
                self.manifest_dataframe[self.manifest_dataframe['center']].collection_id.values

                
                scan_df = lcms_obj.scan_df
                scan_df = scan_df[scan_df.ms_level == ms_level]
                if plot_type == "corrected":
                    # Check that scan_time_aligned is key in scan_df
                    if "scan_time_aligned" not in scan_df.columns:
                        raise ValueError(f"scan_time_aligned not found in scan_df for {lcms_obj.sample_name}")
                    else:
                        ax.plot(scan_df.scan_time_aligned, scan_df.tic, label=lcms_obj.sample_name, c=c, linewidth=0.3)
                elif plot_type == "raw":
                    ax.plot(scan_df.scan_time, scan_df.tic, label=lcms_obj.sample_name, c=c, linewidth=0.3)
            ax.set_xlabel("Retention Time (min," + f" {plot_type})" )
            ax.set_ylabel("TIC")
            if plot_legend:
                ax.legend()
        plt.show()

    def plot_alignments(self, plot_legend=False):
        """Plots the alignment of the LCMS objects in the collection.
        
        Parameters
        -----------
        plot_legend : bool, optional
            If True, plots a legend on the alignment plot that labels each sample. Defaults to False.        
        """
        fig, ax = plt.subplots(figsize=(10, 5))
        colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(self))))

        for lcms_obj in self:
            c = next(colors)
            scan_df = lcms_obj.scan_df
            if "scan_time_aligned" not in scan_df.columns:
                raise ValueError(f"scan_time_aligned not found in scan_df for {lcms_obj.sample_name}")
            scan_df['time_diff'] = scan_df.scan_time - scan_df.scan_time_aligned
            ax.plot(scan_df.scan_time_aligned, scan_df.time_diff, label=lcms_obj.sample_name, c=c, linewidth=0.3)

        ax.set_xlabel("Aligned Retention Time (min)")
        ax.set_ylabel("Time Difference (min)")
        if plot_legend:
            ax.legend()
        plt.show()

    def _drop_isotopologues(self):
        """Drops isotopologues from the mass features in combined_mass_features dataframe."""
        cmb_mf_df = self._combined_mass_features

        # Keep monos or if no monos
        cmb_monos = cmb_mf_df[cmb_mf_df.monoisotopic_mf_id == cmb_mf_df.mf_id]
        cmb_nomonos = cmb_mf_df[cmb_mf_df.monoisotopic_mf_id.isnull()]
        # Keep deconvoluted parent or if no deconvoluted parent
        cmb_decon_parent = cmb_mf_df[cmb_mf_df.mass_spectrum_deconvoluted_parent | cmb_mf_df.monoisotopic_mf_id.isnull()]

        cmb_mf_df2 = pd.concat([cmb_monos, cmb_nomonos, cmb_decon_parent])
        cmb_mf_df2 = cmb_mf_df2[~cmb_mf_df2.index.duplicated(keep='first')]
        self.isotopes_dropped = True
        self._combined_mass_features = cmb_mf_df2
    

    def load_raw_data(self, sample_idx: int, ms_level = 1) -> None:
        """Load raw data for a specific sample index in the collection.
        
        Parameters
        -----------
        sample_idx : int
            The index of the sample in the collection.
        ms_level : int, optional
            The MS level to load raw data for. Defaults to 1.
            
        Raises
        -------
        IndexError
            If the sample index is out of range.
        ValueError
            If raw data for the specified MS level is already loaded for the sample index.
        ValueError
            If the spectra parser is not set for the LCMS object or if the parser type does not support loading raw data.

        Returns
        --------
        None, but updates the LCMS object with the raw data for the specified MS level.
        """
        if sample_idx < 0 or sample_idx >= len(self.samples):
            raise IndexError("Sample index out of range.")

        # Check that the sample does not already have raw data loaded
        if ms_level in self[sample_idx]._ms_unprocessed:
            raise ValueError(f"Raw data for MS{ms_level} already loaded for sample index {sample_idx}. Drop data first if you want to reload it.")

        # Check the parser type of the LCMS object
        if self[sample_idx].spectra_parser is None:
            raise ValueError("Spectra parser is not set for this LCMS object.")

        # Instantiate the parser and load the raw data using the correct method
        parser = self[sample_idx].spectra_parser
        parser_class_name = self[sample_idx].spectra_parser_class.__name__
        scan_df = self[sample_idx].scan_df

        # Get raw data for the specified MS level using the appropriate method
        if parser_class_name == "ImportMassSpectraThermoMSFileReader":
            self[sample_idx]._ms_unprocessed[ms_level] = parser.get_ms_raw(
                spectra=f"ms{ms_level}",
                scan_df=scan_df
            )[f"ms{ms_level}"]

        elif parser_class_name == "MZMLSpectraParser":
            data = parser.load()
            self[sample_idx]._ms_unprocessed[ms_level] = parser.get_ms_raw(
                spectra=f"ms{ms_level}",
                scan_df=scan_df,
                data=data
                )[f"ms{ms_level}"]

        elif parser_class_name == "ReadCoreMSHDFMassSpectra":
            raise ValueError(
                "ReadCoreMSHDFMassSpectra does not have a method to load raw data. Need to instantiate the original parser to access the raw data."
            )

    def drop_raw_data(self, sample_idx: int, ms_level = 1) -> None:
        """Drop raw data for a specific sample index in the collection.

        Parameters
        -----------
        sample_idx : int
            The index of the sample in the collection.
        ms_level : int, optional
            The MS level to drop raw data for. Defaults to 1.

        Raises
        -------
        IndexError
            If the sample index is out of range.
        ValueError
            If raw data for the specified MS level is not loaded for the sample index.

        Returns
        --------
        None
        """
        if sample_idx < 0 or sample_idx >= len(self.samples):
            raise IndexError("Sample index out of range.")

        # Check that the sample has raw data loaded
        if ms_level not in self[sample_idx]._ms_unprocessed:
            raise ValueError(f"No raw data for MS{ms_level} found for sample index {sample_idx}. Load data first if you want to drop it.")

        # Drop the raw data
        del self[sample_idx]._ms_unprocessed[ms_level]

    def update_raw_file_locations(self, new_raw_folder):
        """Update the raw file locations for all LCMS objects in the collection.
        
        This method updates the path to the original raw data files (.raw, .mzML, etc.)
        that were used to create the processed HDF5 files stored in .corems folders.
        
        Parameters
        -----------
        new_raw_folder : str or Path
            The new folder location containing the raw data files (.raw, .mzML, etc.).
            The method will look for raw files with the same base name as each sample.
            
        Raises
        -------
        FileNotFoundError
            If the new raw folder does not exist.
        FileNotFoundError
            If a raw file for a sample is not found in the new folder.
            
        Returns
        --------
        None, but updates the raw_file_location for each LCMS object in the collection.
        
        Examples
        --------
        If raw files were moved from /old/path/ to /new/path/:
        >>> lcms_collection.update_raw_file_locations("/new/path/")
        """
        from pathlib import Path
        
        if isinstance(new_raw_folder, str):
            new_raw_folder = Path(new_raw_folder)
        
        if not new_raw_folder.exists():
            raise FileNotFoundError(f"Raw data folder does not exist: {new_raw_folder}")
        
        # Common raw file extensions
        raw_extensions = ['.raw', '.mzML', '.mzml']
        
        for sample_name in self.samples:
            lcms_obj = self._lcms[sample_name]
            
            # Try to find the raw file with common extensions
            new_raw_file = None
            for ext in raw_extensions:
                candidate = new_raw_folder / f"{sample_name}{ext}"
                if candidate.exists():
                    new_raw_file = candidate
                    break
            
            if new_raw_file is None:
                raise FileNotFoundError(
                    f"Raw file for sample '{sample_name}' not found in {new_raw_folder}. "
                    f"Tried extensions: {', '.join(raw_extensions)}"
                )
            
            # Update the raw file location and set flag that raw files have been relocated
            lcms_obj.raw_file_location = new_raw_file
        self.raw_files_relocated = True

    def collection_pivot_table(self, attribute = 'coll_mf_id', verbose = True):
        """Generate a pivot table of all regular and induced mass features in
        a collection. Default attribute presented is the mass feature ID, also
        prints a list of other available attributes.

        Parameters
        -----------
        attribute : str
            The desired attribute to be presented in the pivot table. Defaults
            to mass feature ID
        verbose : boolean
            Print out all the possible values the fill the pivot table and list
            attributes that are not collected for induced mass features

        Returns
        --------
        pd.DataFrame
            A DataFrame that displays one given attribute across all clusters
            and samples in a collection
        
        """
        
        mf_pivot = self.mass_features_dataframe.copy()
        mf_pivot.reset_index(inplace = True)
        
        # Only include induced mass features if gap-filling has been performed
        if self.induced_mass_features_dataframe is not None:
            imf_pivot = self.induced_mass_features_dataframe.copy()
            imf_pivot.reset_index(inplace = True)
            # Cluster column extracted from mf_id in _prepare_lcms_mass_features_for_combination
            mf_pivot = pd.concat([mf_pivot, imf_pivot], axis = 0)
            mf_pivot.reset_index(drop = True, inplace = True)
        else:
            imf_pivot = None
            
        mf_pivot['cluster'] = mf_pivot['cluster'].astype(int)

        if verbose:
            print(
                'Attributes available for pivot table:\n',
                [x for x in mf_pivot.columns if x not in ['cluster', 'sample_name', 'mf_id', 'partition_idx', 'idx']]
            )
            if imf_pivot is not None:
                print(
                    '\nAttributes that have no value for induced mass features:\n',
                    imf_pivot.columns[imf_pivot.isna().all()].tolist()            
                )
        
        # Create pivot table and reindex to include all samples (even those with no features)
        pivot = mf_pivot.pivot(index = 'cluster', columns = 'sample_name', values = attribute)
        
        # Reindex columns to include all samples in the collection
        all_samples = self.samples
        pivot = pivot.reindex(columns=all_samples)
        
        return pivot

    def cluster_representatives_table(self):
        """Generate a table of representative mass features from each consensus cluster.
        
        This method returns a DataFrame containing all attributes for the
        representative mass feature from each consensus cluster. The representative
        is selected using the same logic as process_consensus_features().

        Returns
        --------
        pd.DataFrame
            A DataFrame with one row per cluster containing all attributes for 
            each cluster's representative mass feature. Includes:
            - cluster: cluster ID (as a column for easy joining)
            - polarity: ionization polarity from the collection
            - n_samples_detected: number of samples where the cluster was detected
            - All other mass feature attributes from the representative
            
        Notes
        -----
        The representative metric used is determined by
        self.parameters.lcms_collection.consensus_representative_metric and
        is the same metric used by process_consensus_features() for consistency.
        Common options include 'intensity' (highest intensity) or 
        'intensity_prefer_ms2' (highest intensity with preference for MS2 data).
        """
        
        mf_df = self.mass_features_dataframe.copy()
        mf_df.reset_index(inplace = True)
        imf_df = self.induced_mass_features_dataframe.copy()
        imf_df.reset_index(inplace = True)
        # Cluster column extracted from mf_id in _prepare_lcms_mass_features_for_combination
        mf_df = pd.concat([mf_df, imf_df], axis = 0)
        mf_df.reset_index(drop = True, inplace = True)
        mf_df['cluster'] = mf_df['cluster'].astype(int)
        
        # Calculate number of samples per cluster
        cluster_sample_counts = mf_df.groupby('cluster')['sample_id'].nunique().to_dict()
        
        # Use the same representative selection logic as process_consensus_features
        # This uses the configured representative_metric from parameters
        representatives = self.get_representative_mass_features_for_all_clusters()
        
        # Get the coll_mf_ids of the representatives
        representative_ids = representatives['coll_mf_id'].tolist()
        
        # Filter mf_df to only include representative features
        consensus_report = mf_df[mf_df.coll_mf_id.isin(representative_ids)].copy()
        
        # Add polarity (get from first sample in collection)
        if len(self) > 0:
            polarity = self[0].polarity
        else:
            polarity = 'unknown'
        consensus_report['polarity'] = polarity
        
        # Add number of samples detected
        consensus_report['n_samples_detected'] = consensus_report['cluster'].map(cluster_sample_counts)
        
        # Reorder columns to put cluster at the front
        cols = consensus_report.columns.tolist()
        if 'cluster' in cols:
            cols.remove('cluster')
            cols = ['cluster'] + cols
            consensus_report = consensus_report[cols]
        
        # Sort by cluster and return with cluster as a regular column
        return consensus_report.sort_values(by='cluster')

    def feature_annotations_table(self, molecular_metadata=None, drop_unannotated=False):
        """Generate a comprehensive annotation table for all loaded mass features across samples.
        
        This method consolidates MS1 molecular formula assignments and MS2 spectral 
        search results for all mass features across all samples in the collection.
        Only includes representative mass features (one per cluster per sample).
        
        Parameters
        ----------
        molecular_metadata : dict, optional
            Dictionary of MolecularMetadata objects, keyed by metabref_mol_id.
            Required for including molecular metadata in MS2 annotations.
            Default is None.
        drop_unannotated : bool, optional
            If True, drops rows where all annotation columns (everything except 
            cluster, MS2 Spectrum, and representative_sample) are NaN.
            Default is False.
        
        Returns
        -------
        pd.DataFrame
            Consolidated annotation report with columns including:
            - cluster: cluster ID
            - sample_name: sample name
            - sample_id: sample ID
            - Mass Feature ID: mass feature ID within the sample
            - Mass feature attributes (mz, scan_time, intensity, etc.)
            - MS1 annotations (if molecular_formula_search was run)
            - MS2 annotations (if ms2_spectral_search was run)
        
        Notes
        -----
        This method uses the standard LCMSMetabolomicsExport.to_report() workflow
        for each sample, then consolidates all results and adds cluster information.
        
        Only mass features that are loaded in each sample's mass_features dict
        are included (typically the representative features if load_representatives
        was used in process_consensus_features).
        """
        from corems.mass_spectra.output.export import LCMSMetabolomicsExport
        
        # Collect reports from all samples
        all_sample_reports = []
        
        for sample_id, lcms_obj in enumerate(self):
            # Skip samples with no loaded mass features
            if not hasattr(lcms_obj, 'mass_features') or len(lcms_obj.mass_features) == 0:
                continue
            
            sample_name = self.samples[sample_id]
            
            # Create exporter and generate report using standard workflow
            exporter = LCMSMetabolomicsExport("temp", lcms_obj)
            sample_report = exporter.to_report(molecular_metadata=molecular_metadata)
            
            # Add sample information
            sample_report['representative_sample'] = sample_name
            sample_report['sample_id'] = sample_id
            
            # Get cluster information from the mass_features_dataframe
            # Build coll_mf_id for each row to look up cluster
            sample_report['coll_mf_id'] = sample_report['sample_id'].astype(str) + "_" + sample_report['Mass Feature ID'].astype(str)
            
            # Get cluster from mass_features_dataframe
            if self.mass_features_dataframe is not None and 'cluster' in self.mass_features_dataframe.columns:
                mf_df = self.mass_features_dataframe.reset_index()
                cluster_lookup = mf_df.set_index('coll_mf_id')['cluster'].to_dict()
                sample_report['cluster'] = sample_report['coll_mf_id'].map(cluster_lookup)
            else:
                sample_report['cluster'] = None
            
            # Drop temporary coll_mf_id column
            sample_report = sample_report.drop(columns=['coll_mf_id'])
            
            all_sample_reports.append(sample_report)
        
        # Combine all sample reports
        if len(all_sample_reports) == 0:
            raise ValueError("No samples with loaded mass features found in collection")
        
        collection_report = pd.concat(all_sample_reports, ignore_index=True)
        
        # Reorder columns to match specified order
        desired_cols = [
            'cluster',
            'Isotopologue Type',
            'Is Largest Ion after Deconvolution',
            'MS2 Spectrum',
            'Calculated m/z',
            'm/z Error (ppm)',
            'm/z Error Score',
            'Isotopologue Similarity',
            'Confidence Score',
            'Ion Formula',
            'Ion Type',
            'Molecular Formula',
            'inchikey',
            'name',
            'ref_ms_id',
            'Entropy Similarity',
            'Library mzs in Query (fraction)',
            'Spectra with Annotation (n)',
            'representative_sample'
        ]
        
        # Include only desired columns that exist, maintaining order
        cols = [col for col in desired_cols if col in collection_report.columns]
        collection_report = collection_report[cols]
        
        # Optionally drop rows without any annotations
        if drop_unannotated:
            # Columns to exclude from the "all NA" check
            exclude_cols = ['cluster', 'MS2 Spectrum', 'representative_sample']
            # Get annotation columns (everything except the excluded ones)
            annot_cols = [col for col in collection_report.columns if col not in exclude_cols]
            # Keep rows where at least one annotation column is not NA
            if len(annot_cols) > 0:
                collection_report = collection_report[collection_report[annot_cols].notna().any(axis=1)]
        
        # Sort by cluster, then by annotation quality
        sort_cols = ['cluster']
        if 'Entropy Similarity' in collection_report.columns:
            sort_cols.extend(['Entropy Similarity', 'Confidence Score'])
            collection_report = collection_report.sort_values(
                by=sort_cols,
                ascending=[True, False, False]
            )
        elif 'Confidence Score' in collection_report.columns:
            sort_cols.append('Confidence Score')
            collection_report = collection_report.sort_values(
                by=sort_cols,
                ascending=[True, False]
            )
        else:
            collection_report = collection_report.sort_values(by=sort_cols)
        
        return collection_report

    @property
    def parameters(self):
        """
        LCMSCollectionParameters : The parameters used for the LCMS collection.
        """
        return self._parameters
    
    @parameters.setter
    def parameters(self, paramsinstance):
        """
        Sets the parameters used for the LCMS analysis collection.

        Parameters
        -----------
        paramsinstance : LCMSCollectionParameters
            The parameters used for the LC-MS analysis.
        """
        self._parameters = paramsinstance
    
    @property
    def mass_features_dataframe(self):
        self._check_mass_features_df()
        return self._combined_mass_features

    @mass_features_dataframe.setter
    def mass_features_dataframe(self, df):
        # Check that the dataframe has the expected columns
        expected_cols = ["sample_name", "sample_id", "mz", "scan_time"]
        if not all([col in df.columns for col in expected_cols]):
            raise ValueError(f"Expected columns not found in dataframe: {expected_cols}")
        
        # Check that coll_mf_id is the index and it is unique
        if df.index.name != "coll_mf_id":
            raise ValueError("coll_mf_id must be the index of the dataframe")
        if not df.index.is_unique:
            raise ValueError("coll_mf_id must be unique")
        self._combined_mass_features = df

    @property
    def induced_mass_features_dataframe(self):
        self._check_mass_features_df(induced_features = True)
        if self._combined_induced_mass_features is not None and len(self._combined_induced_mass_features) > 0:
            # The cluster column is extracted from mf_id in _prepare_lcms_mass_features_for_combination
            # mf_id format for induced features: c{cluster}_{index}_i
            pass
        return self._combined_induced_mass_features

    @induced_mass_features_dataframe.setter
    def induced_mass_features_dataframe(self, df):
        # Check that the dataframe has the expected columns
        expected_cols = ["sample_name", "sample_id", "mz", "scan_time"]
        if not all([col in df.columns for col in expected_cols]):
            raise ValueError(f"Expected columns not found in dataframe: {expected_cols}")
        
        # Check that coll_mf_id is the index and it is unique
        if df.index.name != "coll_mf_id":
            raise ValueError("coll_mf_id must be the index of the dataframe")
        if not df.index.is_unique:
            raise ValueError("coll_mf_id must be unique")
        self._combined_induced_mass_features = df    
    
    @property
    def cluster_summary_dataframe(self):
        return self.summarize_clusters()
    
    @property
    def samples(self):
        manifest_df = self.manifest_dataframe
        # order by batch, then by order
        manifest_df = manifest_df.sort_values(by=['batch', 'order'])
        return manifest_df.index.tolist()
    
    @property
    def manifest(self):
        return self._manifest_dict
    
    @property
    def manifest_dataframe(self):
        return pd.DataFrame(self._manifest_dict).T

    @property
    def raw_files(self):
        """Returns a list of raw files in the collection."""
        return [x.raw_file_location for x in self]
    
    @property
    def rt_alignments(self):
        """Returns a dictionary of retention time alignments for the collection."""
        if self.rt_aligned:
            _rt_alignments = {}
            # Construct a dictionary of aligned retention times (stored on each LCMS object within the collection, not the collection itself)
            for i, lcms_obj in enumerate(self):
                aligned_times = [x for k, x in sorted(lcms_obj._scan_info["scan_time_aligned"].items())]
                _rt_alignments[i] = aligned_times
            return _rt_alignments
        else:
            return None
    
    @property
    def cluster_feature_dictionary(self):
        """Generates a dictionary with clusters for keys and mass feature IDs as entries"""
        df = self.mass_features_dataframe
        cluster_dict = df.groupby('cluster').apply(lambda x: x.index.tolist()).to_dict()
        return cluster_dict
    
    def get_eics_for_cluster(self, cluster_id):
        """
        Retrieve all EICs for mass features in a specific cluster across all samples.
        
        Returns a dictionary mapping sample names to EIC_Data objects for the given cluster.
        Useful for visualizing and comparing chromatographic peaks across samples.
        
        Parameters
        ----------
        cluster_id : int
            The cluster ID to retrieve EICs for
            
        Returns
        -------
        dict
            Dictionary with structure: {sample_name: EIC_Data object}
            Only includes samples where the EIC was loaded.
            
        Examples
        --------
        >>> # Load EICs first
        >>> collection.process_consensus_features(gather_eics=True, ...)
        >>> 
        >>> # Get all EICs for cluster 5
        >>> eics = collection.get_eics_for_cluster(5)
        >>> for sample_name, eic_data in eics.items():
        ...     print(f"{sample_name}: {len(eic_data.scans)} scans")
        
        Notes
        -----
        Requires that EICs have been loaded using gather_eics=True in
        process_consensus_features() or manually loaded via LoadEICsOperation.
        """
        eics_by_sample = {}
        
        # Iterate through all samples
        for sample_id, sample in enumerate(self):
            sample_name = self.samples[sample_id]
            
            # Check if sample has EICs loaded
            if not hasattr(sample, 'eics') or not sample.eics:
                continue
            
            # Find mass features in this cluster for this sample
            # Check both regular and induced mass features
            for mf in list(sample.mass_features.values()) + list(sample.induced_mass_features.values()):
                if hasattr(mf, 'cluster_index') and mf.cluster_index == cluster_id:
                    # Get the EIC for this mass feature's m/z
                    if mf.mz in sample.eics:
                        eics_by_sample[sample_name] = sample.eics[mf.mz]
                        break  # Found the EIC for this sample, move to next sample
        
        return eics_by_sample