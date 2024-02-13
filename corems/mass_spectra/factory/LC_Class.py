"""
Created on Oct 15, 2021
"""
__author__ = "Yuri E. Corilo"
__date__ = "Oct 15, 2021"


from collections.abc import Mapping
from dataclasses import field, dataclass
import logging
from pathlib import Path
from typing import List, Dict, List, Tuple
import sys
import site
import numpy as np
import pandas as pd

from matplotlib import axes
import numpy as np
from corems.encapsulation.factory.parameters import LCMSParameters

from corems.mass_spectra.calc.LC_Calc import LC_Calculations, PH_Calculations
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile

class MassSpectraBase():
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
        A dictionary containing the scan data with columns for scan number, scan time, ms level, precursor m/z, scan text, and scan window (lower and upper). Associated with the property scan_df, which returns a pandas DataFrame or can set this attribute from a pandas DataFrame.
    _ms : dict
        A dictionary containing mass spectra for the dataset, keys of dictionary are scan numbers. Initialized as an empty dictionary.
    _ms_unprocessed: dictionary of pandas.DataFrames or None
        A dictionary of unprocssed mass spectra data, as an (optional) intermediate data product for peak picking.  Key is ms_level, and value is dataframe with columns for scan number, m/z, and intensity. Default is None.
    
    Methods
    --------
    * add_mass_spectra(scan_list, spectrum_mode: str = 'profile', use_parser = True, auto_process=True).
        Add mass spectra (or singlel mass spectrum) to _ms slot, from a list of scans
    """
    def __init__(
            self, 
            file_location, 
            analyzer='Unknown', 
            instrument_label='Unknown', 
            sample_name=None,
            spectra_parser = None
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
            if self.sample_name is not None and self.sample_name != self.spectra_parser.sample_name:
                logging.warning("sample_name provided to MassSpectraBase object does not match sample_name provided to spectra parser object")
            if self.analyzer != self.spectra_parser.analyzer:
                logging.warning("analyzer provided to MassSpectraBase object does not match analyzer provided to spectra parser object")
            if self.instrument_label != self.spectra_parser.instrument_label:
                logging.warning("instrument provided to MassSpectraBase object does not match instrument provided to spectra parser object")
            if self.file_location != self.spectra_parser.file_location:
                logging.warning("file_location provided to MassSpectraBase object does not match file_location provided to spectra parser object")  

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
        if not hasattr(mass_spec, 'scan_number'):
            raise ValueError("Mass spectrum must have a scan_number attribute to be added to the dataset correctly")
        self._ms[mass_spec.scan_number] = mass_spec
    
    
    def add_mass_spectra(self, scan_list, spectrum_mode: str = 'profile', ms_level = 1, use_parser = True, auto_process=True):
        """Add mass spectra to _ms dictionary, from a list of scans or single scan

        Notes
        -----
        The mass spectra will inherit the mass_spectrum, ms_peak, and molecular_search parameters from the LCMSBase object.

        
        Parameters
        -----------
        scan_list : list of ints
            List of scans to use to populate _ms slot
        spectrum_mode : str, optional
            The spectrum mode to use for the mass spectra.  Only profile mode is currently supported. Defaults to 'profile'.
        ms_level : int, optional
            The MS level to use for the mass spectra.  This is used to pass the molecular_search parameters from the LCMS object to the individual MassSpectrum objects. Defaults to 1.
        using_parser : bool
            Whether to use the mass spectra parser to get the mass spectra.  Defaults to True.
        auto_process : bool
            Whether to auto-process the mass spectra.  Defaults to True.

        Raises
        ------
        TypeError
            If scan_list is not a list of ints
        ValueError
            If polarity is not 'positive' or 'negative' or if it is a mix
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
        if self.polarity == set(['negative']):
            polarity = -1
        elif self.polarity == set(['positive']):
            polarity = 1
        else:
            raise ValueError("Polarity not set for dataset, must be a set containing either 'positive' or 'negative'")
        
        # is not using_parser, check that ms1 and ms2 are not None
        if use_parser == False:
            if ms_level not in self._ms_unprocessed.keys():
                raise ValueError("ms_level {} not found in _ms_unprocessed dictionary".format(ms_level))
        
        scan_list = list(set(scan_list))
        scan_list.sort()
        for scan in scan_list: 
            ms = None
            # Instantiate the mass spectrum object using the parser or the unprocessed data
            if use_parser == False:
                if self._ms_unprocessed[ms_level] is not None and scan in self._ms_unprocessed[ms_level].scan.tolist():
                    my_ms_df = self._ms_unprocessed[ms_level][self._ms_unprocessed[ms_level].scan == scan]
                    ms = ms_from_array_profile(my_ms_df.mz, my_ms_df.intensity, self.file_location, polarity=polarity, auto_process=False)
            if use_parser == True:
                ms = self.spectra_parser.get_mass_spectrum_from_scan(scan_number = scan, spectrum_mode = spectrum_mode, polarity = polarity, auto_process=auto_process)
            
            # Set the mass spectrum parameters, auto-process if auto_process is True, and add to the dataset
            if ms is not None:
                ms.parameters.mass_spectrum = self.parameters.mass_spectrum
                ms.parameters.ms_peak = self.parameters.ms_peak
                if ms_level == 1:
                    ms.parameters.molecular_search = self.parameters.ms1_molecular_search
                elif ms_level == 2:
                    ms.parameters.molecular_search = self.parameters.ms2_molecular_search
                else:
                    raise ValueError("ms_level must be 1 or 2")
                ms.scan_number = scan
                if len(ms.mz_exp_profile) > 3:
                    if auto_process:
                        ms.process_mass_spec()
                self.add_mass_spectrum(ms)
    
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
            A pandas DataFrame containing the scan data with columns for scan number, scan time, ms level, precursor m/z, scan text, and scan window (lower and upper).
        """
        self._scan_info = df.to_dict()    

    def __getitem__(self, scan_number):
        
        return self._ms.get(scan_number)


class LCMSBase(MassSpectraBase, LC_Calculations, PH_Calculations):
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
        A dictionary containing extracted ion chromatograms (EICs) for the dataset. Key is the mz of the EIC. Initialized as an empty dictionary.
    mass_features : dictionary of LCMSMassFeature objects
        A dictionary containing mass features for the dataset. Key is mass feature ID. Initialized as an empty dictionary.

    Methods
    --------
    * get_parameters_json(verbose=True).
        Returns the parameters used for the LC-MS analysis in JSON format.
    * add_associated_ms2_dda(add_to_lcmsobj=True, auto_process=True, use_parser=True)
        Adds which MS2 scans are associated with each mass feature to the mass_features dictionary and optionally adds the MS2 spectra to the _ms dictionary.
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
            analyzer='Unknown', 
            instrument_label='Unknown', 
            sample_name=None,
            spectra_parser = None
            ):
        super().__init__(
            file_location,
            analyzer,
            instrument_label,
            sample_name,
            spectra_parser
            )
        self.polarity = ""
        self._parameters = LCMSParameters()
        self._retention_time_list = []
        self._scans_number_list = []
        self._tic_list = []
        self.eics = {}
        self.mass_features = {}

    def get_parameters_json(self, verbose = False):
        """Returns the parameters stored for the LC-MS object in JSON format.

        Parameters
        -----------
        verbose : bool, optional
            If True, prints the JSON string to the console. Default is False.

        Returns
        --------
        str
            The parameters used for the LC-MS analysis in JSON format.
        """
        if verbose:
            print(self.parameters.to_json())
        return self.parameters.to_json()
    
    def add_associated_ms2_dda(self, add_to_lcmsobj = True, auto_process=True, use_parser=True):
        """Add MS2 spectra associated with mass features to the dataset.

        Populates the mass_features ms2_scan_numbers attribute (on mass_features dictionary on LCMSObject)

        Parameters
        -----------
        add_to_lcmsobj : bool, optional
            If True, adds the MS2 spectra to the object's _ms dictionary with the scan as the key. Default is True.
        auto_process : bool, optional
            If True, auto-processes the MS2 spectra before adding it to the object's _ms dictionary. Default is True.
        use_parser : bool, optional
            If True, envoke the spectra parser to get the MS2 spectra. Default is True.
        
        Raises
        ------
        ValueError
            If mass_features is not set, must run find_mass_features() first.
            If no MS2 scans are found in the dataset.
            If no precursor m/z values are found in MS2 scans, not a DDA dataset.
        """
        if self.mass_features is None:
            raise ValueError("mass_features not set, must run find_mass_features() first")
        
        mf_df = self.mass_features_to_df().copy()
        # Find ms2 scans that have a precursor m/z value
        ms2_scans = self.scan_df[self.scan_df.ms_level == 2]
        ms2_scans = ms2_scans[~ms2_scans.precursor_mz.isna()]
        # drop ms2 scans that have no tic
        ms2_scans = ms2_scans[ms2_scans.tic > 0]
        if ms2_scans is None:
            raise ValueError("No DDA scans found in dataset")  

        # set tolerance in rt space (in minutes) and mz space (in daltons)
        time_tol = self.parameters.lc_ms.ms2_dda_rt_tolerance
        mz_tol = self.parameters.lc_ms.ms2_dda_mz_tolerance
       
        # for each mass feature, find the ms2 scans that are within the roi scan time and mz range
        dda_scans = []
        for i, row in mf_df.iterrows():
            ms2_scans_filtered = ms2_scans[ms2_scans.scan_time.between(row.scan_time - time_tol, row.scan_time + time_tol)]
            ms2_scans_filtered = ms2_scans_filtered[ms2_scans_filtered.precursor_mz.between(row.mz - mz_tol, row.mz + mz_tol)]
            dda_scans = dda_scans + ms2_scans_filtered.scan.tolist()
            self.mass_features[i].ms2_scan_numbers = ms2_scans_filtered.scan.tolist()
        if add_to_lcmsobj:
            # add to _ms attribute
            self.add_mass_spectra(scan_list=list(set(dda_scans)), auto_process=auto_process, use_parser=use_parser)
            
            # associate appropriate _ms attribute to appropriate mass feature's ms2_mass_spectra attribute
            for mf_id in self.mass_features:
                if self.mass_features[mf_id].ms2_scan_numbers is not None:
                    for dda_scan in self.mass_features[mf_id].ms2_scan_numbers:
                        if dda_scan in self._ms.keys():
                            self.mass_features[mf_id].ms2_mass_spectra[dda_scan] = self._ms[dda_scan]
    
    def mass_features_to_df(self):
        """Returns a pandas dataframe summarizing the mass features.

        The dataframe contains the following columns: mf_id, mz, apex_scan, scan_time, intensity, persistence, and area.  The index is set to mf_id (mass feature ID).

        Returns
        --------
        pandas.DataFrame
            A pandas dataframe of mass features with the following columns: mf_id, mz, apex_scan, scan_time, intensity, persistence, area.
        """
        cols_in_df = ['mf_id', 'mz', 'apex_scan', 'scan_time', 'intensity', 'persistence', '_area']
        df_mf_list = []
        for mf_id in self.mass_features.keys():
            # Find cols_in_df that are in single_mf
            df_keys = list(set(cols_in_df).intersection(self.mass_features[mf_id].__dir__()))
            dict_mf = {'mf_id': mf_id}
            for key in df_keys:
                dict_mf[key] = getattr(self.mass_features[mf_id], key)
            df_mf_list.append(pd.DataFrame(dict_mf, index=[mf_id]))
        df_mf = pd.concat(df_mf_list).set_index('mf_id', drop=True)
        # rename _area to area
        df_mf = df_mf.rename(columns={'_area': 'area'})
         
        return df_mf
  
    
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

    def set_tic_list_from_data(self, overwrite = False):
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
        if len(self.tic) > 0 and overwrite == False:
            raise ValueError("TIC list already set, use overwrite=True to overwrite")

        self.tic = [self._ms.get(i).tic for i in self.scans_number]
        
    def set_retention_time_from_data(self, overwrite = False):
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
        if len(self.retention_time) > 0 and overwrite == False:
            raise ValueError("Retention time list already set, use overwrite=True to overwrite")

        retention_time_list = []
        for key_ms in sorted(self._ms.keys()):
            retention_time_list.append(self._ms.get(key_ms).retention_time)
        self.retention_time = retention_time_list 

    def set_scans_number_from_data(self, overwrite = False):
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
        if len(self.scans_number) > 0 and overwrite == False:
            raise ValueError("Scan number list already set, use overwrite=True to overwrite")
        
        self.scans_number = sorted(self._ms.keys())

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
    def scans_number(self, l):
        """
        Sets the scan numbers for the dataset.

        Parameters
        -----------
        l : list
            A list of scan numbers for the dataset.
        """
        self._scans_number_list = l

    @property
    def retention_time(self):
        """
        numpy.ndarray : An array of retention times for the dataset.
        """
        return self._retention_time_list

    @retention_time.setter
    def retention_time(self, l):
        """
        Sets the retention times for the dataset.

        Parameters
        -----------
        l : list
            A list of retention times for the dataset.
        """
        self._retention_time_list = np.array(l)

    @property
    def tic(self):
        """
        numpy.ndarray : An array of TIC values for the dataset.
        """
        return self._tic_list
    
    @tic.setter
    def tic(self, l):
        """
        Sets the TIC values for the dataset.

        Parameters
        -----------
        l : list
            A list of TIC values for the dataset.
        """
        self._tic_list = np.array(l)    