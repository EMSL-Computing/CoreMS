from collections import defaultdict
from pathlib import Path
from typing import Optional, Union, List, Tuple

import numpy as np
import pandas as pd
import pymzml
import datetime

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectra.factory.lc_class import LCMSBase, MassSpectraBase
from corems.mass_spectra.input.parserbase import SpectraParserInterface
from corems.mass_spectrum.factory.MassSpectrumClasses import (
    MassSpecCentroid,
    MassSpecProfile,
)


class MZMLSpectraParser(SpectraParserInterface):
    """A class for parsing mzml spectrometry data files into MassSpectraBase or LCMSBase objects

    Parameters
    ----------
    file_location : str or Path
        The path to the RAW file to be parsed.
    analyzer : str, optional
        The type of mass analyzer used in the instrument. Default is "Unknown".
    instrument_label : str, optional
        The name of the instrument used to acquire the data. Default is "Unknown".
    sample_name : str, optional
        The name of the sample being analyzed. If not provided, the stem of the file_location path will be used.

    Attributes
    ----------
    file_location : Path
        The path to the RAW file being parsed.
    analyzer : str
        The type of mass analyzer used in the instrument.
    instrument_label : str
        The name of the instrument used to acquire the data.
    sample_name : str
        The name of the sample being analyzed.

    Methods
    -------
    * load().
        Load mzML file using pymzml.run.Reader and return the data as a numpy array.
    * run(spectra=True).
        Parses the mzml file and returns a dictionary of mass spectra dataframes and a scan metadata dataframe.
    * get_mass_spectrum_from_scan(scan_number, polarity, auto_process=True)
        Parses the mzml file and returns a MassSpecBase object from a single scan.
    * get_mass_spectra_obj().
        Parses the mzml file and instantiates a MassSpectraBase object.
    * get_lcms_obj().
        Parses the mzml file and instantiates an LCMSBase object.
    * get_instrument_info().
        Return instrument information from the mzML file.
    * get_creation_time().
        Return the creation time of the mzML file as a datetime object.

    Inherits from ThermoBaseClass and SpectraParserInterface
    """

    def __init__(
        self,
        file_location,
        analyzer="Unknown",
        instrument_label="Unknown",
        sample_name=None,
    ):
        # implementation details
        if isinstance(file_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)
        if not file_location.exists():
            raise FileExistsError("File does not exist: " + str(file_location))
        self.file_location = file_location
        self.analyzer = analyzer
        self.instrument_label = instrument_label

        if sample_name:
            self.sample_name = sample_name
        else:
            self.sample_name = file_location.stem

    def load(self):
        """
        Load mzML file using pymzml.run.Reader and return the data as a numpy array.

        Returns
        -------
        numpy.ndarray
            The mass spectra data as a numpy array.
        """
        data = pymzml.run.Reader(self.file_location)
        return data

    def get_scan_df(self, data=None, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """
        Return scan data as a pandas DataFrame.

        Parameters
        ----------
        data : pymzml.run.Reader, optional
            The mass spectra data. If None, will load the data.
        time_range : tuple or list of tuples, optional
            Retention time range(s) to filter scans. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, returns all scans.

        Returns
        -------
        pandas.DataFrame
            A pandas DataFrame containing metadata for each scan, including scan number, MS level, polarity, and scan time.
        """
        if data is None:
            data = self.load()
        # Scan dict
        # instatinate scan dict, with empty lists of size of scans
        n_scans = data.get_spectrum_count()
        scan_dict = {
            "scan": np.empty(n_scans, dtype=np.int32),
            "scan_time": np.empty(n_scans, dtype=np.float32),
            "ms_level": [None] * n_scans,
            "polarity": [None] * n_scans,
            "precursor_mz": [None] * n_scans,
            "scan_text": [None] * n_scans,
            "scan_window_lower": np.empty(n_scans, dtype=np.float32),
            "scan_window_upper": np.empty(n_scans, dtype=np.float32),
            "scan_precision": [None] * n_scans,
            "tic": np.empty(n_scans, dtype=np.float32),
            "ms_format": [None] * n_scans,
        }

        # First pass: loop through scans to get scan info
        for i, spec in enumerate(data):
            scan_dict["scan"][i] = spec.ID
            scan_dict["ms_level"][i] = spec.ms_level
            scan_dict["scan_precision"][i] = spec._measured_precision
            scan_dict["tic"][i] = spec.TIC
            if spec.selected_precursors:
                scan_dict["precursor_mz"][i] = spec.selected_precursors[0].get(
                    "mz", None
                )
            if spec["negative scan"] is not None:
                scan_dict["polarity"][i] = "negative"
            if spec["positive scan"] is not None:
                scan_dict["polarity"][i] = "positive"
            if spec["negative scan"] is not None and spec["positive scan"] is not None:
                raise ValueError(
                    "Error: scan {0} has both negative and positive polarity".format(
                        spec.ID
                    )
                )

            scan_dict["scan_time"][i] = spec.get("MS:1000016")
            scan_dict["scan_text"][i] = spec.get("MS:1000512")
            scan_dict["scan_window_lower"][i] = spec.get("MS:1000501")
            scan_dict["scan_window_upper"][i] = spec.get("MS:1000500")
            if spec.get("MS:1000128"):
                scan_dict["ms_format"][i] = "profile"
            elif spec.get("MS:1000127"):
                scan_dict["ms_format"][i] = "centroid"
            else:
                scan_dict["ms_format"][i] = None

        scan_df = pd.DataFrame(scan_dict)
        
        # Apply time range filtering if specified
        if time_range is not None:
            time_ranges = self._normalize_time_range(time_range)
            # Create a mask for scans within any of the time ranges
            mask = np.zeros(len(scan_df), dtype=bool)
            for start_time, end_time in time_ranges:
                mask |= (scan_df["scan_time"] >= start_time) & (scan_df["scan_time"] <= end_time)
            scan_df = scan_df[mask].reset_index(drop=True)

        return scan_df

    def get_ms_raw(self, spectra, scan_df, data=None, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """Return a dictionary of mass spectra data as a pandas DataFrame.

        Parameters
        ----------
        spectra : str
            Which mass spectra data to include in the output.
            Options: None, "ms1", "ms2", "all".
        scan_df : pandas.DataFrame
            Scan dataframe. Output from get_scan_df().
        data : pymzml.run.Reader, optional
            The mass spectra data. If None, will load the data.
        time_range : tuple or list of tuples, optional
            Retention time range(s) to filter scans. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, returns all scans. Note: filtering is typically done at scan_df level.

        Returns
        -------
        dict
            A dictionary containing the mass spectra data as pandas DataFrames, with keys corresponding to the MS level.

        """
        if data is None:
            data = self.load()
        if spectra == "all":
            scan_df_forspec = scan_df
        elif spectra == "ms1":
            scan_df_forspec = scan_df[scan_df.ms_level == 1]
        elif spectra == "ms2":
            scan_df_forspec = scan_df[scan_df.ms_level == 2]
        else:
            raise ValueError("spectra must be 'all', 'ms1', or 'ms2'")

        # Result container
        res = {}

        # Row count container
        counter = {}

        # Column name container
        cols = {}

        # set at float32
        dtype = np.float32

        # First pass: get nrows
        N = defaultdict(lambda: 0)
        for i, spec in enumerate(data):
            if spec.ID in scan_df_forspec.scan.values:
                # Get ms level
                level = "ms{}".format(spec.ms_level)

                # Number of rows
                N[level] += spec.mz.shape[0]

        # Second pass: parse
        for i, spec in enumerate(data):
            if spec.ID in scan_df_forspec.scan.values:
                # Number of rows
                n = spec.mz.shape[0]

                # No measurements
                if n == 0:
                    continue

                # Dimension check
                if len(spec.mz) != len(spec.i):
                    # raise an error if the mz and intensity arrays are not the same length
                    raise ValueError("m/z and intensity array dimension mismatch")

                # Scan/frame info
                id_dict = spec.id_dict

                # Get ms level
                level = "ms{}".format(spec.ms_level)

                # Columns
                cols[level] = list(id_dict.keys()) + ["mz", "intensity"]
                m = len(cols[level])

                # Subarray init
                arr = np.empty((n, m), dtype=dtype)
                inx = 0

                # Populate scan/frame info
                for k, v in id_dict.items():
                    arr[:, inx] = v
                    inx += 1

                # Populate m/z
                arr[:, inx] = spec.mz
                inx += 1

                # Populate intensity
                arr[:, inx] = spec.i
                inx += 1

                # Initialize output container
                if level not in res:
                    res[level] = np.empty((N[level], m), dtype=dtype)
                    counter[level] = 0

                # Insert subarray
                res[level][counter[level] : counter[level] + n, :] = arr
                counter[level] += n

        # Construct ms1 and ms2 mz dataframes
        for level in res.keys():
            res[level] = pd.DataFrame(res[level], columns=cols[level]).drop(
                columns=["controllerType", "controllerNumber"],
                axis=1,
                inplace=False,
            )

        return res

    def run(self, spectra="all", scan_df=None, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """Parse the mzML file and return a dictionary of spectra dataframes and a scan metadata dataframe.

        Parameters
        ----------
        spectra : str, optional
            Which mass spectra data to include in the output. Default is "all".
            Other options: None, "ms1", "ms2".
        scan_df : pandas.DataFrame, optional
            Scan dataframe.  If not provided, the scan dataframe is created from the mzML file.
        time_range : tuple or list of tuples, optional
            Retention time range(s) to load. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, loads all scans.

        Returns
        -------
        tuple
            A tuple containing two elements:
            - A dictionary containing the mass spectra data as numpy arrays, with keys corresponding to the MS level.
            - A pandas DataFrame containing metadata for each scan, including scan number, MS level, polarity, and scan time.
        """

        # Open file
        data = self.load()

        if scan_df is None:
            scan_df = self.get_scan_df(data, time_range=time_range)

        if spectra != "none":
            res = self.get_ms_raw(spectra, scan_df, data)

        else:
            res = None

        return res, scan_df

    def get_mass_spectrum_from_scan(
        self, scan_number, spectrum_mode, auto_process=True
    ):
        """Instatiate a mass spectrum object from the mzML file.

        Parameters
        ----------
        scan_number : int
            The scan number to be parsed.
        spectrum_mode : str
            The type of spectrum to instantiate.  Must be'profile' or 'centroid'.
        polarity : int
            The polarity of the scan.  Must be -1 or 1.
        auto_process : bool, optional
            If True, process the mass spectrum. Default is True.

        Returns
        -------
        MassSpecProfile | MassSpecCentroid
            The MassSpecProfile or MassSpecCentroid object containing the parsed mass spectrum.
        """
        # Use the batch function and return the first result
        result_list = self.get_mass_spectra_from_scan_list(
            [scan_number], spectrum_mode, auto_process
        )
        return result_list[0] if result_list else None

    def get_mass_spectra_from_scan_list(
        self, scan_list, spectrum_mode, auto_process=True
    ):
        """Instatiate mass spectrum objects from the mzML file.

        Parameters
        ----------
        scan_list : list of int
            The scan numbers to be parsed.
        spectrum_mode : str
            The type of spectrum to instantiate.  Must be'profile' or 'centroid'.
        auto_process : bool, optional
            If True, process the mass spectrum. Default is True.

        Returns
        -------
        list of MassSpecProfile | MassSpecCentroid
            List of MassSpecProfile or MassSpecCentroid objects containing the parsed mass spectra.
        """

        def set_metadata(
            scan_number: int,
            polarity: int,
            file_location: str,
            label=Labels.thermo_profile,
        ):
            """
            Set the output parameters for creating a MassSpecProfile or MassSpecCentroid object.

            Parameters
            ----------
            scan_number : int
                The scan number.
            polarity : int
                The polarity of the data.
            file_location : str
                The file location.
            label : str, optional
                The label for the mass spectrum. Default is Labels.thermo_profile.

            Returns
            -------
            dict
                The output parameters ready for creating a MassSpecProfile or MassSpecCentroid object.
            """
            d_params = default_parameters(file_location)
            d_params["label"] = label
            d_params["polarity"] = polarity
            d_params["filename_path"] = file_location
            d_params["scan_number"] = scan_number

            return d_params

        # Open file
        data = self.load()

        mass_spectrum_objects = []
        
        for scan_number in scan_list:
            # Pluck out individual scan mz and intensity
            spec = data[scan_number]

            # Get polarity
            if spec["negative scan"] is not None:
                polarity = -1
            elif spec["positive scan"] is not None:
                polarity = 1

            # Get mass spectrum
            if spectrum_mode == "profile":
                # Check if profile
                if not spec.get("MS:1000128"):
                    raise ValueError("spectrum is not profile")
                data_dict = {
                    Labels.mz: spec.mz,
                    Labels.abundance: spec.i,
                }
                d_params = set_metadata(
                    scan_number,
                    polarity,
                    self.file_location,
                    label=Labels.simulated_profile,
                )
                mass_spectrum_obj = MassSpecProfile(
                    data_dict, d_params, auto_process=auto_process
                )
            elif spectrum_mode == "centroid":
                # Check if centroided
                if not spec.get("MS:1000127"):
                    raise ValueError("spectrum is not centroided")
                data_dict = {
                    Labels.mz: spec.mz,
                    Labels.abundance: spec.i,
                    Labels.rp: [np.nan] * len(spec.mz),
                    Labels.s2n: [np.nan] * len(spec.i),
                }
                d_params = set_metadata(
                    scan_number, polarity, self.file_location, label=Labels.corems_centroid
                )
                mass_spectrum_obj = MassSpecCentroid(
                    data_dict, d_params, auto_process=auto_process
                )
            
            mass_spectrum_objects.append(mass_spectrum_obj)

        return mass_spectrum_objects

    def get_mass_spectra_obj(self, time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """Instatiate a MassSpectraBase object from the mzML file.

        Parameters
        ----------
        time_range : tuple or list of tuples, optional
            Retention time range(s) to load. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, loads all scans.

        Returns
        -------
        MassSpectraBase
            The MassSpectra object containing the parsed mass spectra.
            The object is instatiated with the mzML file, analyzer, instrument, sample name, and scan dataframe.
        """
        _, scan_df = self.run(spectra=False, time_range=time_range)
        mass_spectra_obj = MassSpectraBase(
            self.file_location,
            self.analyzer,
            self.instrument_label,
            self.sample_name,
            self,
        )
        scan_df = scan_df.set_index("scan", drop=False)
        mass_spectra_obj.scan_df = scan_df

        return mass_spectra_obj

    def get_lcms_obj(self, spectra="all", time_range: Optional[Union[Tuple[float, float], List[Tuple[float, float]]]] = None):
        """Instatiates a LCMSBase object from the mzML file.

        Parameters
        ----------
        spectra : str, optional
            Which mass spectra data to include in the output. Default is all.  Other options: none, ms1, ms2.
        time_range : tuple or list of tuples, optional
            Retention time range(s) to load. Can be:
            - Single range: (start_time, end_time) in minutes
            - Multiple ranges: [(start1, end1), (start2, end2), ...] in minutes
            If None, loads all scans. Useful for targeted workflows to improve performance.

        Returns
        -------
        LCMSBase
            LCMS object containing mass spectra data.
            The object is instatiated with the mzML file, analyzer, instrument, sample name, scan dataframe,
            and mz dataframe(s), as well as lists of scan numbers, retention times, and TICs.
        """
        _, scan_df = self.run(spectra="none", time_range=time_range)  # first run it to just get scan info
        if spectra != "none":
            res, scan_df = self.run(
                scan_df=scan_df, spectra=spectra, time_range=time_range
            )  # second run to parse data
        lcms_obj = LCMSBase(
            self.file_location,
            self.analyzer,
            self.instrument_label,
            self.sample_name,
            self,
        )
        if spectra != "none":
            for key in res:
                key_int = int(key.replace("ms", ""))
                res[key] = res[key][res[key].intensity > 0]
                res[key] = res[key].sort_values(by=["scan", "mz"]).reset_index(drop=True)
                lcms_obj._ms_unprocessed[key_int] = res[key]
        lcms_obj.scan_df = scan_df.set_index("scan", drop=False)
        # Check if polarity is mixed
        if len(set(scan_df.polarity)) > 1:
            raise ValueError("Mixed polarities detected in scan data")
        lcms_obj.polarity = scan_df.polarity[0]
        lcms_obj._scans_number_list = list(scan_df.scan)
        lcms_obj._retention_time_list = list(scan_df.scan_time)
        lcms_obj._tic_list = list(scan_df.tic)

        return lcms_obj

    def get_scans_in_time_range(
        self, 
        time_range: Union[Tuple[float, float], List[Tuple[float, float]]],
        ms_level: Optional[int] = None
    ) -> List[int]:
        """
        Return scan numbers within specified retention time range(s).
        
        This method provides efficient filtering of scans by retention time,
        which is particularly useful for targeted workflows where only specific
        time windows are of interest.
        
        Parameters
        ----------
        time_range : tuple or list of tuples
            Retention time range(s) in minutes. Can be:
            - Single range: (start_time, end_time)
            - Multiple ranges: [(start1, end1), (start2, end2), ...]
        ms_level : int, optional
            If specified, only return scans of this MS level (e.g., 1 for MS1, 2 for MS2).
            If None, returns scans of all MS levels.
        
        Returns
        -------
        list of int
            List of scan numbers within the specified time range(s) and MS level.
        
        Examples
        --------
        Get MS1 scans between 1.0 and 2.0 minutes:
        
        >>> scans = parser.get_scans_in_time_range((1.0, 2.0), ms_level=1)
        
        Get scans in multiple time windows:
        
        >>> scans = parser.get_scans_in_time_range([(0.5, 1.5), (3.0, 4.0)])
        """
        # Get scan dataframe filtered by time range
        scan_df = self.get_scan_df(time_range=time_range)
        
        # Further filter by MS level if specified
        if ms_level is not None:
            scan_df = scan_df[scan_df.ms_level == ms_level]
        
        # Return list of scan numbers
        return scan_df.scan.tolist()

    def get_instrument_info(self):
        """
        Return instrument information.

        Returns
        -------
        dict
            A dictionary with the keys 'model' and 'serial_number'.
        """
        # Load the pymzml data
        data = self.load()
        instrument_info = data.info.get('referenceable_param_group_list_element')[0]
        cv_params = instrument_info.findall('{http://psi.hupo.org/ms/mzml}cvParam')

        # Extract details from each cvParam
        params = []
        for param in cv_params:
            accession = param.get('accession')  # Get 'accession' attribute
            name = param.get('name')           # Get 'name' attribute
            value = param.get('value')         # Get 'value' attribute
            params.append({
                'accession': accession,
                'name': name,
                'value': value
            })

        # Loop through params and try to find the relevant information
        instrument_dict = {
            'model': 'Unknown',
            'serial_number': 'Unknown'
        }

        # Assuming there are only two paramters here - one is for the serial number (agnostic to the model) and the other is for the model
        # If there are more than two, we raise an error
        if len(params) < 2:
            raise ValueError("Not enough parameters found in the instrument info, cannot parse.")
        if len(params) > 2:
            raise ValueError("Too many parameters found in the instrument info, cannot parse.")
        for param in params:
            if param['accession'] == 'MS:1000529':
                instrument_dict['serial_number'] = param['value']
            else:
                instrument_dict['model'] = data.OT[param['accession']]     

        return instrument_dict

    def get_creation_time(self) -> datetime.datetime:
        """
        Return the creation time of the mzML file.
        """
        data = self.load()
        write_time = data.info.get('start_time')
        if write_time:
            # Convert the write time to a datetime object
            return datetime.datetime.strptime(write_time, "%Y-%m-%dT%H:%M:%SZ")
        else:
            raise ValueError("Creation time is not available in the mzML file. "
                           "Please ensure the file contains the 'start_time' information.")