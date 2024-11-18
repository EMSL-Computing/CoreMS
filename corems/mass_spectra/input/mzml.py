from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pymzml

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

    def get_scan_df(self, data):
        """
        Return scan data as a pandas DataFrame.

        Parameters
        ----------
        data : pymzml.run.Reader
            The mass spectra data.

        Returns
        -------
        pandas.DataFrame
            A pandas DataFrame containing metadata for each scan, including scan number, MS level, polarity, and scan time.
        """
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

        return scan_df

    def get_ms_raw(self, spectra, scan_df, data):
        """Return a dictionary of mass spectra data as a pandas DataFrame.

        Parameters
        ----------
        spectra : str
            Which mass spectra data to include in the output.
            Options: None, "ms1", "ms2", "all".
        scan_df : pandas.DataFrame
            Scan dataframe. Output from get_scan_df().
        data : pymzml.run.Reader
            The mass spectra data.

        Returns
        -------
        dict
            A dictionary containing the mass spectra data as pandas DataFrames, with keys corresponding to the MS level.

        """
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
            if i in scan_df_forspec.scan:
                # Get ms level
                level = "ms{}".format(spec.ms_level)

                # Number of rows
                N[level] += spec.mz.shape[0]

        # Second pass: parse
        for i, spec in enumerate(data):
            if i in scan_df_forspec.scan:
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

    def run(self, spectra="all", scan_df=None):
        """Parse the mzML file and return a dictionary of spectra dataframes and a scan metadata dataframe.

        Parameters
        ----------
        spectra : str, optional
            Which mass spectra data to include in the output. Default is "all".
            Other options: None, "ms1", "ms2".
        scan_df : pandas.DataFrame, optional
            Scan dataframe.  If not provided, the scan dataframe is created from the mzML file.

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
            scan_df = self.get_scan_df(data)

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
            mass_spectrum_obj = mass_spectrum_obj = MassSpecProfile(
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

        return mass_spectrum_obj

    def get_mass_spectra_obj(self):
        """Instatiate a MassSpectraBase object from the mzML file.


        Returns
        -------
        MassSpectraBase
            The MassSpectra object containing the parsed mass spectra.
            The object is instatiated with the mzML file, analyzer, instrument, sample name, and scan dataframe.
        """
        _, scan_df = self.run(spectra=False)
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

    def get_lcms_obj(self, spectra="all"):
        """Instatiates a LCMSBase object from the mzML file.

        Parameters
        ----------
        spectra : str, optional
            Which mass spectra data to include in the output. Default is all.  Other options: none, ms1, ms2.

        Returns
        -------
        LCMSBase
            LCMS object containing mass spectra data.
            The object is instatiated with the mzML file, analyzer, instrument, sample name, scan dataframe,
            and mz dataframe(s), as well as lists of scan numbers, retention times, and TICs.
        """
        _, scan_df = self.run(spectra="none")  # first run it to just get scan info
        res, scan_df = self.run(
            scan_df=scan_df, spectra=spectra
        )  # second run to parse data
        lcms_obj = LCMSBase(
            self.file_location,
            self.analyzer,
            self.instrument_label,
            self.sample_name,
            self,
        )
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
