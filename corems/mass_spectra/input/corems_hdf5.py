__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"


from threading import Thread
import h5py
import toml
import json
import multiprocessing
from pathlib import Path

import pandas as pd
import warnings

from corems.chroma_peak.factory.chroma_peak_classes import LCMSMassFeature
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_json_parameters_lcms,
    load_and_set_toml_parameters_lcms,
)
from corems.mass_spectra.factory.lc_class import LCMSBase, MassSpectraBase, LCMSCollection
from corems.mass_spectra.factory.chromat_data import EIC_Data
from corems.mass_spectra.input.parserbase import SpectraParserInterface
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum
from corems.molecular_id.factory.spectrum_search_results import SpectrumSearchResults
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.mass_spectra.input.mzml import MZMLSpectraParser


class ReadCoreMSHDFMassSpectra(
    SpectraParserInterface, ReadCoreMSHDF_MassSpectrum, Thread
):
    """Class to read CoreMS HDF5 files and populate a LCMS or MassSpectraBase object.

    Parameters
    ----------
    file_location : str
        The location of the HDF5 file to read, including the suffix.

    Attributes
    ----------
    file_location : str
        The location of the HDF5 file to read.
    h5pydata : h5py.File
        The HDF5 file object.
    scans : list
        A list of the location of individual mass spectra within the HDF5 file.
    scan_number_list : list
        A list of the scan numbers of the mass spectra within the HDF5 file.
    parameters_location : str
        The location of the parameters file (json or toml).

    Methods
    -------
    * import_mass_spectra(mass_spectra).
        Imports all mass spectra from the HDF5 file onto the LCMS or MassSpectraBase object.
    * get_mass_spectrum_from_scan(scan_number).
        Return mass spectrum data object from scan number.
    * load().
        Placeholder method to meet the requirements of the SpectraParserInterface.
    * run(mass_spectra).
        Runs the importer functions to populate a LCMS or MassSpectraBase object.
    * import_scan_info(mass_spectra).
        Imports the scan info from the HDF5 file to populate the _scan_info attribute
        on the LCMS or MassSpectraBase object
    * import_ms_unprocessed(mass_spectra).
        Imports the unprocessed mass spectra from the HDF5 file to populate the
        _ms_unprocessed attribute on the LCMS or MassSpectraBase object
    * import_parameters(mass_spectra).
        Imports the parameters from the HDF5 file to populate the parameters
        attribute on the LCMS or MassSpectraBase object
    * import_mass_features(mass_spectra).
        Imports the mass features from the HDF5 file to populate the mass_features
        attribute on the LCMS or MassSpectraBase object
    * import_eics(mass_spectra).
        Imports the extracted ion chromatograms from the HDF5 file to populate the
        eics attribute on the LCMS or MassSpectraBase object
    * import_spectral_search_results(mass_spectra).
        Imports the spectral search results from the HDF5 file to populate the
        spectral_search_results attribute on the LCMS or MassSpectraBase object
    * get_mass_spectra_obj().
        Return mass spectra data object, populating the _ms list on the LCMS or
        MassSpectraBase object from the HDF5 file
    * get_lcms_obj().
        Return LCMSBase object, populating the majority of the attributes on the
        LCMS object from the HDF5 file

    """

    def __init__(self, file_location: str):
        Thread.__init__(self)
        ReadCoreMSHDF_MassSpectrum.__init__(self, file_location)

        # override the scans attribute on ReadCoreMSHDF_MassSpectrum class to expect a nested location within the HDF5 file
        self.scans = [
            "mass_spectra/" + x for x in list(self.h5pydata["mass_spectra"].keys())
        ]
        self.scan_number_list = sorted(
            [int(float(i)) for i in list(self.h5pydata["mass_spectra"].keys())]
        )

        # set the location of the parameters file (json or toml)
        add_files = [
            x
            for x in self.file_location.parent.glob(
                self.file_location.name.replace(".hdf5", ".*")
            )
            if x.suffix != ".hdf5"
        ]
        if len([x for x in add_files if x.suffix == ".json"]) > 0:
            self.parameters_location = [x for x in add_files if x.suffix == ".json"][0]
        elif len([x for x in add_files if x.suffix == ".toml"]) > 0:
            self.parameters_location = [x for x in add_files if x.suffix == ".toml"][0]
        else:
            self.parameters_location = None

    def get_mass_spectrum_from_scan(self, scan_number):
        """Return mass spectrum data object from scan number."""
        if scan_number in self.scan_number_list:
            mass_spec = self.get_mass_spectrum(scan_number)
            return mass_spec
        else:
            raise Exception("Scan number not found in HDF5 file.")

    def load(self) -> None:
        """ """
        pass

    def get_ms_raw(self, spectra=None, scan_df=None) -> dict:
        """ """
        # Warn if spectra or scan_df are not None that they are not used for CoreMS HDF5 files and should be rerun after instantiation
        if spectra is not None or scan_df is not None:
            SyntaxWarning(
                "get_ms_raw method for CoreMS HDF5 files can only access saved data, consider rerunning after instantiation."
            )
        ms_unprocessed = {}
        dict_group_load = self.h5pydata["ms_unprocessed"]
        dict_group_keys = dict_group_load.keys()
        for k in dict_group_keys:
            ms_up_int = dict_group_load[k][:]
            ms_unprocessed[int(k)] = pd.DataFrame(
                ms_up_int, columns=["scan", "mz", "intensity"]
            )
        return ms_unprocessed

    def get_scan_df(self) -> pd.DataFrame:
        scan_info = {}
        dict_group_load = self.h5pydata["scan_info"]
        dict_group_keys = dict_group_load.keys()
        for k in dict_group_keys:
            scan_info[k] = dict_group_load[k][:]
        scan_df = pd.DataFrame(scan_info)
        scan_df.set_index("scan", inplace=True, drop=False)
        str_df = scan_df.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            scan_df[col] = str_df[col]
        return scan_df

    def run(self, mass_spectra, load_raw=True, load_light=False) -> None:
        """Runs the importer functions to populate a LCMS or MassSpectraBase object.

        Notes
        -----
        The following functions are run in order, if the HDF5 file contains the necessary data:
        1. import_parameters(), which populates the parameters attribute on the LCMS or MassSpectraBase object.
        2. import_mass_spectra(), which populates the _ms list on the LCMS or MassSpectraBase object.
        3. import_scan_info(), which populates the _scan_info on the LCMS or MassSpectraBase object.
        4. import_ms_unprocessed(), which populates the _ms_unprocessed attribute on the LCMS or MassSpectraBase object.
        5. import_mass_features(), which populates the mass_features attribute on the LCMS or MassSpectraBase object.
        6. import_eics(), which populates the eics attribute on the LCMS or MassSpectraBase object.
        7. import_spectral_search_results(), which populates the spectral_search_results attribute on the LCMS or MassSpectraBase object.

        Parameters
        ----------
        mass_spectra : LCMSBase or MassSpectraBase
            The LCMS or MassSpectraBase object to populate with mass spectra, generally instantiated with only the file_location, analyzer, and instrument_label attributes.
        load_raw : bool
            If True, load raw data (unprocessed) from HDF5 files for overall lcms object and individual mass spectra. Default is True.
        load_light : bool
            If True, only load the parameters, mass features, and scan info. Default is False.

        Returns
        -------
        None, but populates several attributes on the LCMS or MassSpectraBase object.

        """
        if self.parameters_location is not None:
            # Populate the parameters attribute on the LCMS object
            self.import_parameters(mass_spectra)

        if "mass_spectra" in self.h5pydata and not load_light:
            # Populate the _ms list on the LCMS object
            self.import_mass_spectra(mass_spectra, load_raw=load_raw)

        if "scan_info" in self.h5pydata:
            # Populate the _scan_info attribute on the LCMS object
            self.import_scan_info(mass_spectra)

        if "ms_unprocessed" in self.h5pydata and load_raw and not load_light:
            # Populate the _ms_unprocessed attribute on the LCMS object
            self.import_ms_unprocessed(mass_spectra)

        if "mass_features" in self.h5pydata:
            # Populate the mass_features attribute on the LCMS object
            self.import_mass_features(mass_spectra)

        if "eics" in self.h5pydata and not load_light:
            # Populate the eics attribute on the LCMS object
            self.import_eics(mass_spectra)

        if "spectral_search_results" in self.h5pydata and not load_light:
            # Populate the spectral_search_results attribute on the LCMS object
            self.import_spectral_search_results(mass_spectra)

    def import_mass_spectra(self, mass_spectra, load_raw=True) -> None:
        """Imports all mass spectra from the HDF5 file.

        Parameters
        ----------
        mass_spectra : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase object to populate with mass spectra.
        load_raw : bool
            If True, load raw data (unprocessed) from HDF5 files for overall lcms object and individual mass spectra. Default

        Returns
        -------
        None, but populates the '_ms' list on the LCMSBase or MassSpectraBase
        object with mass spectra from the HDF5 file.
        """
        for scan_number in self.scan_number_list:
            mass_spec = self.get_mass_spectrum(scan_number, load_raw=load_raw)
            mass_spec.scan_number = scan_number
            mass_spectra.add_mass_spectrum(mass_spec)

    def import_scan_info(self, mass_spectra) -> None:
        """Imports the scan info from the HDF5 file.

        Parameters
        ----------
        lcms : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase objects

        Returns
        -------
        None, but populates the 'scan_df' attribute on the LCMSBase or MassSpectraBase
        object with a pandas DataFrame of the 'scan_info' from the HDF5 file.

        """
        scan_df = self.get_scan_df()
        mass_spectra.scan_df = scan_df

    def import_ms_unprocessed(self, mass_spectra) -> None:
        """Imports the unprocessed mass spectra from the HDF5 file.

        Parameters
        ----------
        lcms : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase objects

        Returns
        -------
        None, but populates the '_ms_unprocessed' attribute on the LCMSBase or MassSpectraBase
        object with a dictionary of the 'ms_unprocessed' from the HDF5 file.

        """
        ms_unprocessed = self.get_ms_raw()
        mass_spectra._ms_unprocessed = ms_unprocessed

    def import_parameters(self, mass_spectra) -> None:
        """Imports the parameters from the HDF5 file.

        Parameters
        ----------
        mass_spectra : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase object to populate with parameters.

        Returns
        -------
        None, but populates the 'parameters' attribute on the LCMS or MassSpectraBase
        object with a dictionary of the 'parameters' from the HDF5 file.

        """
        if ".json" == self.parameters_location.suffix:
            load_and_set_json_parameters_lcms(mass_spectra, self.parameters_location)
        if ".toml" == self.parameters_location.suffix:
            load_and_set_toml_parameters_lcms(mass_spectra, self.parameters_location)
        else:
            raise Exception(
                "Parameters file must be in JSON format, TOML format is not yet supported."
            )

    def import_mass_features(self, mass_spectra) -> None:
        """Imports the mass features from the HDF5 file.

        Parameters
        ----------
        mass_spectra : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase object to populate with mass features.

        Returns
        -------
        None, but populates the 'mass_features' attribute on the LCMSBase or MassSpectraBase
        object with a dictionary of the 'mass_features' from the HDF5 file.

        """
        dict_group_load = self.h5pydata["mass_features"]
        dict_group_keys = dict_group_load.keys()
        for k in dict_group_keys:
            # Instantiate the MassFeature object
            mass_feature = LCMSMassFeature(
                mass_spectra,
                mz=dict_group_load[k].attrs["_mz_exp"],
                retention_time=dict_group_load[k].attrs["_retention_time"],
                intensity=dict_group_load[k].attrs["_intensity"],
                apex_scan=dict_group_load[k].attrs["_apex_scan"],
                persistence=dict_group_load[k].attrs["_persistence"],
                id=int(k),
            )

            # Populate additional attributes on the MassFeature object
            for key in dict_group_load[k].attrs.keys() - {
                "_mz_exp",
                "_mz_cal",
                "_retention_time",
                "_intensity",
                "_apex_scan",
                "_persistence",
            }:
                setattr(mass_feature, key, dict_group_load[k].attrs[key])

            # Populate attributes on MassFeature object that are lists
            for key in dict_group_load[k].keys():
                setattr(mass_feature, key, dict_group_load[k][key][:])

            mass_spectra.mass_features[int(k)] = mass_feature

        # Associate mass features with ms1 and ms2 spectra, if available
        for mf_id in mass_spectra.mass_features.keys():
            if mass_spectra.mass_features[mf_id].apex_scan in mass_spectra._ms.keys():
                mass_spectra.mass_features[mf_id].mass_spectrum = mass_spectra._ms[
                    mass_spectra.mass_features[mf_id].apex_scan
                ]
            if mass_spectra.mass_features[mf_id].ms2_scan_numbers is not None:
                for ms2_scan in mass_spectra.mass_features[mf_id].ms2_scan_numbers:
                    if ms2_scan in mass_spectra._ms.keys():
                        mass_spectra.mass_features[mf_id].ms2_mass_spectra[ms2_scan] = (
                            mass_spectra._ms[ms2_scan]
                        )

    def import_eics(self, mass_spectra):
        """Imports the extracted ion chromatograms from the HDF5 file.

        Parameters
        ----------
        mass_spectra : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase object to populate with extracted ion chromatograms.

        Returns
        -------
        None, but populates the 'eics' attribute on the LCMSBase or MassSpectraBase
        object with a dictionary of the 'eics' from the HDF5 file.

        """
        dict_group_load = self.h5pydata["eics"]
        dict_group_keys = dict_group_load.keys()
        for k in dict_group_keys:
            my_eic = EIC_Data(
                scans=dict_group_load[k]["scans"][:],
                time=dict_group_load[k]["time"][:],
                eic=dict_group_load[k]["eic"][:],
            )
            for key in dict_group_load[k].keys():
                if key not in ["scans", "time", "eic"]:
                    setattr(my_eic, key, dict_group_load[k][key][:])
                    # if key is apexes, convert to a tuple of a list
                    if key == "apexes" and len(my_eic.apexes) > 0:
                        my_eic.apexes = [tuple(x) for x in my_eic.apexes]
            # Add to mass_spectra object
            mass_spectra.eics[dict_group_load[k].attrs["mz"]] = my_eic

        # Add to mass features
        for idx in mass_spectra.mass_features.keys():
            mz = mass_spectra.mass_features[idx].mz
            if mz in mass_spectra.eics.keys():
                mass_spectra.mass_features[idx]._eic_data = mass_spectra.eics[mz]

    def import_spectral_search_results(self, mass_spectra):
        """Imports the spectral search results from the HDF5 file.

        Parameters
        ----------
        mass_spectra : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase object to populate with spectral search results.

        Returns
        -------
        None, but populates the 'spectral_search_results' attribute on the LCMSBase or MassSpectraBase
        object with a dictionary of the 'spectral_search_results' from the HDF5 file.

        """
        overall_results_dict = {}
        ms2_results_load = self.h5pydata["spectral_search_results"]
        for k in ms2_results_load.keys():
            overall_results_dict[int(k)] = {}
            for k2 in ms2_results_load[k].keys():
                ms2_search_res = SpectrumSearchResults(
                    query_spectrum=mass_spectra._ms[int(k)],
                    precursor_mz=ms2_results_load[k][k2].attrs["precursor_mz"],
                    spectral_similarity_search_results={},
                )

                for key in ms2_results_load[k][k2].keys() - {"precursor_mz"}:
                    setattr(ms2_search_res, key, list(ms2_results_load[k][k2][key][:]))
                overall_results_dict[int(k)][
                    ms2_results_load[k][k2].attrs["precursor_mz"]
                ] = ms2_search_res

        # add to mass_spectra
        mass_spectra.spectral_search_results.update(overall_results_dict)

        # If there are mass features, associate the results with each mass feature
        if len(mass_spectra.mass_features) > 0:
            for mass_feature_id, mass_feature in mass_spectra.mass_features.items():
                scan_ids = mass_feature.ms2_scan_numbers
                for ms2_scan_id in scan_ids:
                    precursor_mz = mass_feature.mz
                    try:
                        mass_spectra.spectral_search_results[ms2_scan_id][precursor_mz]
                    except KeyError:
                        pass
                    else:
                        mass_spectra.mass_features[
                            mass_feature_id
                        ].ms2_similarity_results.append(
                            mass_spectra.spectral_search_results[ms2_scan_id][
                                precursor_mz
                            ]
                        )

    def get_mass_spectra_obj(self, load_raw=True, load_light=False) -> MassSpectraBase:
        """
        Return mass spectra data object, populating the _ms list on MassSpectraBase object from the HDF5 file.

        Parameters
        ----------
        load_raw : bool
            If True, load raw data (unprocessed) from HDF5 files for overall spectra object and individual mass spectra. Default is True.
        load_light : bool
            If True, only load the parameters, mass features, and scan info. Default is False.

        """
        # Instantiate the LCMS object
        spectra_obj = MassSpectraBase(
            file_location=self.file_location,
            analyzer=self.analyzer,
            instrument_label=self.instrument_label,
            sample_name=self.sample_name,
        )

        # This will populate the _ms list on the LCMS or MassSpectraBase object
        self.run(spectra_obj, load_raw=load_raw, load_light=load_light)

        return spectra_obj

    def get_lcms_obj(
        self, load_raw=True, load_light=False, use_original_parser=True, raw_file_path=None
    ) -> LCMSBase:
        """
        Return LCMSBase object, populating attributes on the LCMSBase object from the HDF5 file.

        Parameters
        ----------
        load_raw : bool
            If True, load raw data (unprocessed) from HDF5 files for overall lcms object and individual mass spectra. Default is True.
        load_light : bool
            If True, only load the parameters, mass features, and scan info. Default is False.
        use_original_parser : bool
            If True, use the original parser to populate the LCMS object. Default is True.
        raw_file_path : str
            The location of the raw file to parse if attempting to use original parser.
            Default is None, which attempts to get the raw file path from the HDF5 file.
            If the original file path has moved, this parameter can be used to specify the new location.
        """
        # Instantiate the LCMS object
        lcms_obj = LCMSBase(
            file_location=self.file_location,
            analyzer=self.analyzer,
            instrument_label=self.instrument_label,
            sample_name=self.sample_name,
        )

        # This will populate the majority of the attributes on the LCMS object
        self.run(lcms_obj, load_raw=load_raw, load_light=load_light)

        # Set final attributes of the LCMS object
        lcms_obj.polarity = self.h5pydata.attrs["polarity"]
        lcms_obj._scans_number_list = list(lcms_obj.scan_df.scan)
        lcms_obj._retention_time_list = list(lcms_obj.scan_df.scan_time)
        lcms_obj._tic_list = list(lcms_obj.scan_df.tic)

        # If use_original_parser is True, instantiate the original parser and populate the LCMS object
        if use_original_parser:
            lcms_obj = self.add_original_parser(lcms_obj, raw_file_path=raw_file_path)
        else:
            lcms_obj.spectra_parser_class = self.__class__

        return lcms_obj

    def get_raw_file_location(self):
        """
        Get the raw file location from the HDF5 file attributes.

        Returns
        -------
        str
            The raw file location.
        """
        if "original_file_location" in self.h5pydata.attrs:
            return self.h5pydata.attrs["original_file_location"]
        else:
            return None
    
    def add_original_parser(self, mass_spectra, raw_file_path=None):
        """
        Add the original parser to the mass spectra object.

        Parameters
        ----------
        mass_spectra : MassSpectraBase | LCMSBase
            The MassSpectraBase or LCMSBase object to add the original parser to.
        raw_file_path : str
            The location of the raw file to parse. Default is None, which attempts to get the raw file path from the HDF5 file.
        """
        # Get the original parser type
        og_parser_type = self.h5pydata.attrs["parser_type"]

        # If raw_file_path is None, get it from the HDF5 file attributes
        if raw_file_path is None:
            raw_file_path = self.get_raw_file_location()
            if raw_file_path is None:
                raise ValueError(
                    "Raw file path not found in HDF5 file attributes, cannot instantiate original parser."
                )
            
        # Set the raw file path on the mass_spectra object so the parser knows where to find the raw file
        mass_spectra.raw_file_location = raw_file_path

        if og_parser_type == "ImportMassSpectraThermoMSFileReader":
            # Check that the parser can be instantiated with the raw file path
            parser = ImportMassSpectraThermoMSFileReader(raw_file_path)
        elif og_parser_type == "MZMLSpectraParser":
            # Check that the parser can be instantiated with the raw file path
            parser = MZMLSpectraParser(raw_file_path)

        # Set the spectra parser class on the mass_spectra object so the spectra_parser property can be used with the original parser
        mass_spectra.spectra_parser_class = parser.__class__

        return mass_spectra
    
    def get_creation_time(self):
        """
        Raise a NotImplemented Warning, as creation time is not available in CoreMS HDF5 files and returning None.
        """
        warnings.warn(
            "Creation time is not available in CoreMS HDF5 files, returning None." \
            "This should be accessed through the original parser.",
        )
        return None
    
    def get_instrument_info(self):
        """
        Raise a NotImplemented Warning, as instrument info is not available in CoreMS HDF5 files and returning None.
        """
        warnings.warn(
            "Instrument info is not available in CoreMS HDF5 files, returning None." \
            "This should be accessed through the original parser.",
        )
        return None


class ReadCoreMSHDFMassSpectraCollection:
    """Class to read a collection of CoreMS HDF5 files and populate a LCMSCollection object.
    
    Parameters
    ----------
    folder_location : str
        The location of the folder containing the CoreMS HDF5 files.
    manifest_file : str
        The location of the manifest file containing the sample names, order, and batch.
        This must be a csv with the following columns: 'sample_name', 'order', 'batch'.
        Other fields can be included in the manifest file, but these are required.
    cores : int
        The number of cores to use for multiprocessing. Default is 1.

    Attributes
    ----------
    folder_location : str
        The location of the folder containing the CoreMS HDF5 files.
    manifest_filepath : str
        The location of the manifest file containing the sample names, order, and batch.
    """
    def __init__(
            self, 
            folder_location: str, 
            manifest_file: str, 
            cores: int = 1
            ):
        # Check for folder location and manifest file
        if not folder_location.exists():
            raise FileNotFoundError(f"Folder location {folder_location} not found.")
        if not manifest_file.exists():
            raise FileNotFoundError(f"Manifest file {manifest_file} not found.")

        # Check if the manifest file is a CSV
        if manifest_file.suffix != ".csv":
            raise ValueError("Manifest file must be a CSV.")

        self.folder_location = folder_location
        self._manifest_dict = None
        self._parse_manifest(manifest_file)
        self._validate_manifest()
        self._validate_parameters()
        self._validate_cores(cores)
    
    def _validate_cores(self, cores):
        # Check if the cores parameter is an integer greater than 0 and less than the number of cores available
        if not isinstance(cores, int) or cores < 1:
            raise ValueError("Cores must be an integer greater than 0.")
        if cores > multiprocessing.cpu_count():
            raise ValueError(
                f"Cores must be less than or equal to the number of cores available ({multiprocessing.cpu_count()})."
            )
        self._cores = cores

    def _parse_manifest(self, manifest_file):
        """Parse the manifest file and set the manifest dictionary."""
        self.manifest_filepath = manifest_file
        manifest = pd.read_csv(manifest_file)
        # Check if the following columns exisit in the manifest file
        if not all(
            col in manifest.columns for col in ["sample_name", "order", "batch"]
        ):
            raise ValueError(
                "Manifest file must contain the following columns: 'sample_name', 'order', 'batch'."
            )
        # Set index to the 'sample_name' column
        manifest.set_index("sample_name", inplace=True)
        self._manifest_dict = manifest.to_dict(orient="index")

    def _validate_manifest(self):
        """Validate the manifest dictionary against the CoreMS folder location."""
        # Check if the folder location contains HDF5 files for each sample
        for sample_name in self._manifest_dict.keys():
            corems_dir = self.folder_location / f"{sample_name}.corems"
            if not corems_dir.exists():
                raise FileNotFoundError(f"CoreMS folder for {sample_name} not found.")
            hdf5_file = corems_dir / f"{sample_name}.hdf5"
            if not hdf5_file.exists():
                raise FileNotFoundError(f"HDF5 file for {sample_name} not found.")

    def _validate_parameters(self):
        """Validate that the parameters used for all samples within a batch are the same."""
        # Check if parameters files are saved as JSON or TOML
        if self.parameters_files[0].suffix == ".json":
            importer = json
            suffix = ".json"

        elif self.parameters_files[0].suffix == ".toml":
            importer = toml
            suffix = ".toml"

        manfiest_df = self.manifest_dataframe

        # Split up samples by batch
        batches = manfiest_df["batch"].unique()

        for batch in batches:
            samples = manfiest_df[manfiest_df["batch"] == batch].index
            # check if self.parameters_files end with .json or .toml
            batch_param_files = [
                self.folder_location / f"{sample_name}.corems/{sample_name}{suffix}"
                for sample_name in self._manifest_dict.keys()
                if sample_name in samples
            ]
            with open(
                batch_param_files[0],
                "r",
                encoding="utf8",
            ) as stream:
                first_parameters = importer.load(stream)
            for parameters_file in batch_param_files[1:]:
                with open(
                    parameters_file,
                    "r",
                    encoding="utf8",
                ) as stream:
                    parameters = importer.load(stream)
                if parameters != first_parameters:
                    raise ValueError(
                        f"Parameters files for samples in batch {batch} are not equal."
                    )        
    
    def get_lcms_obj(self, sample_name: str, load_raw=False, load_light=True, use_original_parser=True, raw_file_path=None) -> LCMSBase:
        """Return a LCMSBase object for a given sample name within the collection.
        
        Parameters
        ----------
        sample_name : str
            The sample name to retrieve the LCMS object for.
        load_raw : bool
            If True, load raw data from HDF5 files. Default is False. 
        load_light : bool
            If True, only load the parameters, mass features, and scan info are initially loaded for each lcms object. Default is True.   
        """
        hdf5_file = self.folder_location / f"{sample_name}.corems/{sample_name}.hdf5"
        parser = ReadCoreMSHDFMassSpectra(hdf5_file)
        lcms_obj = parser.get_lcms_obj(load_raw=load_raw, load_light=load_light, use_original_parser=use_original_parser, raw_file_path=raw_file_path)
        if load_light:
            mf_df = lcms_obj.mass_features_to_df()
            lcms_obj.mass_features = {}
            lcms_obj.light_mf_df = mf_df
        return lcms_obj
    
    def get_lcms_collection(self, load_raw = False, load_light = True) -> LCMSCollection:
        """Return a LCMSCollection object
        
        Parameters
        ----------
        load_raw : bool
            If True, load raw data from HDF5 files. Default is False. 
        load_light : bool
            If True, only load the parameters, mass features, and scan info are initially loaded for each lcms object. 
            After concatenating the mass_features, remove the mass_features attribute from the individual LCMS objects for memory efficiency. Default is True.
            Default is True.   
        """
        # Instantiate the LCMSCollection object
        lcms_coll = LCMSCollection(
            collection_location=self.folder_location,
            manifest=self.manifest,
            collection_parser=self
        )

        # Set the number of cores on the LCMSCollection object from the ReadCoreMSHDFMassSpectraCollection object
        lcms_coll.parameters.lcms_collection.cores = self._cores

        # Add LCMS objects to the collection
        samples = self._manifest_dict.keys()

        # Initialize the LCMS object dictionary
        if self._cores > 1:
            if self._cores > len(samples):
                ncores = len(samples)
            else:
                ncores = self._cores
            # Create a pool of workers (one for each core or sample, whichever is smaller)
            pool = multiprocessing.Pool(ncores)
            # Load the LCMS objects in parallel - do not instantiate the original parser by default
            use_original_parser = True
            args = [(sample, load_raw, load_light, use_original_parser) for sample in samples]
            lcms_objs = pool.starmap(self.get_lcms_obj, args)
            for sample_name, lcms_obj in zip(samples, lcms_objs):
                lcms_coll._lcms[sample_name] = lcms_obj

        elif self._cores == 1:
            # Load the LCMS objects sequentially - do not instantiate the original parser by default
            for sample_name in samples:
                lcms_coll._lcms[sample_name] = self.get_lcms_obj(sample_name, load_raw=load_raw, load_light=load_light, use_original_parser=False)

        else:
            raise ValueError("Number of cores must be greater than 0 and set on the ReadCoreMSHDFMassSpectraCollection object.")

        # Check that all LCMS objects have the same polarity
        if len(set([x.polarity for k, x in lcms_coll._lcms.items()])) != 1:
            raise ValueError("All samples must have the same polarity.")
        
        # Set ids on the LCMS objects in the manifest
        i = 0
        for sample in lcms_coll.samples:
            lcms_coll._manifest_dict[sample]["collection_id"] = i
            i += 1
        
        # Reorder the LCMS objects
        lcms_coll._reorder_lcms_objects()

        # Collect the mass features from the LCMS objects and combine them into a single dataframe for the collection
        lcms_coll._combine_mass_features()
        
        # If load_light, remove the mass_feature attribute from the individual LCMS objects
        if load_light:
            for sample_name in lcms_coll.samples:
                lcms_coll._lcms[sample_name].mass_features = {}
                # Remove the light_mf_df attribute from the individual LCMS objects
                del lcms_coll._lcms[sample_name].light_mf_df


        return lcms_coll

    @property
    def manifest(self):
        return self._manifest_dict

    @property
    def manifest_dataframe(self):
        return pd.DataFrame(self._manifest_dict).T

    @property
    def hdf5_files(self):
        return [
            self.folder_location / f"{sample_name}.corems/{sample_name}.hdf5"
            for sample_name in self._manifest_dict.keys()
        ]

    @property
    def parameters_files(self):
        # Check if parameters files are saved as JSON or TOML
        json_files = [
            self.folder_location / f"{sample_name}.corems/{sample_name}.json"
            for sample_name in self._manifest_dict.keys()
        ]
        toml_files = [
            self.folder_location / f"{sample_name}.corems/{sample_name}.toml"
            for sample_name in self._manifest_dict.keys()
        ]
        if all([x.exists() for x in json_files]):
            return json_files
        elif all([x.exists() for x in toml_files]):
            return toml_files
        else:
            raise ValueError("Parameters files are not saved for all samples.")

class ReadSavedLCMSCollection(ReadCoreMSHDFMassSpectraCollection):
    """
    Subclass to read and re-instantiate a LCMSCollection from a saved HDF5 file.
    
    
    Parameters
    ----------
    collection_hdf5_path : str or Path
        Path to the saved LCMSCollection HDF5 file.
    cores : int, optional
        Number of cores for processing. Default is 1.
    """
    
    def __init__(
        self, 
        collection_hdf5_path: str, 
        cores: int = 1
    ):
        # Convert to Path objects
        self.collection_hdf5_path = Path(collection_hdf5_path)
        
        # Validate the collection file exists
        if not self.collection_hdf5_path.exists():
            raise FileNotFoundError(f"Collection HDF5 file {self.collection_hdf5_path} not found.")
        
        # Validate cores
        self._validate_cores(cores)

        # Set data
        self.h5pydata = h5py.File(self.collection_hdf5_path, "r")

        # Load metadata from saved collection
        self._load_collection_metadata()

        if not self.folder_location.exists():
            raise FileNotFoundError(f"Folder location {self.folder_location} not found.")

        # Load the mass spectra data
        self._validate_manifest()

    def _load_collection_metadata(self):
        """Load metadata and manifest from the saved collection HDF5 file."""
        with h5py.File(self.collection_hdf5_path, 'r') as f:
            self.folder_location = Path(f.attrs.get('lcms_objects_folder', ''))
            manifest_json = f.attrs.get('manifest', '{}')
            if isinstance(manifest_json, bytes):
                manifest_json = manifest_json.decode('utf-8')
            self._manifest_dict = json.loads(manifest_json)
