__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"


from threading import Thread
from pathlib import Path

import pandas as pd

from corems.chroma_peak.factory.chroma_peak_classes import LCMSMassFeature
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_json_parameters_lcms,
    load_and_set_toml_parameters_lcms,
)
from corems.mass_spectra.factory.lc_class import LCMSBase, MassSpectraBase
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

    def run(self, mass_spectra, load_raw=True) -> None:
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
        Returns
        -------
        None, but populates several attributes on the LCMS or MassSpectraBase object.

        """
        if self.parameters_location is not None:
            # Populate the parameters attribute on the LCMS object
            self.import_parameters(mass_spectra)

        if "mass_spectra" in self.h5pydata:
            # Populate the _ms list on the LCMS object
            self.import_mass_spectra(mass_spectra, load_raw=load_raw)

        if "scan_info" in self.h5pydata:
            # Populate the _scan_info attribute on the LCMS object
            self.import_scan_info(mass_spectra)

        if "ms_unprocessed" in self.h5pydata and load_raw:
            # Populate the _ms_unprocessed attribute on the LCMS object
            self.import_ms_unprocessed(mass_spectra)

        if "mass_features" in self.h5pydata:
            # Populate the mass_features attribute on the LCMS object
            self.import_mass_features(mass_spectra)

        if "eics" in self.h5pydata:
            # Populate the eics attribute on the LCMS object
            self.import_eics(mass_spectra)

        if "spectral_search_results" in self.h5pydata:
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

    def get_mass_spectra_obj(self, load_raw=True) -> MassSpectraBase:
        """
        Return mass spectra data object, populating the _ms list on MassSpectraBase object from the HDF5 file.

        Parameters
        ----------
        load_raw : bool
            If True, load raw data (unprocessed) from HDF5 files for overall spectra object and individual mass spectra. Default is True.

        """
        # Instantiate the LCMS object
        spectra_obj = MassSpectraBase(
            file_location=self.file_location,
            analyzer=self.analyzer,
            instrument_label=self.instrument_label,
            sample_name=self.sample_name,
        )

        # This will populate the _ms list on the LCMS or MassSpectraBase object
        self.run(spectra_obj, load_raw=load_raw)

        return spectra_obj

    def get_lcms_obj(
        self, load_raw=True, use_original_parser=True, raw_file_path=None
    ) -> LCMSBase:
        """
        Return LCMSBase object, populating attributes on the LCMSBase object from the HDF5 file.

        Parameters
        ----------
        load_raw : bool
            If True, load raw data (unprocessed) from HDF5 files for overall lcms object and individual mass spectra. Default is True.
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
        self.run(lcms_obj, load_raw=load_raw)

        # Set final attributes of the LCMS object
        lcms_obj.polarity = self.h5pydata.attrs["polarity"]
        lcms_obj._scans_number_list = list(lcms_obj.scan_df.scan)
        lcms_obj._retention_time_list = list(lcms_obj.scan_df.scan_time)
        lcms_obj._tic_list = list(lcms_obj.scan_df.tic)

        # If use_original_parser is True, instantiate the original parser and populate the LCMS object
        if use_original_parser:
            lcms_obj = self.add_original_parser(lcms_obj, raw_file_path=raw_file_path)

        return lcms_obj

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
        # Try to get the raw file path from the HDF5 file
        if raw_file_path is None:
            raw_file_path = self.h5pydata.attrs["original_file_location"]
            # Check if og_file_location exists, if not raise an error
            raw_file_path = self.h5pydata.attrs["original_file_location"]

        raw_file_path = Path(raw_file_path)
        if not raw_file_path.exists():
            raise FileExistsError(
                "File does not exist: " + str(raw_file_path),
                ". Cannot use original parser for instatiating the lcms_obj.",
            )

        # Get the original parser type
        og_parser_type = self.h5pydata.attrs["parser_type"]

        if og_parser_type == "ImportMassSpectraThermoMSFileReader":
            parser = ImportMassSpectraThermoMSFileReader(raw_file_path)
        elif og_parser_type == "MZMLSpectraParser":
            parser = MZMLSpectraParser(raw_file_path)

        mass_spectra.spectra_parser_class = parser.__class__
        mass_spectra.spectra_parser = parser

        return mass_spectra
