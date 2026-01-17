__author__ = "Yuri E. Corilo"
__date__ = "Oct 29, 2019"


from threading import Thread
import h5py
import toml
import json
import multiprocessing
from pathlib import Path
import datetime

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


def create_manifest_from_folder(
    folder_path: Path,
    output_path: Path = None,
    batch_time_threshold_hours: float = 12.0,
    center_name: str = None,
    overwrite: bool = False
) -> Path:
    """
    Create a manifest CSV file for ReadCoreMSHDFMassSpectraCollection from CoreMS HDF5 files.
    
    Scans a folder for .corems subdirectories and generates a manifest with columns:
    sample_name, batch, order, center, time. Files are batched by creation time, and
    one sample is designated as the retention time alignment center.
    
    Parameters
    ----------
    folder_path : Path
        Path to folder containing .corems subdirectories with HDF5 files.
    output_path : Path, optional
        Output manifest CSV path. Default: folder_path/manifest.csv.
    batch_time_threshold_hours : float, optional
        Time gap in hours for batch separation. Default: 12.0.
    center_name : str, optional
        Sample name to designate as RT alignment center (must exist in samples).
        If None, the middle sample (by creation time) is used.
    overwrite : bool, optional
        Whether to overwrite existing manifest. Default: False.
        
    Returns
    -------
    Path
        Path to created manifest file.
        
    Raises
    ------
    FileNotFoundError
        If folder_path doesn't exist or contains no .corems subdirectories.
    FileExistsError
        If output file exists and overwrite is False.
    ValueError
        If no HDF5 files found, or center_name doesn't match any sample.
    """
    if not folder_path.exists():
        raise FileNotFoundError(f"Folder {folder_path} does not exist.")
    
    # Set default output path if not provided
    if output_path is None:
        output_path = folder_path / "manifest.csv"
    
    # Check if output file exists
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Manifest file {output_path} already exists. "
            "Set overwrite=True to replace it."
        )
    
    # Find all .corems subdirectories
    corems_dirs = sorted([d for d in folder_path.iterdir() if d.is_dir() and d.suffix == ".corems"])
    
    if not corems_dirs:
        raise FileNotFoundError(
            f"No .corems subdirectories found in {folder_path}. "
            "Ensure the folder contains processed CoreMS data."
        )
    
    # Collect sample information
    sample_data = []
    
    for corems_dir in corems_dirs:
        sample_name = corems_dir.stem  # Remove .corems extension
        hdf5_file = corems_dir / f"{sample_name}.hdf5"
        
        if not hdf5_file.exists():
            print(f"Warning: HDF5 file not found for {sample_name}, skipping.")
            continue
        
        # Get creation time using the ReadCoreMSHDFMassSpectra method
        try:
            # Use context manager to ensure file is properly closed
            with ReadCoreMSHDFMassSpectra(str(hdf5_file)) as parser:
                # Use the get_original_creation_time() method which checks HDF5 attrs first,
                # then falls back to original parser if needed
                creation_time = parser.get_original_creation_time()
                
                # Skip sample if creation time unavailable
                if creation_time is None:
                    print(f"Warning: Could not get original creation time for {sample_name}, skipping.")
                    continue
            
        except Exception as e:
            print(f"Warning: Error getting creation time for {sample_name}: {e}, skipping.")
            continue
        
        sample_data.append({
            'sample_name': sample_name,
            'creation_time': creation_time,
            'hdf5_path': hdf5_file
        })
    
    if not sample_data:
        raise ValueError(
            f"No valid HDF5 files found in {folder_path}. "
            "Ensure .corems subdirectories contain .hdf5 files."
        )
    
    # Sort by creation time
    sample_data.sort(key=lambda x: x['creation_time'])
    
    # Assign batches based on time threshold
    batch_assignments = []
    current_batch = 1
    
    for i, sample in enumerate(sample_data):
        if i == 0:
            batch_assignments.append(current_batch)
        else:
            time_diff = sample['creation_time'] - sample_data[i-1]['creation_time']
            time_diff_hours = time_diff.total_seconds() / 3600
            
            if time_diff_hours > batch_time_threshold_hours:
                current_batch += 1
            
            batch_assignments.append(current_batch)
    
    # Determine which sample should be the center for retention time alignment
    sample_names = [s['sample_name'] for s in sample_data]
    
    if center_name is not None:
        # Validate that center_name is in the discovered samples
        if center_name not in sample_names:
            raise ValueError(
                f"Specified center_name '{center_name}' not found in discovered samples. "
                f"Available samples: {', '.join(sample_names)}"
            )
        center_sample = center_name
    else:
        # Use the middle sample (by creation time) as center
        middle_idx = len(sample_data) // 2
        center_sample = sample_data[middle_idx]['sample_name']
        print(f"Auto-selected center sample: {center_sample} (index {middle_idx} of {len(sample_data)}, middle by creation time)")
    
    # Create manifest dataframe with center column as TRUE/FALSE
    manifest_df = pd.DataFrame({
        'sample_name': sample_names,
        'batch': batch_assignments,
        'order': list(range(1, len(sample_data) + 1)),
        'center': ['TRUE' if name == center_sample else 'FALSE' for name in sample_names],
        'time': [s['creation_time'].strftime('%Y-%m-%dT%H:%M:%SZ') for s in sample_data]
    })
    
    # Sort manifest by time before saving to ensure proper order
    manifest_df = manifest_df.sort_values('time').reset_index(drop=True)
    # Update order column to reflect sorted order
    manifest_df['order'] = list(range(1, len(manifest_df) + 1))
    
    # Save manifest
    manifest_df.to_csv(output_path, index=False)
    
    print(f"Manifest created successfully at {output_path}")
    print(f"Total samples: {len(sample_data)}")
    print(f"Number of batches: {current_batch}")
    print(f"Batch assignments: {dict(zip(range(1, current_batch + 1), [batch_assignments.count(b) for b in range(1, current_batch + 1)]))}")
    
    return output_path


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
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - closes the HDF5 file."""
        if hasattr(self, 'h5pydata') and self.h5pydata is not None:
            self.h5pydata.close()
        return False
    
    def close(self):
        """Explicitly close the HDF5 file."""
        if hasattr(self, 'h5pydata') and self.h5pydata is not None:
            self.h5pydata.close()

    def get_mass_spectrum_from_scan(self, scan_number):
        """Return mass spectrum data object from scan number."""
        if scan_number in self.scan_number_list:
            mass_spec = self.get_mass_spectrum(scan_number)
            return mass_spec
        else:
            raise Exception("Scan number not found in HDF5 file.")

    def get_mass_spectra_from_scan_list(
        self, scan_list, spectrum_mode, auto_process=True
    ):
        """Return a list of mass spectrum data objects from a list of scan numbers.

        Parameters
        ----------
        scan_list : list
            A list of scan numbers to retrieve mass spectra for.
        spectrum_mode : str
            The spectrum mode to use when retrieving the mass spectra.
            Note that this parameter is not used for CoreMS HDF5 files, as the spectra are already processed and only
            centroided spectra are saved.
        auto_process : bool
            If True, automatically process the mass spectra when retrieving them.
            Note that this parameter is not used for CoreMS HDF5 files, as the spectra are already processed and only
            centroided spectra are saved.

        Returns
        -------
        list
            A list of mass spectrum data objects corresponding to the provided scan numbers.
        """
        mass_spectra_list = []
        for scan_number in scan_list:
            if scan_number in self.scan_number_list:
                mass_spec = self.get_mass_spectrum_from_scan(scan_number)
                mass_spectra_list.append(mass_spec)
            else:
                warnings.warn(f"Scan number {scan_number} not found in HDF5 file.")
        return mass_spectra_list

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

    def import_mass_features(self, mass_spectra, mf_ids=None) -> None:
        """Imports the mass features from the HDF5 file.

        Parameters
        ----------
        mass_spectra : LCMSBase | MassSpectraBase
            The MassSpectraBase or LCMSBase object to populate with mass features.
        mf_ids : list, optional
            A list of mass feature IDs to import. If None, all mass features are imported.

        Returns
        -------
        None, but populates the 'mass_features' attribute on the LCMSBase or MassSpectraBase
        object with a dictionary of the 'mass_features' from the HDF5 file.

        """
        dict_group_load = self.h5pydata["mass_features"]
        dict_group_keys = dict_group_load.keys()
        for k in dict_group_keys:
            if mf_ids is not None and int(k) not in mf_ids:
                continue
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
                # Convert _noise_score from array to tuple
                if key == "_noise_score":
                    mass_feature._noise_score = tuple(mass_feature._noise_score)
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
            parser_class = ImportMassSpectraThermoMSFileReader
        elif og_parser_type == "MZMLSpectraParser":
            # Check that the parser can be instantiated with the raw file path
            parser_class = MZMLSpectraParser

        # Set the spectra parser class on the mass_spectra object so the spectra_parser property can be used with the original parser
        mass_spectra.spectra_parser_class = parser_class

        return mass_spectra

    def get_original_creation_time(self):
        """
        Get the creation time of the original raw data file.
        
        First checks if creation_time is saved in the HDF5 file attributes.
        If not found, attempts to instantiate the original parser and get the creation time.
        
        Returns
        -------
        datetime
            The creation time of the original raw data file, or None if not available.
        """
        # Check if creation_time is saved in HDF5 attributes
        if "creation_time" in self.h5pydata.attrs:
            from datetime import datetime
            return datetime.fromisoformat(self.h5pydata.attrs["creation_time"])
        
        # Fall back to using original parser to get creation time
        try:
            # Get the original parser type and raw file path
            og_parser_type = self.h5pydata.attrs.get("parser_type")
            raw_file_path = self.get_raw_file_location()
            
            if og_parser_type is None or raw_file_path is None:
                warnings.warn(
                    "Cannot retrieve creation time: parser_type or original_file_location not found in HDF5 attributes."
                )
                return None
            
            # Check if raw file exists
            from pathlib import Path
            if not Path(raw_file_path).exists():
                warnings.warn(
                    f"Cannot retrieve creation time: original raw file not found at {raw_file_path}"
                )
                return None
            
            # Instantiate the original parser
            if og_parser_type == "ImportMassSpectraThermoMSFileReader":
                parser = ImportMassSpectraThermoMSFileReader(raw_file_path)
            elif og_parser_type == "MZMLSpectraParser":
                parser = MZMLSpectraParser(raw_file_path)
            else:
                warnings.warn(
                    f"Unknown parser type: {og_parser_type}, cannot retrieve creation time."
                )
                return None
            
            # Get creation time from parser
            return parser.get_creation_time()
            
        except Exception as e:
            warnings.warn(
                f"Failed to retrieve creation time from original parser: {e}"
            )
            return None
    
    def get_creation_time(self):
        """
        Get the creation time of the original raw data file.
        
        This is an alias for get_original_creation_time() for backward compatibility.
        
        Returns
        -------
        datetime
            The creation time of the original raw data file, or None if not available.
        """
        return self.get_original_creation_time()

    def get_instrument_info(self):
        """
        Raise a NotImplemented Warning, as instrument info is not available in CoreMS HDF5 files and returning None.
        """
        warnings.warn(
            "Instrument info is not available in CoreMS HDF5 files, returning None."
            "This should be accessed through the original parser.",
        )
        return None


class ReadCoreMSHDFMassSpectraCollection:
    """Read a collection of CoreMS HDF5 files and populate an LCMSCollection object.
    
    Parameters
    ----------
    folder_location : Path
        Folder containing .corems subdirectories with HDF5 files.
    manifest_file : Path, optional
        Manifest CSV with columns: sample_name, order, batch, center, time.
        One sample must have center='TRUE' for RT alignment.
        If None, auto-generates from folder contents. Default: None.
    cores : int, optional
        Number of cores for multiprocessing. Default: 1.
    auto_manifest_batch_threshold_hours : float, optional
        Time gap (hours) for auto-generated batch separation. Default: 12.0.
    auto_manifest_center_name : str, optional
        Sample name for RT alignment center when auto-generating.
        Must match a discovered sample. If None, uses middle sample. Default: None.

    Attributes
    ----------
    folder_location : Path
        Folder containing CoreMS HDF5 files.
    manifest_filepath : Path
        Path to manifest file.
    manifest : dict
        Manifest data indexed by sample_name.
    """
    def __init__(
            self, 
            folder_location: Path, 
            manifest_file: Path = None, 
            cores: int = 1,
            auto_manifest_batch_threshold_hours: float = 12.0,
            auto_manifest_center_name: str = None
            ):
        # Check for folder location
        folder_location = Path(folder_location)
        if not folder_location.exists():
            raise FileNotFoundError(f"Folder location {folder_location} not found.")
        
        # Auto-generate manifest if not provided
        if manifest_file is None:
            print(f"No manifest file provided. Auto-generating manifest from {folder_location}")
            manifest_file = create_manifest_from_folder(
                folder_path=folder_location,
                output_path=folder_location / "manifest_auto.csv",
                batch_time_threshold_hours=auto_manifest_batch_threshold_hours,
                center_name=auto_manifest_center_name,
                overwrite=True
            )
        else:
            manifest_file = Path(manifest_file)
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
        
        # Check that at least one sample has center='TRUE' for retention time alignment
        center_values = [sample_data.get('center') for sample_data in self._manifest_dict.values()]
        if not any(center_val == 'TRUE' or center_val == True for center_val in center_values):
            raise ValueError(
                "Manifest must contain at least one sample with center='TRUE' for retention time alignment. "
                "None of the samples in the manifest have center='TRUE'."
            )

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
        with ReadCoreMSHDFMassSpectra(hdf5_file) as parser:
            lcms_obj = parser.get_lcms_obj(load_raw=load_raw, load_light=load_light, use_original_parser=use_original_parser, raw_file_path=raw_file_path)
            if load_light:
                mf_df = lcms_obj.mass_features_to_df()
                lcms_obj.mass_features = {}
                lcms_obj.light_mf_df = mf_df
        return lcms_obj

    def get_lcms_collection(self, load_raw = False, load_light = True, use_original_parser = True) -> LCMSCollection:
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
            args = [(sample, load_raw, load_light, use_original_parser) for sample in samples]
            lcms_objs = pool.starmap(self.get_lcms_obj, args)
            for sample_name, lcms_obj in zip(samples, lcms_objs):
                lcms_coll._lcms[sample_name] = lcms_obj

        elif self._cores == 1:
            # Load the LCMS objects sequentially
            for sample_name in samples:
                lcms_coll._lcms[sample_name] = self.get_lcms_obj(sample_name, load_raw=load_raw, load_light=load_light, use_original_parser=use_original_parser)

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

        # Load metadata from saved collection
        self._load_collection_metadata()

        if not self.folder_location.exists():
            raise FileNotFoundError(f"Folder location {self.folder_location} not found.")

        # Load the mass spectra data
        self._validate_manifest()
        
        # Set the parameters file location
        self.parameters_location = self._get_parameters_location()

    def _get_parameters_location(self):
        """Find the parameters file (JSON or TOML) associated with the collection HDF5 file."""
        # Check for TOML file first (preferred)
        toml_path = self.collection_hdf5_path.with_suffix('.toml')
        if toml_path.exists():
            return toml_path
        
        # Check for JSON file
        json_path = self.collection_hdf5_path.with_suffix('.json')
        if json_path.exists():
            return json_path
        
        # No parameters file found
        return None

    def _load_collection_metadata(self):
        """Load metadata and manifest from the saved collection HDF5 file."""
        with h5py.File(self.collection_hdf5_path, 'r') as f:
            self.folder_location = Path(f.attrs.get('lcms_objects_folder', ''))
            self.missing_mass_features_searched = f.attrs.get('missing_mass_features_searched', False)

            # Call the _load_manifest function to process the manifest
            self._manifest_dict = self._load_manifest(f)

    def _load_manifest(self, hdf_handle):
        """Load and clean the manifest from the HDF5 file."""
        manifest_json = hdf_handle.attrs.get('manifest', '{}')
        if isinstance(manifest_json, bytes):
            manifest_json = manifest_json.decode('utf-8')
        loaded_manifest = json.loads(manifest_json)

        # Convert integer values for 'use_rt_alignment' back to booleans
        def convert_back_to_bool(data):
            if isinstance(data, dict):
                # Process each key-value pair recursively
                return {k: (bool(v) if k == 'use_rt_alignment' and isinstance(v, int) else convert_back_to_bool(v)) for k, v in data.items()}
            elif isinstance(data, list):
                # Recursively process lists
                return [convert_back_to_bool(item) for item in data]
            else:
                # Return non-dict/list types unchanged
                return data

        # Clean the loaded manifest
        return convert_back_to_bool(loaded_manifest)
    
    def _load_rt_alignments(self, lcms_collection):
        """Load retention time alignments from the saved collection HDF5 file."""
        with h5py.File(self.collection_hdf5_path, 'r') as f:
            if "rt_alignments" in f:
                # Set the lcms_collection 
                lcms_collection.rt_aligned = True
                # Iterate over the group `rt_alignments` containing datasets and add to the corresponding lcms object
                rt_alignments_group = f["rt_alignments"]
                for sample_idx, lcms_obj in zip(rt_alignments_group.keys(), lcms_collection):
                    alignment_data = rt_alignments_group[sample_idx][:]
                    scan_df = lcms_obj.scan_df
                    scan_df["scan_time_aligned"] = alignment_data
                    lcms_obj.scan_df = scan_df

    def _load_cluster_assignments(self, lcms_collection):
        """Load cluster assignments from the saved collection HDF5 file."""
        with h5py.File(self.collection_hdf5_path, 'r') as f:
            if "cluster_assignments" in f:
                # Access the group containing cluster assignments
                cluster_grp = f["cluster_assignments"]
                
                # Reload index and cluster data
                index = cluster_grp["index"][:]  # Extract index
                index = [idx.decode('utf-8') for idx in index]  # Convert byte strings back to regular strings
                cluster_data = cluster_grp["cluster"][:]  # Extract cluster column
                
                # Reassemble the DataFrame
                cluster_df = pd.DataFrame({"cluster": cluster_data}, index=index)

                # Assign cluster data back to lcms_collection.mass_features_dataframe
                lcms_collection.mass_features_dataframe = lcms_collection.mass_features_dataframe.join(cluster_df, how='left')

                # Drop rows with NaN cluster values
                lcms_collection.mass_features_dataframe.dropna(subset=['cluster'], inplace=True)

    def get_lcms_collection(self, load_raw=False, load_light=False):
        """Get the LCMS collection from the saved HDF5 file."""
        # First load the LCMSCollection object exactly as in the parent class
        lcms_collection = super().get_lcms_collection(load_raw=load_raw, load_light=load_light)
        
        # Set the missing_mass_features_searched flag from saved metadata
        lcms_collection.missing_mass_features_searched = self.missing_mass_features_searched

        # Load parameters if a parameters file exists
        if self.parameters_location:
            self._load_parameters(lcms_collection)

        # Add retention time alignments if they exist
        self._load_rt_alignments(lcms_collection)

        # Add cluster assignments if they exist
        self._load_cluster_assignments(lcms_collection)
        
        # Load induced mass features if they exist
        self._load_induced_mass_features(lcms_collection)
        
        # Combine induced mass features into the collection-level dataframe if any were loaded
        if lcms_collection.missing_mass_features_searched:
            lcms_collection._combine_mass_features(induced_features=True)

        return lcms_collection
    
    def _load_parameters(self, lcms_collection):
        """Load collection-level parameters from the saved parameters file."""
        from corems.encapsulation.input.parameter_from_json import (
            load_and_set_json_parameters_lcms_collection,
            load_and_set_toml_parameters_lcms_collection,
        )
        
        if self.parameters_location.suffix == ".json":
            load_and_set_json_parameters_lcms_collection(lcms_collection, self.parameters_location)
        elif self.parameters_location.suffix == ".toml":
            load_and_set_toml_parameters_lcms_collection(lcms_collection, self.parameters_location)
        else:
            warnings.warn(f"Unknown parameter file format: {self.parameters_location.suffix}. Skipping parameter loading.")
    
    def _load_induced_mass_features(self, lcms_collection):
        """Load induced mass features from the saved collection HDF5 file.
        
        Induced mass features are gap-filled features that exist at the collection level.
        This method loads them from the collection HDF5 file with all their attributes
        and datasets, and distributes them to individual LCMS objects.
        
        Parameters
        ----------
        lcms_collection : LCMSCollection
            The LCMS collection object to populate with induced mass features.
        """
        with h5py.File(self.collection_hdf5_path, 'r') as f:
            if "induced_mass_features" not in f:
                return
            
            # Access the top-level induced mass features group
            imf_group = f["induced_mass_features"]
            
            # Iterate through each sample's induced mass features
            for sample_idx in imf_group.keys():
                lcms_obj = lcms_collection[int(sample_idx)]
                sample_group = imf_group[sample_idx]
                
                # Load each mass feature for this sample
                for mf_id_str in sample_group.keys():
                    mf_group = sample_group[mf_id_str]
                    
                    # The mf_id in HDF5 is stored as the collection ID (e.g., 'c10006_422_i' or '0_c10006_422_i')
                    # Extract the integer ID - it's the second-to-last part when split by '_'
                    # Format: sample_id_cCluster_mf_id_i
                    parts = mf_id_str.split('_')
                    # Find the part that's a number (should be second-to-last before 'i')
                    mf_id = int(parts[-2]) if len(parts) > 1 else int(mf_id_str)
                    
                    # Instantiate the LCMSMassFeature object with required attributes
                    mass_feature = LCMSMassFeature(
                        lcms_obj,
                        mz=mf_group.attrs["_mz_exp"],
                        retention_time=mf_group.attrs["_retention_time"],
                        intensity=mf_group.attrs["_intensity"],
                        apex_scan=mf_group.attrs["_apex_scan"],
                        persistence=mf_group.attrs.get("_persistence", 0),
                        id=mf_id,
                    )
                    
                    # Populate additional attributes from HDF5 attributes
                    for key in mf_group.attrs.keys() - {
                        "_mz_exp",
                        "_mz_cal",
                        "_retention_time",
                        "_intensity",
                        "_apex_scan",
                        "_persistence",
                    }:
                        setattr(mass_feature, key, mf_group.attrs[key])
                    
                    # Populate attributes from HDF5 datasets (arrays)
                    for key in mf_group.keys():
                        setattr(mass_feature, key, mf_group[key][:])
                        # Convert _noise_score from array to tuple
                        if key == "_noise_score":
                            mass_feature._noise_score = tuple(mass_feature._noise_score)
                    
                    # Add to the LCMS object's induced_mass_features dictionary
                    lcms_obj.induced_mass_features[mf_id] = mass_feature

