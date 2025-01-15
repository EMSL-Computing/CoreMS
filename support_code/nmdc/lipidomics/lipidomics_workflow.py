"""
Notes
--------
Assumes that ms1 are collected in profile mode, persistent homology not applicable for centroided data.
"""

import sys

sys.path.append("./")
from multiprocessing import Pool
from pathlib import Path
import datetime
import toml
import warnings

import pandas as pd
import time

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.mass_spectra.output.export import LipidomicsExport
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas, SearchMolecularFormulasLC
from corems.molecular_id.search.database_interfaces import MetabRefLCInterface
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_toml_parameters_lcms,
)


def instantiate_lcms_obj(file_in):
    """Instantiate a corems LCMS object from a binary file.  Pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)

    Parameters
    ----------
    file_in : str or Path
        Path to binary file
    verbose : bool
        Whether to print verbose output

    Returns
    -------
    myLCMSobj : corems LCMS object
        LCMS object with ms1 spectra in dataframe
    """
    # Instantiate parser based on binary file type
    if ".raw" in str(file_in):
        parser = ImportMassSpectraThermoMSFileReader(file_in)

    if ".mzML" in str(file_in):
        parser = MZMLSpectraParser(file_in)

    # Instantiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="ms1")

    return myLCMSobj


def set_params_on_lcms_obj(myLCMSobj, params_toml, verbose):
    """Set parameters on the LCMS object

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to set parameters on
    params_toml : str or Path
        Path to toml file with parameters

    Returns
    -------
    None, sets parameters on the LCMS object
    """
    # Load parameters from toml file
    load_and_set_toml_parameters_lcms(myLCMSobj, params_toml)

    # If myLCMSobj is a positive mode, remove Cl from atoms used in molecular search
    # This cuts down on the number of molecular formulas searched hugely
    if myLCMSobj.polarity == "positive":
        myLCMSobj.parameters.mass_spectrum["ms1"].molecular_search.usedAtoms.pop("Cl")
    elif myLCMSobj.polarity == "negative":
        myLCMSobj.parameters.mass_spectrum["ms1"].molecular_search.usedAtoms.pop("Na")

    if verbose:
        print("Parameters set on LCMS object")


def load_scan_translator(scan_translator=None):
    """Translate scans using a scan translator

    Parameters
    ----------
    scan_translator : str or Path
        Path to scan translator yaml file

    Returns
    -------
    scan_dict : dict
        Dict with keys as parameter keys and values as lists of scans
    """
    # Convert the scan translator to a dictionary
    if scan_translator is None:
        scan_translator_dict = {"ms2": {"scan_filter": "", "resolution": "high"}}
    else:
        # Convert the scan translator to a dictionary
        if isinstance(scan_translator, str):
            scan_translator = Path(scan_translator)
        # read in the scan translator from toml
        with open(scan_translator, "r") as f:
            scan_translator_dict = toml.load(f)
    for param_key in scan_translator_dict.keys():
        if scan_translator_dict[param_key]["scan_filter"] == "":
            scan_translator_dict[param_key]["scan_filter"] = None
    return scan_translator_dict


def check_scan_translator(myLCMSobj, scan_translator):
    """Check if scan translator is provided and that it maps correctly to scans and parameters"""
    scan_translator_dict = load_scan_translator(scan_translator)
    # Check that the scan translator maps correctly to scans and parameters
    scan_df = myLCMSobj.scan_df
    scans_pulled_out = []
    for param_key in scan_translator_dict.keys():
        assert param_key in myLCMSobj.parameters.mass_spectrum.keys()
        assert "scan_filter" in scan_translator_dict[param_key].keys()
        assert "resolution" in scan_translator_dict[param_key].keys()
        # Pull out scans that match the scan filter
        scan_df_sub = scan_df[
            scan_df.scan_text.str.contains(
                scan_translator_dict[param_key]["scan_filter"]
            )
        ]
        scans_pulled_out.extend(scan_df_sub.scan.tolist())
        if len(scan_df_sub) == 0:
            raise ValueError(
                "No scans pulled out by scan translator for parameter key: ",
                param_key,
                " and scan filter: ",
                scan_translator_dict[param_key]["scan_filter"],
            )

    # Check that the scans pulled out by the scan translator are not overlapping and assert error if they are
    if len(set(scans_pulled_out)) != len(scans_pulled_out):
        raise ValueError("Overlapping scans pulled out by scan translator")


def add_mass_features(myLCMSobj, scan_translator):
    """Process ms1 spectra and perform molecular search

    This includes peak picking, adding and processing associated ms1 spectra,
    integration of mass features, annotation of c13 mass features, deconvolution of ms1 mass features,
    and adding of peak shape metrics of mass features to the mass feature dataframe.

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    scan_translator : str or Path
        Path to scan translator yaml file

    Returns
    -------
    None, processes the LCMS object
    """
    myLCMSobj.find_mass_features()
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )
    myLCMSobj.integrate_mass_features(drop_if_fail=True)
    # Count and report how many mass features are left after integration
    print("Number of mass features after integration: ", len(myLCMSobj.mass_features))
    myLCMSobj.find_c13_mass_features()
    myLCMSobj.deconvolute_ms1_mass_features()
    myLCMSobj.add_peak_metrics()

    scan_dictionary = load_scan_translator(scan_translator=scan_translator)
    for param_key in scan_dictionary.keys():
        scan_filter = scan_dictionary[param_key]["scan_filter"]
        if scan_filter == "":
            scan_filter = None
        myLCMSobj.add_associated_ms2_dda(
            spectrum_mode="centroid", ms_params_key=param_key, scan_filter=scan_filter
        )


def molecular_formula_search(myLCMSobj):
    """Perform molecular search on ms1 spectra

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process

    Returns
    -------
    None, processes the LCMS object
    """
    mol_search = SearchMolecularFormulasLC(myLCMSobj)
    mol_search.run_mass_feature_search()
    print("Finished molecular search")


def export_results(myLCMSobj, out_path, molecular_metadata=None, final=False):
    """Export results to hdf5 and csv as a lipid report

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    out_path : str or Path
        Path to output file
    molecular_metadata : dict
        Dict with molecular metadata
    final : bool
        Whether to export final results

    Returns
    -------
    None, exports results to hdf5 and csv as a lipid report
    """
    exporter = LipidomicsExport(out_path, myLCMSobj)
    exporter.to_hdf(overwrite=True)
    if final:
        # Do not show warnings, these are expected
        exporter.report_to_csv(molecular_metadata=molecular_metadata)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exporter.report_to_csv()


def save_times(myLCMSobj, time_start, out_path, time_end=None):
    """Get times for processing steps

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    time_start : float
        Start time of processing
    out_path : str or Path
        Path to output file
    time_end : float
        End time of processing

    Returns
    -------
    None, writes out times to a file within the output directory
    """
    # Check if out_path (with .corems) exisits
    out_dir = Path(str(out_path) + ".corems/")
    if not out_dir.exists():
        print("Output directory does not exist")

    time_toml_path = out_dir / "times.toml"
    if not time_toml_path.exists():
        raw_data_creation_time = myLCMSobj.spectra_parser.get_creation_time().strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        )
        processed_data_creation_time = time_start.strftime("%Y-%m-%dT%H:%M:%SZ")
        time_dict = {
            "raw_data_creation_time": raw_data_creation_time,
            "metabolomics_workflow_start_time": processed_data_creation_time,
        }
        # save as a toml file
        toml_string = toml.dumps(time_dict)  # Output to a string
        with open(time_toml_path, "w") as f:
            f.write(toml_string)
    elif time_toml_path.exists() and time_end is not None:
        time_dict = toml.load(time_toml_path)
        time_dict["metabolomics_workflow_end_time"] = time_end.strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        )
        toml_string = toml.dumps(time_dict)
        with open(time_toml_path, "w") as f:
            f.write(toml_string)


def process_ms2(myLCMSobj, metadata, scan_translator):
    """Process ms2 spectra and perform molecular search

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    metadata : dict
        Dict with keys "mzs", "fe", and "molecular_metadata" with values of dicts of precursor mzs (negative and positive), flash entropy search databases (negative and positive), and molecular metadata, respectively

    Returns
    -------
    None, processes the LCMS object
    """
    # Perform molecular search on ms2 spectra
    # Grab fe from metatdata associated with polarity (this is inherently high resolution as its from a in-silico high res library)
    fe_search = metadata["fe"][myLCMSobj.polarity]

    scan_dictionary = load_scan_translator(scan_translator)
    ms2_scan_df = myLCMSobj.scan_df[myLCMSobj.scan_df.ms_level == 2]

    # Process high resolution MS2 scans
    # Collect all high resolution MS2 scans using the scan translator
    for param_key in scan_dictionary.keys():
        ms2_scans_oi_hr = []
        if scan_dictionary[param_key]["resolution"] == "high":
            scan_filter = scan_dictionary[param_key]["scan_filter"]
            if scan_filter is not None:
                ms2_scan_df_hr = ms2_scan_df[
                    ms2_scan_df.scan_text.str.contains(scan_filter)
                ]
            else:
                ms2_scan_df_hr = ms2_scan_df
            ms2_scans_oi_hr_i = [
                x for x in ms2_scan_df_hr.scan.tolist() if x in myLCMSobj._ms.keys()
            ]
            ms2_scans_oi_hr.extend(ms2_scans_oi_hr_i)
    # Perform search on high res scans
    if len(ms2_scans_oi_hr) > 0:
        myLCMSobj.fe_search(
            scan_list=ms2_scans_oi_hr, fe_lib=fe_search, peak_sep_da=0.01
        )

    # Process low resolution MS2 scans
    # Collect all low resolution MS2 scans using the scan translator
    for param_key in scan_dictionary.keys():
        ms2_scans_oi_lr = []
        if scan_dictionary[param_key]["resolution"] == "low":
            scan_filter = scan_dictionary[param_key]["scan_filter"]
            if scan_filter is not None:
                ms2_scan_df_lr = ms2_scan_df[
                    ms2_scan_df.scan_text.str.contains(scan_filter)
                ]
            else:
                ms2_scan_df_lr = ms2_scan_df
            ms2_scans_oi_lri = [
                x for x in ms2_scan_df_lr.scan.tolist() if x in myLCMSobj._ms.keys()
            ]
            ms2_scans_oi_lr.extend(ms2_scans_oi_lri)
    # Perform search on low res scans
    if len(ms2_scans_oi_lr) > 0:
        # Recast the flashentropy search database to low resolution
        metabref = MetabRefLCInterface()
        fe_search_lr = metabref._to_flashentropy(
            metabref_lib=fe_search,
            normalize=True,
            fe_kwargs={
                "normalize_intensity": True,
                "min_ms2_difference_in_da": 0.4,
                "max_ms2_tolerance_in_da": 0.2,
                "max_indexed_mz": 3000,
                "precursor_ions_removal_da": None,
                "noise_threshold": 0,
            },
        )
        myLCMSobj.fe_search(
            scan_list=ms2_scans_oi_lr, fe_lib=fe_search_lr, peak_sep_da=0.3
        )


def run_lipid_sp_ms1(
    file_in,
    out_path,
    params_toml,
    scan_translator=None,
    verbose=True,
    return_mzs=True,
    ms1_molecular_search=True,
):
    """Run signal processing, get associated ms1, add associated ms2, do ms1 molecular search, and export intermediate results

    Parameters
    ----------
    file_in : str or Path
        Path to binary file
    out_path : str or Path
        Path to output file
    params_toml : str or Path
        Path to toml file with parameters
    verbose : bool
        Whether to print verbose output
    return_mzs : bool
        Whether to return precursor mzs

    Returns
    -------
    mz_dict : dict
        Dict with keys "positive" and "negative" and values of lists of precursor mzs
    """
    time_start = datetime.datetime.now()
    myLCMSobj = instantiate_lcms_obj(file_in)
    set_params_on_lcms_obj(myLCMSobj, params_toml, verbose)
    check_scan_translator(myLCMSobj, scan_translator)
    add_mass_features(myLCMSobj, scan_translator)
    myLCMSobj.remove_unprocessed_data()
    #myLCMSobj.parameters.mass_spectrum['ms1'].molecular_search.verbose_processing = False
    if ms1_molecular_search:
        molecular_formula_search(myLCMSobj)
    export_results(myLCMSobj, out_path=out_path, final=False)
    save_times(myLCMSobj, time_start, out_path)
    if return_mzs:
        precursor_mz_list = list(
            set(
                [
                    v.mz
                    for k, v in myLCMSobj.mass_features.items()
                    if len(v.ms2_scan_numbers) > 0 and v.isotopologue_type is None
                ]
            )
        )
        mz_dict = {myLCMSobj.polarity: precursor_mz_list}
        return mz_dict


def prep_metadata(mz_dicts, out_dir):
    """Prepare metadata for ms2 spectral search

    Parameters
    ----------
    mz_dicts : list of dicts
        List of dicts with keys "positive" and "negative" and values of lists of precursor mzs
    out_dir : Path
        Path to output directory

    Returns
    -------
    metadata : dict
        Dict with keys "mzs", "fe", and "molecular_metadata" with values of dicts of precursor mzs (negative and positive), flash entropy search databases (negative and positive), and molecular metadata, respectively

    Notes
    -------
    Also writes out files for the flash entropy search databases and molecular metadata
    """
    metadata = {
        "mzs": {"positive": None, "negative": None},
        "fe": {"positive": None, "negative": None},
        "molecular_metadata": {},
    }
    for d in mz_dicts:
        metadata["mzs"].update(d)

    metabref = MetabRefLCInterface()

    print("Preparing positive lipid library")
    if metadata["mzs"]["positive"] is not None:
        metabref_positive, lipidmetadata_positive = metabref.get_lipid_library(
            mz_list=metadata["mzs"]["positive"],
            polarity="positive",
            mz_tol_ppm=5,
            format="flashentropy",
            normalize=True,
            fe_kwargs={
                "normalize_intensity": True,
                "min_ms2_difference_in_da": 0.02,  # for cleaning spectra
                "max_ms2_tolerance_in_da": 0.01,  # for setting search space
                "max_indexed_mz": 3000,
                "precursor_ions_removal_da": None,
                "noise_threshold": 0,
            },
        )
        metadata["fe"]["positive"] = metabref_positive
        metadata["molecular_metadata"].update(lipidmetadata_positive)
        fe_positive_df = pd.DataFrame.from_dict(
            {k: v for k, v in enumerate(metadata["fe"]["positive"])}, orient="index"
        )
        fe_positive_df.to_csv(out_dir / "ms2_db_positive.csv")

    print("Preparing negative lipid library")
    if metadata["mzs"]["negative"] is not None:
        metabref_negative, lipidmetadata_negative = metabref.get_lipid_library(
            mz_list=metadata["mzs"]["negative"],
            polarity="negative",
            mz_tol_ppm=5,
            mz_tol_da_api=0.01,
            format="flashentropy",
            normalize=True,
            fe_kwargs={
                "normalize_intensity": True,
                "min_ms2_difference_in_da": 0.02,  # for cleaning spectra
                "max_ms2_tolerance_in_da": 0.01,  # for setting search space
                "max_indexed_mz": 3000,
                "precursor_ions_removal_da": None,
                "noise_threshold": 0,
            },
        )
        metadata["fe"]["negative"] = metabref_negative
        metadata["molecular_metadata"].update(lipidmetadata_negative)
        fe_negative_df = pd.DataFrame.from_dict(
            {k: v for k, v in enumerate(metadata["fe"]["negative"])}, orient="index"
        )
        fe_negative_df.to_csv(out_dir / "ms2_db_negative.csv")

    mol_metadata_df = pd.concat(
        [
            pd.DataFrame.from_dict(v.__dict__, orient="index").transpose()
            for k, v in metadata["molecular_metadata"].items()
        ],
        ignore_index=True,
    )
    mol_metadata_df.to_csv(out_dir / "molecular_metadata.csv")
    return metadata


def run_lipid_ms2(out_path, metadata, scan_translator=None):
    """Run ms2 spectral search and export final results

    Parameters
    ----------
    out_path : str or Path
        Path to output file
    metadata : dict
        Dict with keys "mzs", "fe", and "molecular_metadata" with values of dicts of precursor mzs (negative and positive), flash entropy search databases (negative and positive), and molecular metadata, respectively

    Returns
    -------
    None, runs ms2 spectral search and exports final results
    """
    # Read in the intermediate results
    out_path_hdf5 = str(out_path) + ".corems/" + out_path.stem + ".hdf5"
    parser = ReadCoreMSHDFMassSpectra(out_path_hdf5)
    myLCMSobj = parser.get_lcms_obj()

    # Process ms2 spectra, perform spectral search, and export final results
    process_ms2(myLCMSobj, metadata, scan_translator=scan_translator)
    export_results(myLCMSobj, str(out_path), metadata["molecular_metadata"], final=True)
    time_end = datetime.datetime.now()
    save_times(myLCMSobj, time_start=None, out_path=out_path, time_end=time_end)


def run_lipid_workflow(
    file_dir,
    out_dir,
    params_toml,
    scan_translator=None,
    verbose=True,
    ms1_molecular_search=True,
    cores=1,
):
    """Run lipidomics workflow

    Parameters
    ----------
    file_dir : str or Path
        Path to directory with raw or mzml lipid files
    out_dir : str or Path
        Path to output directory
    params_toml : str or Path
        Path to toml file with parameters
    verbose : bool
        Whether to print verbose output
    cores : int
        Number of cores to use

    Returns
    -------
    None, runs lipidomics workflow and exports final results
    """
    # Make output dir and get list of files to process
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = list(file_dir.glob("*.raw"))
    out_paths_list = [out_dir / f.stem for f in files_list]

    # Run signal processing, get associated ms1, add associated ms2, do ms1 molecular search, and export temp results
    if cores == 1 or len(files_list) == 1:
        mz_dicts = []
        for file_in, file_out in list(zip(files_list, out_paths_list)):
            mz_dict = run_lipid_sp_ms1(
                file_in=str(file_in),
                out_path=str(file_out),
                params_toml=params_toml,
                scan_translator=scan_translator,
                verbose=verbose,
                ms1_molecular_search=ms1_molecular_search,
            )
            mz_dicts.append(mz_dict)
    elif cores > 1:
        pool = Pool(cores)
        args = [
            (
                str(file_in),
                str(file_out),
                params_toml,
                scan_translator,
                verbose,
                ms1_molecular_search,
            )
            for file_in, file_out in list(zip(files_list, out_paths_list))
        ]
        mz_dicts = pool.starmap(run_lipid_sp_ms1, args)
        pool.close()
        pool.join()
    # Prepare ms2 spectral search space
    metadata = prep_metadata(mz_dicts, out_dir)

    # Run ms2 spectral search and export final results
    if cores == 1 or len(files_list) == 1:
        for file_out in out_paths_list:
            mz_dicts = run_lipid_ms2(
                file_out, metadata, scan_translator=scan_translator
            )
    elif cores > 1:
        pool = Pool(cores)
        args = [(file_out, metadata, scan_translator) for file_out in out_paths_list]
        mz_dicts = pool.starmap(run_lipid_ms2, args)
        pool.close()
        pool.join()
    print("Finished processing, data are written in " + str(out_dir))


if __name__ == "__main__":
    # Set input variables to run
    cores = 1
    file_dir = Path("/Users/heal742/LOCAL/corems_dev/corems/tmp_data/thermo_raw_mini")
    out_dir = Path("tmp_data/_test_250115")
    params_toml = Path("/Users/heal742/LOCAL/05_NMDC/02_MetaMS/data_processing/configurations/emsl_lipidomics_corems_params.toml")
    verbose = True
    scan_translator = Path("tmp_data/thermo_raw_collection/scan_translator.toml")

    # Set up output directory
    out_dir.mkdir(parents=True, exist_ok=True)

    # if cores > 1, don't use verbose output
    if cores > 2:
        verbose = False

    run_lipid_workflow(
        file_dir=file_dir,
        out_dir=out_dir,
        params_toml=params_toml,
        scan_translator=scan_translator,
        verbose=verbose,
        cores=cores,
    )
