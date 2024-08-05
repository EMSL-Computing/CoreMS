"""
Notes
--------
Assumes that ms1 are collected in profile mode, persistent homology not applicable for centroided data.
"""

import sys

sys.path.append("./")
import cProfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.mass_spectra.output.export import LipidomicsExport
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.search.database_interfaces import MetabRefLCInterface
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_toml_parameters_lcms,
)

def instantiate_lcms_obj(file_in, verbose):
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
    myLCMSobj = parser.get_lcms_obj(spectra="ms1", verbose=verbose)

    return myLCMSobj


def set_params_on_lcms_obj(myLCMSobj, params_toml):
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


def signal_processing_lcms(myLCMSobj, verbose):
    """Signal processing for LCMS object.

    This includes peak picking, peak grouping, peak integration, annotation of c13 mass features.

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    verbose : bool
        Whether to print verbose output

    Returns
    -------
    None, processes the LCMS object
    """
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.
    myLCMSobj.find_mass_features(verbose=verbose)
    myLCMSobj.integrate_mass_features(drop_if_fail=True)
    myLCMSobj.find_c13_mass_features(verbose=verbose)


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
    i = 1
    # get df of mass features
    mf_df = myLCMSobj.mass_features_to_df()

    # get unique scans to search
    unique_scans = mf_df.apex_scan.unique()

    # search molecular formulas for each mass feature
    for scan in unique_scans:
        if i > 3:  # only search first 3 scans for testing
            break
        print("searching mz for scan: ", str(i), " of ", str(len(unique_scans)))
        # gather mass features for this scan
        mf_df_scan = mf_df[mf_df.apex_scan == scan]
        peaks_to_search = [
            myLCMSobj.mass_features[x].ms1_peak for x in mf_df_scan.index.tolist()
        ]
        SearchMolecularFormulas(
            myLCMSobj._ms[scan],
            first_hit=False,
            find_isotopologues=True,
        ).run_worker_ms_peaks(peaks_to_search)
        i += 1

    print("Finished molecular search")


def process_ms1(myLCMSobj, ms1_molecular_search=True):
    """Process ms1 spectra and perform molecular search

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process
    ms1_molecular_search : bool
        Whether to perform molecular search on ms1 spectra

    Returns
    -------
    None, processes the LCMS object
    """
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )
    myLCMSobj.deconvolute_ms1_mass_features()
    # myLCMSobj.mass_features_to_df()
    if ms1_molecular_search:
        molecular_formula_search(myLCMSobj)


def add_ms2(myLCMSobj):
    """Add ms2 spectra to LCMS object

    Parameters
    ----------
    myLCMSobj : corems LCMS object
        LCMS object to process

    Returns
    -------
    None, processes the LCMS object
    """
    myLCMSobj.add_associated_ms2_dda(spectrum_mode="centroid")


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
    metabref.set_token("tmp_data/thermo_raw_NMDC/metabref.token")

    print("Preparing positive lipid library")
    if metadata["mzs"]["positive"] is not None:
        metabref_positive, lipidmetadata_positive = metabref.get_lipid_library(
            mz_list=metadata["mzs"]["positive"],
            polarity="positive",
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


def process_ms2(myLCMSobj, metadata):
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
    # Grab fe from metatdata associated with polarity
    fe_search = metadata["fe"][myLCMSobj.polarity]

    # separate high and low res ms2 for searching
    hcd_ms2_scan_df = myLCMSobj.scan_df[
        myLCMSobj.scan_df.scan_text.str.contains("hcd")
        & (myLCMSobj.scan_df.ms_level == 2)
    ]
    cid_ms2_scan_df = myLCMSobj.scan_df[
        myLCMSobj.scan_df.scan_text.str.contains("cid")
        & (myLCMSobj.scan_df.ms_level == 2)
    ]

    if len(hcd_ms2_scan_df) > 0:
        # Perform search on high res scans
        ms2_scans_oi_hr = [
            x for x in hcd_ms2_scan_df.scan.tolist() if x in myLCMSobj._ms.keys()
        ]
        myLCMSobj.fe_search(
            scan_list=ms2_scans_oi_hr, fe_lib=fe_search, peak_sep_da=0.01
        )

    if len(cid_ms2_scan_df):
        # Prep low resolution database and search on low res scans
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
        ms2_scans_oi_lr = [
            x for x in cid_ms2_scan_df.scan.tolist() if x in myLCMSobj._ms.keys()
        ]
        myLCMSobj.fe_search(
            scan_list=ms2_scans_oi_lr, fe_lib=fe_search_lr, peak_sep_da=0.3
        )


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
        exporter.report_to_csv(molecular_metadata=molecular_metadata)


def run_lipid_sp_ms1(
    file_in, out_path, params_toml, verbose=True, return_mzs=True
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
    myLCMSobj = instantiate_lcms_obj(file_in, verbose)
    set_params_on_lcms_obj(myLCMSobj, params_toml)
    signal_processing_lcms(myLCMSobj, verbose)
    process_ms1(myLCMSobj, ms1_molecular_search=False) #TODO: change to True when ready
    add_ms2(myLCMSobj)
    export_results(myLCMSobj, out_path=out_path, final=False)
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


def run_lipid_ms2(out_path, metadata):
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
    process_ms2(myLCMSobj, metadata)
    export_results(myLCMSobj, str(out_path), metadata["molecular_metadata"], final=True)


def run_lipid_workflow(file_dir, out_dir, params_toml, verbose, cores=1):
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
                str(file_in), str(file_out), params_toml, verbose
            )
            mz_dicts.append(mz_dict)
    elif cores > 1:
        pool = Pool(cores)
        args = [
            (str(file_in), str(file_out), params_toml, verbose)
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
            mz_dicts = run_lipid_ms2(file_out, metadata)
    elif cores > 1:
        pool = Pool(cores)
        args = [(file_out, metadata) for file_out in out_paths_list]
        mz_dicts = pool.starmap(run_lipid_ms2, args)
        pool.close()
        pool.join()

    print("Finished processing, data are written in " + str(out_dir))


if __name__ == "__main__":
    # Set input variables to run
    cores = 1
    file_dir = Path("tmp_data/thermo_raw_NMDC_mini")
    out_dir = Path("tmp_data/NMDC_processed_0730b")
    params_toml = Path("tmp_data/thermo_raw_NMDC_mini/nmdc_lipid_params.toml")

    verbose = True

    # Set up output directory
    out_dir.mkdir(parents=True, exist_ok=True)

    # if cores > 1, don't use verbose output
    if cores > 1:
        verbose = False

    profile = False
    if profile:
        with cProfile.Profile() as pr:
            run_lipid_workflow(file_dir, out_dir, params_toml, verbose, cores)

        df = pd.DataFrame(
            pr.getstats(),
            columns=["func", "ncalls", "ccalls", "tottime", "cumtime", "callers"],
        )
        df.sort_values("cumtime", ascending=False, inplace=True)
        df.drop("callers", axis=1, inplace=True)
        df["func"] = df["func"].apply(lambda x: str(x))
        df_small = df.head(1000)
        df_small.to_csv(out_dir / "profile.csv", index=False)
        pr.dump_stats(out_dir / "profile.prof")
    else:
        run_lipid_workflow(file_dir, out_dir, params_toml, verbose, cores)
