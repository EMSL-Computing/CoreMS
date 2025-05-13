# TODO KRH: This is a work in progress. It is not yet complete but serves as the development testbed for the
# LCMS metabolomics workflow.
import os
from pathlib import Path
from multiprocessing import Pool
import warnings

from corems.molecular_id.search.database_interfaces import MSPInterface
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.output.export import LCMSMetabolomicsExport

from support_code.nmdc.lipidomics.lipidomics_workflow import (
    instantiate_lcms_obj,
    set_params_on_lcms_obj,
    check_scan_translator,
    add_mass_features,
    molecular_formula_search,
    process_ms2,
)


def prepare_metadata(msp_file_path):
    print("Parsing MSP file...")
    my_msp = MSPInterface(file_path=msp_file_path)
    print("Parsing MSP file complete.")
    metadata = {
        "fe": {"positive": None, "negative": None},
        "molecular_metadata": {},
    }
    msp_positive, metabolite_metadata_positive = (
        my_msp.get_metabolomics_spectra_library(
            polarity="positive",
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
    )
    metadata["fe"]["positive"] = msp_positive
    metadata["molecular_metadata"] = metabolite_metadata_positive

    msp_negative, metabolite_metadata_negative = (
        my_msp.get_metabolomics_spectra_library(
            polarity="negative",
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
    )
    metadata["fe"]["negative"] = msp_negative
    metadata["molecular_metadata"].update(metabolite_metadata_negative)

    return metadata


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
    exporter = LCMSMetabolomicsExport(out_path, myLCMSobj)
    exporter.to_hdf(overwrite=True)
    if final:
        # Do not show warnings, these are expected
        exporter.report_to_csv(molecular_metadata=molecular_metadata)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exporter.report_to_csv()


def run_ms2_search(out_path, metadata, scan_translator=None):
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
    out_path = Path(out_path)
    out_path_hdf5 = str(out_path) + ".corems/" + out_path.stem + ".hdf5"
    parser = ReadCoreMSHDFMassSpectra(out_path_hdf5)
    myLCMSobj = parser.get_lcms_obj()
    process_ms2(myLCMSobj, metadata, scan_translator=scan_translator)
    export_results(
        myLCMSobj,
        out_path=str(out_path),
        molecular_metadata=metadata["molecular_metadata"],
        final=True,
    )


def run_lcms_metabolomics_workflow(
    file_in,
    out_path,
    params_toml,
    scan_translator,
    metadata,
):
    myLCMSobj = instantiate_lcms_obj(file_in)
    set_params_on_lcms_obj(myLCMSobj, params_toml, verbose)
    
    ms1_scan_df = myLCMSobj.scan_df[myLCMSobj.scan_df.ms_level == 1]
    if all(x == "centroid" for x in ms1_scan_df.ms_format.to_list()):
        # Switch peak picking method to centroided persistent homology
        myLCMSobj.parameters.lc_ms.peak_picking_method = "centroided_persistent_homology"
        myLCMSobj.parameters.mass_spectrum[
            "ms1"
        ].mass_spectrum.noise_threshold_method = "relative_abundance"

    myLCMSobj.parameters.mass_spectrum[
        "ms1"
    ].mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    check_scan_translator(myLCMSobj, scan_translator)
    add_mass_features(myLCMSobj, scan_translator)
    myLCMSobj.remove_unprocessed_data()
    molecular_formula_search(myLCMSobj)
    # Just for testing, not needed for workflow
    export_results(
        myLCMSobj,
        out_path=str(out_path),
        molecular_metadata=metadata["molecular_metadata"],
        final=False,
    )
    process_ms2(myLCMSobj, metadata, scan_translator=scan_translator)
    export_results(
        myLCMSobj,
        out_path=str(out_path),
        molecular_metadata=metadata["molecular_metadata"],
        final=True,
    )


def run_lcms_metabolomics_workflow_batch(
    file_dir,
    out_dir,
    params_toml,
    msp_file_path,
    scan_translator=None,
    cores=1,
):
    """Run the LCMS metabolomics workflow on a batch of files

    Parameters
    ----------
    file_dir : str or Path
        Path to directory with raw files
    out_dir : str or Path
        Path to output directory
    params_toml : str or Path
        Path to toml file with parameters for all samples
    msp_file_path : str or Path
        Path to msp file with spectral library for MS2 search
    scan_translator : str or Path
        Path to toml file with scan translator for MS2 search
    cores : int
        Number of cores to use for processing
    """
    # Make output dir and get list of files to process
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = [
        f for f in file_dir.iterdir() if f.suffix.lower() in {".raw", ".mzml"}
    ]
    out_paths_list = [out_dir / f.stem for f in files_list]

    # Prepare search databases for ms2 search
    my_msp_FE = prepare_metadata(msp_file_path)

    # Run signal processing, get associated ms1, add associated ms2, do ms1 molecular search, and export temp results
    # Note that this is exactly the same as the lipidomics workflow
    if cores == 1 or len(files_list) == 1:
        for file_in, file_out in list(zip(files_list, out_paths_list)):
            print(f"Processing {file_in}")
            run_lcms_metabolomics_workflow(
                file_in=file_in,
                out_path=file_out,
                params_toml=params_toml,
                scan_translator=scan_translator,
                metadata=my_msp_FE,
            )
    elif cores > 1:
        raise ValueError(
            "Parallel processing is not yet supported for LCMS metabolomics workflow."
        )


if __name__ == "__main__":
    # Set input variables to run
    msp_file_path = "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/data/databases/20250407_gnps_curated.msp"
    # msp_file_path = "/Users/heal742/LOCAL/corems_dev/corems/tests/tests_data/lcms/test_db.msp"
    file_dir = Path(
        "/Users/heal742/Library/CloudStorage/OneDrive-PNNL/Documents/_DMS_data/_NMDC/_lcms_metab_test_data/centroid"
    )
    params_toml = Path(
        "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/data_processing/configurations/emsl_lcms_metabolomics_corems_params.toml"
    )
    scan_translator = Path(
        "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/data_processing/configurations/emsl_lcms_metabolomics_scan_translator.toml"
    )
    out_dir = Path("tmp_data/__lcms_metab_test_data_250509")
    verbose = True
    cores = 1

    # Set up output directory
    out_dir.mkdir(parents=True, exist_ok=True)

    run_lcms_metabolomics_workflow_batch(
        file_dir=file_dir,
        out_dir=out_dir,
        params_toml=params_toml,
        msp_file_path=msp_file_path,
        scan_translator=scan_translator,
        cores=cores,
    )
