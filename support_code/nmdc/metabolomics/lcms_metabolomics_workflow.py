# TODO KRH: This is a work in progress. It is not yet complete but serves as the development testbed for the
# LCMS metabolomics workflow.
import os
from pathlib import Path
from multiprocessing import Pool
from corems.molecular_id.search.database_interfaces import MSPInterface

from support_code.nmdc.lipidomics.lipidomics_workflow import (
    instantiate_lcms_obj,
    set_params_on_lcms_obj,
    check_scan_translator,
    add_mass_features,
    molecular_formula_search,
    export_results,
    run_lipid_sp_ms1,
)


def run_lcms_metabolomics_workflow(
    file_dir,
    out_dir,
    params_toml,
    msp_file_path,
    scan_translator=None,
    verbose=True,
    cores=1,
):
    # Make output dir and get list of files to process
    out_dir.mkdir(parents=True, exist_ok=True)
    files_list = list(file_dir.glob("*.raw"))
    out_paths_list = [out_dir / f.stem for f in files_list]

    # Prepare search databases for ms2 search
    my_msp_FE = prepare_metadata(msp_file_path)

    # Run signal processing, get associated ms1, add associated ms2, do ms1 molecular search, and export temp results
    # Note that this is exactly the same as the lipidomics workflow
    if cores == 1 or len(files_list) == 1:
        for file_in, file_out in list(zip(files_list, out_paths_list)):
            print(f"Processing {file_in}")
            run_lipid_sp_ms1(
                file_in=str(file_in),
                out_path=str(file_out),
                params_toml=params_toml,
                scan_translator=scan_translator,
                verbose=verbose,
                return_mzs=False,
            )
    elif cores > 1:
        raise ValueError(
            "Parallel processing is not yet supported for LCMS metabolomics workflow."
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
            format="df",
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
    metadata["molecular_metadata"]["positive"] = metabolite_metadata_positive

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


if __name__ == "__main__":
    # Set input variables to run
    msp_file_path = "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/data/databases/20250407_gnps_curated.msp"
    file_dir = Path(
        "/Users/heal742/Library/CloudStorage/OneDrive-PNNL/Documents/_DMS_data/_NMDC/_lcms_metab_test_data"
    )
    params_toml = Path(
        "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/data_processing/configurations/emsl_lcms_metabolomics_corems_params.toml"
    )
    scan_translator = Path(
        "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/data_processing/configurations/emsl_lcms_metabolomics_scan_translator.toml"
    )
    out_dir = Path("tmp_data/__lcms_metab_test_data")
    verbose = True
    cores = 1

    # Set up output directory
    out_dir.mkdir(parents=True, exist_ok=True)

    # if cores > 1, don't use verbose output
    if cores > 2:
        verbose = False

    run_lcms_metabolomics_workflow(
        file_dir=file_dir,
        out_dir=out_dir,
        params_toml=params_toml,
        msp_file_path=msp_file_path,
        scan_translator=scan_translator,
        verbose=verbose,
        cores=cores,
    )
