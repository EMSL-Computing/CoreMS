#!/usr/bin/env python3
"""Run the LC-MS collection workflow on a folder of Thermo raw files.

Same steps as ``examples/notebooks/LCMS_Collection_Tutorial.ipynb``, on the
files in ``DATA_FOLDER``, with more than one core:

1. Peak-pick each raw file and write a CoreMS HDF5 file.
2. Load the files as one collection.
3. Align retention times.
4. Build consensus features.
5. Gap-fill features that a sample is missing.
6. Reload the representative feature for each consensus feature.

Formula search and MS2 library search are left off, so no database is
required.

Set ``DATA_FOLDER`` below, then run this file. In an editor, use Run or
Debug with no extra arguments. From the repository root:

    python examples/spawn_pool_collection.py

On Linux, set ``PYTHONNET_RUNTIME=coreclr`` in the environment before
starting Python, the same way the CoreMS container does. Leave that
variable unset on macOS.

The parent process opens one raw file before any process pool. Each step
prints when it starts and when it finishes. A step that never prints
"finished" is the hang this script is here to catch.
"""

import os
import shutil
import sys
import time
from pathlib import Path

# This file lives in examples/. Running it puts that directory on sys.path,
# and an editable install can point at a different checkout. Import the
# CoreMS next to this file. Spawn workers reimport this module, so the
# insert has to happen at import time.
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# =============================================================================
# Configuration. Edit these, then run this file. There are no arguments.
# =============================================================================

# Folder of Thermo .raw files. At least two files. Leave None until you set it.
DATA_FOLDER = None

# Where the .corems folders and manifest.csv are written.
# None uses DATA_FOLDER / "corems_collection_out".
OUTPUT_FOLDER = None

# Must be 2 or more. One core never starts a process pool.
CORES = 2

# True deletes OUTPUT_FOLDER and peak-picks every raw file again.
REPROCESS = True


def _stop(message):
    print(message, flush=True)
    raise SystemExit(1)


def _checked_folders():
    """Return the raw-data folder and the output folder."""
    if not isinstance(CORES, int) or CORES < 2:
        _stop("Set CORES to an integer of 2 or more. One core does not start a pool.")
    if not isinstance(REPROCESS, bool):
        _stop("Set REPROCESS to True or False.")
    if DATA_FOLDER is None or str(DATA_FOLDER).strip() == "":
        _stop(
            "Set DATA_FOLDER at the top of examples/spawn_pool_collection.py "
            "to a folder that contains at least two Thermo .raw files."
        )

    data_folder = Path(DATA_FOLDER).expanduser()
    if not data_folder.is_dir():
        _stop(f"DATA_FOLDER is not a folder: {data_folder}")

    if OUTPUT_FOLDER is None or str(OUTPUT_FOLDER).strip() == "":
        output_folder = data_folder / "corems_collection_out"
    else:
        output_folder = Path(OUTPUT_FOLDER).expanduser()

    data_resolved = data_folder.resolve()
    output_resolved = output_folder.resolve()
    if output_resolved == data_resolved or data_resolved.is_relative_to(output_resolved):
        _stop(
            "OUTPUT_FOLDER would include the raw files. "
            "Pick a folder that is not the raw-data folder or a parent of it."
        )
    return data_folder, output_folder


def _raw_files(data_folder):
    """Return Thermo raw files in ``data_folder``, one path per file name."""
    found = list(data_folder.glob("*.raw")) + list(data_folder.glob("*.RAW"))
    unique = {path.resolve(): path for path in found}
    raw_files = sorted(unique.values(), key=lambda path: path.name.lower())
    if len(raw_files) < 2:
        _stop(
            f"Need at least two Thermo .raw files in {data_folder}. "
            f"Found {len(raw_files)}. Files in subfolders are not included."
        )
    stems = [path.stem for path in raw_files]
    if len(stems) != len(set(stems)):
        _stop("Two raw files have the same file name. Each file needs a unique name.")
    available = os.cpu_count() or 1
    if CORES > available:
        _stop(f"CORES is {CORES}, and this machine reports {available} cores.")
    return raw_files


def _hdf_path(output_folder, sample_name):
    return output_folder / f"{sample_name}.corems" / f"{sample_name}.hdf5"


def _configure_peak_picking(lcms_obj):
    """Use the fast persistent-homology settings from the collection notebook."""
    import numpy as np
    from corems.encapsulation.factory.parameters import LCMSParameters, hush_output

    hush_output()
    lcms_obj.parameters = LCMSParameters(use_defaults=True)
    lcms_obj.parameters.lc_ms.peak_picking_method = "persistent homology"
    lcms_obj.parameters.lc_ms.ph_inten_min_rel = 0.001
    lcms_obj.parameters.lc_ms.ph_persis_min_rel = 0.05
    lcms_obj.parameters.lc_ms.ph_smooth_it = 0
    lcms_obj.parameters.lc_ms.ms1_scans_to_average = 3

    ms1_params = lcms_obj.parameters.mass_spectrum["ms1"]
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz = 0
    ms1_params.mass_spectrum.min_picking_mz = 0
    ms1_params.mass_spectrum.noise_max_mz = np.inf
    ms1_params.mass_spectrum.max_picking_mz = np.inf

    ms2_params = lcms_obj.parameters.mass_spectrum["ms2"]
    ms2_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms2_params.mass_spectrum.noise_threshold_min_relative_abundance = 1.0


def process_one_raw(raw_file, output_folder):
    """Peak-pick one Thermo raw file and write its CoreMS HDF5 file.

    A spawn worker reimports this module and calls this function. Keep it at
    module level, and pass paths rather than an open reader.
    """
    from corems.mass_spectra.input.rawFileReader import (
        ImportMassSpectraThermoMSFileReader,
    )
    from corems.mass_spectra.output.export import LCMSMetabolomicsExport

    raw_file = Path(raw_file)
    output_folder = Path(output_folder)
    print(f"Peak picking {raw_file.name}", flush=True)
    parser = ImportMassSpectraThermoMSFileReader(str(raw_file))
    try:
        lcms_obj = parser.get_lcms_obj(spectra="ms1")
        _configure_peak_picking(lcms_obj)
        lcms_obj.find_mass_features(assign_ms2_scans=True)
        lcms_obj.integrate_mass_features(drop_if_fail=True)
        n_features = len(lcms_obj.mass_features)
        sample_name = raw_file.stem
        exporter = LCMSMetabolomicsExport(str(output_folder / sample_name), lcms_obj)
        exporter.to_hdf(overwrite=True)
    finally:
        parser.close_file()
    print(f"Finished {raw_file.name}: {n_features} mass features", flush=True)
    return sample_name, n_features


def _open_raw_in_parent(raw_file):
    """Start the .NET runtime in this process, then close the file."""
    from corems.mass_spectra.input.rawFileReader import (
        ImportMassSpectraThermoMSFileReader,
    )

    print(f"Opening {raw_file.name} in the parent process, before any pool", flush=True)
    parser = ImportMassSpectraThermoMSFileReader(str(raw_file))
    parser.close_file()
    print("Parent closed that raw file", flush=True)


def _peak_pick(raw_files, output_folder):
    from corems.encapsulation.spawn_pool import spawn_pool

    if output_folder.exists():
        shutil.rmtree(output_folder)
    output_folder.mkdir(parents=True)
    ncores = min(CORES, len(raw_files))
    args = [(str(path), str(output_folder)) for path in raw_files]
    print(
        f"Peak picking {len(raw_files)} files with {ncores} cores",
        flush=True,
    )
    with spawn_pool(ncores) as pool:
        results = pool.starmap(process_one_raw, args)
    for sample_name, _n_features in results:
        hdf_path = _hdf_path(output_folder, sample_name)
        if not hdf_path.is_file():
            _stop(f"Peak picking did not write {hdf_path}")
    return results


def _write_manifest(output_folder, sample_names):
    import csv

    manifest_path = output_folder / "manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_name", "batch", "order", "center"])
        for order, sample_name in enumerate(sample_names, start=1):
            writer.writerow([sample_name, 1, order, order == 1])
    print(f"Wrote {manifest_path}", flush=True)
    print(f"Center sample for alignment: {sample_names[0]}", flush=True)
    return manifest_path


def _step(name):
    print(f"\n=== {name} ===", flush=True)
    return time.perf_counter()


def _finished(name, started):
    elapsed = time.perf_counter() - started
    print(f"=== {name} finished in {elapsed:.1f} s ===", flush=True)


def main():
    data_folder, output_folder = _checked_folders()
    raw_files = _raw_files(data_folder)
    sample_names = [path.stem for path in raw_files]
    print(f"Raw files ({len(raw_files)}):", flush=True)
    for path in raw_files:
        print(f"  {path.name}", flush=True)
    print(f"Output folder: {output_folder}", flush=True)
    print(f"Cores: {CORES}", flush=True)

    # Importing the Thermo reader starts pythonnet. Do that, and open one
    # file, before the first pool. That is the order that deadlocks on Linux
    # when a pool uses fork.
    _open_raw_in_parent(raw_files[0])

    started = _step("Peak picking")
    if REPROCESS:
        _peak_pick(raw_files, output_folder)
    else:
        missing = [
            name
            for name in sample_names
            if not _hdf_path(output_folder, name).is_file()
        ]
        if missing:
            _stop(
                "REPROCESS is False, and these samples have no HDF5 file yet: "
                + ", ".join(missing)
            )
        print("Using the HDF5 files already in the output folder", flush=True)
    _finished("Peak picking", started)

    manifest_path = _write_manifest(output_folder, sample_names)

    started = _step("Load collection")
    from corems.encapsulation.factory.parameters import hush_output
    from corems.mass_spectra.input.corems_hdf5 import (
        ReadCoreMSHDFMassSpectraCollection,
    )

    hush_output()
    reader = ReadCoreMSHDFMassSpectraCollection(
        folder_location=output_folder,
        manifest_file=manifest_path,
        cores=CORES,
    )
    lcms_collection = reader.get_lcms_collection()
    n_features = len(lcms_collection.mass_features_dataframe)
    print(
        f"Loaded {len(lcms_collection)} samples, {n_features} mass features",
        flush=True,
    )
    _finished("Load collection", started)

    # Alignment settings from the collection notebook.
    collection_settings = lcms_collection.parameters.lcms_collection
    collection_settings.mass_feature_anchor_technique = ["relative_intensity"]
    collection_settings.mass_feature_anchor_relative_intensity_threshold = 0.0
    collection_settings.alignment_mz_tol_ppm = 5
    collection_settings.alignment_rt_tol = 0.2

    started = _step("Align retention times")
    lcms_collection.align_lcms_objects()
    print(f"Alignment attempted: {lcms_collection.rt_alignment_attempted}", flush=True)
    _finished("Align retention times", started)

    started = _step("Consensus features")
    lcms_collection.add_consensus_mass_features()
    features = lcms_collection.mass_features_dataframe
    n_clusters = features["cluster"].nunique()
    print(
        f"{len(features)} mass features in {n_clusters} consensus features",
        flush=True,
    )
    _finished("Consensus features", started)

    started = _step("Gap fill")
    gap_fill = lcms_collection.process_consensus_features(
        load_representatives=False,
        perform_gap_filling=True,
        add_ms1=False,
        add_ms2=False,
        molecular_formula_search=False,
        ms2_spectral_search=False,
        spectral_lib=None,
        molecular_metadata=None,
        gather_eics=False,
        keep_raw_data=False,
        show_progress=True,
    )
    induced = lcms_collection.induced_mass_features_dataframe
    n_induced = 0 if induced is None else len(induced)
    print(f"Gap-filled features: {n_induced}", flush=True)
    print(f"Gap-fill result keys: {sorted(gap_fill)}", flush=True)
    _finished("Gap fill", started)

    started = _step("Reload representative features")
    reloaded = lcms_collection.process_consensus_features(
        load_representatives=True,
        perform_gap_filling=False,
        add_ms1=False,
        add_ms2=False,
        molecular_formula_search=False,
        ms2_spectral_search=False,
        spectral_lib=None,
        molecular_metadata=None,
        gather_eics=False,
        keep_raw_data=False,
        show_progress=True,
    )
    print(f"Reload result keys: {sorted(reloaded)}", flush=True)
    _finished("Reload representative features", started)
    print("\nCollection run finished.", flush=True)


if __name__ == "__main__":
    main()
