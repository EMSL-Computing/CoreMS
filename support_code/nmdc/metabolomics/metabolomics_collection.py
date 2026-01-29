from pathlib import Path
import time
import pandas as pd
import numpy as np
from multiprocessing import Pool
import shutil

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection
from corems.mass_spectra.output.export import LCMSMetabolomicsExport
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.encapsulation.factory.parameters import LCMSParameters
from corems.molecular_id.search.database_interfaces import MSPInterface
from corems.encapsulation.factory.parameters import hush_output
from corems.mass_spectra.output.export import LCMSCollectionExport

"""
Example showing the new pipeline-based sample processing approach.

The new approach combines multiple sample-level operations (gap-filling, 
feature reloading, MS1/MS2 searches, etc.) into a single parallelized pass,
which is more efficient than processing samples multiple times.

Two usage patterns:
1. High-level convenience method: process_consensus_features()
2. Advanced pipeline builder: process_samples_pipeline() with custom operations
"""

def summarize_processing_results(lcms_collection):
    """
    Summarize the processing state of the LCMS collection.
    
    Reports on completed processing steps by inspecting the collection
    and sample objects directly. Useful for verifying which operations
    were performed during process_consensus_features().
    
    Parameters
    ----------
    lcms_collection : LCMSCollection
        The LCMS collection to summarize
    """
    print("\n" + "="*60)
    print("LCMS Collection Processing Summary")
    print("="*60)
    
    # Basic collection info
    n_samples = len(lcms_collection)
    
    # Collection-level feature count (from dataframe - all features)
    collection_mf_count = len(lcms_collection.mass_features_dataframe) if lcms_collection.mass_features_dataframe is not None else 0
    
    # Sample-level loaded feature count (from mass_features dict - loaded for processing)
    loaded_mf_count = sum(len(lcms_obj.mass_features) for lcms_obj in lcms_collection)
    
    # Cluster count
    if 'cluster' in lcms_collection.mass_features_dataframe.columns:
        n_clusters = lcms_collection.mass_features_dataframe['cluster'].nunique()
    else:
        n_clusters = 0
    
    print(f"\nSamples: {n_samples}")
    print(f"Total features in collection: {collection_mf_count}")
    print(f"Representative features loaded: {loaded_mf_count}")
    if n_clusters > 0:
        print(f"Consensus clusters: {n_clusters}")
    
    # Gap filling - check for induced mass features
    induced_counts = [len(lcms_obj.induced_mass_features) for lcms_obj in lcms_collection]
    total_induced = sum(induced_counts)
    samples_with_induced = sum(1 for c in induced_counts if c > 0)
    if total_induced > 0:
        print(f"\nGap Filling: ✓ Complete")
        print(f"  {samples_with_induced}/{n_samples} samples have induced features ({total_induced} total)")
    
    # Feature loading - check if mass features have MS1/MS2 spectra
    mf_with_ms1 = 0
    mf_with_ms2 = 0
    mf_with_ms2_scans = 0
    total_ms2_spectra = 0
    total_ms2_scan_numbers = 0
    
    for lcms_obj in lcms_collection:
        for mf in lcms_obj.mass_features.values():
            if hasattr(mf, 'mass_spectrum') and mf.mass_spectrum is not None:
                mf_with_ms1 += 1
            if hasattr(mf, 'ms2_mass_spectra') and mf.ms2_mass_spectra:
                mf_with_ms2 += 1
                total_ms2_spectra += len(mf.ms2_mass_spectra)
            if hasattr(mf, 'ms2_scan_numbers') and mf.ms2_scan_numbers is not None and len(mf.ms2_scan_numbers) > 0:
                mf_with_ms2_scans += 1
                total_ms2_scan_numbers += len(mf.ms2_scan_numbers)
    
    if mf_with_ms1 > 0 or mf_with_ms2 > 0:
        print(f"\nMS Data Association: ✓ Complete")
        if mf_with_ms1 > 0:
            print(f"  MS1: {mf_with_ms1}/{loaded_mf_count} loaded features ({mf_with_ms1/loaded_mf_count*100:.1f}%)")
        if mf_with_ms2 > 0:
            print(f"  MS2: {mf_with_ms2}/{loaded_mf_count} loaded features ({total_ms2_spectra} spectra)")
        if mf_with_ms2_scans > 0:
            print(f"  MS2 scan numbers: {mf_with_ms2_scans}/{loaded_mf_count} loaded features ({total_ms2_scan_numbers} scans)")
    
    # Molecular formula search
    mf_with_formulas = 0
    total_formulas = 0
    
    for lcms_obj in lcms_collection:
        for mf in lcms_obj.mass_features.values():
            if hasattr(mf, 'mass_spectrum') and mf.mass_spectrum is not None:
                try:
                    ms1_peak = mf.ms1_peak
                    if hasattr(ms1_peak, 'molecular_formulas') and ms1_peak.molecular_formulas:
                        mf_with_formulas += 1
                        total_formulas += len(ms1_peak.molecular_formulas)
                except (AttributeError, IndexError):
                    pass
    
    if mf_with_formulas > 0:
        print(f"\nMolecular Formula Search: ✓ Complete")
        print(f"  {mf_with_formulas}/{loaded_mf_count} loaded features assigned ({total_formulas} total formulas)")
        print(f"  Average {total_formulas/mf_with_formulas:.1f} formulas per feature")
    
    # MS2 spectral search
    mf_with_spectral_matches = 0
    total_spectral_matches = 0
    scans_searched = 0
    
    for lcms_obj in lcms_collection:
        if hasattr(lcms_obj, 'spectral_search_results') and lcms_obj.spectral_search_results:
            scans_searched += len(lcms_obj.spectral_search_results)
        
        for mf in lcms_obj.mass_features.values():
            if hasattr(mf, 'ms2_similarity_results') and mf.ms2_similarity_results:
                mf_with_spectral_matches += 1
                total_spectral_matches += len(mf.ms2_similarity_results)
    
    if mf_with_spectral_matches > 0:
        print(f"\nMS2 Spectral Search: ✓ Complete")
        print(f"  {scans_searched} MS2 scans with library search results")
        print(f"  {mf_with_spectral_matches}/{loaded_mf_count} loaded features matched ({total_spectral_matches} total matches)")
        if hasattr(lcms_collection, 'spectral_search_molecular_metadata'):
            print(f"  Library size: {len(lcms_collection.spectral_search_molecular_metadata)} entries")
    
    # Check for loaded EICs
    total_eics_loaded = 0
    samples_with_eics = 0
    
    for lcms_obj in lcms_collection:
        if hasattr(lcms_obj, 'eics') and lcms_obj.eics:
            samples_with_eics += 1
            total_eics_loaded += len(lcms_obj.eics)
    
    if total_eics_loaded > 0:
        print(f"\nEIC Loading: ✓ Complete")
        print(f"  {total_eics_loaded} EICs loaded across {samples_with_eics}/{n_samples} samples")
        print(f"  Average {total_eics_loaded/samples_with_eics:.1f} EICs per sample")
    
    # Memory management check
    raw_data_present = any(1 in lcms_obj._ms_unprocessed and not lcms_obj._ms_unprocessed[1].empty 
                          for lcms_obj in lcms_collection)
    if not raw_data_present:
        print(f"\nMemory: ✓ Raw MS1 data cleaned")
    
    print("\n" + "="*60)


def validate_save_load(lcms_collection, collection_save_path, ncores=1):
    """
    Validate that the LCMS collection can be saved and reloaded correctly.
    
    This function tests the save/load functionality by reloading the collection
    from HDF5 and comparing various attributes to ensure data integrity.
    
    Parameters
    ----------
    lcms_collection : LCMSCollection
        The original LCMS collection to validate against
    collection_save_path : Path
        Path to the saved collection HDF5 file (without extension)
    ncores : int, optional
        Number of cores to use for loading. Default is 1.
    
    Returns
    -------
    LCMSCollection
        The reloaded LCMS collection for further testing
    """
    print("\n" + "="*60)
    print("SAVE/LOAD VALIDATION TEST")
    print("="*60)
    
    # Reload the collection from HDF5
    from corems.mass_spectra.input.corems_hdf5 import ReadSavedLCMSCollection
    reader = ReadSavedLCMSCollection(
        collection_hdf5_path=str(collection_save_path.with_suffix('.hdf5')),
        cores=ncores)
    lcms_collection2 = reader.get_lcms_collection(
        load_raw=False, load_light=True,
        load_representatives=True, load_eics=True,
        load_ms1=True, load_ms2=True)

    # Check that len of mass features matches
    mf_count_1 = len(lcms_collection.mass_features_dataframe)
    mf_count_2 = len(lcms_collection2.mass_features_dataframe)
    if mf_count_1 == mf_count_2:
        print(f"✓ Mass feature count matches after reload: {mf_count_1}")
    else:
        print(f"✗ Mass feature count mismatch after reload! Original: {mf_count_1}, Reloaded: {mf_count_2}")
    
    # Check that the len of induced mass features matches
    induced_count_1 = len(lcms_collection.induced_mass_features_dataframe)
    induced_count_2 = len(lcms_collection2.induced_mass_features_dataframe)
    if induced_count_1 == induced_count_2:
        print(f"✓ Induced mass feature count matches after reload: {induced_count_1}")
    else:
        print(f"✗ Induced mass feature count mismatch after reload! Original: {induced_count_1}, Reloaded: {induced_count_2}")
    
    # Check that the len of loaded EICs matches
    total_eics_1 = sum(len(lcms_obj.eics) for lcms_obj in lcms_collection)
    total_eics_2 = sum(len(lcms_obj.eics) for lcms_obj in lcms_collection2)
    if total_eics_1 == total_eics_2:
        print(f"✓ Total loaded EIC count matches after reload: {total_eics_1}")
    else:
        print(f"✗ Total loaded EIC count mismatch after reload! Original: {total_eics_1}, Reloaded: {total_eics_2}")
    
    # Check that the _ms dictionary matches (mass spectra loaded)
    total_ms_1 = sum(len(lcms_obj._ms) for lcms_obj in lcms_collection)
    total_ms_2 = sum(len(lcms_obj._ms) for lcms_obj in lcms_collection2)
    if total_ms_1 == total_ms_2:
        print(f"✓ Total loaded mass spectra (_ms) count matches after reload: {total_ms_1}")
    else:
        print(f"✗ Total loaded mass spectra (_ms) count mismatch after reload! Original: {total_ms_1}, Reloaded: {total_ms_2}")
    
    # Check that we can replot the first cluster from the reloaded collection
    if len(lcms_collection2.cluster_summary_dataframe) > 0:
        first_cluster_id = lcms_collection2.cluster_summary_dataframe.index[0]
        print(f"Re-plotting cluster {first_cluster_id} from reloaded collection")
        lcms_collection2.plot_cluster(
            cluster_id=first_cluster_id,
            to_plot=["EIC", "MS1", "MS2"],
            plot_smoothed_eic=False,
            plot_eic_datapoints=False
        )

    # Check that scan numbers in _ms match for first sample
    if len(lcms_collection) > 0 and len(lcms_collection2) > 0:
        ms_scans_1 = set(lcms_collection[0]._ms.keys())
        ms_scans_2 = set(lcms_collection2[0]._ms.keys())
        if ms_scans_1 == ms_scans_2:
            print(f"✓ Mass spectra scan numbers match for first sample: {len(ms_scans_1)} scans")
        else:
            print(f"✗ Mass spectra scan numbers mismatch for first sample!")
            print(f"  Original: {len(ms_scans_1)} scans, Reloaded: {len(ms_scans_2)} scans")
            only_in_1 = ms_scans_1 - ms_scans_2
            only_in_2 = ms_scans_2 - ms_scans_1
            if only_in_1:
                print(f"  Only in original: {list(sorted(only_in_1))[:10]}...")
            if only_in_2:
                print(f"  Only in reloaded: {list(sorted(only_in_2))[:10]}...")
    
    print("\n" + "="*60)
    
    return lcms_collection2


def preprocess_raw_samples(raw_data_path, processed_folder, ncores=1, reprocess=False):
    """
    Preprocess raw LCMS sample files into HDF5 format.
    
    Parameters
    ----------
    raw_data_path : Path
        Path to folder containing raw data files
    processed_folder : Path
        Path to folder where processed HDF5 files will be saved
    ncores : int, optional
        Number of cores to use for parallel processing. Default is 1.
    reprocess : bool, optional
        If True, deletes existing processed folder and reprocesses all files.
        If False, skips preprocessing. Default is False.
    
    Returns
    -------
    list or None
        List of processed HDF5 file paths if reprocess=True, None otherwise
    """
    if not reprocess:
        print("\n=== Skipping sample preprocessing (using existing processed data) ===")
        return None
    
    # Delete existing processed dir if reprocessing
    if processed_folder.exists():
        shutil.rmtree(processed_folder)
    
    # Create processed folder
    processed_folder.mkdir(parents=True, exist_ok=True)
    
    # Find all raw files (adjust extension based on your data format)
    raw_files = list(raw_data_path.glob("*.raw")) + list(raw_data_path.glob("*.mzML"))
    
    if not raw_files:
        raise ValueError(f"No raw files found in {raw_data_path}")
    
    print(f"\n=== Preprocessing {len(raw_files)} samples in parallel using {ncores} cores ===")
    start_time = time.time()
    
    # Get configured parameters once (will be shared across all workers)
    params = get_configured_lcms_parameters()
    
    # Prepare arguments for parallel processing
    process_args = [(raw_file, processed_folder, params) for raw_file in raw_files]
    
    # Process samples in parallel
    with Pool(processes=ncores) as pool:
        processed_files = pool.map(process_single_sample, process_args)
    
    print(f"Preprocessing complete: {time.time() - start_time:.1f} seconds using {ncores} cores")
    print(f"Processed {len(processed_files)} samples\n")
    
    return processed_files


def get_configured_lcms_parameters():
    """
    Create and configure LCMSParameters for sample processing.
    
    Returns
    -------
    LCMSParameters
        Configured parameters object with all processing settings
    """
    # Suppress verbose output before creating parameters
    hush_output()
    
    # Create parameters (use_defaults=False respects hush_output)
    params = LCMSParameters()
    
    # Persistent homology parameters
    params.lc_ms.peak_picking_method = "persistent homology"
    params.lc_ms.ph_inten_min_rel = 0.0005
    params.lc_ms.ph_persis_min_rel = 0.01
    params.lc_ms.ph_smooth_it = 0
    params.lc_ms.ms2_min_fe_score = 0.3
    params.lc_ms.ms1_scans_to_average = 1
    
    # MSParameters for ms1 mass spectra
    ms1_params = params.mass_spectrum['ms1']
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz = 0
    ms1_params.mass_spectrum.min_picking_mz = 0
    ms1_params.mass_spectrum.noise_max_mz = np.inf
    ms1_params.mass_spectrum.max_picking_mz = np.inf
    ms1_params.ms_peak.legacy_resolving_power = False
    ms1_params.molecular_search.url_database = ""
    ms1_params.molecular_search.usedAtoms = {
        'C': (5, 30),
        'H': (18, 200),
        'O': (1, 23),
        'N': (0, 3),
        'P': (0, 1),
        'S': (0, 1),
    }
    
    # Settings for ms2 data (HCD scans)
    ms2_params_hcd = ms1_params.copy()
    params.mass_spectrum['ms2'] = ms2_params_hcd
    
    # Reporting settings
    params.lc_ms.export_eics = True
    params.lc_ms.export_profile_spectra = False
    
    # Peak metrics filtering settings
    params.lc_ms.remove_mass_features_by_peak_metrics = True
    params.lc_ms.mass_feature_attribute_filter_dict = {
        'dispersity_index': {'value': 0.5, 'operator': '<'}
    }
    
    return params


def process_single_sample(args):
    """
    Process a single LCMS sample file.
    
    Parameters, params)
    
    Returns
    -------
    str
        Path to the processed HDF5 file
    """
    raw_file_path, processed_folder, params = args
    
    # Import the raw data
    print(f"Processing {raw_file_path.name}...\n")
    parser = ImportMassSpectraThermoMSFileReader(str(raw_file_path))
    lcms_obj = parser.get_lcms_obj(spectra="ms1")
    
    # Use the pre-configured parameters
    lcms_obj.parameters = params
    
    # Get configured parameters
    lcms_obj.parameters = get_configured_lcms_parameters()

    # Use persistent homology to find mass features in the lc-ms data and integrate
    # Assign MS2 scan numbers during peak picking for choosing representatives with MS2
    lcms_obj.find_mass_features(assign_ms2_scans=True, ms2_scan_filter=None)
    lcms_obj.integrate_mass_features(drop_if_fail=True)
    
    # Add peak metrics and filter mass features based on the new parameters
    lcms_obj.add_peak_metrics(remove_by_metrics=True)
    
    # Save to HDF5
    output_name = raw_file_path.stem
    exporter = LCMSMetabolomicsExport(str(processed_folder / output_name), lcms_obj)
    exporter.to_hdf(overwrite=True)
    
    print(f"✓ Completed {raw_file_path.name} - {len(lcms_obj.mass_features)} mass features")
    return str(processed_folder / f"{output_name}.hdf5")

if __name__ == "__main__":
    # =============================================================================
    # Configuration
    # =============================================================================
    ncores = 1
    reprocess_samples = False  # Set to True to reprocess raw data
    perform_ms2_search = True  # Set to True to perform MS2 spectral library search

    # Paths
    base_path = Path("/Volumes/LaCie/nmdc_data/collection_testing/dev_test/")
    collection_save_path = base_path / "collection"
    raw_data_path = base_path / "raw"
    processed_folder = base_path / "processed2"
    msp_file_location = Path("/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/test_data/test_lcms_metab_data/20250407_database.msp")
    new_raw_data_path = raw_data_path  # Update raw file paths in collection if moved after the first step
    
    # =============================================================================
    # Step 1: Preprocess Individual Samples (Optional)
    # =============================================================================
    if reprocess_samples:
        print("\n=== Preprocessing Raw Samples and Doing Initial Peak Picking===")
        preprocess_raw_samples(
            raw_data_path=raw_data_path,
            processed_folder=processed_folder,
            ncores=ncores,
            reprocess=reprocess_samples
        )
   
    # =============================================================================
    # Step 2: Load LCMS Collection
    # =============================================================================
    print("\n=== Loading LCMS Collection ===")
    parser = ReadCoreMSHDFMassSpectraCollection(
        folder_location=processed_folder,
        cores=ncores
    )
    print(f"Found {len(parser.manifest)} samples")
    
    # Load collection (light loading for efficiency)
    start_time = time.time()
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    print(f"Loaded in {time.time() - start_time:.1f} seconds")
    print(f"Total mass features: {len(lcms_collection.mass_features_dataframe)}")
    
    # Update raw file locations
    lcms_collection.update_raw_file_locations(new_raw_folder=str(new_raw_data_path))
    
    # =============================================================================
    # Step 3: Align Retention Times Across Samples
    # =============================================================================
    print("\n=== Aligning Retention Times ===")
    start_time = time.time()
    lcms_collection.align_lcms_objects()
    print(f"Alignment complete: {time.time() - start_time:.1f} seconds")
    
    # =============================================================================
    # Step 4: Generate Consensus Mass Features
    # =============================================================================
    print("\n=== Generating Consensus Mass Features ===")
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    print(f"Generated {len(lcms_collection.cluster_summary_dataframe)} consensus clusters")
    print(f"Consensus generation: {time.time() - start_time:.1f} seconds")
    
    # =============================================================================
    # Step 5: Prepare MS2 Spectral Library
    # =============================================================================
    if perform_ms2_search:
        print("\n=== Preparing MS2 Spectral Library ===")
        my_msp = MSPInterface(file_path=msp_file_location)
        spectral_lib, molecular_metadata = my_msp.get_metabolomics_spectra_library(
            polarity="positive",
            format="flashentropy",
            normalize=True,
            fe_kwargs={
                "normalize_intensity": True,
                "min_ms2_difference_in_da": 0.02,
                "max_ms2_tolerance_in_da": 0.01,
                "max_indexed_mz": 3000,
                "precursor_ions_removal_da": None,
                "noise_threshold": 0,
            },
        )
        print(f"Loaded spectral library: {len(molecular_metadata)} entries")
    else: 
        spectral_lib = None
        molecular_metadata = None
        print("Skipping MS2 spectral library preparation")

    # =============================================================================
    # Step 6: Process Consensus Features with Integrated Pipeline
    # =============================================================================
    print("\n=== Processing Consensus Features ===")
    start_time = time.time()
    pipeline_results = lcms_collection.process_consensus_features(
        load_representatives=True,
        perform_gap_filling=True,
        add_ms1=True,  
        add_ms2=True,
        molecular_formula_search=False,
        ms2_spectral_search=perform_ms2_search,
        spectral_lib=spectral_lib,
        molecular_metadata=molecular_metadata,
        gather_eics=True,
        keep_raw_data=False
    )
    print(f"Pipeline complete: {time.time() - start_time:.1f} seconds using {ncores} cores")
    
    # =============================================================================
    # Step 6.5: Test Consensus Cluster Plotting
    # =============================================================================
    print("\n=== Testing Consensus Cluster Plotting ===")
    
    # Plot the first cluster as a test
    if len(lcms_collection.cluster_summary_dataframe) > 0:
        first_cluster_id = lcms_collection.cluster_summary_dataframe.index[0]
        # 3116 is a good one to look at :)
        print(f"Plotting cluster {first_cluster_id}")
        lcms_collection.plot_cluster(
            cluster_id=first_cluster_id,
            to_plot=["EIC", "MS1", "MS2"],
            plot_smoothed_eic=False,
            plot_eic_datapoints=False
        )
    
    # =============================================================================
    # Step 7: Summarize Processing Results
    # =============================================================================
    summarize_processing_results(lcms_collection)
    

    # =============================================================================
    # Step 8: Save and Export Results
    # =============================================================================
    print("\n=== Exporting LCMS Collection ===")
    # Create pivot tables summarizing the collection across samples
    pivot_table_intensity = lcms_collection.collection_pivot_table(attribute='intensity', verbose=False)
    pivot_table_ids = lcms_collection.collection_pivot_table(verbose=False)
    pivot_table_intensity.to_csv("example_collection_pivot_intensity.csv")

    # Describe each cluster with its representative mass feature
    cluster_reps = lcms_collection.cluster_representatives_table()
    cluster_reps.to_csv("example_cluster_representatives.csv", index=False)
    print(f"Cluster representatives table: {len(cluster_reps)} clusters")
    
    # Summarize the annotations for each cluster
    feature_annotations = lcms_collection.feature_annotations_table(
        molecular_metadata=molecular_metadata,
        drop_unannotated=True
    )
    print(f"Feature annotations table: {len(feature_annotations)} rows across {feature_annotations['cluster'].nunique()} clusters")

    # Save the feature annotations table to CSV for inspection
    feature_annotations.to_csv("example_feature_annotations.csv", index=False)

    # Save the entire collection to HDF5
    exporter = LCMSCollectionExport(
        out_file_path=str(collection_save_path),
        mass_spectra_collection=lcms_collection)
    exporter.export_to_hdf5(overwrite=True, save_parameters=True, parameter_format="toml")
    
    # =============================================================================
    # Step 9: Validate Save/Load Functionality
    # =============================================================================
    lcms_collection2 = validate_save_load(lcms_collection, collection_save_path, ncores=ncores)


    """
    # Make some more plots
    lcms_collection.plot_mz_features_across_samples()
    lcms_collection.plot_mz_features_per_cluster()
    lcms_collection.plot_consensus_mz_features() ## zoomed out
    lcms_collection.plot_consensus_mz_features(xb = 10, xt = 15, yb = 500, yt = 600) ## zoomed in 
    lcms_collection.cluster_inspection_plot(0)
       
    dim_list = [
        'mz', 
        'scan_time_aligned', 
        'half_height_width', 
        'tailing_factor', 
        'dispersity_index', 
        'intensity', 
        'persistence'
    ]
    lcms_collection.plot_cluster_outlier_frequency(dim_list, clu_size_thresh = 0.25)
    
    # Create pivot tables and reports
    results = lcms_collection.collection_pivot_table()
    results1 = lcms_collection.collection_pivot_table(attribute = 'intensity')
    results2 = lcms_collection.collection_consensus_report(how = 'intensity')
    results3 = lcms_collection.collection_consensus_report(how = 'median')
    
    #TODO KRH: Add visualization of matched spectrum with consensus mass feature
    """