from pathlib import Path
import time
import pandas as pd
import numpy as np
from multiprocessing import Pool
import shutil

from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection, ReadSavedLCMSCollection
from corems.mass_spectra.output.export import LCMSCollectionExport, LCMSMetabolomicsExport
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.encapsulation.factory.parameters import LCMSParameters
from corems.molecular_id.search.database_interfaces import MSPInterface

"""
Example showing the new pipeline-based sample processing approach.

The new approach combines multiple sample-level operations (gap-filling, 
feature reloading, MS1/MS2 searches, etc.) into a single parallelized pass,
which is more efficient than processing samples multiple times.

Two usage patterns:
1. High-level convenience method: process_consensus_features()
2. Advanced pipeline builder: process_samples_pipeline() with custom operations
"""
def process_single_sample(args):
    """
    Process a single LCMS sample file.
    
    Parameters
    ----------
    args : tuple
        (raw_file_path, processed_folder)
    
    Returns
    -------
    str
        Path to the processed HDF5 file
    """
    raw_file_path, processed_folder = args
    
    # Import the raw data
    print(f"Processing {raw_file_path.name}...")
    parser = ImportMassSpectraThermoMSFileReader(str(raw_file_path))
    lcms_obj = parser.get_lcms_obj(spectra="ms1")
    
    # Set parameters to the defaults for reproducible testing
    lcms_obj.parameters = LCMSParameters(use_defaults=True)

    # Set parameters on the LCMS object that are reasonable for testing
    ## persistent homology parameters
    lcms_obj.parameters.lc_ms.peak_picking_method = "persistent homology"
    lcms_obj.parameters.lc_ms.ph_inten_min_rel = 0.0005
    lcms_obj.parameters.lc_ms.ph_persis_min_rel = 0.01
    lcms_obj.parameters.lc_ms.ph_smooth_it = 0
    lcms_obj.parameters.lc_ms.ms2_min_fe_score = 0.3
    lcms_obj.parameters.lc_ms.ms1_scans_to_average = 1

    ## MSParameters for ms1 mass spectra
    ms1_params = lcms_obj.parameters.mass_spectrum['ms1']
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz, ms1_params.mass_spectrum.min_picking_mz = 0, 0
    ms1_params.mass_spectrum.noise_max_mz, ms1_params.mass_spectrum.max_picking_mz = np.inf, np.inf
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

    ## settings for ms2 data (HCD scans)
    ms2_params_hcd = ms1_params.copy()
    lcms_obj.parameters.mass_spectrum['ms2'] = ms2_params_hcd

    ## reporting settings
    lcms_obj.parameters.lc_ms.export_eics = True
    lcms_obj.parameters.lc_ms.export_profile_spectra = True

    ## peak metrics filtering settings
    lcms_obj.parameters.lc_ms.remove_mass_features_by_peak_metrics = True
    lcms_obj.parameters.lc_ms.mass_feature_attribute_filter_dict = {
        'dispersity_index': {'value': 0.5, 'operator': '<'}
    }

    # Use persistent homology to find mass features in the lc-ms data
    lcms_obj.find_mass_features()
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
    ncores = 1
    reprocess_samples = False # Set to True to reprocess raw data, False to use existing processed data

    # Set paths
    raw_data_path = Path("/Volumes/LaCie/nmdc_data/collection_testing/dev_test/raw")
    processed_folder = Path("/Volumes/LaCie/nmdc_data/collection_testing/dev_test/processed")
    msp_file_location = Path("/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/test_data/test_lcms_metab_data/20250407_database.msp")
    
    if reprocess_samples:
        # Delete existing processed dir if reprocessing
        if processed_folder.exists():
            shutil.rmtree(processed_folder)

        # Create processed folder if it doesn't exist
        processed_folder.mkdir(parents=True, exist_ok=True)
        
        # Find all raw files (adjust extension based on your data format)
        raw_files = list(raw_data_path.glob("*.raw")) + list(raw_data_path.glob("*.mzML"))
        
        if not raw_files:
            raise ValueError(f"No raw files found in {raw_data_path}")
        
        print(f"\n=== Preprocessing {len(raw_files)} samples in parallel using {ncores} cores ===")
        start_time = time.time()
        
        # Prepare arguments for parallel processing
        process_args = [(raw_file, processed_folder) for raw_file in raw_files]
        
        # Process samples in parallel
        with Pool(processes=ncores) as pool:
            processed_files = pool.map(process_single_sample, process_args)
        
        print(f"Time to preprocess all samples: {time.time() - start_time:.1f} seconds using {ncores} cores")
        print(f"Processed {len(processed_files)} samples to {processed_folder}\n")
   
    # Set the path to the collection of LCMS runs (previously processed)
    collection_path = processed_folder
    
    # Instantiate the parser (manifest will be auto-generated if it doesn't exist)
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            cores = ncores
            )
    print("\n=== Loading LCMS collection with", len(parser.manifest), "samples using", ncores, "cores ===")

    # Load the LCMS collection (minimally load the data)
    start_time = time.time()
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    print("Time to load LCMS collection ", time.time() - start_time, "seconds -", len(lcms_collection), " LCMS runs and ", ncores, " cores")

    # Update raw file locations to point to the raw data folder
    lcms_collection.update_raw_file_locations(
        new_raw_folder = str(raw_data_path)
        )   
    print("Number of total mass features: ", len(lcms_collection.mass_features_dataframe))

    # Align the LCMS runs between each other
    # For now, adjusting this parameter to force alignment for testing
    lcms_collection.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold = -1
    lcms_collection.parameters.lcms_collection.alignment_acceptance_technique = ['fraction_improved']
    print("Aligning LCMS collection")
    start_time = time.time()
    assert not lcms_collection.rt_aligned, "LCMS collection should not be marked as retention time aligned yet."
    assert lcms_collection.rt_alignments is None, "LCMS collection should not have rt_alignments yet."
    lcms_collection.align_lcms_objects()
    assert lcms_collection.rt_aligned, "LCMS collection should be marked as retention time aligned."
    assert lcms_collection.rt_alignments is not None, "LCMS collection should have rt_alignments now."
    print("Time to align LCMS collection: ", time.time() - start_time, "seconds") 

    # Make consensus mass features from the consolidated mass features
    print("Generating consensus mass features across the LCMS collection")
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    print("Time to generate consensus mass features: ", time.time() - start_time, "seconds -", len(lcms_collection.mass_features_dataframe), " total mass features", ncores, " cores")   

    # Tell the user how many clusters were generated
    print(f"Total clusters formed: {len(lcms_collection.cluster_summary_dataframe)}")

    # Prepare spectral library for MS2 search (mimicking test_lcms_metabolomics)
    print("\n=== Preparing spectral library for MS2 search ===")
    
    # Check if MSP file exists before attempting to load
    if msp_file_location.exists():
        my_msp = MSPInterface(file_path=msp_file_location)
        spectral_lib, molecular_metadata = my_msp.get_metabolomics_spectra_library(
            polarity="negative",  # Change to match your data polarity
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
        print(f"Loaded spectral library with {len(molecular_metadata)} entries")
        enable_spectral_search = True
    else:
        raise FileNotFoundError(f"MSP file not found at {msp_file_location}. Cannot perform spectral search.")
    
    # PIPELINE APPROACH: Gap fill, reload, add MS1/MS2, run molecular formula and spectral search in a single pass
    print("\n=== Testing new pipeline approach with MS1, MS2, molecular formula, and spectral search ===")
    start_time = time.time()
    pipeline_results = lcms_collection.process_consensus_features(
        perform_gap_filling=True,
        reload_representatives=True,
        add_ms1=True,  
        add_ms2=True,
        molecular_formula_search=True,
        ms2_spectral_search=True,
        spectral_lib=spectral_lib,
        molecular_metadata=molecular_metadata,
        keep_raw_data=False
    )
    print("Time for combined reload, MS1, MS2, MF search, and spectral search: ", time.time() - start_time, "seconds, using", ncores, " cores")
    print("Gap-filled features in", len([s for s in pipeline_results.get('gap_fill', {}).values() if s]), "samples")
    print("Reloaded features in", len([s for s in pipeline_results.get('reload', {}).values() if s]), "samples")
    print("Molecular formula search completed on", len([s for s in pipeline_results.get('mf_search', {}).values() if s]), "samples")
    if enable_spectral_search:
        print("MS2 spectral search completed on", len([s for s in pipeline_results.get('ms2_search', {}).values() if s and s > 0]), "samples")
        total_ms2_searched = sum([s for s in pipeline_results.get('ms2_search', {}).values() if s])
        print(f"Total MS2 spectra searched: {total_ms2_searched}")
    
    # Verify that mass features were reloaded
    total_mf_reloaded = sum([len(lcms_obj.mass_features) for lcms_obj in lcms_collection])
    print(f"Total mass features reloaded: {total_mf_reloaded}")
    assert total_mf_reloaded > 0, "Should have reloaded some mass features"
    
    # Check for MS1 associations
    total_ms1_with_spectra = 0
    total_mf_checked = 0
    for lcms_obj in lcms_collection:
        for mf_id, mf in lcms_obj.mass_features.items():
            total_mf_checked += 1
            if hasattr(mf, 'mass_spectrum') and mf.mass_spectrum is not None:
                total_ms1_with_spectra += 1
    
    print(f"Total mass features with MS1 spectra: {total_ms1_with_spectra} out of {total_mf_checked}")
    if total_ms1_with_spectra > 0:
        print(f"✓ MS1 spectra successfully associated with {total_ms1_with_spectra/total_mf_checked*100:.1f}% of mass features")
        assert total_ms1_with_spectra > 0, "Should have MS1 spectra associated with mass features"
    else:
        print("⚠ No MS1 spectra associated")
    
    # Check for MS2 associations
    total_ms2 = 0
    for lcms_obj in lcms_collection:
        for mf_id, mf in lcms_obj.mass_features.items():
            if hasattr(mf, 'ms2_mass_spectra') and mf.ms2_mass_spectra:
                total_ms2 += len(mf.ms2_mass_spectra)
    print(f"Total MS2 spectra associated: {total_ms2}")
    if total_ms2 > 0:
        print("✓ MS2 spectra successfully associated with mass features")
    else:
        print("⚠ No MS2 spectra associated (this may be expected if no MS2 data exists)")
    
    # Check for molecular formula assignments
    total_mf_with_formulas = 0
    total_formula_assignments = 0
    for lcms_obj in lcms_collection:
        for mf_id, mf in lcms_obj.mass_features.items():
            if hasattr(mf, 'mass_spectrum') and mf.mass_spectrum is not None:
                # Check if the ms1_peak has molecular formula assignments
                try:
                    ms1_peak = mf.ms1_peak
                    if hasattr(ms1_peak, 'molecular_formulas') and ms1_peak.molecular_formulas:
                        total_mf_with_formulas += 1
                        total_formula_assignments += len(ms1_peak.molecular_formulas)
                except (AttributeError, IndexError):
                    # Skip if ms1_peak can't be determined
                    pass
    
    print(f"Total mass features with molecular formula assignments: {total_mf_with_formulas} out of {total_mf_checked}")
    print(f"Total molecular formula assignments: {total_formula_assignments}")
    if total_mf_with_formulas > 0:
        print(f"✓ Molecular formula search successfully assigned formulas to {total_mf_with_formulas/total_mf_checked*100:.1f}% of mass features")
        print(f"  Average {total_formula_assignments/total_mf_with_formulas:.1f} formulas per assigned feature")
    else:
        print("⚠ No molecular formula assignments (check search parameters)")
    
    # Verify raw data was cleaned up (unless keep_raw_data=True)
    raw_data_present = any(1 in lcms_obj._ms_unprocessed and not lcms_obj._ms_unprocessed[1].empty 
                          for lcms_obj in lcms_collection)
    if not raw_data_present:
        print("✓ Raw MS1 data successfully cleaned up after pipeline")
    else:
        print("⚠ Raw MS1 data still present (expected if keep_raw_data=True)")
    
    # Check for spectral match results (if spectral search was performed)
    if enable_spectral_search:
        total_spectral_matches = 0
        total_mf_with_matches = 0
        for lcms_obj in lcms_collection:
            if hasattr(lcms_obj, 'spectral_search_results') and lcms_obj.spectral_search_results:
                total_spectral_matches += len(lcms_obj.spectral_search_results)
                total_mf_with_matches += 1
                break  # Count each mass feature only once
        
        print(f"\nSpectral Search Results:")
        print(f"Total mass features with spectral matches: {total_mf_with_matches}")
        print(f"Total spectral matches: {total_spectral_matches}")
        if total_mf_with_matches > 0:
            print(f"✓ MS2 spectral search successfully found matches")
            print(f"  Average {total_spectral_matches/total_mf_with_matches:.1f} matches per feature")
            # Access the molecular metadata if needed for export
            if hasattr(lcms_collection, 'spectral_search_molecular_metadata'):
                print(f"  Molecular metadata available with {len(lcms_collection.spectral_search_molecular_metadata)} entries")
        else:
            print("⚠ No spectral matches found (check library and search parameters)")

    """
    # OLD APPROACH (commented out - replaced by pipeline above):
    # Gap fill missing cluster features BEFORE saving
    start_time = time.time()
    lcms_collection.fill_missing_cluster_features()
    print("Time to gap fill missing cluster features: ", time.time() - start_time, "seconds, using", ncores, " cores")

    # Reload representative mass features with MS2 data associated
    sample_mf_map = lcms_collection.reload_representative_mass_features(
            add_ms2=True,
            auto_process_ms2=True,
            ms2_spectrum_mode=None,
            ms2_scan_filter=None
        )
    """
    
    """
    # ADVANCED PIPELINE APPROACH (for custom workflows):
    # Build a custom pipeline with full control over operations
    from corems.mass_spectra.calc.lc_calc_operations import (
        GapFillOperation, 
        ReloadFeaturesOperation,
        CustomOperation
    )
    
    # Define custom operation
    def my_custom_processing(sample_id, collection, **params):
        sample = collection[sample_id]
        # Do custom processing here
        # e.g., normalization, quality checks, etc.
        return None
    
    # Build pipeline
    ops = [
        GapFillOperation('gap_fill', expand_on_miss=True),
        ReloadFeaturesOperation('reload', add_ms2=True, auto_process_ms2=True),
        CustomOperation('custom', func=my_custom_processing)
    ]
    
    # Execute
    results = lcms_collection.process_samples_pipeline(ops, description="Custom workflow")
    """

    """
    # Check save and load functionality for LCMSCollection
    print("Saving and re-loading LCMS collection to test save/load functionality")
    print(f"Before saving: missing_mass_features_searched = {lcms_collection.missing_mass_features_searched}")
    exporter = LCMSCollectionExport(
        out_file_path="/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/collection",
        mass_spectra_collection=lcms_collection)
    exporter.export_to_hdf5(overwrite=True, save_parameters=True, parameter_format="toml")
    
    # Reload the collection
    reader = ReadSavedLCMSCollection(
        collection_hdf5_path="/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/collection.hdf5",
        cores=ncores)
    lcms_collection2 = reader.get_lcms_collection(load_raw=False, load_light=True)
    
    # Verify parameters, induced features, and missing mass features searched flag upon reload
    assert lcms_collection2.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold == -1, \
        f"Re-loaded threshold parameter should be -1, got {lcms_collection2.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold}"
    assert lcms_collection2.missing_mass_features_searched, "Re-loaded collection should have missing_mass_features_searched=True"
    original_induced_count = sum([len(lcms_obj.induced_mass_features) for lcms_obj in lcms_collection])
    reloaded_induced_count = sum([len(lcms_obj.induced_mass_features) for lcms_obj in lcms_collection2])
    assert reloaded_induced_count == original_induced_count, \
        f"Re-loaded induced mass features count mismatch: {reloaded_induced_count} vs {original_induced_count}"
    assert reloaded_induced_count > 0, "Collection should have some induced mass features after gap filling"

    # Verify the pivot table outputs are the same before and after save/load
    original_pivot = lcms_collection.collection_pivot_table()
    reloaded_pivot = lcms_collection2.collection_pivot_table()
    pd.testing.assert_frame_equal(original_pivot, reloaded_pivot, check_dtype=False)

    print('Test completed successfully! LCMSCollection save and load functionality works as expected.')
    
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
    
    #TODO KRH: Add visualization of a consensus mass feature
    #TODO KRH: Add visualization of matched spectrum with consensus mass feature
    #TODO KRH: Add code to deal with annotations of features in the collection context
    """