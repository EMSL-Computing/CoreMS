from pathlib import Path
import time
import pandas as pd
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection, ReadSavedLCMSCollection
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

if __name__ == "__main__":
   
    # Set the path to the collection of LCMS runs (previously processed)
    collection_path = Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2")
    # Path to manifest file
    manifest_file = Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2/manifest_tiny.csv")

    # Set the number of cores to use for loading the data (the parser is parallelized)
    ncores = 6

    # Instantiate the parser
    # Note, these samples have not had the DDA MS2 scans associated or MS1 data loaded 
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            manifest_file = manifest_file,
            cores = ncores
            )
    print("Loading LCMS collection with", len(parser.manifest), "samples using", ncores, "cores")

    # Load the LCMS collection (minimally load the data)
    start_time = time.time()
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    print("Time to load LCMS collection ", time.time() - start_time, "seconds -", len(lcms_collection), " LCMS runs and ", ncores, " cores")

    # Update raw file locations (optionally, but common)
    og_file_location = lcms_collection[0].raw_file_location
    lcms_collection.update_raw_file_locations(
        new_raw_folder = "/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test2"
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

    """# Make some plots 
    lcms_collection.plot_tics(type="both")
    lcms_collection.plot_alignments()    
    """

    # Make consensus mass features from the consolidated mass features
    print("Generating consensus mass features across the LCMS collection")
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    print("Time to generate consensus mass features: ", time.time() - start_time, "seconds -", len(lcms_collection.mass_features_dataframe), " total mass features", ncores, " cores")   

    # NEW PIPELINE APPROACH: Gap fill, reload, and add MS1/MS2 in a single pass
    print("\n=== Testing new pipeline approach with MS1 and MS2 ===")
    start_time = time.time()
    pipeline_results = lcms_collection.process_consensus_features(
        perform_gap_filling=True,
        reload_representatives=True,
        add_ms1=True,
        add_ms2=True,
        auto_process_ms2=True,
        ms2_scan_filter=None,
        keep_raw_data=False
    )
    print("Time for combined gap-fill, reload, MS1, and MS2: ", time.time() - start_time, "seconds, using", ncores, " cores")
    print("Gap-filled features in", len([s for s in pipeline_results.get('gap_fill', {}).values() if s]), "samples")
    print("Reloaded features in", len([s for s in pipeline_results.get('reload', {}).values() if s]), "samples")
    
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
    
    # Verify raw data was cleaned up (unless keep_raw_data=True)
    raw_data_present = any(1 in lcms_obj._ms_unprocessed and not lcms_obj._ms_unprocessed[1].empty 
                          for lcms_obj in lcms_collection)
    if not raw_data_present:
        print("✓ Raw MS1 data successfully cleaned up after pipeline")
    else:
        print("⚠ Raw MS1 data still present (expected if keep_raw_data=True)")

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