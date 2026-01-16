from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection, ReadSavedLCMSCollection
from corems.mass_spectra.output.export import LCMSCollectionExport

if __name__ == "__main__":
   
    # Set the path to the collection of LCMS runs (previously processed)
    collection_path = Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2")
    # Path to manifest file
    manifest_file = Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2/manifest.csv")

    # Set the number of cores to use for loading the data (the parser is parallelized)
    ncores = 6

    # Instantiate the parser
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

    # Check save and load functionality for LCMSCollection
    print("Saving and re-loading LCMS collection to test save/load functionality")
    exporter = LCMSCollectionExport(
        out_file_path="/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/collection",
        mass_spectra_collection=lcms_collection)
    exporter.export_to_hdf5(save_parameters=True, parameter_format="toml")
    
    # Reload the collection
    reader = ReadSavedLCMSCollection(
        collection_hdf5_path="/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/collection.hdf5",
        cores=ncores)
    lcms_collection2 = reader.get_lcms_collection(load_raw=False, load_light=True)
    
    # Verify basic collection structure
    assert len(lcms_collection2) == len(lcms_collection), "Re-loaded LCMS collection should have the same number of samples."
    assert len(lcms_collection2.mass_features_dataframe) == len(lcms_collection.mass_features_dataframe), "Re-loaded LCMS collection should have the same number of mass features."
    assert str(lcms_collection2[0].raw_file_location) == str(lcms_collection[0].raw_file_location), "Re-loaded LCMS collection should have the same raw file locations."
    assert lcms_collection2.rt_aligned == lcms_collection.rt_aligned and lcms_collection.rt_aligned, "Both LCMS collections should be marked as retention time aligned."
    assert "cluster" in lcms_collection2.mass_features_dataframe.columns, "Re-loaded LCMS collection mass features should have cluster assignments."
    
    # Verify parameters were saved and loaded correctly (using the modified parameters from earlier)
    assert lcms_collection2.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold == -1, \
        f"Re-loaded threshold parameter should be -1, got {lcms_collection2.parameters.lcms_collection.alignment_acceptance_fraction_improved_threshold}"
    assert lcms_collection2.parameters.lcms_collection.alignment_acceptance_technique == ['fraction_improved'], \
        f"Re-loaded technique parameter should be ['fraction_improved'], got {lcms_collection2.parameters.lcms_collection.alignment_acceptance_technique}"
    
    """# Make some more plots
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
    """
    
    # Gap fill missing cluster features
    start_time = time.time()
    lcms_collection.fill_missing_cluster_features()
    print("Time to gap fill missing cluster features: ", time.time() - start_time, "seconds, using", ncores, " cores")

    results = lcms_collection.collection_pivot_table()
    results1 = lcms_collection.collection_pivot_table(attribute = 'intensity')
    results2 = lcms_collection.collection_consensus_report(how = 'intensity')
    results3 = lcms_collection.collection_consensus_report(how = 'median')
        
    #TODO: Add code to save and load collection to HDF5 file
    #TODO: Add code to deal with annotations of features in the collection context
