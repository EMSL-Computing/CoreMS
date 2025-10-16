from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    # Set the path to the collection of LCMS runs (previously processed)
    collection_path = Path("/Users/heal742/Library/CloudStorage/OneDrive-PNNL/Documents/_DMS_data/_NMDC/_blanchard_lipidomics/mini_collection_test_out")
    # Path to manifest file
    manifest_file = Path("/Users/heal742/Library/CloudStorage/OneDrive-PNNL/Documents/_DMS_data/_NMDC/_blanchard_lipidomics/mini_collection_test_out/manifest.csv")

    # Set the number of cores to use for loading the data (the parser is parallelized)
    ncores = 5

    # Instantiate the parser
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            manifest_file = manifest_file,
            cores = ncores
            )
    print("Loading LCMS collection with", len(parser.manifest), "samples using", ncores, " cores")

    # Load the LCMS collection (minimally load the data)
    start_time = time.time()
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    print("Time to load LCMS collection ", time.time() - start_time, "seconds -", len(lcms_collection), " LCMS runs and ", ncores, " cores")
    #10s for 7 samples, 10 cores; 162s for 70 samples, 10 cores

    # Check and demonstrate the parsers' ability to load raw data
    lcms_collection.load_raw_data(sample_idx=0, ms_level=1)
    assert lcms_collection[0]._ms_unprocessed[1] is not None, "Raw data for MS1 should be loaded successfully."
    lcms_collection.drop_raw_data(sample_idx=0, ms_level=1)

    # Set flag to call _drop_isotopologue() when running _check_mass_features_df()
    lcms_collection.parameters.lcms_collection.drop_isotopologues = True
    print("Number of total mass features: ", len(lcms_collection.mass_features_dataframe))

    # Align the LCMS runs between each other
    print("Aligning LCMS collection")
    start_time = time.time()
    lcms_collection.align_lcms_objects()
    print("Time to align LCMS collection: ", time.time() - start_time, "seconds") 
    #1.5s for 7 samples; 15s for 70 samples

    # Make some plots 
    lcms_collection.plot_tics(type="both")
    lcms_collection.plot_alignments()    
    # TODO: Think about other plots that would be useful to have here for assessing the quality of the data and alignment

    # Make consensus mass features from the consolidated mass features
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    # THIS FUNCTION NEEDS WORK AND NEEDS MECHANISM TO EVALUATE THE RESULTS (plots? reports?)
    print("Time to roll up consensus mass features: ", time.time() - start_time, "seconds -", len(lcms_collection.mass_features_dataframe), " total mass features", ncores, " cores")   

    lcms_collection.plot_mz_features_across_samples()
    lcms_collection.plot_mz_features_per_cluster()
    lcms_collection.plot_consensus_mz_features() ## zoomed out
    lcms_collection.plot_consensus_mz_features(xb = 10, xt = 15, yb = 500, yt = 600) ## zoomed in 
    lcms_collection.cluster_inspection_plot(11391)
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

    lcms_collection.search_for_missing_mass_features_in_collection()
    lcms_collection.collection_pivot_table(attribute = 'mz')
    lcms_collection.collection_consensus_report(how = 'intensity')
    lcms_collection.collection_consensus_report(how = 'mean')
    lcms_collection.collection_consensus_report(how = 'median')
    
    #TODO: Add code to load and save information about chromatographic settings
    #TODO: Add code to save and load collection to HDF5 file
    #TODO: Generate report of summarize_clusters as a table with both regular and induced mass features