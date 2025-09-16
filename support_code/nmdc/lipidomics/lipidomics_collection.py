from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection, ReadSavedLCMSCollection
from corems.mass_spectra.output.export import LCMSCollectionExporter

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
    # Check and demonstrate the parsers' ability to load raw data
    lcms_collection.load_raw_data(sample_idx=0, ms_level=1)
    assert lcms_collection[0]._ms_unprocessed[1] is not None, "Raw data for MS1 should be loaded successfully."
    lcms_collection.drop_raw_data(sample_idx=0, ms_level=1)
    #10s for 7 samples, 10 cores; 162s for 70 samples, 10 cores


    # Set flag to call _drop_isotopologue() when running _check_mass_features_df()
    lcms_collection.parameters.lcms_collection.drop_isotopologues = True
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

    # Make some plots 
    lcms_collection.plot_tics(type="both")
    lcms_collection.plot_alignments()  
    #1.5s for 7 samples; 15s for 70 samples
    

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
    # BROKEN - needs fixing
    #lcms_collection.plot_cluster_outlier_frequency(dim_list, clu_size_thresh = 0.25)

    # Save the LCMS collection to a new location
    print("Saving LCMS collection to test_lcms_collection_out.hdf5")
    exporter = LCMSCollectionExporter(
        out_file_path="test_lcms_collection_out",
        mass_spectra_collection=lcms_collection
    )
    exporter.export_to_hdf5(overwrite=True)

    # Reload the LCMS collection from the saved location and check that we can load raw data
    print("Reloading LCMS collection from test_lcms_collection_out.hdf5 and checking functionality")
    parser2 = ReadSavedLCMSCollection(
        collection_hdf5_path=Path("test_lcms_collection_out.hdf5"),
        cores=ncores
    )
    lcms_collection2 = parser2.get_lcms_collection(load_raw=False, load_light=True)
    lcms_collection2.plot_alignments() 
    lcms_collection2.plot_consensus_mz_features() ## zoomed in 
    lcms_collection2.plot_consensus_mz_features(xb = 10, xt = 15, yb = 500, yt = 600) ## zoomed in 
    lcms_collection2.cluster_inspection_plot(11391)
    #lcms_collection2.plot_cluster_outlier_frequency(dim_list, clu_size_thresh = 0.25)
    assert "cluster" in lcms_collection2.mass_features_dataframe.columns, "Reloaded LCMS collection should have cluster assignments in mass_features_dataframe."
    assert lcms_collection2.rt_aligned, "Reloaded LCMS collection should be marked as retention time aligned." 
    lcms_collection2.load_raw_data(sample_idx=0, ms_level=1)
    assert lcms_collection2[0]._ms_unprocessed[1] is not None, "Raw data for MS1 should be loaded successfully."
    lcms_collection2.drop_raw_data(sample_idx=0, ms_level=1)
    del parser2, lcms_collection2
    print("Reloaded LCMS collection successfully and checked functionality, then deleted temporary objects.")


    ## WORK IN PROGRESS: temporary code for testing
    ## want to adjust function to iterate throught samples by index, not name
    ## want to be able to do that in parallel/multiprocess
    samplename = 'Blanch_Nat_Lip_H_11_AB_M_13_POS_23Jan18_Brandi-WCSH5801'    
    lcms_collection.search_for_missing_mass_features_in_one_sample(samplename)
    print(lcms_collection._lcms[samplename].mass_features_to_df(induced_features = True))

    #TODO: Add code to plot a consensus mass feature