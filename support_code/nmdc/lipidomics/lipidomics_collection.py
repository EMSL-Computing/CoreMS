from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    # Set the path to the collection of LCMS runs (previously processed)
    collection_path = Path("/Users/heal742/Library/CloudStorage/OneDrive-PNNL/Documents/_DMS_data/_NMDC/_blanchard_lipidomics/mini_collection_test_out")
    # Path to manifest file
    manifest_file = Path("/Users/heal742/Library/CloudStorage/OneDrive-PNNL/Documents/_DMS_data/_NMDC/_blanchard_lipidomics/mini_collection_test_out/manifest.csv")
    # This file will need to be created by the user or helper script?
    chromatography_file = Path("/Users/heal742/LOCAL/10_lcms_collection_testing/UDN_neg/processed_data/long_lipid_gradient_chroma.csv") 

    # Set the number of cores to use for loading the data (the parser is parallelized)
    ncores = 5

    # Instantiate the parser
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            manifest_file = manifest_file,
            chromatography_file=chromatography_file,
            cores = ncores
            )
    print("Loading LCMS collection with", len(parser.manifest), "samples using", ncores, " cores")

    # Load the LCMS collection (minimally load the data)
    start_time = time.time()
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    print("Time to load LCMS collection ", time.time() - start_time, "seconds -", len(lcms_collection), " LCMS runs and ", ncores, " cores")
    #10s for 7 samples, 10 cores; 162s for 70 samples, 10 cores

    # Check parsers and the ability to load raw data
    #TODO KRH: Make this into two methods on the collection - lcms_collection.load_raw_data(sample_idx);   lcms_collection.drop_raw_data(sample_idx)
    # Add error handling into those methods; make sure it works with both thermo and mzml data
    lcms_collection[0]._ms_unprocessed[1] = lcms_collection[0].spectra_parser.get_ms_raw(spectra="ms1", scan_df=lcms_collection[0].scan_df)['ms1']

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
    
    #TODO: Add code to load and save information about chromatographic settings
    #TODO: Add code to save and load collection to HDF5 file
    #TODO: Add code to plot a consensus mass feature