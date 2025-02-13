from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    collection_path = Path("/Users/heal742/LOCAL/10_lcms_collection_testing/UDN_neg/processed_data")
    # Will need a bunch of already processed data for this to work (KRH to raw data, process using the lipidomics_workflow script, no need to do MS2 matching)
    manifest_file = collection_path / "manifest_very_small.csv"
    # This file will need to be created by the user or helper script?  For now, we have a small and large manifest file for testing (KRH to provide an example)
    chromatography_file = collection_path / "long_lipid_gradient_chroma.csv"
    # This file will need to be created by the user or helper script?  
    # For now, setting to None should also work, it isn't used at the moment (or KRH to provide an example if setting to None doesn't work)

    # Set the number of cores to use for loading the data (the parser is parallelized)
    ncores = 10

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

    # Honestly, I can't quite remember what this does - I think it is to remove isotopologues from the mass features for future steps
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
    # TODO: Think about other plots that would be useful to have here

    # Consolidate the mass features from the LCMS runs into a single dataframe
    mass_feature_df = lcms_collection.mass_features_to_df()

    # Make consensus mass features from the consolidated mass features
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    # THIS FUNCTION NEEDS WORK AND NEEDS MECHANISM TO EVALUATE THE RESULTS (plots? reports?)
    print("Time to roll up consensus mass features: ", time.time() - start_time, "seconds -", len(lcms_collection.mass_features_dataframe), " total mass features", ncores, " cores")   
    
    #TODO KRH: Add code to load and save information about chromatographic settings
    #TODO KRH: Add code to save and load collection to HDF5 file
    #TODO KRH: Add code to plot a consensus mass feature