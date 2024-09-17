from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    collection_path = Path("/Users/heal742/LOCAL/10_lcms_collection_testing/KidsFirst_T-ALL_neg/processed_files")
    manifest_file = collection_path / "manifest_small.csv"
    
    '''
    # Read in manifest file
    import pandas as pd
    manifest_df = pd.read_csv(manifest_file)
    # Remove any rows with missing file paths
    for index, row in manifest_df.iterrows():
        sample_folder = Path(row["sample_name"] + ".corems")
        # Check if sample_folder exists in collection_path
        if not (collection_path / sample_folder).exists():
            manifest_df.drop(index, inplace=True)

    manifest_df.to_csv(manifest_file, index=False)
    '''
    ncores = 10
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            manifest_file = manifest_file,
            cores = ncores
            )
    
    print("Loading LCMS collection with", len(parser.manifest), "samples using", ncores, " cores")
    start_time = time.time()
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    print("Time to load LCMS collection ", time.time() - start_time, "seconds -", len(lcms_collection), " LCMS runs and ", ncores, " cores") 
    #10s for 7 samples, 10 cores; 162s for 70 samples, 10 cores
    lcms_collection.mass_features_dataframe


    print("Aligning LCMS collection")
    start_time = time.time()
    lcms_collection.align_lcms_objects()
    print("Time to align LCMS collection: ", time.time() - start_time, "seconds") 
    #1.5s for 7 samples; 15s for 70 samples

    #lcms_collection.plot_tics(type="both")
    #lcms_collection.plot_alignments()
    #mass_feature_df = lcms_collection.mass_features_to_df()
    #print("Adding consensus mass features")
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    print("Time to calculate distance matrices: ", time.time() - start_time, "seconds -", len(lcms_collection.mass_features_dataframe), " total mass features", ncores, " cores")
    # 33 seconds for 22K mass features (7 samples) - 10 cores; 290 seconds for 250K mass features (70 samples) - 10 cores
    lcms_collection.mass_features_dataframe.to_csv(collection_path / "collection_mass_features_ward.csv", index=True)
    print("Here")

    #lcms_collection.