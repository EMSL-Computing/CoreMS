from pathlib import Path
import time
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    collection_path = Path("/Users/heal742/LOCAL/10_lcms_collection_testing/UDN_neg/processed_data")
    manifest_file = collection_path / "manifest_working.csv"
    """
    # Read in manifest file
    import pandas as pd
    manifest_df = pd.read_csv(manifest_file)
    missing_files = []
    # Remove any rows with missing file paths
    for index, row in manifest_df.iterrows():
        sample_folder = Path(row["sample_name"] + ".corems")
        # Check if sample_folder exists in collection_path
        if not (collection_path / sample_folder).exists():
            missing_files.append(row["sample_name"])
            manifest_df.drop(index, inplace=True)


    manifest_df.to_csv(collection_path / "manifest_working.csv", index=False)
    # save missing files to a text file
    with open(collection_path / "missing_files.txt", "w") as f:
        for item in missing_files:
            f.write("%s\n" % item)
    """
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
    lcms_collection.parameters.lcms_collection.drop_isotopologues = True
    print("Number of total mass features: ", len(lcms_collection.mass_features_dataframe))

    print("Aligning LCMS collection")
    start_time = time.time()
    lcms_collection.align_lcms_objects()
    print("Time to align LCMS collection: ", time.time() - start_time, "seconds") 
    #1.5s for 7 samples; 15s for 70 samples

    #lcms_collection.plot_tics(type="both")
    lcms_collection.plot_alignments()
    #mass_feature_df = lcms_collection.mass_features_to_df()
    #print("Adding consensus mass features")
    start_time = time.time()    
    lcms_collection.add_consensus_mass_features()
    print("Time to calculate distance matrices: ", time.time() - start_time, "seconds -", len(lcms_collection.mass_features_dataframe), " total mass features", ncores, " cores")
    # 33 seconds for 22K mass features (7 samples) - 10 cores; 290 seconds for 250K mass features (70 samples) - 10 cores
    

    # Prepare mass features dataframe for export for examination
    collection_mass_features = lcms_collection.mass_features_dataframe

    # Pivot the mass features dataframe, sample_id as columns and cluster_id as index, then remove the index name
    mass_feature_pivot = collection_mass_features.pivot(index='cluster', columns='sample_name', values='area').reset_index()
    mass_feature_pivot.columns.name = None

    # Get average mz and rt for each cluster
    cluster_avg_mz_rt = collection_mass_features.groupby('cluster').agg({'mz': 'median', 'scan_time_aligned': 'median', 'sample_id': 'nunique', 'area': 'max'}).reset_index()
    # Rename sample_id to n_samples
    cluster_avg_mz_rt.rename(columns={'sample_id': 'n_samples', 'area': 'max_area'}, inplace=True)

    # Add average mz and rt to mass_feature_pivot
    mass_feature_pivot = cluster_avg_mz_rt.merge(mass_feature_pivot, on='cluster', how='left')

    # Reorder by mz and reset index
    mass_feature_pivot.sort_values('mz', inplace=True)

    # Save mass_feature_pivot to csv
    mass_feature_pivot.to_csv(collection_path / "collection_mass_features_pivot_rollup_area.csv", index=True)

    print("here")