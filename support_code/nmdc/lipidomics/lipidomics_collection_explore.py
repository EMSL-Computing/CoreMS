import pandas as pd
from pathlib import Path

collection_path = Path("/Users/heal742/LOCAL/10_lcms_collection_testing/KidsFirst_T-ALL_neg/processed_files")
collection_mass_features = pd.read_csv(collection_path / "collection_mass_features_ward.csv")

# Pivot the mass features dataframe, sample_id as columns and cluster_id as index
mass_feature_pivot = collection_mass_features.pivot(index='cluster', columns='sample_id', values='intensity')

# Get average mz and rt for each cluster
cluster_avg_mz_rt = collection_mass_features.groupby('cluster').agg({'mz': 'mean', 'scan_time_aligned': 'mean', 'sample_id': 'nunique', 'intensity': 'mean'}).reset_index()

# Add average mz and rt to mass_feature_pivot
mass_feature_pivot = mass_feature_pivot.join(cluster_avg_mz_rt, on='cluster')

# Rearrange columns
mass_feature_pivot = mass_feature_pivot[['mz', 'scan_time_aligned', 'sample_id', 'intensity'] + list(mass_feature_pivot.columns[2:])]

# Save mass_feature_pivot to csv
mass_feature_pivot.to_csv(collection_path / "collection_mass_features_pivot_ward.csv", index=True)

print("here")

