"""
Simple examples demonstrating manifest creation utility usage.
"""

from pathlib import Path
from corems.mass_spectra.input.corems_hdf5 import create_manifest_from_folder

# Example 1: Basic usage - auto-select middle sample as center
print("Example 1: Basic usage (middle sample as center)")
print("-" * 60)
manifest_path = create_manifest_from_folder(
    folder_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2"),
    overwrite=True
)
print()

# Example 2: Custom batch threshold
print("Example 2: Custom batch threshold (24 hours)")
print("-" * 60)
manifest_path = create_manifest_from_folder(
    folder_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2"),
    output_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2/manifest_24hr.csv"),
    batch_time_threshold_hours=24.0,
    overwrite=True
)
print()

# Example 3: Specify a sample as center by name
print("Example 3: Specify center sample by name")
print("-" * 60)
manifest_path = create_manifest_from_folder(
    folder_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2"),
    output_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2/manifest_custom_center.csv"),
    batch_time_threshold_hours=12.0,
    center_name="Blanch_Nat_Lip_C_7_AB_M_07_POS_23Jan18_Brandi-WCSH5801",
    overwrite=True
)
print()

# Example 4: Strict batch threshold
print("Example 4: Strict batch threshold (1 hour)")
print("-" * 60)
manifest_path = create_manifest_from_folder(
    folder_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2"),
    output_path=Path("/Volumes/LaCie/nmdc_data/collection_testing/blanchard_lipid/mini_collection_test_out2/manifest_1hr.csv"),
    batch_time_threshold_hours=1.0,
    overwrite=True
)
print()
