from pathlib import Path
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    collection_path = Path("tmp_data/NMDC_processed_collection_0813")
    manifest_file = collection_path / "manifest.csv"
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            manifest_file = manifest_file
            )
    lcms_collection = parser.get_lcms_collection(load_raw=False)
    lcms_collection[0].mass_features_to_df()
    lcms_collection.align_lcms_objects()
    print("Here")