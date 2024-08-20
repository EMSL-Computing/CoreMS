from pathlib import Path
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectraCollection


if __name__ == "__main__":
    collection_path = Path("tmp_data/NMDC_processed_collection_0819")
    manifest_file = collection_path / "manifest.csv"
    parser = ReadCoreMSHDFMassSpectraCollection(
            folder_location = collection_path,
            manifest_file = manifest_file
            )
    lcms_collection = parser.get_lcms_collection(load_raw=False, load_light=True)
    #lcms_collection.plot_tics()
    lcms_collection.align_lcms_objects()
    #lcms_collection.plot_tics(type="both")
    #mass_feature_df = lcms_collection.mass_features_to_df()
    lcms_collection.add_consensus_mass_features()
    print("Here")

    #lcms_collection.