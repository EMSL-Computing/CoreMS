from pathlib import Path
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra


if __name__ == "__main__":

    out_path = Path("tmp_data/NMDC_processed_0812/EMSL_49991_Brodie_123_Lipids_Neg_12Aug19_Lola-WCSH417820")
    out_path_hdf5 = str(out_path) + ".corems/" + out_path.stem + ".hdf5"
    parser = ReadCoreMSHDFMassSpectra(out_path_hdf5)
    myLCMSobj = parser.get_lcms_obj()
    myLCMSobj.mass_features_to_df()