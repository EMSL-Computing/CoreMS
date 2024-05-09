# %% Import libs
import sys

sys.path.append("./")
from pathlib import Path

from corems.mass_spectra.input.coremsHDF5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from corems.mass_spectra.output.export import LCMSExport


def test_lcms_peakpick():
    # Instantiate parser based on binary file type
    file_raw = (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.raw"
    )
    parser = ImportMassSpectraThermoMSFileReader(file_raw)

    # Instatiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="ms1", verbose=False)

    # Set parameters on the LCMS object that are reasonable for testing
    myLCMSobj.parameters.lc_ms.ph_inten_min_rel = 0.0005
    myLCMSobj.parameters.lc_ms.ph_persis_min_rel = 0.05
    myLCMSobj.parameters.lc_ms.ph_smooth_it = 0
    myLCMSobj.parameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    myLCMSobj.parameters.mass_spectrum.noise_threshold_min_relative_abundance = 1
    myLCMSobj.parameters.mass_spectrum.noise_min_mz = 0
    myLCMSobj.parameters.mass_spectrum.noise_max_mz = 2500
    myLCMSobj.parameters.mass_spectrum.min_picking_mz = 0
    myLCMSobj.parameters.mass_spectrum.max_picking_mz = 2500

    # Use persistent homology to find mass features in the lc-ms data
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.
    myLCMSobj.find_mass_features(verbose=False)
    myLCMSobj.cluster_mass_features(verbose=False)
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )
    myLCMSobj.integrate_mass_features(drop_if_fail=True)
    myLCMSobj.find_c13_mass_features(verbose=False)
    
    assert len(myLCMSobj.mass_features) == 130

    # Add ms2 data to lcms object
    myLCMSobj.add_associated_ms2_dda(spectrum_mode="centroid")

    # Export the mass features to a pandas dataframe
    df = myLCMSobj.mass_features_to_df()
    assert df.shape == (130, 9)

    # Plot a mass feature
    myLCMSobj.mass_features[1].plot(return_fig=False)

    # Export the lcms object to an hdf5 file
    exporter = LCMSExport(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801", myLCMSobj
    )
    exporter.to_hdf(overwrite=True)

    # Import the hdf5 file, assert that its df is same as above and that we can plot a mass feature
    parser = ReadCoreMSHDFMassSpectra(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems/Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.hdf5"
    )
    myLCMSobj2 = parser.get_lcms_obj()
    assert myLCMSobj2.mass_features_to_df().shape == (130, 9)
    myLCMSobj2.mass_features[1].plot(return_fig=False)

if __name__ == "__main__":
    test_lcms_peakpick()
