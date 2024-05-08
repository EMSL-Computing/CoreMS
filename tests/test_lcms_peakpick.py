# %% Import libs
import sys

sys.path.append("./")
from corems.mass_spectra.input.rawFileReader import ImportMassSpectraThermoMSFileReader
from pathlib import Path


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
    myLCMSobj.parameters.lc_ms.ph_inten_min = (
        myLCMSobj._ms_unprocessed[1].intensity.max() * 0.05 * 0.01
    )
    myLCMSobj.parameters.lc_ms.ph_persis_min = (
        myLCMSobj._ms_unprocessed[1].intensity.max() * 0.05
    )
    myLCMSobj.parameters.lc_ms.ph_smooth_it = 0

    # Use persistent homology to find mass features in the lc-ms data
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.
    myLCMSobj.find_mass_features(verbose=False)
    myLCMSobj.cluster_mass_features(verbose=False)
    myLCMSobj.integrate_mass_features(drop_if_fail=True)
    myLCMSobj.find_c13_mass_features(verbose=False)

    assert len(myLCMSobj.mass_features) == 130

    # Plot a mass feature
    myLCMSobj.mass_features[1].plot(return_fig=False)


if __name__ == "__main__":
    test_lcms_peakpick()
