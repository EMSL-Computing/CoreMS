import sys
from pathlib import Path
import pytest

sys.path.append(".")

from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.AutoRecalibration import HighResRecalibration
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.transient.input.brukerSolarix import ReadBrukerSolarix


@pytest.fixture
def mass_spectrum_ftms():
    """Creates a mass spectrum object to be used in the tests"""

    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d/"
    bruker_reader = ReadBrukerSolarix(file_location)
    bruker_transient = bruker_reader.get_transient()

    # Instantiate the mass spectrum object
    mass_spectrum = bruker_transient.get_mass_spectrum(
        plot_result=False, auto_process=False, keep_profile=True
    )

    mass_spectrum.settings.noise_threshold_method = "log"
    mass_spectrum.settings.noise_threshold_log_nsigma = 12
    mass_spectrum.mspeaks_settings.peak_min_prominence_percent = 0.01

    # Process the mass spectrum
    mass_spectrum.process_mass_spec()

    return mass_spectrum

@pytest.fixture
def centroid_mass_spectrum():
    """Creates a centroid mass spectrum object to be used in the tests"""
    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA_UnCal_Unassign.csv"

    MSParameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(
        file_location, analyzer="ICR", instrument_label="12T"
    )

    mass_spectrum = mass_list_reader.get_mass_spectrum(loadSettings=False)

    return mass_spectrum


def test_mz_domain_calibration(mass_spectrum_ftms):
    print("test_mz_domain_calibration")
    mass_spectrum_ftms.settings.min_calib_ppm_error = -10
    mass_spectrum_ftms.settings.max_calib_ppm_error = 10
    mass_spectrum_ftms.filter_by_noise_threshold()

    ref_file_location = Path.cwd() / "tests/tests_data/ftms/SRFA.ref"

    # Check that mass_spectrum_ftms has not been calibrated
    assert mass_spectrum_ftms.calibration_RMS is None
    MzDomainCalibration(mass_spectrum_ftms, ref_file_location).run()

    # Check that the calibration was successful
    assert mass_spectrum_ftms.calibration_RMS < 0.6

