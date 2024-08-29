import pytest
from pathlib import Path

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
def ref_file_location():
    """Returns the location of the reference file for calibration for the tests"""
    return Path.cwd() / "tests/tests_data/ftms/SRFA.ref"
