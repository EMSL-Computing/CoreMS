
import pytest
import sys
from pathlib import Path
sys.path.append(".")

from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak
from corems.transient.input.brukerSolarix import ReadBrukerSolarix

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

def test_mspeaks_fit():

    from matplotlib import pyplot
    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d/"

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)
    i = 0

    mass_spectrum_obj[3].plot_simulation()
    mass_spectrum_obj[3].plot_simulation(sim_type="gaussian", oversample_multiplier=10)
    mass_spectrum_obj[3].plot_simulation()

    # pyplot.show()

def test_mspeak_calculations():

    kendrick_base = {'C': 1, 'H': 2}

    ion_charge = +1
    mz_exp = 212.1234
    abundance = 200
    resolving_power = 1000000
    signal_to_noise = 200
    massspec_index = (300, 300, 300)
    index = 1
    mspeak = ICRMassPeak(ion_charge, mz_exp, abundance,
                         resolving_power, signal_to_noise, massspec_index, index)

    mspeak.change_kendrick_base(kendrick_base)

    mspeak._calc_kmd(kendrick_base)
    mspeak.calc_area()

    mspeak.plot()

    assert round(mspeak.kendrick_mass, 3) == 211.887
    print(round(mspeak.kmd * 100, 0)) == -89
    assert mspeak.knm == 211

    mspeak.set_calc_resolving_power(50, 3)

if __name__ == '__main__':

    # test_mspeaks_fit()
    test_mspeak_calculations()
