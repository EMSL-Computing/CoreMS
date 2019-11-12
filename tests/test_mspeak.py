
import pytest
import sys
from pathlib import Path
sys.path.append(".")

from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak
__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

def simulate_peak(mspeak):

    mz_lo, abund_lo = mspeak.lorentz_pdf()
    mz_gaus, abund_gaus = mspeak.gaussian_pdf()

    return round(mz_lo[0], 3), round(mz_gaus[0], 3), round(abund_lo[0], 3), round(abund_gaus[0], 3)
   

def test_mspeak_calculations():

    kendrick_base = {'C': 1, 'H': 2}

    ion_charge = +1
    mz_exp = 212.1234
    abundance = 200
    resolving_power = 1000000
    signal_to_noise = 200
    massspec_index = 300
    index = 1
    mspeak = ICRMassPeak(ion_charge, mz_exp, abundance,
                         resolving_power, signal_to_noise, massspec_index, index)

    mspeak.change_kendrick_base(kendrick_base)

    assert round(mspeak.kendrick_mass, 3) == 211.887
    assert round(mspeak.kmd, 3) == -89.0
    assert mspeak.knm == 211

    assert simulate_peak(mspeak) == (212.123, 212.123, 3.077, 0.0)

    mspeak.set_threoretical_resolving_power(50, 3)

    simulate_peak(mspeak)


if __name__ == '__main__':

    test_mspeak_calculations()
