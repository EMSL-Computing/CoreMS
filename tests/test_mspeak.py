
import pytest
import sys
from pathlib import Path
sys.path.append(".")

from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak
__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

def simulate_peak(mspeak):

    mz_lo, abund_lo = mspeak.lorentz_pdf()
    mz_gauss, abund_gauss = mspeak.gaussian_pdf()

    return round(mz_lo[0], 3), round(mz_gauss[0], 3), round(abund_lo[0], 3), round(abund_gauss[0], 3)
 
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
    
    mspeak._calc_kdm(kendrick_base)
    
    mspeak.lorentz_pdf()
    mspeak.gaussian_pdf()

    assert round(mspeak.kendrick_mass, 3) == 211.887
    print( round(mspeak.kmd*100, 0)) == -89
    assert mspeak.knm == 211

    print(simulate_peak(mspeak))
    assert simulate_peak(mspeak) == (211.9, 211.9, 0.0, 0.0)

    mspeak.set_theoretical_resolving_power(50, 3)

    simulate_peak(mspeak)


if __name__ == '__main__':

    test_mspeak_calculations()
