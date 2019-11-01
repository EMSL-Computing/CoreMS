import sys
from pathlib import Path
sys.path.append(".")

import pytest
from matplotlib import pyplot
from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak


__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

def simulate_peak(mspeak):

    mz_lo, abund_lo = mspeak.lorentz_pdf()
    mz_gaus, abund_gaus = mspeak.gaussian_pdf()
    
    pyplot.plot(mz_lo, abund_lo, 'r', label='Lorentz')
    pyplot.plot(mz_gaus, abund_gaus, 'g', label='Gauss')
    pyplot.legend()
    pyplot.show()

def test_mspeak_calculations():

    kendrick_base = {'C': 1, 'H': 2}
    
    mspeak = ICRMassPeak(+1, 212.1234, 200, 1000000, 200, 300, 1)
    
    mspeak.change_kendrick_base(kendrick_base)

    assert round(mspeak.kendrick_mass, 3) == 211.887
    assert round(mspeak.kmd, 3) == -89.0
    assert mspeak.knm == 211

    simulate_peak(mspeak)

    mspeak.set_threoretical_resolving_power(50, 3)

    simulate_peak(mspeak)

if __name__ == '__main__':
    test_mspeak_calculations()