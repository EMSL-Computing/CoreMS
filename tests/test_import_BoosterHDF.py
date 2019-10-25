
__author__ = "Yuri E. Corilo"
__date__ = "Oct 23, 2019"

import os, sys
sys.path.append(".")
from corems.mass_spectrum.input.boosterHDF import ReadHDF_Booster
import pytest

def test_import_booster_hdf():

    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
 
    #polariy need to be set or read from the file
    polariy = -1

    booster_reader = ReadHDF_Booster(file_location, polariy)

    mass_spectrum = booster_reader.get_mass_spectrum(auto_process=True, auto_noise=False)

    #mass_spectrum.plot_mz_domain_profile()

    print(
        "number_average_molecular_weight",
        mass_spectrum.number_average_molecular_weight(),
    )
    print(
        "weight_average_molecular_weight",
        mass_spectrum.weight_average_molecular_weight(),
    )

    mass_spectrum.filter_by_s2n(100)
    
    assert round(mass_spectrum[0].mz_exp,3) == 248.961

if __name__ == '__main__':
    test_import_booster_hdf()    

    
 