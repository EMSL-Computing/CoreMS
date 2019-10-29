
__author__ = "Yuri E. Corilo"
__date__ = "Oct 23, 2019"

import sys
from pathlib import Path
import pytest

sys.path.append(".")
from corems.mass_spectrum.input.boosterHDF import ReadHDF_BoosterMassSpectrum
from corems.mass_spectra.input.boosterHDF import ReadHDF_BoosterMassSpectra

def test_import_booster_mass_spectrum_hdf():

    file_path = Path.cwd() / "tests/tests_data/" / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    
    if file_path.exists:
        #polariy need to be set or read from the file
        polariy = -1

        booster_reader = ReadHDF_BoosterMassSpectrum(file_path, polariy)

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
    
    else:
        
        FileNotFoundError(file_path)



def test_import_booster_mass_spectra_hdf():

    file_path = Path.cwd() / "tests/tests_data/" / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    
    if file_path.exists:
        #polariy need to be set or read from the file
        polariy = -1

        booster_reader = ReadHDF_BoosterMassSpectra(file_path, polariy)

        booster_reader.start()
        booster_reader.join()

        #lcms = booster_reader.get_lcms_obj()


if __name__ == '__main__':

    #test_import_booster_mass_spectrum_hdf()        
    test_import_booster_mass_spectra_hdf()