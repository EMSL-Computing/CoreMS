
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

        mass_spectrum = booster_reader.get_mass_spectrum(auto_process=True, auto_noise=True)

        #mass_spectrum.plot_mz_domain_profile()
        
        print(
            "number_average_molecular_weight",
            mass_spectrum.number_average_molecular_weight(),
        )
        print(
            "weight_average_molecular_weight",
            mass_spectrum.weight_average_molecular_weight(),
        )

        assert round(mass_spectrum[0].mz_exp,3) == 220.147
    
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

    test_import_booster_mass_spectrum_hdf()        
    #test_import_booster_mass_spectra_hdf()
    #import h5py
    #file_path = Path.cwd() / "tests/tests_data/" / "ESFA_100k_9767-13548_chB.A_re_pc_Averaged1min.h5"
    #hdf_obj =  h5py.File(file_path, 'r')
    #x = hdf_obj['2']
    #for key, metadata in x.attrs.items():
    #        print(key, ":", metadata)

    #for key, value in hdf_obj.items():
    #    print(key, value[0])
    #    for keyse, metadata in value.attrs.items():
    #        print(keyse, ":", metadata)


    