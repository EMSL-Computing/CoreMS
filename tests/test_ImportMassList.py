

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


import sys
from pathlib import Path
sys.path.append(".")

import pytest

from corems.mass_spectrum.input.massList import ReadMassList, ReadCoremsMasslist
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum



def test_import_corems_hdf5():

    file_location = Path.cwd() / "tests/tests_data/" / "neg_esi_srfa_1ppm_test.hdf5"
    
    #polariy need to be set or read from the file
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoreMSHDF_MassSpectrum(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum()

    for mspeak in mass_spectrum:
        
        if mspeak:
            
            for mf in mspeak:
                
                print(mf.to_string)

def test_import_corems_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_SRFA_COREMS.csv"
    
    #polariy need to be set or read from the file
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(file_location, delimiter=",")

    mass_spectrum = mass_list_reader.get_mass_spectrum()

    for mspeak in mass_spectrum:
        if mspeak:
            for mf in mspeak:
                print(mf.to_string)

def test_import_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_ESFA.ascii"
    
    #polariy need to be set or read from the file
    polarity = -1

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadMassList(file_location, delimiter="  ")

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity, auto_process=True)

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
    
    mass_list_reader = ReadMassList(file_location, isCentroid=False, delimiter="  ")

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity,auto_process=True)

if __name__ == '__main__':
    test_import_corems_hdf5()
    test_import_corems_mass_list()
    test_import_mass_list()