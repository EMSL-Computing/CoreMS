import sys
from pathlib import Path
sys.path.append(".")
from corems.mass_spectrum.input.textMassList import ReadMassList, ReadCoremsMasslist
import pytest

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


def xtest_import_corems_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_SRFA_COREMS.csv"
    
    #polariy need to be set or read from the file
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum()

def test_import_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_ESFA.ascii"
    
    #polariy need to be set or read from the file
    polariy = -1

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadMassList(file_location, polariy, delimiter="  ")

    mass_spectrum = mass_list_reader.get_mass_spectrum(auto_process=True)

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
    
    mass_list_reader = ReadMassList(file_location, polariy, isCentroid=False, delimiter="  ")

    mass_spectrum = mass_list_reader.get_mass_spectrum(auto_process=True)


if __name__ == '__main__':
    test_import_mass_list()
    #test_import_corems_mass_list()