import os, sys
import pathlib
sys.path.append(".")
from enviroms.mass_spectrum.input.textMassList import Read_MassList
import pytest

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"

def test_import_mass_list():

    file_location = os.path.join(os.getcwd(), "tests/tests_data/") + "ESI_NEG_ESFA.ascii"
    
    #polariy need to be set or read from the file
    polariy = -1

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = Read_MassList(file_location, polariy, delimiter="  ")

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
    
    print( mass_spectrum[0])

if __name__ == '__main__':
    test_import_mass_list()