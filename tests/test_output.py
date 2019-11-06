
import sys
from pathlib import Path
sys.path.append(".")

import pytest

from corems.mass_spectrum.output.export import MassSpecExport
from corems.mass_spectrum.input.massList import ReadCoremsMasslist


def import_corems_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_SRFA_COREMS.csv"
    
    #polariy need to be set or read from the file
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(file_location, delimiter=",")

    mass_spectrum = mass_list_reader.get_mass_spectrum()

    return mass_spectrum

def test_export():

    mass_spectrum = import_corems_mass_list()

    exportMS= MassSpecExport('neg_esi_srfa_1ppm_test', mass_spectrum)

    exportMS.to_pandas()
    exportMS.to_excel()
    exportMS.to_csv()
    exportMS.to_hdf()

if __name__ == "__main__":
                        
    test_export()

    #data_list_dict = exportMS.get_list_dict_data()
    import h5py
    file1 = h5py.File('neg_esi_srfa_1ppm_test.hdf5', 'r')
    
    print(file1['0'].attrs.keys())
    #print(file1['0'].attrs['test'])
    #for f in file1['0']:
    #    print(type(f), f)