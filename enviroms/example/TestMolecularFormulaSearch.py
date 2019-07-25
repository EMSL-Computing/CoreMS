__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"

import time, sys, os
sys.path.append(".")
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms
from enviroms.emsl.yec.molecular_id.calc.MolecularLookupTable import  MolecularCombinations
from enviroms.emsl.yec.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupTableSetting
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList
from TestMolecularTableLookup import create_lookup_table

def creat_mass_spectrum():

    directory = os.path.join(os.getcwd(), "data/")
    
    file_name = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"
    
    file_location = directory + file_name
   
    #polariy need to be set or read from the file
    polariy = -1

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = Read_MassList(file_location, polariy, delimiter="  ")

    mass_spectrum_obj  = mass_list_reader.get_mass_spectrum(auto_process=True)
    
    return mass_spectrum_obj

if __name__ == "__main__":
    
    time0 = time.time()
    
    dict_molecular_lookup_table = create_lookup_table()
    
    time1 = time.time()
    
    print("create the molecular lookup table took %.2f seconds", time1-time0)
    
    time2 = time.time()
    
    mass_spectrum_obj = creat_mass_spectrum()

    print("Loading the data took %.2f seconds", time2-time1)

    for mspeaks_index in mass_spectrum_obj:
        
        print(mspeaks_index)

