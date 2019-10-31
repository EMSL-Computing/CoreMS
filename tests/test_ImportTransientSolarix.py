import pickle, sys
from pathlib import Path
sys.path.append(".")

import pytest

from corems.transient.input.BrukerSolarix import ReadBrukerSolarix


__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

def test_import_transient():
    
    # from corems.structure.input.MidasDatFile import ReadMidasDatFile
    
    file_location = Path.cwd() / "tests/tests_data/" / "20190709_WK_CADY_Auto_SRFA_QC_O1_1_01_32.d"

    #setting for signal processing
    apodization_method = "Hanning"
    number_of_truncations = 0
    number_of_zero_fills = 1

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    bruker_transient.set_processing_parameter(
        apodization_method, number_of_truncations, number_of_zero_fills
    )

    mass_spec = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    print(mass_spec.mspeaks[0].mz_exp, mass_spec.mspeaks[-1].mz_exp)
   
    with open("test.pkl", "wb") as file:
        pickle.dump(bruker_transient, file, protocol=pickle.HIGHEST_PROTOCOL)

def test_read_scan_xml():
    
    from bs4 import BeautifulSoup
    
    file_location = Path.cwd() / "tests/tests_data/" / "20190709_WK_CADY_Auto_SRFA_QC_O1_1_01_32.d" / "scan.xml"

    soup = BeautifulSoup(file_location.open(),'xml')

    for s in soup.find_all('scan'):
        
        rt = s.find_all('minutes')[0].text
        tic = s.find_all('tic')[0].text
        
        print(rt,tic )

if __name__ == "__main__":
    #test_import_transient()
    test_read_scan_xml()