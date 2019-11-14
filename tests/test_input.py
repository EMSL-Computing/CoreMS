

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


import sys
from pathlib import Path
sys.path.append(".")

import pytest

from corems.transient.input.BrukerSolarix import ReadBrukerSolarix

from corems.mass_spectra.input.brukerSolarix import ReadBruker_SolarixTransientMassSpectra
from corems.mass_spectra.input.boosterHDF5 import ReadHDF_BoosterMassSpectra
from corems.mass_spectra.input.coremsHDF5 import ReadCoreMSHDF_MassSpectra
from corems.mass_spectra.input.massList import ReadCoremsMassSpectraText


from corems.mass_spectrum.input.massList import ReadMassList, ReadCoremsMasslist
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum
from corems.mass_spectrum.input.boosterHDF5 import ReadHDF_BoosterMassSpectrum



def test_import_booster_mass_spectrum_hdf():

    file_path = Path.cwd() / "tests/tests_data/" / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    
    if file_path.exists:
        
        #polariy need to be set or read from the file
        
        booster_reader = ReadHDF_BoosterMassSpectrum(file_path, isCentroid=False)

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
        
def test_import_lcms_from_transient():

    file_location = Path.cwd() / "tests/tests_data/" / "NEG_ESI_SRFA_Auto.d"#"SOM_LC_PeatMix_2p8_0p6_2_30AUG19_GIMLI_ZORBAX-1186_1_01_259.d"

    read_lcms = ReadBruker_SolarixTransientMassSpectra(file_location)

    read_lcms.start()
    read_lcms.join()

    lcms = read_lcms.get_lcms_obj()
    lcms.find_nearest_scan(1)
    lcms.get_scans_number()
    lcms.set_retention_time_from_data()
    lcms.set_tic_list_from_data()
    lcms.get_retention_time()
    lcms.get_tic()
    lcms[0]
    
    for ms in lcms:
        #assign mf
        for mspeak in ms:
            #mspeak.mz_exp,mspeak.mz_abund 
            for mf in mspeak:
                mf.to_string, mf.mz_theor, mf.is_isotopologue    
                pass


def test_import_transient():
    
    # from corems.structure.input.MidasDatFile import ReadMidasDatFile
    
    file_location = Path.cwd() / "ESI_NEG_SRFA.d"

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

def test_import_corems_hdf5():

    file_location = Path.cwd() / "tests/tests_data/" / "NEG_ESI_SRFA_CoreMS.hdf5"
    
    #polariy need to be set or read from the file
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoreMSHDF_MassSpectrum(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum()

    for mspeak in mass_spectrum:
        
        if mspeak:
            
            for mf in mspeak:
                
                print('mass_spectrum', mf.to_string)
    
    read_lc_ms = ReadCoreMSHDF_MassSpectra(file_location)

    read_lc_ms.start()
    read_lc_ms.join()
    
    mass_spectra = read_lc_ms.get_lcms_obj()

    for mspeak in mass_spectra[0]:
        
        if mspeak:
            
            for mf in mspeak:
                
                print('mass_spectra', mf.to_string)
    

def test_import_corems_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_SRFA_COREMS.csv"
    
    #polariy need to be set or read from the file
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(file_location,  analyzer='ICR', instrument_label='12T')

    mass_spectrum = mass_list_reader.get_mass_spectrum()

    for mspeak in mass_spectrum:
        
        if mspeak:
            
            for mf in mspeak:
                print(mf.to_string)

    file_location = Path.cwd() / "tests/tests_data/" /  "NEG_ESI_SRFA_CoreMS.corems"

    read_lc_ms = ReadCoremsMassSpectraText(file_location)

    read_lc_ms.start()
    read_lc_ms.join()
    
    mass_spectra = read_lc_ms.get_lcms_obj()

    for mspeak in mass_spectra[0]:
        
        if mspeak:
            
            for mf in mspeak:
                
                print('mass_spectra', mf.to_string)                

def test_import_mass_list():

    file_location = Path.cwd() / "tests/tests_data/" / "NEG_ESI_SRFA_CoreMS.xlsx"
    
    mass_list_reader = ReadMassList(file_location)

    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_ESFA.ascii"
    
    mass_list_reader = ReadMassList(file_location)

    #polariy need to be set or read from the file
    polarity = -1
     

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadMassList(file_location)

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
    
    mass_list_reader = ReadMassList(file_location, isCentroid=False, )

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity,auto_process=True)

if __name__ == '__main__':
    
    test_import_booster_mass_spectrum_hdf()
    test_import_booster_mass_spectra_hdf()
    test_import_lcms_from_transient()
    test_import_transient()
    test_import_corems_hdf5()
    test_import_corems_mass_list()
    test_import_mass_list()
    