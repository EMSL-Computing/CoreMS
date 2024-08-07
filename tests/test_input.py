__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


import sys

sys.path.append(".")
from pathlib import Path


import pytest
from matplotlib import pyplot


from corems.mass_spectra.input.boosterHDF5 import ReadHDF_BoosterMassSpectra
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.mass_spectra.input.brukerSolarix import ReadBruker_SolarixTransientMassSpectra
from corems.mass_spectra.input.coremsHDF5 import ReadCoreMSHDF_MassSpectra
from corems.mass_spectra.input.massList import ReadCoremsMassSpectraText
from corems.mass_spectrum.input.boosterHDF5 import ReadHDF_BoosterMassSpectrum
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum
from corems.mass_spectrum.input.massList import ReadCoremsMasslist, ReadMassList
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectra.input import rawFileReader


def test_andi_netcdf_gcms():

    file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"

    reader_gcms = ReadAndiNetCDF(file_path)
	
    reader_gcms.run()
    
def test_import_booster_mass_spectrum_hdf():

    file_path = Path.cwd() / "tests/tests_data/ftms/" / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    
    if file_path.exists():
        
        #polarity need to be set or read from the file
        
        booster_reader = ReadHDF_BoosterMassSpectrum(file_path, isCentroid=False)

        mass_spectrum = booster_reader.get_mass_spectrum(auto_process=True)

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

    file_path = Path.cwd() / "tests/tests_data/ftms/" / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    
    if file_path.exists():
        #polarity need to be set or read from the file
        polarity = -1

        booster_reader = ReadHDF_BoosterMassSpectra(file_path, polarity)

        booster_reader.start()
        booster_reader.join()
        #lcms = booster_reader.get_lcms_obj()
        
def test_import_lcms_from_transient():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_Auto.d"#"SOM_LC_PeatMix_2p8_0p6_2_30AUG19_GIMLI_ZORBAX-1186_1_01_259.d"

    MSParameters.mass_spectrum.noise_threshold_method = 'log'
    MSParameters.mass_spectrum.noise_threshold_log_nsigma = 20
    MSParameters.ms_peak.peak_min_prominence_percent = 1
    
    read_lcms = ReadBruker_SolarixTransientMassSpectra(file_location)

    read_lcms.start()
    read_lcms.join()

    lcms = read_lcms.get_lcms_obj()
    lcms.find_nearest_scan(0)
    lcms.scans_number
    lcms.set_retention_time_from_data()
    lcms.set_tic_list_from_data()
    lcms.retention_time
    lcms.tic
    lcms[0]
    
    for ms in lcms:
        #assign mf
        
        for mspeak in ms:
            #mspeak.mz_exp,mspeak.mz_abund 
            for mf in mspeak:
                mf.string, mf.mz_calc, mf.is_isotopologue
                pass

def test_import_transient():
    
    # from corems.structure.input.MidasDatFile import ReadMidasDatFile
    # file_location = Path.cwd() / "tests/tests_data/ftms/SRFAII_20ppm_14Jul2020_IATp08_After_WebEx_1_01_54136.d/"
    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA.d"
    
    MSParameters.transient.apodization_method = "Hanning"
    MSParameters.transient.number_of_truncations = 0
    MSParameters.transient.number_of_zero_fills = 1

    with ReadBrukerSolarix(file_location) as bruker_transient:
        
        #MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
        #MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 1

        #MSParameters.mass_spectrum.noise_threshold_method = 'signal_noise'
        #MSParameters.mass_spectrum.noise_threshold_min_s2n = 50

        MSParameters.mass_spectrum.noise_threshold_method = 'log'
        MSParameters.mass_spectrum.noise_threshold_log_nsigma = 20
        MSParameters.ms_peak.peak_min_prominence_percent = 1
    
        mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)
        #from corems.encapsulation.constant import Labels
        #from corems.mass_spectrum.input import numpyArray
        
        #mass_spectrum_test = numpyArray.ms_from_array_profile(mz=mass_spectrum_obj.mz_exp_profile,
                                                    # abundance=mass_spectrum_obj.abundance_profile,
                                                    # dataname='test',
                                                    # polarity=-1,
                                                    # data_type=Labels.booster_profile,
                                                    # )

        #mass_spectrum_test.plot_mz_domain_profile()

        mass_spectrum_obj.plot_profile_and_noise_threshold()
        
        #pyplot.show()
        
        #mass_spectrum_test.plot_profile_and_noise_threshold()
        
        #mass_spectrum_obj.filter_by_noise_threshold()

        #print(mass_spectrum_obj.get_noise_threshold())     
        
        # pyplot.show()

        #print(len(mass_spectrum_obj))
    
        #print(mass_spectrum_obj.mspeaks[0].mz_exp, mass_spectrum_obj.mspeaks[-1].mz_exp)

def test_import_corems_hdf5():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_CoreMS.hdf5"
    
    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoreMSHDF_MassSpectrum(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum()
    mass_spectrum.to_dataframe()

    for mspeak in mass_spectrum:
        
        if mspeak:
            
            for mf in mspeak:
                
                print('mass_spectrum', mf.string)
    
    read_lc_ms = ReadCoreMSHDF_MassSpectra(file_location)

    read_lc_ms.start()
    read_lc_ms.join()
    
    mass_spectra = read_lc_ms.get_lcms_obj()

    for mspeak in mass_spectra[0]:
        
        if mspeak:
            
            for mf in mspeak:
                
                print('mass_spectra', mf.string)
 
def test_import_corems_mass_list():

    file_location = Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA_COREMS_withdupes.csv"
    
    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(file_location,  analyzer='ICR', instrument_label='12T')

    mass_spectrum = mass_list_reader.get_mass_spectrum(loadSettings=False)

    for mspeak in mass_spectrum:
        
        if mspeak:
            
            for mf in mspeak:
                print(mf.string)

    file_location = Path.cwd() / "tests/tests_data/ftms/" /  "NEG_ESI_SRFA_CoreMS.corems"

    read_lc_ms = ReadCoremsMassSpectraText(file_location)

    read_lc_ms.start()
    read_lc_ms.join()
    
    
    mass_spectra = read_lc_ms.get_lcms_obj()

    for mspeak in mass_spectra[0]:
        
        if mspeak:
            
            for mf in mspeak:
                
                print('mass_spectra', mf.string)                

def test_import_thermo_profile_mass_list():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "Thermo_Profile_MassList.txt" 
    
    mass_list_reader = ReadMassList(file_location, header_lines=7, isCentroid=False, isThermoProfile=True)

    polarity = +1

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity, auto_process=True, loadSettings=False)
    
    #mass_spectrum.plot_profile_and_noise_threshold()
    
    from corems.encapsulation.constant import Labels
    from corems.mass_spectrum.input import numpyArray
    mass_spectrum_test = numpyArray.ms_from_array_profile(mz=mass_spectrum.mz_exp_profile,
                                                abundance=mass_spectrum.abundance_profile,
                                                dataname='test',
                                                polarity=-1,
                                                data_type=Labels.booster_profile)

    mass_spectrum_test.plot_mz_domain_profile()

    # pyplot.show()

def test_import_maglab_pks():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "SRFA.pks"
    
    ref_file_location = Path.cwd() / "tests/tests_data/ftms/SRFA.ref"

    mass_list_reader = ReadMassList(file_location)

    polarity = -1

    #MSParameters.mass_spectrum.min_calib_ppm_error = 3
    #MSParameters.mass_spectrum.max_calib_ppm_error = 4

    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity)

    #MzDomainCalibration(mass_spectrum, ref_file_location).run()

def test_import_mass_list():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_CoreMS.xlsx"
    
    mass_list_reader = ReadMassList(file_location)

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "ESI_NEG_ESFA.ascii"
    
    mass_list_reader = ReadMassList(file_location)

    #polarity need to be set or read from the file
    polarity = -1

    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 1

    # MSParameters.mass_spectrum.noise_threshold_method = 'signal_noise'
    # MSParameters.mass_spectrum.noise_threshold_min_s2n = 100

    #MSParameters.mass_spectrum.noise_threshold_method = 'log'
    #MSParameters.mass_spectrum.noise_threshold_min_std = 32

    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadMassList(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity, auto_process=True)
    
    print(mass_spectrum.baseline_noise, mass_spectrum.baseline_noise_std)
    mass_spectrum.filter_by_noise_threshold()
    print(len(mass_spectrum))
    #mass_spectrum.plot_mz_domain_profile()
    mass_spectrum.plot_profile_and_noise_threshold()
    # pyplot.show()
    print(
        "number_average_molecular_weight",
        mass_spectrum.number_average_molecular_weight(),
    )
    print(
        "weight_average_molecular_weight",
        mass_spectrum.weight_average_molecular_weight(),
    )

    mass_spectrum.filter_by_s2n(100)
    
    # mass_list_reader = ReadMassList(file_location, isCentroid=False,)

    # mass_spectrum = mass_list_reader.get_mass_spectrum(polarity,auto_process=True)

def test_import_thermo_average():
    
    file_location = Path.cwd() / "tests/tests_data/ftms/" / "SRFA_NEG_ESI_ORB.raw"

        # change parameters here
    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 1

    # creates the parser obj
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

    # sums all the mass spectra

    parser.chromatogram_settings.scans = (-1, -1)
    mass_spectrum = parser.get_average_mass_spectrum(spectrum_mode='profile')

    # sums scans in selected range
    parser.chromatogram_settings.scans = (1, 1)
    mass_spectrum = parser.get_average_mass_spectrum(spectrum_mode='profile')

    parser.chromatogram_settings.scans = [1]

    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum(spectrum_mode='profile')

    mass_spectrum.plot_mz_domain_profile()
    mass_spectrum.plot_profile_and_noise_threshold()

    #print("polarity", mass_spectrum.polarity)
    #pyplot.show()

    
if __name__ == '__main__':
    
    pass
    # test_import_booster_mass_spectrum_hdf()
    # test_import_booster_mass_spectra_hdf()
    #test_import_lcms_from_transient()
    #test_import_thermo_profile_mass_list()
    # test_import_transient()
    test_import_corems_hdf5()
    #test_import_corems_mass_list()
    #test_import_mass_list()
    #test_import_maglab_pks()
    #test_andi_netcdf_gcms()
    #test_import_corems_mass_list()
    #test_import_thermo_average()

