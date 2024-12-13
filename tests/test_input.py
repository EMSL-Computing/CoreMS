__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


import sys

from pathlib import Path

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import MSParameters, reset_ms_parameters
from corems.mass_spectra.input import rawFileReader
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.mass_spectra.input.boosterHDF5 import ReadHDF_BoosterMassSpectra
from corems.mass_spectra.input.brukerSolarix import (
    ReadBruker_SolarixTransientMassSpectra,
)
from corems.mass_spectra.input.massList import ReadCoremsMassSpectraText
from corems.mass_spectrum.input.boosterHDF5 import ReadHDF_BoosterMassSpectrum
from corems.mass_spectrum.input.coremsHDF5 import ReadCoreMSHDF_MassSpectrum
from corems.mass_spectrum.input.massList import ReadCoremsMasslist, ReadMassList
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile

def test_andi_netcdf_gcms():
    file_path = (
        Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"
    )

    reader_gcms = ReadAndiNetCDF(file_path)

    reader_gcms.run()

    gcms = reader_gcms.get_gcms_obj()

    assert len(gcms.tic) > 0


def test_import_booster_mass_spectrum_hdf():
    file_path = (
        Path.cwd()
        / "tests/tests_data/ftms/"
        / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    )

    booster_reader = ReadHDF_BoosterMassSpectrum(file_path, isCentroid=False)

    mass_spectrum = booster_reader.get_mass_spectrum(auto_process=False)
    mass_spectrum.parameters = MSParameters(use_defaults=True)
    mass_spectrum.process_mass_spec()

    assert len(mass_spectrum) > 0
    assert mass_spectrum.number_average_molecular_weight() > 0
    assert mass_spectrum.weight_average_molecular_weight() > 0
    assert round(mass_spectrum[0].mz_exp, 3) == 220.147


def test_import_booster_mass_spectra_hdf():
    file_path = (
        Path.cwd()
        / "tests/tests_data/ftms/"
        / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"
    )

    polarity = -1

    booster_reader = ReadHDF_BoosterMassSpectra(file_path, polarity)

    booster_reader.start()
    booster_reader.join()
    mass_spectra = booster_reader.get_lcms_obj()
    assert len(mass_spectra) == 1


def test_import_lcms_from_transient():
    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_Auto.d"

    MSParameters.mass_spectrum.noise_threshold_method = "log"
    MSParameters.mass_spectrum.noise_threshold_log_nsigma = 20
    MSParameters.ms_peak.peak_min_prominence_percent = 1

    read_lcms = ReadBruker_SolarixTransientMassSpectra(file_location)

    read_lcms.start()
    read_lcms.join()

    lcms = read_lcms.get_lcms_obj()
    lcms.find_nearest_scan(0)
    lcms.scans_number
    lcms.set_retention_time_from_data(overwrite=True)
    lcms.set_tic_list_from_data(overwrite=True)
    assert lcms.retention_time[0] > 0
    assert len(lcms.tic) > 0
    assert len(lcms) > 0

    # Return the MSParameters to the default values
    reset_ms_parameters()

def test_import_transient(mass_spectrum_ftms):
    # This test is using the fixture mass_spectrum_ftms
    mass_spectrum_ftms.plot_profile_and_noise_threshold()
    assert len(mass_spectrum_ftms) > 0


def test_import_corems_hdf5():
    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_CoreMS.hdf5"
    
    #load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoreMSHDF_MassSpectrum(file_location)

    # Import processed mass spectrum, check that the mass spectrum is loaded correctly
    mass_spectrum = mass_list_reader.get_mass_spectrum()
    mass_spectrum.to_dataframe()
    assert round(mass_spectrum[0].mz_exp,3) == 576.075
    assert mass_spectrum[0][0].string == 'C25 H20 O16'
    assert mass_spectrum.to_dataframe().shape == (20, 26)
    assert mass_spectrum.settings.noise_threshold_method == 'log'
    assert len(mass_spectrum) == 20

    # Import unprocessed mass spectrum, check that the mass spectrum is loaded correctly
    mass_spectrum2 = mass_list_reader.get_mass_spectrum(
        load_settings=False, 
        auto_process=False,
        load_molecular_formula=False
    )
    mass_spectrum2.parameters.mass_spectrum.noise_threshold_method = 'relative_abundance'

    assert mass_spectrum2.settings.noise_threshold_method == 'relative_abundance' 
    assert len(mass_spectrum2) == 0
 
def test_import_corems_mass_list():
    file_location = (
        Path.cwd() / "tests/tests_data/ftms/ESI_NEG_SRFA_COREMS_withdupes.csv"
    )

    MSParameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(
        file_location, analyzer="ICR", instrument_label="12T"
    )

    mass_spectrum = mass_list_reader.get_mass_spectrum(loadSettings=False)
    assert mass_spectrum.to_dataframe().shape[1] == 26
    assert mass_spectrum.to_dataframe().shape[0] > 0
    assert round(mass_spectrum[0].mz_exp, 0) == 576
    assert mass_spectrum[0][0].string == "C25 H20 O16"

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_CoreMS.corems"

    read_lc_ms = ReadCoremsMassSpectraText(file_location)

    read_lc_ms.start()
    read_lc_ms.join()

    mass_spectra = read_lc_ms.get_lcms_obj()
    assert len(mass_spectra) > 0
    assert mass_spectra[0].to_dataframe().shape[0] > 0
    assert round(mass_spectra[0][0].mz_exp, 0) == 227

    # Return the MSParameters to the default values
    reset_ms_parameters()


def test_import_thermo_profile_mass_list():
    file_location = (
        Path.cwd() / "tests/tests_data/ftms/" / "Thermo_Profile_MassList.txt"
    )

    mass_list_reader = ReadMassList(
        file_location, header_lines=7, isCentroid=False, isThermoProfile=True
    )

    polarity = +1

    mass_spectrum = mass_list_reader.get_mass_spectrum(
        polarity, auto_process=False, loadSettings=False
    )
    mass_spectrum.parameters = MSParameters(use_defaults=True)
    mass_spectrum.process_mass_spec()

    assert mass_spectrum.to_dataframe().shape[0] > 0
    assert round(mass_spectrum[0].mz_exp, 0) == 59


def test_import_numpy_array_profile(mass_spectrum_ftms):
    mass_spectrum_new = ms_from_array_profile(
        mz=mass_spectrum_ftms.mz_exp_profile,
        abundance=mass_spectrum_ftms.abundance_profile,
        dataname="test",
        polarity=-1,
        data_type=Labels.booster_profile,
        auto_process=False
    )
    mass_spectrum_new.parameters = mass_spectrum_ftms.parameters
    mass_spectrum_new.process_mass_spec()
    
    assert mass_spectrum_new.to_dataframe().shape == mass_spectrum_ftms.to_dataframe().shape
    assert round(mass_spectrum_new[0].mz_exp, 0) == round(mass_spectrum_ftms[0].mz_exp, 0)
    assert not mass_spectrum_new.is_centroid

    mass_spectrum_new.plot_mz_domain_profile()


def test_import_maglab_pks():
    file_location = Path.cwd() / "tests/tests_data/ftms/" / "SRFA.pks"

    mass_list_reader = ReadMassList(file_location)

    polarity = -1

    MSParameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity)

    assert mass_spectrum.to_dataframe().shape[0] > 0
    assert round(mass_spectrum[0].mz_exp, 0) == 131

    # Return the MSParameters to the default values
    reset_ms_parameters()


def test_import_xml_mass_list():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "srfa_neg_xml_example.xml"

    mass_list_reader = ReadMassList(file_location, isCentroid=True, isThermoProfile=False)
    polarity = -1

    MSParameters.mass_spectrum.noise_threshold_method = 'absolute_abundance' 
    MSParameters.mass_spectrum.noise_threshold_absolute_abundance = 1000 

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity, auto_process=True, loadSettings=False)
    # check there are lots of peaks (should be ~36k)
    assert len(mass_spectrum)>30_000
    # check the 100th peak is as expected 
    assert round(mass_spectrum.mz_exp[100],3) == 118.049
    
    # Return the MSParameters to the default values
    reset_ms_parameters()


def test_import_xml_mass_list():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "srfa_neg_xml_example.xml"

    mass_list_reader = ReadMassList(file_location, isCentroid=True, isThermoProfile=False)
    polarity = -1

    MSParameters.mass_spectrum.noise_threshold_method = 'absolute_abundance' 
    MSParameters.mass_spectrum.noise_threshold_absolute_abundance = 1000 

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity, auto_process=True, loadSettings=False)
    # check there are lots of peaks (should be ~36k)
    assert len(mass_spectrum)>30_000
    # check the 100th peak is as expected 
    assert round(mass_spectrum.mz_exp[100],3) == 118.049

    # Return the MSParameters to the default values
    reset_ms_parameters()


def test_import_mass_list():
    file_location = Path.cwd() / "tests/tests_data/ftms/" / "NEG_ESI_SRFA_CoreMS.xlsx"

    mass_list_reader = ReadMassList(file_location)

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "ESI_NEG_ESFA.ascii"

    mass_list_reader = ReadMassList(file_location)

    # polarity need to be set or read from the file
    polarity = -1

    mass_list_reader = ReadMassList(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum(polarity, auto_process=False)
    mass_spectrum.parameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    mass_spectrum.parameters.mass_spectrum.noise_threshold_min_relative_abundance = 1
    mass_spectrum.process_mass_spec()

    assert mass_spectrum.baseline_noise > 10000
    assert mass_spectrum.baseline_noise_std > 10000
    mass_spectrum.filter_by_noise_threshold()
    assert mass_spectrum.to_dataframe().shape[0] > 0
    assert len(mass_spectrum) > 0
    assert round(mass_spectrum.number_average_molecular_weight()) > 200
    assert round(mass_spectrum.weight_average_molecular_weight()) > 200

    mass_spectrum.plot_profile_and_noise_threshold()
    og_len = len(mass_spectrum)

    mass_spectrum.filter_by_s2n(100)
    assert len(mass_spectrum) < og_len


def test_import_thermo_average():
    file_location = Path.cwd() / "tests/tests_data/ftms/" / "SRFA_NEG_ESI_ORB.raw"

    # creates the parser obj
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

    # sums all the mass spectra
    parser.chromatogram_settings.scans = (-1, -1)
    mass_spectrum = parser.get_average_mass_spectrum(spectrum_mode="profile", auto_process=False)
    mass_spectrum.parameters = MSParameters(use_defaults=True)
    mass_spectrum.parameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    mass_spectrum.parameters.mass_spectrum.noise_threshold_min_relative_abundance = 1
    mass_spectrum.process_mass_spec()
    assert len(mass_spectrum) == 762

    # sums scans in selected range
    parser.chromatogram_settings.scans = (1, 1)
    mass_spectrum = parser.get_average_mass_spectrum(spectrum_mode="profile")
    mass_spectrum.parameters = MSParameters(use_defaults=True)
    mass_spectrum.parameters.mass_spectrum.noise_threshold_method = "relative_abundance"
    mass_spectrum.parameters.mass_spectrum.noise_threshold_min_relative_abundance = 1
    mass_spectrum.process_mass_spec()
    assert len(mass_spectrum) == 953

    parser.chromatogram_settings.scans = [1]

    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum(spectrum_mode="profile")

    mass_spectrum.plot_mz_domain_profile()
    mass_spectrum.plot_profile_and_noise_threshold()

    assert mass_spectrum.to_dataframe().shape[0] == 1518
    assert round(mass_spectrum[0].mz_exp, 0) == 100
