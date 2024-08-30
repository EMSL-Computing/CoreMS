from corems.encapsulation.factory.parameters  import MSParameters
from corems.encapsulation.factory.processingSetting  import TransientSetting
from corems.mass_spectra.input.brukerSolarix import ReadBrukerSolarix


def xtest_hamming_sn(ftms_file_location):
    """Test the creation and processing of a mass spectrum object with the Hamming apodization method and signal to noise noise threshold method"""
    # MSParameters.transient.apodization_method = 'Hamming' 
    # ABOVE DOES NOT WORK TO SET THE APDIZATION METHOD

    # TransientSetting.apodization_method = 'Hamming'
    # ABOVE DOES NOT WORK TO SET THE APDIZATION METHOD

    #bruker_transient.parameters.apodization_method = 'Hamming'
    # ABOVE DOES NOT WORK TO SET THE APDIZATION METHOD
    bruker_reader = ReadBrukerSolarix(ftms_file_location)
    bruker_transient = bruker_reader.get_transient()

    mass_spectrum_obj = bruker_transient.get_mass_spectrum( plot_result=False, auto_process=False)

    #mass_spectrum_obj.settings.noise_threshold_method = "signal_noise"
    #mass_spectrum_obj.settings.noise_threshold_min_s2n = 15
    #mass_spectrum_obj.process_mass_spec()
    #assert mass_spectrum_obj.transient_settings.apodization_method == "Hamming"
    #assert len(mass_spectrum_obj) > 0 #504?

def xtest_blackman_relative_abundance(bruker_transient):
    """Test the creation and processing of a mass spectrum object with the Blackman apodization method and relative abundance noise threshold method"""
    mass_spectrum_obj = bruker_transient.get_mass_spectrum( plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "relative_abundance"
    mass_spectrum_obj.settings.noise_threshold_min_relative_abundance = 20
    mass_spectrum_obj.transient_settings.apodization_method = "Blackman"
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.transient_settings.apodization_method == "Blackman"
    assert len(mass_spectrum_obj) > 0 #212?

def xtest_fullsine_absolute_abundance(bruker_transient):
    """Test the creation and processing of a mass spectrum object with the full-sine apodization method and absolute abundance noise threshold method"""
    
    mass_spectrum_obj = bruker_transient.get_mass_spectrum( plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 100
    mass_spectrum_obj.transient_settings.apodization_method = "Full-Sine"
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.transient_settings.apodization_method == "Full-Sine"
    assert len(mass_spectrum_obj) > 0

def xtest_kaiser_minima(bruker_transient):
    """Test the creation and processing of a mass spectrum object with the Kaiser apodization method and signal to noise noise threshold method"""
    mass_spectrum_obj = bruker_transient.get_mass_spectrum( plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "minima"
    mass_spectrum_obj.settings.noise_threshold_min_std = 10
    mass_spectrum_obj.transient_settings.apodization_method = "Kaiser"
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.transient_settings.apodization_method == "Kaiser"
    assert len(mass_spectrum_obj) > 0

def xtest_mass_spectrum_filtering(mass_spectrum_ftms):
    """Test the filtering methods of the mass spectrum object"""
    og_peaks = len(mass_spectrum_ftms)

    # Test the filtering methods and check that the number of peaks has decreased each time
    mass_spectrum_ftms.filter_by_mz(200, 1000)
    assert len(mass_spectrum_ftms) < og_peaks
    mass_spectrum_ftms.reset_indexes()
    assert len(mass_spectrum_ftms) == og_peaks
    mass_spectrum_ftms.filter_by_abundance(0, 1000)
    assert len(mass_spectrum_ftms) < og_peaks
    mass_spectrum_ftms.reset_indexes()
    mass_spectrum_ftms.filter_by_max_resolving_power(2, 3)
    assert len(mass_spectrum_ftms) < og_peaks
    mass_spectrum_ftms.reset_indexes()
    mass_spectrum_ftms.filter_by_min_resolving_power(15, 3)
    assert len(mass_spectrum_ftms) < og_peaks
    mass_spectrum_ftms.reset_indexes()
    mass_spectrum_ftms.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_ftms.settings.noise_threshold_absolute_abundance = 10000000
    mass_spectrum_ftms.filter_by_noise_threshold()
    assert len(mass_spectrum_ftms) < og_peaks

def xtest_mass_spectrum_properties(mass_spectrum_ftms):
    """Test the properties of the mass spectrum object"""
    res = mass_spectrum_ftms.get_mz_and_abundance_peaks_tuples()
    assert len(res) == len(mass_spectrum_ftms)
    mass_spectrum_ftms.get_masses_count_by_nominal_mass()
    mass_spectrum_ftms.resolving_power_calc(12, 1)


