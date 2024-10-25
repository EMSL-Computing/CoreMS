# Tests for adpodization methods
def test_hamming(bruker_transient):
    """Test the creation of a mass spectrum object with the Hamming apodization method"""
    bruker_transient.set_processing_parameter(apodization_method='Hamming',
                                              number_of_truncations=0,
                                              number_of_zero_fills=1)        
    assert bruker_transient.parameters.apodization_method == "Hamming"                      
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.plot_mz_domain_profile()

def test_blackman(bruker_transient):
    """Test the creation of a mass spectrum object with the Blackman apodization method"""
    bruker_transient.set_processing_parameter(apodization_method='Blackman',
                                                number_of_truncations=0,
                                                number_of_zero_fills=1)
    assert bruker_transient.parameters.apodization_method == "Blackman"
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.plot_mz_domain_profile()

def xtest_fullsine(bruker_transient):
    """Test the creation of a mass spectrum object with the Full-Sine apodization method"""
    # This test is disabled because the Full-Sine apodization method is behaving strangely, see issue #163
    bruker_transient.set_processing_parameter(apodization_method='Full-Sine',
                                                number_of_truncations=0,
                                                number_of_zero_fills=1)
    assert bruker_transient.parameters.apodization_method == "Full-Sine"
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.plot_mz_domain_profile()

def test_kaiser(bruker_transient):
    """Test the creation of a mass spectrum object with the Kaiser apodization method"""
    bruker_transient.set_processing_parameter(apodization_method='Kaiser',
                                                number_of_truncations=0,
                                                number_of_zero_fills=1)
    assert bruker_transient.parameters.apodization_method == "Kaiser"
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.plot_mz_domain_profile()

# Tests for noise threshold methods (note that mass_spectrum_ftms is processed by log in the fixture, no need to test it)
def test_relative_abundance(bruker_transient, mass_spectrum_ftms):
    """Test the creation of a mass spectrum object with the relative abundance noise threshold method"""
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "relative_abundance"
    mass_spectrum_obj.settings.noise_threshold_relative_abundance = 0.01
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.settings.noise_threshold_method == "relative_abundance"
    assert len(mass_spectrum_obj) != len(mass_spectrum_ftms)

def test_absolute_abundance(bruker_transient, mass_spectrum_ftms):
    """Test the creation of a mass spectrum object with the absolute abundance noise threshold method"""
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 20000000
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.settings.noise_threshold_method == "absolute_abundance"
    assert len(mass_spectrum_obj) != len(mass_spectrum_ftms)

def test_signal_to_noise(bruker_transient, mass_spectrum_ftms):
    """Test the creation of a mass spectrum object with the signal to noise noise threshold method"""
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "signal_noise"
    mass_spectrum_obj.settings.noise_threshold_min_s2n = 4
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.settings.noise_threshold_method == "signal_noise"
    assert len(mass_spectrum_obj) != len(mass_spectrum_ftms)

def test_minima(bruker_transient, mass_spectrum_ftms):
    """Test the creation of a mass spectrum object with the minima noise threshold method"""
    mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
    mass_spectrum_obj.settings.noise_threshold_method = "minima"
    mass_spectrum_obj.settings.noise_threshold_min_std = 10
    mass_spectrum_obj.process_mass_spec()
    assert mass_spectrum_obj.settings.noise_threshold_method == "minima"
    assert len(mass_spectrum_obj) != len(mass_spectrum_ftms)

# Tests for peak filtering methods
def test_mass_spectrum_filtering(mass_spectrum_ftms):
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

# Tests for mass spectrum properties
def test_mass_spectrum_properties(mass_spectrum_ftms):
    """Test the properties of the mass spectrum object"""
    res = mass_spectrum_ftms.get_mz_and_abundance_peaks_tuples()
    assert len(res) == len(mass_spectrum_ftms)
    mass_spectrum_ftms.get_masses_count_by_nominal_mass()
    mass_spectrum_ftms.resolving_power_calc(12, 1)