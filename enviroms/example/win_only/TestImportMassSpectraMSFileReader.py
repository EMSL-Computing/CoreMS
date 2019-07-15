import sys
sys.path.append(".")
from enviroms.emsl.yec.encapsulation.settings.ProcessingSetting import MassSpectrumSetting
from enviroms.emsl.yec.mass_spectra.input.win_only.ThermoMSFileReader import ImportLCMSThermoMSFileReader

if __name__ == "__main__":

    file_location = "C:\\Users\\eber373\\Desktop\\data\\WK_ps_lignin_190301112616.raw"

    lcms_reader = ImportLCMSThermoMSFileReader(file_location)

    #setting the threshould method
    MassSpectrumSetting.threshold_method = "signal_noise"
    
    all_scans = lcms_reader.get_scans_numbers()

    print("There are a total of %i scans" % all_scans)

    # a.initial_scan_number = 100
    # a.final_scan_number = 103

    lcms = lcms_reader.get_mass_spectra(auto_process=False)

    """to use the thread
    lcms_reader.start()
    do something
    lcms_reader.join()
    lcms = lcms_reader.get_lcms"""

    mass_spec = lcms.get_mass_spec_by_scan_number(200)
    mass_spec.plot_mz_domain_profile()
    mass_spec.plot_mz_domain_profile_and_noise_threshold()

