__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"

import sys
sys.path.append(".")
from corems.encapsulation.settings.input.ProcessingSetting import MassSpectrumSetting
from corems.mass_spectra.input.win_only.ThermoMSFileReader import ImportLCMSThermoMSFileReader


if __name__ == "__main__":

    file_location = "C:\\Users\\eber373\\Desktop\\data\\WK_ps_lignin_190301112616.raw"

    lcms_reader = ImportLCMSThermoMSFileReader(file_location)

    #setting the threshould method
    MassSpectrumSetting.threshold_method = "signal_noise"
    
    all_scans = lcms_reader.get_scans_numbers()

    print("There are a total of %i scans" % all_scans)

    # a.initial_scan_number = 100
    # a.final_scan_number = 103

    lcms = lcms_reader.get_mass_spectra(auto_process=True)

    """to use the thread
    lcms_reader.start()
    do something
    lcms_reader.join()
    lcms = lcms_reader.get_lcms"""
    kendrick_base = {"C" : 1, "H" : 0, "O" :1}
    for mass_spec in lcms:
       
       # should throw a Exception error because auto_process is set to False
       #print(mass_spec.number_average_molecular_weight())
       #mass_spec.change_kendrick_base_all_mspeaks(kendrick_base)
       for ms_peak in mass_spec:
           
           print(ms_peak.mz_exp)
           print(ms_peak.abundance)
           
           print(ms_peak.kendrick_mass)
           kendrick_base = {"C" : 1, "H" : 0, "O" :1}
           ms_peak.change_kendrick_base(kendrick_base)
           print(ms_peak.kendrick_mass)
           
           for molecular_formula in ms_peak:
               print(molecular_formula.to_string)
               print(molecular_formula.mass_error)

    mass_spec = lcms.get_mass_spec_by_scan_number(200)
    mass_spec.plot_mz_domain_profile()
    mass_spec.plot_mz_domain_profile_and_noise_threshold()

