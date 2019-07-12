import sys

sys.path.append(".")
from enviroms.emsl.yec.encapsulation.settings.ProcessingSetting import (
    MassSpectrumSetting,
)
from enviroms.emsl.yec.mass_spectra.input.win_only.ThermoMSFileReader import (
    ImportLCMSThermoMSFileReader,
)


if __name__ == "__main__":
    a = ImportLCMSThermoMSFileReader(
        "C:\\Users\\eber373\\Desktop\\data\\WK_ps_lignin_190301112616.raw"
    )
    MassSpectrumSetting.threshold_method = "signal_noise"
    all_scans = a.get_scans_numbers()

    print("There are a total of %i scans" % all_scans)

    # a.initial_scan_number = 100
    # a.final_scan_number = 103

    a.start()
    '''can do something else here while the code runs
    usefull for progress bars not for speed'''
    a.join()
    lcms = a.LCMS
    mass_spec = lcms.get_mass_spec_by_scan_number(200)
    mass_spec.plot_mz_domain_profile()
    mass_spec.plot_mz_domain_profile_and_noise_threshold()
