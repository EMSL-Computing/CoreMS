
import sys
sys.path.append(".")

from enviroms.emsl.yec.lcms.input.win_only.BrukerCompassXtract import ImportLCMSBrukerCompassXtract




if __name__ == '__main__':
    file_location = "C:\\Users\\eber373\\Desktop\\Data\\20190205_WK_SRFA_opt_000001.d"
    a = ImportLCMSBrukerCompassXtract("C:\\Users\\eber373\\Desktop\\Data\\20190205_WK_SRFA_opt_000001.d")

    all_scans = a.get_scans_numbers()
    print("There are a total of %i scans" % all_scans)
    #a.initial_scan_number = 100
    #a.final_scan_number = 103

    a.start()
    a.join()
    lcms = a.LCMS
    mass_spec = lcms.get_mass_spec_by_scan_number(1)
    mass_spec.plot_mz_domain_profile()
    mass_spec.plot_mz_domain_profile_and_noise_threshold()
