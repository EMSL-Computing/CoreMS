import sys

sys.path.append(".")

from enviroms.emsl.yec.mass_spectra.input.win_only.BrukerCompassXtract import (
    ImportLCMSBrukerCompassXtract,
)

if __name__ == "__main__":
    file_location = "C:\\Users\\eber373\\Desktop\\Data\\20190205_WK_SRFA_opt_000001.d"

    lcms_reader = ImportLCMSBrukerCompassXtract(file_location)

    all_scans = lcms_reader.get_scans_numbers()

    print("There are a total of %i scans" % all_scans)
    # a.initial_scan_number = 100
    # a.final_scan_number = 103

    lcms = lcms_reader.get_mass_spectra(auto_process=True)
    """to use the thread
    lcmc_reader.start()
    do something 
    lcmc_reader.join()
    lcms = lcmc_reader.get_lcms"""

    mass_spec = lcms.get_mass_spec_by_scan_number(1)
    mass_spec.plot_mz_domain_profile()
    mass_spec.plot_mz_domain_profile_and_noise_threshold()
