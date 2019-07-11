

import sys
sys.path.append(".")
from enviroms.emsl.yec.lcms.input.win_only.ThermoMSFileReader import ImportLCMSThermoMSFileReader
from enviroms.emsl.yec.encapsulation.settings.ProcessingSetting import MassSpectrumSetting



a = ImportLCMSThermoMSFileReader(
    "C:\\Users\\eber373\\Desktop\\data\\WK_ps_lignin_190301112616.raw")
MassSpectrumSetting.threshold_method = "relative_abudance"
all_scans = a.get_scans_numbers()

print("There are a total of %i scans" % all_scans)

#a.initial_scan_number = 100
#a.final_scan_number = 103

a.start()
"can do something else here in like run a progress bar etc"
a.join()
lcms = a.LCMS
mass_spec = lcms.get_mass_spec_by_scan_number(1)
mass_spec.plot_mz_domain_profile()

