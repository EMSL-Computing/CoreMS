from pathlib import Path

from matplotlib import pyplot

from corems.encapsulation.factory.parameters import MSParameters
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

file_location =  "tests/tests_data/ftms/ESI_NEG_SRFA.d"

MSParameters.transient.apodization_method = "Hanning"
MSParameters.transient.number_of_truncations = 0
MSParameters.transient.number_of_zero_fills = 1

MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 1

#MSParameters.mass_spectrum.noise_threshold_method = 'signal_noise'
#MSParameters.mass_spectrum.noise_threshold_min_s2n = 50

#MSParameters.mass_spectrum.noise_threshold_method = 'log'
#MSParameters.mass_spectrum.noise_threshold_min_std = 32

MSParameters.ms_peak.peak_min_prominence_percent = 1
        
def import_transient():
    
    with ReadBrukerSolarix(file_location) as bruker_transient:

        mass_spectrum = bruker_transient.get_mass_spectrum(plot_result=True, auto_process=True)

        mass_spectrum.plot_profile_and_noise_threshold()

        print("m/z count", len(mass_spectrum))

        print('first m/z', mass_spectrum.mspeaks[0].mz_exp, 'final m/z', mass_spectrum.mspeaks[-1].mz_exp)
    
    return mass_spectrum

mass_spectrum = import_transient()


MSParameters.mass_spectrum.min_calib_ppm_error = -5
MSParameters.mass_spectrum.max_calib_ppm_error = 5

mass_spectrum = import_transient()

ref_file_location = 'tests/tests_data/ftms/SRFA.ref'

MzDomainCalibration(mass_spectrum, ref_file_location).run()

from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.factory.classification import HeteroatomsClassification

mass_spectrum.molecular_search_settings.url_database = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"

mass_spectrum.molecular_search_settings.error_method = 'None'
mass_spectrum.molecular_search_settings.min_ppm_error  = -1
mass_spectrum.molecular_search_settings.max_ppm_error = 1

mass_spectrum.molecular_search_settings.min_dbe = 0
mass_spectrum.molecular_search_settings.max_dbe = 50

mass_spectrum.molecular_search_settings.isProtonated = True 
mass_spectrum.molecular_search_settings.isRadical= False 
mass_spectrum.molecular_search_settings.isadduct = True 

mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1,20)
mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,0)
mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,0)

SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
mass_spectrum.percentile_assigned(report_error=True)

mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)
mass_spectrum_by_classes.plot_ms_assigned_unassigned()

for mspeaks in mass_spectrum.sort_by_abundance():
   if mspeaks: #or just if mspeak:
        for mf in mspeaks:
            print(mf.mz_calc, mf.dbe, mf.class_label, mf.string_formated)

#exporting data
mass_spectrum.to_csv("test")

# save pandas Datarame to pickle
mass_spectrum.to_pandas("test")

# get pandas Dataframe
df = mass_spectrum.to_dataframe()

df.head()