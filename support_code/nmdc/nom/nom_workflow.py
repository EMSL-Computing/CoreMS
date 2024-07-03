import warnings


warnings.filterwarnings("ignore")

import sys
sys.path.append("./")


from pathlib import Path
import cProfile
import json
import pstats
from multiprocessing import Pool, Process

from pandas import DataFrame
from matplotlib import pyplot as plt

from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

from corems.mass_spectrum.input.massList import ReadMassList
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems import SuppressPrints, get_filename, get_dirnames
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.mass_spectra.input import rawFileReader
from corems.encapsulation.factory.parameters import MSParameters
from support_code.nmdc.metadata.nmdc_registration import DMS_Mapping, NMDC_Metadata

def run_bruker(file_location):

    with ReadBrukerSolarix(file_location) as transient:

        MSParameters.mass_spectrum.noise_threshold_methodmethod = 'log'
        MSParameters.mass_spectrum.noise_threshold_min_s2n = 6

        mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)
        # mass_spectrum.plot_profile_and_noise_threshold()
        # plt.show()

        return mass_spectrum, transient.transient_time

def run_thermo(file_location):

    MSParameters.mass_spectrum.noise_threshold_method = 'log'
    MSParameters.mass_spectrum.noise_threshold_min_std = 3
    MSParameters.ms_peak.peak_min_prominence_percent = 0.2

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

    # mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)

    mass_spectrum = parser.get_average_mass_spectrum()

    return mass_spectrum, 3


def calspec(msobj, refmasslist, order=2):

    calfn = MzDomainCalibration(msobj, refmasslist)
    ref_mass_list_fmt = calfn.load_ref_mass_list()

    imzmeas, mzrefs = calfn.find_calibration_points(ref_mass_list_fmt,
                                                    calib_ppm_error_threshold=(-1.0, 1.0),
                                                    calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-1.5, 1.5),
                                                        calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-3, 3),
                                                        calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-5, 5),
                                                        calib_snr_threshold=4)    

    if len(mzrefs) < 5:

        imzmeas, mzrefs = calfn.find_calibration_points(ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-7, 7),
                                                        calib_snr_threshold=4) 

    if len(mzrefs) < 5:

        imzmeas, mzrefs = calfn.find_calibration_points(ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-10, 10),
                                                        calib_snr_threshold=4)

    calfn.recalibrate_mass_spectrum(imzmeas, mzrefs, order=order)

def set_parameters(mass_spectrum, field_strength=12, pos=False):

    if field_strength == 12:

        mass_spectrum.settings.max_calib_ppm_error = 10
        mass_spectrum.settings.min_calib_ppm_error = -10

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -1
        mass_spectrum.molecular_search_settings.max_ppm_error = 1

        mass_spectrum.settings.calib_sn_threshold = 4

    elif field_strength == 15:

        mass_spectrum.settings.max_calib_ppm_error = 3
        mass_spectrum.settings.min_calib_ppm_error = -3

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -1
        mass_spectrum.molecular_search_settings.max_ppm_error = 1

        mass_spectrum.settings.calib_sn_threshold = 4

    elif field_strength == 7:

        mass_spectrum.settings.max_calib_ppm_error = 3
        mass_spectrum.settings.min_calib_ppm_error = -3

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -1
        mass_spectrum.molecular_search_settings.max_ppm_error = 1

        mass_spectrum.settings.calib_sn_threshold = 4
    # 21T Data
    else:

        mass_spectrum.settings.max_calib_ppm_error = 1
        mass_spectrum.settings.min_calib_ppm_error = -1

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -0.5
        mass_spectrum.molecular_search_settings.max_ppm_error = 0.5

    mass_spectrum.molecular_search_settings.url_database = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 40

    if pos:

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 200)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 12)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 3)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)

    else:

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 200)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 22)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)

    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 0)

    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical = False
    mass_spectrum.molecular_search_settings.isAdduct = False

def run_nmdc_workflow(args):
    # mass_spectrum = get_masslist(file_location)
    try:
        file_location, ref_calibration_file, field_strength = args

        if field_strength == 21:

            # return "21T", None
            # print("{}   {}".format("21T", file_location))
            print("{} {}  {}".format("processing", field_strength, file_location))
            mass_spectrum, transient_time = run_thermo(file_location)

        else:

            print("{} {}  {}".format("processing", field_strength, file_location))
            mass_spectrum, transient_time = run_bruker(file_location)
            # return "not 21T", None

        is_pos = True if mass_spectrum.polarity > 0 else False
            
        if is_pos:
            
            return "positive", None

        if len(mass_spectrum) < 30:

            print("{}   {}".format("too few peaks", file_location))
            return "too few peaks", None

        set_parameters(mass_spectrum, field_strength=field_strength, pos=is_pos)

        if ref_calibration_file:

            calspec(mass_spectrum, ref_calibration_file)
            # MzDomainCalibration(mass_spectrum, ref_calibration_file).run()

        # mass_spectrum.filter_by_max_resolving_power(field_strength, transient_time)

        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)
        mass_spectrum.molecular_search_settings.score_method = "prob_score"
        mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

        return "all_good", mass_spectrum
    
    except Exception as inst:

        print(type(inst))    # the exception instance
        print(inst.args)     # arguments stored in .args
        print(inst)
        return 'error', None

def monitor(target):

    import psutil
    import time

    worker_process = Process(target=target)
    worker_process.start()
    p = psutil.Process(worker_process.pid)

    # log cpu usage of `worker_process` every 10 ms
    cpu_percents = []
    while worker_process.is_alive():
        cpu_percents.append(p.cpu_percent())
        time.sleep(0.01)

    worker_process.join()
    return cpu_percents

def worker(file_location):

    cProfile.runctx('run_assignment(file_location)', globals(), locals(), 'di-fticr-di.prof')
    # stats = pstats.Stats("topics.prof")
    # stats.strip_dirs().sort_stats("time").print_stats()

def run_nom_nmdc_data_processing():

    file_ext = '.raw'  # '.d' 
    data_dir = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/")
    dms_file_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/SPRUCE_FTICR_Peat.xlsx")

    results_dir = Path("results/")
    registration_path = results_dir / "spruce_ftms_nom_data_products.json"
    failed_files = results_dir / "nom_failed_files.json"
    pos_files = results_dir / "pos_files.json"

    field_strength = 21
    cores = 4
    ref_calibration_path = False

    # file_paths = get_dirnames()
    ref_calibration_path = Path("db/Hawkes_neg.ref")

    dms_mapping = DMS_Mapping(dms_file_path)
    selected_files = dms_mapping.get_selected_sample_list()

    failed_list = []
    pos_list = []
    
    for file_name, field_strength in selected_files:

        in_file_path = data_dir / file_name / file_name.with_suffix(file_ext)    

        try:

            issue, ms = run_nmdc_workflow((in_file_path, ref_calibration_path, field_strength))

            if ms:

                file_name = Path(in_file_path.name)

                output_file_path = results_dir / file_name.with_suffix('.csv')
                ms.to_csv(output_file_path, write_metadata=False)

                nmdc = NMDC_Metadata(in_file_path, ref_calibration_path, output_file_path, dms_file_path)
                nmdc.create_nmdc_ftms_metadata(ms, registration_path)

            else:
                print(issue)

        except Exception as inst:

            print(type(inst))    # the exception instance
            print(inst.args)     # arguments stored in .args
            print(inst)
            failed_list.append(str(in_file_path))

    with failed_files.open('w') as json_file:

        json_file.write(json.dumps(failed_list, indent=1))


if __name__ == "__main__":

    # run_multiprocess()
    # cpu_percents = monitor(target=run_multiprocess)
    # print(cpu_percents)
    run_nom_nmdc_data_processing()
    # file_location = get_filename()
    # if file_location:
    #    run_assignment(file_location)
