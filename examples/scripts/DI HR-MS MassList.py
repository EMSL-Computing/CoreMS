import warnings

from pandas.core.frame import DataFrame
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path
import cProfile
import json
import pstats
from multiprocessing import Pool, Process

import pandas as pd
from matplotlib import pyplot as plt

from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

from corems.mass_spectrum.input.massList import ReadMassList
from corems.molecular_id.factory.classification import HeteroatomsClassification
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems import SuppressPrints, get_filename, get_filenames
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.mass_spectra.input import rawFileReader

from corems.encapsulation.constant import Atoms
from corems.encapsulation.factory.parameters import MSParameters

def mzdomain_calibration(mass_spectrum):

    mass_spectrum.settings.min_calib_ppm_error = 0
    mass_spectrum.settings.max_calib_ppm_error = 1

    #file_location = Path.cwd() / "tests/tests_data/ESI_NEG_SRFA.d/"

def run_bruker(file_location):

    with ReadBrukerSolarix(file_location) as transient:

        MSParameters.mass_spectrum.noise_threshold_method = 'log'
        MSParameters.mass_spectrum.noise_threshold_min_s2n = 6

        mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)
        # mass_spectrum.plot_profile_and_noise_threshold()
        # plt.show()
        # find_formula_thread = FindOxygenPeaks(mass_spectrum)
        # find_formula_thread.run()

        # mspeaks_results = find_formula_thread.get_list_found_peaks()
        
        # mass_spectrum.clear_molecular_formulas()

        return mass_spectrum, transient.transient_time

def run_thermo(file_location):

    MSParameters.mass_spectrum.noise_threshold_method = 'log'
    MSParameters.mass_spectrum.noise_threshold_min_s2n = 6

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

    # mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)

    mass_spectrum = parser.get_average_mass_spectrum()

    return mass_spectrum, 3

def get_masslist(file_location):

    return(ReadMassList(file_location).get_mass_spectrum(polarity=-1))

def calspec(msobj, refmasslist, order=2):

    calfn = MzDomainCalibration(msobj, refmasslist)
    ref_mass_list_fmt = calfn.load_ref_mass_list(refmasslist)

    imzmeas, mzrefs = calfn.find_calibration_points(msobj, ref_mass_list_fmt,
                                                    calib_ppm_error_threshold=(-1.0, 1.0),
                                                    calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(msobj, ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-1.5, 1.5),
                                                        calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(msobj, ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-3, 3),
                                                        calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(msobj, ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-5, 5),
                                                        calib_snr_threshold=4)

    if len(mzrefs) < 5:
        imzmeas, mzrefs = calfn.find_calibration_points(msobj, ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-7, 7),
                                                        calib_snr_threshold=4)

    if len(mzrefs) < 5:

        imzmeas, mzrefs = calfn.find_calibration_points(msobj, ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(-10, 10),
                                                        calib_snr_threshold=4)

    calfn.recalibrate_mass_spectrum(msobj, imzmeas, mzrefs, order=order)

def set_parameters(mass_spectrum, field_strength=12, pos=False):

    if field_strength == 12:

        mass_spectrum.settings.max_calib_ppm_error = 5
        mass_spectrum.settings.min_calib_ppm_error = -5

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -1
        mass_spectrum.molecular_search_settings.max_ppm_error = 1

        mass_spectrum.settings.calib_sn_threshold = 2

    elif field_strength == 15:

        mass_spectrum.settings.max_calib_ppm_error = 3
        mass_spectrum.settings.min_calib_ppm_error = -3

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -1
        mass_spectrum.molecular_search_settings.max_ppm_error = 1

        mass_spectrum.settings.calib_sn_threshold = 2

    else:

        mass_spectrum.settings.max_calib_ppm_error = 1
        mass_spectrum.settings.min_calib_ppm_error = -1

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -0.5
        mass_spectrum.molecular_search_settings.max_ppm_error = 0.5

    mass_spectrum.molecular_search_settings.url_database = None
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

def merge_files(file_paths: list, variable='Peak Height'):

    master_data_dict = []
    list_filenames = []
    for filepath in file_paths:
        
        filepath = Path(filepath)
        
        with filepath.open('r') as f:
            
            data = json.loads(json.load(f))
            
            df = DataFrame(data)
            idx = df.groupby(['Molecular Formula'])['Confidence Score'].transform(max) == df['Confidence Score']

            df = df[idx]
            df.fillna(0, inplace=True)

            name_column = "{} ({})".format(variable, filepath.stem)
            df.rename({variable: name_column}, inplace=True, axis=1)

            list_filenames.append(name_column)
            master_data_dict.extend(df.to_dict('records'))

    formula_dict = {}
    for record in master_data_dict:
        molecular_formula = record.get('Molecular Formula')
        
        if molecular_formula in formula_dict.keys():
            formula_dict[molecular_formula].append(record)
        else:
            formula_dict[molecular_formula] = [record]
    
    def dict_mean(dict_list, average_keys):
        mean_dict = {}
        
        for key in average_keys:
            
            mean_dict[key] = sum(d[key] for d in dict_list) / len(dict_list)
        
        return mean_dict
        
    average_records = []
    
    average_keys = ['m/z', 'Calibrated m/z', 'Calculated m/z', 'Peak Area', 'Resolving Power', 'S/N', 'm/z Error (ppm)', 'm/z Error Score', 
                    'Isotopologue Similarity', 'Mono Isotopic Index', 'Confidence Score']
    average_keys.extend(list_filenames)

    for formula, records in formula_dict.items():
        
        #mean_dict = dict_mean(records, average_keys)
        mean_dict = {}
        for record in  records: 
            #get the selected variable
            for filename in list_filenames:
                if filename in record.keys():
                    mean_dict[filename] = record[filename]
            
        for record in  records: 
            #than get the rest of the data
            for key in record.keys():
                if key not in average_keys:
                    mean_dict[key] = record[key]
    
        average_records.append(mean_dict)                
    master_df = pd.DataFrame(average_records)
    
    master_df.set_index('Molecular Formula', inplace=True)
    print(master_df)

    master_df.to_csv('{}.csv')
    #grouped = master_df.groupby(["Molecular Formula", "Sample Name", "Peak Height"])
    
   

def run_assignment(file_location, field_strength=12):
    #mass_spectrum = get_masslist(file_location)

    mass_spectrum, transient_time = run_bruker(file_location)
    set_parameters(mass_spectrum, field_strength=field_strength, pos=False)
    #mass_spectrum.filter_by_max_resolving_power(field_strength, transient_time)

    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

    mass_spectrum.percentile_assigned(report_error=True)
    
    mass_spectrum.to_csv(mass_spectrum.sample_name, write_metadata=False)
    
    mass_spectrum.molecular_search_settings.score_method = "prob_score"
    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"
    data_table = mass_spectrum.to_json()

    with open(mass_spectrum.sample_name + '.json', 'w') as outfile:
        json.dump(data_table, outfile)

    # export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")

    # mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)
    # mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    # mass_spectrum_by_classes.plot_mz_error()

    # mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    # plt.show()
    # mass_spectrum_by_classes.plot_mz_error()
    # plt.show()
    # mass_spectrum_by_classes.plot_ms_class()
    # plt.show()
    # dataframe = mass_spectrum_by_classes.to_dataframe()
    # return (mass_spectrum, mass_spectrum_by_classes)

    # class_plot(dataframe)

def get_all_used_atoms_in_order(mass_spectrum):

    atoms_in_order = Atoms.atoms_order
    all_used_atoms = set()
    if mass_spectrum:
        for ms_peak in mass_spectrum:
            if ms_peak:
                for m_formula in ms_peak:
                    for atom in m_formula.atoms:
                        all_used_atoms.add(atom)

    def sort_method(atom):
        return [atoms_in_order.index(atom)]

    return sorted(all_used_atoms, key=sort_method)

def export_calc_isotopologues(mass_spectrum, out_filename):

    columns_label = ["Mono Isotopic Index", "Calculated m/z", "Calculated Peak Height", 'Heteroatom Class', "Molecular Formula"]

    atoms_order_list = get_all_used_atoms_in_order(mass_spectrum)

    column_labels = columns_label + atoms_order_list

    dict_data_list = []

    for index, ms_peak in enumerate(mass_spectrum):

        if ms_peak:
            for m_formula in ms_peak:
                if not m_formula.is_isotopologue:
                    for imf in m_formula.expected_isotopologues:

                        formula_dict = imf.to_dict()
                        dict_result = {"Mono Isotopic Index": index,
                                       "Calculated m/z": imf.mz_calc,
                                       "Calculated Peak Height": imf.abundance_calc,
                                       'Heteroatom Class': imf.class_label,
                                       'H/C': imf.H_C,
                                       'O/C': imf.O_C,
                                       'Ion Type': imf.ion_type.lower(),
                                       }

                        for atom in atoms_order_list:
                            if atom in formula_dict.keys():
                                dict_result[atom] = formula_dict.get(atom)

                        dict_data_list.append(dict_result)

    df = DataFrame(dict_data_list, columns=column_labels)
    df.to_csv(out_filename + ".csv", index=False)

def monitor(target):
    ''' psutil is not installed by default, use the requirement_dev.txt to install non essential packages'''
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

def run_multiprocess():

    cores = 4
    # file_location = get_dirname()
    file_location = get_filename()
    p = Pool(cores)
    args = [(file_path) for file_path in [file_location] * 1]
    ms_collection = p.map(worker, args)
    p.close()
    p.join()

    for ms in ms_collection:
        ms[0].to_hdf('test')

if __name__ == "__main__":

    # run_multiprocess()
    # cpu_percents = monitor(target=run_multiprocess)
    # print(cpu_percents)
    file_location = get_filenames()
    if file_location:
        merge_files(file_location)
        #run_assignment(file_location)
