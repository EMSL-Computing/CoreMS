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

from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
from corems.encapsulation.constant import Atoms
from corems.encapsulation.factory.parameters import MSParameters


from doc.examples.nmdc.NMDC_Metadata import DMS_Mapping, NMDC_Metadata

def run_bruker(file_location):

    with ReadBrukerSolarix(file_location) as transient:

        MSParameters.mass_spectrum.threshold_method = 'auto'
        MSParameters.mass_spectrum.s2n_threshold = 6
        
        mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)
        #mass_spectrum.plot_profile_and_noise_threshold()
        #plt.show()
        # find_formula_thread = FindOxygenPeaks(mass_spectrum)
        # find_formula_thread.run()

        # mspeaks_results = find_formula_thread.get_list_found_peaks()
        # calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
        # calibrate.ledford_calibration()

        # mass_spectrum.clear_molecular_formulas()

        return mass_spectrum, transient.transient_time

def run_thermo(file_location):

    MSParameters.mass_spectrum.threshold_method = 'auto'
    MSParameters.mass_spectrum.s2n_threshold = 6

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

    #mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)

    mass_spectrum = parser.get_average_mass_spectrum_in_scan_range()

    return mass_spectrum, 3

def get_masslist(file_location):

    return(ReadMassList(file_location).get_mass_spectrum(polarity=-1))

def calspec(msobj,refmasslist,order=2):
    
    calfn = MzDomainCalibration(msobj, refmasslist)
    ref_mass_list_fmt = calfn.load_ref_mass_list(refmasslist)
    
    imzmeas, mzrefs = calfn.find_calibration_points(msobj,ref_mass_list_fmt,
                                  calib_ppm_error_threshold=(-1.0,1.0),
                                  calib_snr_threshold=4)

    
    if len(mzrefs)<5:
            imzmeas, mzrefs = calfn.find_calibration_points(msobj,ref_mass_list_fmt,
                                  calib_ppm_error_threshold=(-1.5,1.5),
                                  calib_snr_threshold=4)
    

    if len(mzrefs)<5:
            imzmeas, mzrefs = calfn.find_calibration_points(msobj,ref_mass_list_fmt,
                                  calib_ppm_error_threshold=(-3,3),
                                  calib_snr_threshold=4)
    
    if len(mzrefs)<5:
            imzmeas, mzrefs = calfn.find_calibration_points(msobj,ref_mass_list_fmt,
                                  calib_ppm_error_threshold=(-5,5),
                                  calib_snr_threshold=4)    
    
    if len(mzrefs)<5:
            
            imzmeas, mzrefs = calfn.find_calibration_points(msobj,ref_mass_list_fmt,
                                  calib_ppm_error_threshold=(-7,7),
                                  calib_snr_threshold=4) 

    if len(mzrefs)<5:
            
            imzmeas, mzrefs = calfn.find_calibration_points(msobj,ref_mass_list_fmt,
                                  calib_ppm_error_threshold=(-10,10),
                                  calib_snr_threshold=4)

    calfn.recalibrate_mass_spectrum(msobj, imzmeas, mzrefs,order=order)

def set_parameters(mass_spectrum, field_strength=12, pos=False):

    if field_strength == 12:

        mass_spectrum.settings.max_calib_ppm_error = 10
        mass_spectrum.settings.min_calib_ppm_error = -10

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

    # mass_spectrum.molecular_search_settings.url_database = "postgres://coremsdb:coremsmolform@localhost:5432/molformula"
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
    file_location, ref_calibration_file, field_strength = args

    if field_strength == 21:

        #return "21T", None
        #print("{}   {}".format("21T", file_location))
        print("{} {}  {}".format("processing", field_strength,  file_location))
        mass_spectrum, transient_time = run_thermo(file_location)
        #TODO add_thermo_processing_here
    
    else:    
        
        return "not 21T", None
        print("{} {}  {}".format("processing", field_strength,  file_location))
        mass_spectrum, transient_time = run_bruker(file_location)
        # return "not 21T", None
        
    if mass_spectrum.polarity > 0:
        
        print("{}   {}".format("positive mode", file_location))
        is_pos = True
        return "positive mode", None

    else:
        
        is_pos = False

    if len(mass_spectrum) < 30:

        print("{}   {}".format("too few peaks", file_location))
        return "too few peaks", None

    
    
    set_parameters(mass_spectrum, field_strength=field_strength, pos=is_pos)

    if ref_calibration_file:

        calspec(mass_spectrum, ref_calibration_file)
        #MzDomainCalibration(mass_spectrum, ref_calibration_file).run()

    #mass_spectrum.filter_by_max_resolving_power(field_strength, transient_time)

    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
    mass_spectrum.percentile_assigned(report_error=True)
    mass_spectrum.molecular_search_settings.score_method = "prob_score"
    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

    return "all_good", mass_spectrum

    # export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")

    # mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)
    # mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    # mass_spectrum_by_classes.plot_mz_error()

    plt.show()
    # dataframe = mass_spectrum_by_classes.to_dataframe()
    # return (mass_spectrum, mass_spectrum_by_classes)

    # class_plot(dataframe)

def run_assignment(file_location, field_strength=12):
    # mass_spectrum = get_masslist(file_location)

    mass_spectrum, transient_time = run_bruker(file_location)
    set_parameters(mass_spectrum, field_strength=field_strength, pos=False)
    mass_spectrum.filter_by_max_resolving_power(field_strength, transient_time)

    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
    mass_spectrum.percentile_assigned(report_error=True)
    mass_spectrum.molecular_search_settings.score_method = "prob_score"
    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

    mass_spectrum.to_csv(mass_spectrum.sample_name, write_metadata=False)
    # export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")

    # mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)
    # mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    # mass_spectrum_by_classes.plot_mz_error()

    plt.show()
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

    columns_label = ["Mono Isotopic Index", "Calculated m/z", "Calculated Peak Height", 'Heteroatom Class', "Molecular Formula" ]

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

    import psutil, time

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


def run_nmdc_data_processing():

    data_dir = Path("C:\\Users\\eber373\\Desktop\\nmdc ftms pilot data")
    dms_file_path = Path("db/NOM Data to Process.xlsx")
    
    results_dir = Path("results")
    registration_path = results_dir / "ftms_nom_data_products.json"
    failed_files = results_dir / "nom_failed_files.json"
    pos_files = results_dir / "pos_files.json"

    field_strength = 12
    cores = 4
    ref_calibration_path = False

    #file_paths = get_dirnames()
    ref_calibration_path = Path("db/Hawkes_neg.ref")

    # args = [(file_path, ref_calibration_path, field_strength) for file_path in file_paths]
    #print(args)
    #p = Pool(cores)
    #ms_collection = [run_nmdc_workflow(args[0])]
    # ms_collection = p.map(run_nmdc_workflow, args)
    dms_mapping = DMS_Mapping(dms_file_path)
    selected_files = dms_mapping.get_selected_sample_list()
    
    failed_list = []
    pos_list = []
    
    for file_name, field_strength in selected_files:
        
        if str(file_name)[0:6] == 'Brodie':
            
            in_file_path = data_dir / 'Brodie' / file_name / file_name.with_suffix('.raw')    
        
        else:
            
            in_file_path = data_dir / 'Stegen' / file_name / file_name.with_suffix('.raw')    

        try:
            
            issue, ms = run_nmdc_workflow((in_file_path, ref_calibration_path, field_strength))
            
            if ms:
                
                file_name = Path(in_file_path.name)

                output_file_path = results_dir / file_name.with_suffix('.csv')
                ms.to_csv(output_file_path, write_metadata=False)

                nmdc = NMDC_Metadata(in_file_path, ref_calibration_path, output_file_path, dms_file_path)
                nmdc.create_nmdc_ftms_metadata(ms, registration_path)
                
            else:

                if issue == 'positive mode':
                    pos_list.append(str(in_file_path))
                else:
                    pass
        except Exception as inst:
            
            print(type(inst))    # the exception instance
            print(inst.args)     # arguments stored in .args
            print(inst)  
            failed_list.append(str(in_file_path))
            
    with failed_files.open('w') as json_file:

        json_file.write(json.dumps(failed_list, indent=1))

    with pos_files.open('w') as json_file:

        json_file.write(json.dumps(pos_list, indent=1))    

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
    run_nmdc_data_processing()
    #file_location = get_filename()
    #if file_location:
    #    run_assignment(file_location)
