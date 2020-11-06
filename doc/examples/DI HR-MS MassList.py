import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path
import cProfile
import pstats
from multiprocessing import Pool, Process

from pandas import DataFrame
from matplotlib import pyplot as plt

from corems.mass_spectrum.input.massList import ReadMassList
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems import SuppressPrints, get_filename, get_dirname
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.mass_spectrum.calc.CalibrationCalc import FreqDomain_Calibration
from corems.encapsulation.constant import Atoms

def run_bruker(file_location):
    
    with ReadBrukerSolarix(file_location) as transient:

        mass_spectrum = transient.get_mass_spectrum(plot_result=False, auto_process=True)

        # find_formula_thread = FindOxygenPeaks(mass_spectrum)
        # find_formula_thread.run()
        
        # mspeaks_results = find_formula_thread.get_list_found_peaks()
        # calibrate = FreqDomain_Calibration(mass_spectrum, mspeaks_results)
        # calibrate.ledford_calibration()
        
        # mass_spectrum.clear_molecular_formulas()

        return mass_spectrum
        
def get_masslist(file_location):

    return(ReadMassList(file_location).get_mass_spectrum(polarity=-1))

def run_assignment(file_location):
    
    #mass_spectrum = run_bruker(file_location)
    mass_spectrum = get_masslist(file_location)

    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_ppm_error  = -5
    mass_spectrum.molecular_search_settings.max_ppm_error = 5

    # mass_spectrum.molecular_search_settings.url_database = "postgres://coremsdb:coremsmolform@localhost:5432/molformula"
    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 50

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,100)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1,22)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0,0)
    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical= False
    mass_spectrum.molecular_search_settings.isAdduct = False
    
    #mass_spectrum.filter_by_max_resolving_power(15, 2)
    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
    mass_spectrum.percentile_assigned(report_error=True)
    mass_spectrum.molecular_search_settings.score_method = "prob_score"
    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"
    
    mass_spectrum.to_csv("15T_Neg_ESI_SRFA")
    
    #export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")
    
    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)
    #mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    mass_spectrum_by_classes.plot_mz_error()

    plt.show()
    # dataframe = mass_spectrum_by_classes.to_dataframe()
    return (mass_spectrum, mass_spectrum_by_classes)
    
    #class_plot(dataframe)

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
                        dict_result = { "Mono Isotopic Index" : index, 
                                        "Calculated m/z": imf.mz_calc,
                                        "Calculated Peak Height" : imf.abundance_calc,
                                        'Heteroatom Class':  imf.class_label,
                                        'H/C':  imf.H_C,
                                        'O/C':  imf.O_C,
                                        'Ion Type': imf.ion_type.lower(),  
                                        }
                        
                        for atom in atoms_order_list:
                            if atom in formula_dict.keys():
                                dict_result[atom] = formula_dict.get(atom)
                        
                        dict_data_list.append(dict_result)  
    
    df = DataFrame(dict_data_list, columns=column_labels)
    df.to_csv(out_filename+".csv", index=False)

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
    #stats = pstats.Stats("topics.prof")
    #stats.strip_dirs().sort_stats("time").print_stats() 

def run_multiprocess():
    
    cores = 4
    #file_location = get_dirname()
    file_location = get_filename()
    p = Pool(cores)
    args = [(file_path) for file_path in [file_location]*1]
    ms_collection = p.map(worker, args)
    p.close()
    p.join()
    
    for ms in ms_collection:
        ms[0].to_hdf('test')

if __name__ == "__main__":

    #run_multiprocess()
    #cpu_percents = monitor(target=run_multiprocess)
    #print(cpu_percents)
    file_location = get_filename()
    if file_location:
        run_assignment(file_location)
    
