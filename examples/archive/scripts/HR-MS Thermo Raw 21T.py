__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

def run_thermo(file_location):

    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 0.1

    # MSParameters.mass_spectrum.noise_threshold_method = 'log'
    # MSParameters.mass_spectrum.noise_threshold_min_s2n = 6

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

    parser.chromatogram_settings.scans = (-1, -1)
    
    tic_data, ax = parser.get_tic(ms_type='MS', plot=True)

    plt.show()

    print(parser.get_all_filters())

    transient_time_list = parser.get_icr_transient_times()

    print(transient_time_list)

    # sums all the mass spectra
    mass_spectrum = parser.get_average_mass_spectrum()

    # sums scans in selected range
    parser.chromatogram_settings.scans = (1, 10)
    
    mass_spectrum = parser.get_average_mass_spectrum()
    
    parser.chromatogram_settings.scans = [1]
    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum()

    return mass_spectrum

def run_assignment(file_location):

    # mass_spectrum = run_bruker(file_location)
    # mass_spectrum = get_masslist(file_location)
    mass_spectrum = run_thermo(file_location)

    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_ppm_error = -5
    mass_spectrum.molecular_search_settings.max_ppm_error = 5

    mass_spectrum.molecular_search_settings.url_database = None
    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 50

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 100)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 30)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 0)

    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical = False
    mass_spectrum.molecular_search_settings.isAdduct = False

    # mass_spectrum.filter_by_max_resolving_power(15, 2)
    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

    mass_spectrum.percentile_assigned(report_error=True)
    mass_spectrum.molecular_search_settings.score_method = "prob_score"
    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

    # export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)

    mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    plt.show()
    mass_spectrum_by_classes.plot_mz_error()
    plt.show()
    mass_spectrum_by_classes.plot_ms_class("O2")
    plt.show()
    # dataframe = mass_spectrum_by_classes.to_dataframe()
    return mass_spectrum

    # class_plot(dataframe)

if __name__ == "__main__":

    # app = QApplication(sys.argv)
    # file_dialog = QFileDialog()
    # file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    # file_location = file_dialog.getOpenFileName()[0]
    # app.quit()

    file_location = "tests/tests_data/ftms/Exploris_SRFA_Example.raw"
    # change parameters here

    mass_spectrum = run_thermo(file_location)

    print("acquisition_time:",  mass_spectrum.acquisition_time )
    ax = mass_spectrum.plot_mz_domain_profile()

    #for mspeak in mass_spectrum:
    #    mspeak.plot(ax=ax)
    #plt.show()
    #plt.savefig("test.png")

    #mass_spectrum.plot_profile_and_noise_threshold()
    #plt.show()
    #plt.savefig("test.png")

    mass_spectrum = run_assignment(file_location)

    mass_spectrum.to_csv(mass_spectrum.sample_name)
    # print("polarity", mass_spectrum.polarity)
