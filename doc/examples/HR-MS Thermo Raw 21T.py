import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt
from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location = file_dialog.getOpenFileName()[0]
    app.quit()
    
    # change parameters here
    MSParameters.mass_spectrum.threshold_method = 'relative_abundance'
    # MSParameters.mass_spectrum.relative_abundance_threshold = 10

    # creates the parser obj
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)
    
    # sums all the mass spectra
    mass_spectrum = parser.get_average_mass_spectrum_in_scan_range()
    
    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum_in_scan_range(first_scan=1, last_scan=5)
    
    scans_list = [1,4,6,9]
    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans_list)
    
    mass_spectrum.plot_mz_domain_profile()
    mass_spectrum.plot_profile_and_noise_threshold()
    
    #print("polarity", mass_spectrum.polarity)
    plt.show()
   