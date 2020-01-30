import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path

from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.mass_spectrum.input.massList import ReadMassList
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment


if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)

    file_location = file_dialog.getOpenFileName()[0]
    app.quit()

    mass_spectrum = ReadMassList(file_location).get_mass_spectrum(polarity=-1)
    print(mass_spectrum)

    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_mz_error = -1
    mass_spectrum.molecular_search_settings.max_mz_error = 1

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0,20)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,3)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0,0)

    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical= False
    mass_spectrum.molecular_search_settings.isAdduct = True

    OxygenPriorityAssignment(mass_spectrum).run()

    mass_spectrum.percentile_assigned()