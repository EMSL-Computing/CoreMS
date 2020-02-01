import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path

import matplotlib.pyplot as plt

from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.mass_spectrum.input.massList import ReadMassList
from corems.mass_spectrum.factory.classification import HeteroatomsClassification, Labels 

from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment

def class_plot(df):

    import seaborn as sns

    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    sns.pairplot(df, vars=[ 'mz','abundance'],  hue="class")

    plt.show()

def dbe_c_number_plot(carbon_number, dbe, abundance ):

    import seaborn as sns

    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    sns.jointplot(carbon_number, dbe, kind="hex", color="#4CB391")

    plt.show()

if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location = file_dialog.getOpenFileName()[0]
    app.quit()

    mass_spectrum = ReadMassList(file_location).get_mass_spectrum(polarity=-1)
    mass_spectrum.plot_centroid()
    
    plt.show()

    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_mz_error = -0.5
    mass_spectrum.molecular_search_settings.max_mz_error = 0.5
    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (2,20)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0,0)
    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical= False
    mass_spectrum.molecular_search_settings.isAdduct = False

    OxygenPriorityAssignment(mass_spectrum).run()

    mass_spectrum.percentile_assigned()

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum)

    mass_spectrum_by_classes.plot_ms_assigned_unassigned()

    dataframe = mass_spectrum_by_classes.to_dataframe()

    
    for classe in mass_spectrum_by_classes.get_classes(1):

        if classe != Labels.unassigned:

            print(classe)

            carbon_number = mass_spectrum_by_classes.carbon_number(classe)
            dbe = mass_spectrum_by_classes.dbe(classe)
            abundance = mass_spectrum_by_classes.abundance(classe)

            dbe_c_number_plot(carbon_number, dbe, abundance)

    plt.show()