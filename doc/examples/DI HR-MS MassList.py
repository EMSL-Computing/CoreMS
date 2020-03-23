import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt
from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.mass_spectrum.input.massList import ReadMassList
from corems.mass_spectrum.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems import SuppressPrints

def class_plot(df):

    import seaborn as sns

    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    #sns.pairplot(df, vars=[ 'mz','abundance'],  hue="class")
    g = sns.PairGrid(df, vars=[ 'mz','abundance'], hue="class")
    g = g.map_upper(sns.kdeplot)
    g = g.map_upper(sns.kdeplot)
    g = g.map_diag(sns.kdeplot, lw=2)
    plt.show()

if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location = file_dialog.getOpenFileName()[0]
    app.quit()

    mass_spectrum = ReadMassList(file_location).get_mass_spectrum(polarity=-1)
    mass_spectrum.plot_centroid()
    
    #plt.show()

    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_ppm_error  = -0.5
    mass_spectrum.molecular_search_settings.max_ppm_error = 0.5

    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 50

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0,22)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,1)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0,1)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0,0)
    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical= False
    mass_spectrum.molecular_search_settings.isAdduct = True

    
    mass_spectrum.filter_by_max_resolving_power(15, 2)

    #plt.show()

    #with SuppressPrints():
    OxygenPriorityAssignment(mass_spectrum).run()

    mass_spectrum.percentile_assigned()

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum)

    mass_spectrum_by_classes.plot_ms_assigned_unassigned()

    plt.show()

    #dataframe = mass_spectrum_by_classes.to_dataframe()
    
    #class_plot(dataframe)

    all_classes = 0
    
    colors = ["r","blue","g","purple","black","orange",]
    classes = ["O7","O9","O12","O15","O18","O21",]
    
    color_dictionary = dict(zip(classes, colors))

    for classe in mass_spectrum_by_classes.get_classes(threshold_perc=0, isotopologue=False):
        
#    for index, classe in enumerate(classes):
        
        #plt.subplot(2, 3, index+1)
        mass_spectrum_by_classes.plot_dbe_vs_carbon_number(classe, color='g')

        #mass_spectrum_by_classes.plot_ms_class(classe, color=color_dictionary.get(classe)) 
        

    plt.show()

    print("Sum Relative Abundance = %.2f" % all_classes)