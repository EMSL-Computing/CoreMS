import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt


from corems.mass_spectrum.input.massList import ReadMassList
from corems.mass_spectrum.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems import SuppressPrints, get_filename


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
    
    file_location = get_filename()
    
    mass_spectrum = ReadMassList(file_location).get_mass_spectrum(polarity=-1)
    mass_spectrum.plot_centroid()
    
    #plt.show()

    mass_spectrum.molform_search_settings.error_method = 'None'
    mass_spectrum.molform_search_settings.min_ppm_error  = -1
    mass_spectrum.molform_search_settings.max_ppm_error = 1

    mass_spectrum.molform_search_settings.min_dbe = 0
    mass_spectrum.molform_search_settings.max_dbe = 50

    mass_spectrum.molform_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molform_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molform_search_settings.usedAtoms['O'] = (0,22)
    mass_spectrum.molform_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molform_search_settings.usedAtoms['S'] = (0,0)
    mass_spectrum.molform_search_settings.usedAtoms['Cl'] = (0,0)
    mass_spectrum.molform_search_settings.usedAtoms['P'] = (0,0)
    mass_spectrum.molform_search_settings.usedAtoms['Na'] = (0,0)
    mass_spectrum.molform_search_settings.isProtonated = True
    mass_spectrum.molform_search_settings.isRadical= False
    mass_spectrum.molform_search_settings.isAdduct = False
    
    mass_spectrum.filter_by_max_resolving_power(15, 2)

    #plt.show()

    #with SuppressPrints():
    SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
    #OxygenPriorityAssignment(mass_spectrum).run()
    mass_spectrum.percentile_assigned()

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum)
    #plt.show()
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
        mass_spectrum_by_classes.plot_dbe_vs_carbon_number(classe, color='PuBu_r')

        #mass_spectrum_by_classes.plot_ms_class(classe, color=color_dictionary.get(classe)) 
        

    plt.show()

    print("Sum Relative Abundance = %.2f" % all_classes)