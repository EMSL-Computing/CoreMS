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
    
    # MSParameters.mass_spectrum.threshold_method = 'relative_abundance'
    # MSParameters.mass_spectrum.relative_abundance_threshold = 10

    mass_spectrum = rawFileReader.ImportLCMSThermoMSFileReader(file_location).get_summed_mass_spectrum(2,8)
    #print(mass_spectrum)
    # mass_spectrum.plot_mz_domain_profile()
    # mass_spectrum.plot_profile_and_noise_threshold()
    #print("polarity", mass_spectrum.polarity)
    # plt.show()
    '''
    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_ppm_error  = -0.4
    mass_spectrum.molecular_search_settings.max_ppm_error = 0.0

    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 30

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1,20)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,3)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,3)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['K'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Fe'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Se'] = (0,0)
    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical= False
    mass_spectrum.molecular_search_settings.isAdduct = True
    
    
    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
    
    mass_spectrum.percentile_assigned()

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=False)

    mass_spectrum_by_classes.plot_ms_assigned_unassigned()

    # plt.show()
    all_classes = 0
    for classe in mass_spectrum_by_classes.get_classes(threshold_perc=0, isotopologue=True):
        
        mass_spectrum_by_classes.plot_dbe_vs_carbon_number(classe)
        #plt.show()

        mass_spectrum_by_classes.plot_ms_class(classe)
        plt.show()

        mass_spectrum_by_classes.plot_mz_error_class(classe)

        mz_calc_l = mass_spectrum_by_classes.mz_calc(classe)

        mf_l = mass_spectrum_by_classes.molecular_formula(classe)

        print("# Name; m/z value; charge; ion formula; collision cross section [A^2]")
        

        for index, mz_calc in enumerate(mz_calc_l):

            #print("%s : %s : %s : %.6f : 1+ : %s" %(classe, mf_l[index].class_label, mf_l[index].string, mz_calc, mf_l[index].string))
            
            print("%s %.6f 1+ %s" %(mf_l[index].string.replace(' ', ''), mz_calc, mf_l[index].string.replace(' ', '')))
        
        plt.show()

    print("Sum Relative Abundance = %.2f" % all_classes)
    '''