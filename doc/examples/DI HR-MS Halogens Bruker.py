import sys
sys.path.append("./")

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import matplotlib.pyplot as plt
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QApplication, QFileDialog

from corems import SuppressPrints
from corems.mass_spectrum.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


def get_mass_spectrum(file_location):

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    return bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location = file_dialog.getExistingDirectory()
    if file_dialog:

        mass_spectrum = get_mass_spectrum(file_location)

        mass_spectrum.molecular_search_settings.min_ppm_error  = 0
        mass_spectrum.molecular_search_settings.max_ppm_error = 2

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (5,5)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,0)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (1,1)
        mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (2,2)

        #enables the search for Br compounds 
        mass_spectrum.molecular_search_settings.adduct_atoms_neg = ['Cl', 'F']
    
        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical= False
        mass_spectrum.molecular_search_settings.isAdduct = False

        mass_spectrum.molecular_search_settings.use_min_peaks_filter = False

        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
        
        ax = mass_spectrum.plot_mz_domain_profile()
        
        for mspeak in mass_spectrum:
            
            if mspeak: 
                
                for mf in mspeak: 
                
                    ax.plot(mspeak.mz_exp,mspeak.abundance, color='black', linewidth=0, marker='v')
                
                    print( mf.to_string, mf.mz_calc, mf.mz_error) 
                    
                    ax.annotate(mf.to_string, (mspeak.mz_exp,mspeak.abundance))

                    for imf in mf.expected_isotopologues: 
                        
                        print(imf.to_string, imf.mz_calc, imf.prob_ratio)

                        ax.plot(imf.mz_calc, imf.abundance_calc, color='red', linewidth=0, marker='v')

                        ax.annotate(imf.to_string, (imf.mz_calc, imf.abundance_calc))
                

        plt.show()
        
        # Bromothymol Blue – C27H28Br2O5S
        # Bromocresol Green – C21H14Br4O5S
        # Chlorophenol Red – C19H12Cl2O5S
        # Bromphenol Blue – C19H9Br4O5SNa
        
        