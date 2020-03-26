import sys
sys.path.append("./")

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path

import matplotlib.pyplot as plt
from numpy import unique

from corems import SuppressPrints, get_dirnames
from corems.mass_spectrum.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


def get_mass_spectrum(file_location):

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    return bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

def export(list_dict, export_name):
    
    import pandas as pd
    
    columns = [ "Sample Name", 
                "Calculated m/z",
                "Experimental m/z",
                "Error m/z", 
                "Peak Height",
                "Signal Noise Ratio",
                "Peak Datapoint Count",
                "Resolving Power",
                "Molecular Formula",
                "Peak Area",
                "Area Ratio Error",
                "Probability Ratio",
                "Abundance Error"
               ]
    
    df = pd.DataFrame(list_dict, columns=columns)
    df.to_pickle(export_name+".pkl")
    df.to_csv(export_name+".csv", index=False)


def set_settings_for_bromothymol_blue(mass_spectrum):

    mass_spectrum.molform_search_settings.min_ppm_error  = -2
    mass_spectrum.molform_search_settings.max_ppm_error = 5

    mass_spectrum.molform_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molform_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molform_search_settings.usedAtoms['O'] = (5,5)
    mass_spectrum.molform_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molform_search_settings.usedAtoms['S'] = (1,1)
    mass_spectrum.molform_search_settings.usedAtoms['Br'] = (2,2)

    #enables the search for Br compounds 
    mass_spectrum.molform_search_settings.adduct_atoms_neg = ['Cl', 'F']

    mass_spectrum.molform_search_settings.isProtonated = True
    mass_spectrum.molform_search_settings.isRadical= False
    mass_spectrum.molform_search_settings.isAdduct = False

    mass_spectrum.molform_search_settings.use_min_peaks_filter = False

def set_settings_for_chlorophenol_red(mass_spectrum):
    #Chlorophenol Red –  C19H12Cl2O5S -  421.978241

    mass_spectrum.molform_search_settings.min_ppm_error  = -2
    mass_spectrum.molform_search_settings.max_ppm_error = 5

    mass_spectrum.molform_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molform_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molform_search_settings.usedAtoms['O'] = (5,5)
    mass_spectrum.molform_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molform_search_settings.usedAtoms['S'] = (1,1)
    mass_spectrum.molform_search_settings.usedAtoms['Cl'] = (2,2)

    #enables the search for Cl compounds 
    mass_spectrum.molform_search_settings.adduct_atoms_neg = ['Br', 'F']

    mass_spectrum.molform_search_settings.isProtonated = True
    mass_spectrum.molform_search_settings.isRadical= False
    mass_spectrum.molform_search_settings.isAdduct = False

    mass_spectrum.molform_search_settings.use_min_peaks_filter = False

if __name__ == "__main__":
    
    list_dict = []
    
    dirnames = get_dirnames()
    
    if dirnames:
        
        for file_location in dirnames:
            
            print(file_location)
            
            mass_spectrum = get_mass_spectrum(file_location)

            set_settings_for_bromothymol_blue(mass_spectrum)
            #set_settings_for_chlorophenol_red(mass_spectrum)

            SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
            
            ax = mass_spectrum.plot_mz_domain_profile()
            #plt.show()
            
            for mspeak in mass_spectrum:
                
                if mspeak:
                    
                    for mf in mspeak: 
                    
                        ax.plot(mspeak.mz_exp,mspeak.abundance, color='black', linewidth=0, marker='s', markersize = 8, label='Assigned')
                    
                        if mf.is_isotopologue:
                            
                            ax.annotate(round(mf.area_error,1), (mspeak.mz_exp+0.2, mspeak.abundance), label="Error (%)")
                            list_dict.append( {
                                
                                "Sample Name" : mass_spectrum.sample_name, 
                                "Resolving Power" : mspeak.resolving_power, 
                                "Peak Height" : mspeak.abundance, 
                                "Signal Noise Ratio" : mspeak.signal_to_noise,
                                "Calculated m/z" : mf.mz_calc,
                                "Experimental m/z" : mspeak.mz_exp, 
                                "Error m/z" : mf.mz_error, 
                                "Peak Area": mspeak.area, 
                                "Molecular Formula" : mf.to_string,
                                "Probability Ratio": mf.prob_ratio,
                                "Area Ratio Error": mf.area_error,
                                "Abundance Error": mf.abundance_error,
                                "Peak Datapoint Count" : mspeak.final_index-mspeak.start_index
                            })
                        
                        else:
                            
                            print( mf.to_string, mf.mz_calc, mf.mz_error) 
                            
                            list_dict.append( {
                                
                                "Sample Name" : mass_spectrum.sample_name, 
                                "Resolving Power" : mspeak.resolving_power,
                                "Peak Height" : mspeak.abundance, 
                                "Signal Noise Ratio" : mspeak.signal_to_noise,
                                "Calculated m/z" : mf.mz_calc,
                                "Experimental m/z" : mspeak.mz_exp, 
                                "Error m/z" : mf.mz_error, 
                                "Peak Area": mspeak.area,
                                "Molecular Formula" : mf.to_string,
                                "Peak Datapoint Count" : mspeak.final_index-mspeak.start_index
                            })
                        #ax.annotate(mspeak.mz_exp, (mspeak.mz_exp, mspeak.abundance))

                        for imf in mf.expected_isotopologues: 
                            
                            ax.plot(imf.mz_calc, imf.abundance_calc, color='red', linewidth=0, marker='v', markersize = 8, label='Calculated')

                            #ax.annotate(imf.to_string, (imf.mz_calc, imf.abundance_calc))
            
            hand, labl = ax.get_legend_handles_labels()
            handout=[]
            lablout=[]
            for h,l in zip(hand,labl):
                if l not in lablout:
                        lablout.append(l)
                        handout.append(h)
            
            plt.legend(handout, lablout)
            plt.xlim(620,630)
            plt.show()
            #plt.savefig(mass_spectrum.sample_name+"_area"+".png")
            plt.cla()
            
        export(list_dict, "Chlorophenol Red C19 H12 Cl2 O5 S1")
        #export(list_dict, "Chlorophenol Red C19 H12 Cl2 O5 S1")    
        
        # Bromothymol Blue –  C27H28Br2O5S -  622.002380
        # Bromocresol Green – C21H14Br4O5S - 693.729492
        # Chlorophenol Red –  C19H12Cl2O5S -  421.978241
        # Bromphenol Blue –   C19H10Br4O5S(Na) -  665.698242
        
        