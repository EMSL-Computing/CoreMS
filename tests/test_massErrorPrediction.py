import sys 
sys.path.append('.')

from pathlib import Path
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.mass_spectrum.calc.MassErrorPrediction import MassErrorPrediction
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from test_molecularFormulaSearch import create_mass_spectrum


def test_error_prediction():

    mass_spectrum = create_mass_spectrum()

    mass_error_prediction = MassErrorPrediction(mass_spectrum)
    
    mass_error_prediction.get_results()

    return mass_spectrum

if __name__ == "__main__":

    mass_spectrum = test_error_prediction()
    
    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_ppm_error  = -1
    mass_spectrum.molecular_search_settings.max_ppm_error = 0

    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 50

    mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1,90)
    mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4,200)
    mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0,22)
    mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0,0)
    mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0,0)
    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical= False
    mass_spectrum.molecular_search_settings.isAdduct = True

    #plt.show()

    #with SuppressPrints():
    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

    mass_spectrum.percentile_assigned()

    for mspeak in mass_spectrum:

         for mf in mspeak:
            mf.confidence_score
