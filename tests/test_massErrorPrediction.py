import sys 
sys.path.append('.')

from pathlib import Path
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.mass_spectrum.calc.MassErrorPrediction import MassErrorPrediction
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from test_molecularFormulaSearch import create_mass_spectrum


def x_test_error_prediction():
    'This function will be removed in CoreMS 2.0. adding x to skip test'
    mass_spectrum = create_mass_spectrum()

    mass_error_prediction = MassErrorPrediction(mass_spectrum)
    
    mass_error_prediction.get_results()

    return mass_spectrum

if __name__ == "__main__":
    x_test_error_prediction()
