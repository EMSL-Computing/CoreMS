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
