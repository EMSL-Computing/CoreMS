import pickle, sys
from pathlib import Path
sys.path.append(".")

import pytest

from corems.transient.input.BrukerSolarix import ReadBrukerSolarix

__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

def test_import_transient():
    
    # from corems.structure.input.MidasDatFile import ReadMidasDatFile
    
    file_location = Path.cwd() / "tests/tests_data/" / "ESI_NEG_SRFA.d"

    #setting for signal processing
    apodization_method = "Hanning"
    number_of_truncations = 0
    number_of_zero_fills = 1

    bruker_reader = ReadBrukerSolarix(file_location)

    bruker_transient = bruker_reader.get_transient()

    bruker_transient.set_processing_parameter(
        apodization_method, number_of_truncations, number_of_zero_fills
    )

    mass_spec = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)

    print(mass_spec.mspeaks[0].mz_exp, mass_spec.mspeaks[-1].mz_exp)
   
    with open("test.pkl", "wb") as file:
        pickle.dump(bruker_transient, file, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    test_import_transient()
