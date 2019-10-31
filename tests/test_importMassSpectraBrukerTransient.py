
import sys
from pathlib import Path
sys.path.append(".")

from corems.transient.input.BrukerSolarix import ReadBrukerSolarix
from corems.mass_spectra.input.brukerSolarix import ReadBruker_SolarixTransientMassSpectra


def test_import_lcms_from_transient():

    file_location = Path.cwd() / "tests/tests_data/" / "NEG_ESI_SRFA_Auto.d"#"SOM_LC_PeatMix_2p8_0p6_2_30AUG19_GIMLI_ZORBAX-1186_1_01_259.d"

    read_lcms = ReadBruker_SolarixTransientMassSpectra(file_location)

    read_lcms.start()
    read_lcms.join()

if __name__ == "__main__":
    
    #test_import_lcms_from_transient()
    import timeit
    t = timeit.Timer("test_import_lcms_from_transient()", setup="from __main__ import test_import_lcms_from_transient")

    print(t.timeit(1))

    #timeit.timeit("test_import_lcms_from_transient()", setup="from __main__ import test_import_lcms_from_transient"))
