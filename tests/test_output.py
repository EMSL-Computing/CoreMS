
import sys
from pathlib import Path
sys.path.append(".")

import pytest

from corems.mass_spectra.output.export import HighResMassSpectraExport
from corems.mass_spectrum.output.export import HighResMassSpecExport
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.mass_spectra.input.boosterHDF5 import ReadHDF_BoosterMassSpectra
from corems.encapsulation.factory.parameters import MSParameters

def import_corems_mass_list():

    file_location = Path.cwd() / "tests/tests_data/ftms/" / "ESI_NEG_SRFA_COREMS.csv"

    # polarity need to be set or read from the file
    MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 1
    
    # load any type of mass list file, change the delimeter to read another type of file, i.e : "," for csv, "\t" for tabulated mass list, etc
    mass_list_reader = ReadCoremsMasslist(file_location)

    mass_spectrum = mass_list_reader.get_mass_spectrum()

    return mass_spectrum

def import_booster_mass_spectra_hdf():

    file_path = Path.cwd() / "tests/tests_data/ftms/" / "ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5"

    if file_path.exists():
        # polarity need to be set or read from the file
        booster_reader = ReadHDF_BoosterMassSpectra(file_path)

        booster_reader.start()
        booster_reader.join()

    return booster_reader.get_lcms_obj()


def test_export_mass_spectra():

    mass_spectra = import_booster_mass_spectra_hdf()

    exportMS = HighResMassSpectraExport('NEG_ESI_SRFA_CoreMS', mass_spectra)

    exportMS.get_mass_spectra_attrs(mass_spectra)
    exportMS.get_pandas_df()
    exportMS.to_pandas()
    exportMS.to_excel()
    exportMS.to_csv()
    exportMS.to_hdf()


def test_export_mass_spectrum():

    mass_spectrum = import_corems_mass_list()

    exportMS = HighResMassSpecExport('NEG_ESI_SRFA_CoreMS', mass_spectrum)

    # exportMS.to_pandas()
    # exportMS.to_excel()
    # exportMS.to_csv()
    # exportMS.to_hdf()

    exportMS._output_type = 'excel'
    exportMS.save()
    exportMS._output_type = 'csv'
    exportMS.save()
    exportMS._output_type = 'pandas'
    exportMS.save()
    exportMS._output_type = 'hdf5'
    exportMS.save()
    exportMS.get_pandas_df()    
    exportMS.to_json() 

    mass_spectrum.to_excel('NEG_ESI_SRFA_CoreMS')
    mass_spectrum.to_dataframe()

    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"
    mass_spectrum.to_csv('NEG_ESI_SRFA_CoreMS_prob_score')

    mass_spectrum.to_json()
    mass_spectrum.to_pandas('NEG_ESI_SRFA_CoreMS')

if __name__ == "__main__":
                  
    test_export_mass_spectra()
    #test_export_mass_spectrum()
