import sys
sys.path.append("./ext_lib")
sys.path.append(".")
from corems.encapsulation.settings.input import InputSetting
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.encapsulation.constant import Labels

import clr
clr.AddReference("ChemstationMSFileReader")
import ChemstationMSFileReader 

file_loc = 'tests/tests_data/DATA.MS'
clsChemstation = ChemstationMSFileReader.clsChemstationDataMSFileReader(file_loc)

header = clsChemstation.ReadHeaders(file_loc)

clsChemstation.GetSpectrum()

#System.Int32,ChemstationMSFileReader.clsSpectralRecord@,System.Int32@)