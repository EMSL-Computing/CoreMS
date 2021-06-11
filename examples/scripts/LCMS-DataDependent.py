import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt

from corems.mass_spectra.input import rawFileReader
from corems import get_filename

def run_thermo(file_location):

    parser = rawFileReader.ImportDataDependentThermoMSFileReader(file_location)

    parser.get_tic(ms_type='MS', plot=False)

    parser.get_tic(ms_type='MS2', plot=False)

    plt.show()

    # print(len(parser.selected_mzs))

    parser.get_eics(parser.selected_mzs[0:100],
                    ppm_tolerance=1,
                    plot=True)
    plt.show()
    # print(parser.get_all_filters())

if __name__ == "__main__":
    dirpath = "C:\\Users\\eber373\\Desktop\\Data\\LCMS\\RAW Files\\HILIC"
    filepath = "\\NEG\\LCMS_5191_CapDev_HILIC_Mix1_NEG_30Apr2021.raw"
    # run_multiprocess()
    # cpu_percents = monitor(target=run_multiprocess)
    # print(cpu_percents)
    # file_location = get_filename()
    file_location = dirpath + filepath
    if file_location:
        run_thermo(file_location)
