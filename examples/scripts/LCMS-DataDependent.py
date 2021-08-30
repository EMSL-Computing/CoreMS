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
    
    print(file_location)
    
    parser = rawFileReader.ImportDataDependentThermoMSFileReader(file_location)

    tic_data, ax_tic = parser.get_tic(ms_type='MS', plot=True)

    ms2_tic, ax_ms2_tic = parser.get_tic(ms_type='MS2', plot=False)

    centroid_tic = parser.centroid_detector(tic_data['TIC'], tic_data['Time'])
    
    for i in centroid_tic:
        print(i)
    #print(data)

    # get selected data dependent mzs 
    target_mzs = parser.selected_mzs

    data = parser.get_eics(target_mzs[0:1000],
                    ppm_tolerance=1,
                    plot=True,
                    ax=ax_tic)
    
    plt.show()
    #print(parser.get_all_filters())

if __name__ == "__main__":
    
    dirpath = "/Users/eber373/OneDrive - PNNL/Documents/Data/LCMS/RAW Files/C18/1st run/NEG/"
    
    filepath = "LCMS_5191_CapDev_C18_Mix1_NEG_28Apr2021.raw"
    # run_multiprocess()
    # cpu_percents = monitor(target=run_multiprocess)
    # print(cpu_percents)
    # file_location = get_filename()
    file_location = dirpath + filepath
    if file_location:
        run_thermo(file_location)
