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

    parser.chromatogram_settings.start_scan = -1
    parser.chromatogram_settings.end_scan = -1

    parser.chromatogram_settings.smooth_window = 5
    
    parser.chromatogram_settings.min_peak_datapoints = 3
    parser.chromatogram_settings.peak_height_min_percent = 0.1

    parser.chromatogram_settings.eic_signal_threshold = 0.1
    parser.chromatogram_settings.eic_tolerance_ppm = 5
    parser.chromatogram_settings.enforce_target_ms2 = True

    tic_data, ax_tic = parser.get_tic(ms_type='MS', peak_detection=True, 
                                      smooth=True, plot=False)

    ms2_tic, ax_ms2_tic = parser.get_tic(ms_type='MS2', peak_detection=False, plot=False)

    #print(data)

    # get selected data dependent mzs 
    target_mzs = parser.selected_mzs

    eics_data, ax_eic = parser.get_eics(target_mzs,
                                        tic_data,
                                        smooth=True,
                                        plot=True,
                                        legend=False,
                                        peak_detection=True,
                                        ax=ax_tic)
    
    ax_eic.plot(tic_data['Time'], tic_data['TIC'], c='black')
    #print(parser.get_all_filters())

    plt.show()
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
