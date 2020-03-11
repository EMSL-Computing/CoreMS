import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("C:\\Users\\eber373\\Desenvolvimento\\Projects-Python\\CoreMS")

from pathlib import Path

from matplotlib import pyplot
from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch

def sql_database(file_location):
    
    print(file_location)

    sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

    min_max_rt = (18.037, 19.037)

    min_max_ri = (1637.30, 1737.30)

    sqlLite_obj.query_min_max_ri((1637.30, 1638.30)) 
    sqlLite_obj.query_min_max_rt((17.111, 18.111))            
    sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111)) 


def andi_netcdf_gcms():

    file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"

    ref_file_path = Path.cwd() / "tests/tests_data/gcms/" / "FAMES_REF.MSL"

    reader_gcms = ReadAndiNetCDF(file_path)
	
    reader_gcms.run()
    
    gcms = reader_gcms.get_gcms_obj()

    #gcms.plot_chromatogram()

    #gcms.plot_smoothed_chromatogram()

    gcms.process_chromatogram()

    gcms.plot_processed_chromatogram()

    gcms.plot_gc_peaks()

    #gcms.plot_baseline_subtraction()

    #gcms.plot_detected_baseline()
    lowResSearch = LowResMassSpectralMatch(gcms, ref_file_path)

    lowResSearch.run()

    

#app = QApplication(sys.argv)
#file_dialog = QFileDialog()
#file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
#file_location = file_dialog.getOpenFileName()[0]
#app.quit()

andi_netcdf_gcms()
pyplot.show()
