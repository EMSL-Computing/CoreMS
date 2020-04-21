import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path
from multiprocessing import Pool
import cProfile, pstats

from numpy import array, polyfit, poly1d
from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch
from corems.mass_spectra.calc.GC_RI_Calibration import get_rt_ri_pairs

def start_sql_from_file():
    
    ref_lib_path = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"
    sql_obj = ReadNistMSI(ref_lib_path).get_sqlLite_obj()
    return sql_obj

def sql_database(file_location):
    
    sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

    min_max_rt = (18.037, 18.037)

    min_max_ri = (1637.30, 1737.30)

    sqlLite_obj.query_min_max_ri((1637.30, 1638.30))
    sqlLite_obj.query_min_max_rt((17.111, 18.111))
    sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111))

def get_gcms(file_path):
    
    reader_gcms = ReadAndiNetCDF(file_path)
	
    reader_gcms.run()
    
    gcms = reader_gcms.get_gcms_obj()

    gcms.process_chromatogram()

    return gcms

def get_reference_dict():

    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_path = file_dialog.getOpenFileName(None, "FAMES REF FILE", filter="*.cdf")[0]
    file_dialog.close()
    app.exit()
    
    if not file_path: return  None
    
    else:
        
        gcms_ref_obj = get_gcms(file_path)

        rt_ri_pairs = get_rt_ri_pairs(gcms_ref_obj)

        return rt_ri_pairs
        
def run(args):
    
    file_path, ref_dict = args
    
    gcms = get_gcms(file_path)
    
    gcms.calibrate_ri(ref_dict)
    
    # sql_obj = start_sql_from_file()
    # lowResSearch = LowResMassSpectralMatch(gcms, sql_obj=sql_obj)
    # !!!!!! READ !!!!! use the previous two lines if db/pnnl_lowres_gcms_compounds.sqlite does not exist
    # and comment the next line
    lowResSearch = LowResMassSpectralMatch(gcms)
    lowResSearch.run()

    return gcms

def calibrate_and_search(out_put_file_name, cores):
    
    import csv
    
    ref_dict = get_reference_dict()
    
    if ref_dict:
        
        file_dialog = QFileDialog()
        file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
        
        if file_dialog:
            
            file_locations = file_dialog.getOpenFileNames(None, "Standard Compounds Files", filter="*.cdf")
            file_dialog.close()
            
            # run in multiprocessing mode
            pool = Pool(cores)
            args = [(file_path, ref_dict) for file_path in file_locations[0]]
            gcmss = pool.map(run, args)
            pool.close()
            pool.join()
            for gcms in gcmss:
                
                gcms.to_csv(out_put_file_name, highest_score=False)
                #gcms.to_excel(out_put_file_name, highest_score=False)
                #gcms.to_pandas(out_put_file_name)
                
                #df = gcms.get_dataframe()
                #json_data = gcms.to_json()
                
                #print(json_data)

                #gcms.plot_processed_chromatogram()
                
                #gcms.plot_gc_peaks()

                #gcms.plot_chromatogram()

                #gcms.plot_smoothed_chromatogram()

                #gcms.plot_baseline_subtraction()

                #gcms.plot_detected_baseline()

                #matplotlib.pyplot.show()

def worker(args):

    cProfile.runctx('run(args)', globals(), locals(), 'gc-ms.prof')
      
if __name__ == '__main__':                           
    import matplotlib
    matplotlib.use('TkAgg')

    cores = 4
    out_put_file_name = 'Group 6_Standards'
    calibrate_and_search(out_put_file_name, cores)


