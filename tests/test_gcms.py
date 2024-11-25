
__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"
import warnings
warnings.filterwarnings("ignore")

import sys

import pytest
from pathlib import Path

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch
from corems.mass_spectra.calc.GC_RI_Calibration import get_rt_ri_pairs

def start_sql_from_file():
    
    ref_lib_path = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"
    sql_obj = ReadNistMSI(ref_lib_path).get_sqlLite_obj()
    return sql_obj

def test_sql_database():
    
    file_location = Path.cwd() / "tests/tests_data/gcms/" / "FAMES_REF.MSL"

    sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

    min_max_rt = (18.037, 19.037)

    min_max_ri = (1637.30, 1737.30)

    sqlLite_obj.query_min_max_ri((1637.30, 1638.30))
    sqlLite_obj.query_min_max_rt((17.111, 18.111))
    sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111))

def get_gcms(file_path):
    
    reader_gcms = ReadAndiNetCDF(file_path)
	
    reader_gcms.run()
    
    gcms = reader_gcms.get_gcms_obj()

    #gcms.process_chromatogram()

    return gcms

def get_reference_dict():

    file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"
    
    gcms = get_gcms(file_path)

    sql_obj = start_sql_from_file()

    if not file_path: return  None
    
    else:
        
        gcms_ref_obj = get_gcms(file_path)

        rt_ri_pairs = get_rt_ri_pairs(gcms_ref_obj, sql_obj)

        return rt_ri_pairs, file_path

def run(args):
    
    file_path, ref_dict, cal_file_path = args
    
    gcms = get_gcms(file_path)
    
    gcms.calibrate_ri(ref_dict, cal_file_path)
    
    gcms.peaks_rt_tic()

    gcms.peaks_rt_tic(json_string=True)

    sql_obj = start_sql_from_file()
    lowResSearch = LowResMassSpectralMatch(gcms, sql_obj=sql_obj)
    # !!!!!! READ !!!!! use the previous two lines if db/pnnl_lowres_gcms_compounds.sqlite does not exist
    # and comment the next line
    #lowResSearch = LowResMassSpectralMatch(gcms)
    lowResSearch.run()

    return gcms

def calibrate_and_search(out_put_file_name):
    
    import csv
    
    ref_dict, cal_file_path = get_reference_dict()
    
    if ref_dict:
        
        file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"
        gcms = run((file_path, ref_dict, cal_file_path))
        
        gcms.to_csv(out_put_file_name)
        gcms.to_excel(out_put_file_name)
        gcms.to_pandas(out_put_file_name)
        
        df = gcms.to_dataframe()
        json_data = gcms.to_json()
        
        #print(json_data)

        gcms.plot_processed_chromatogram()

        gcms.plot_gc_peaks()

        gcms.plot_chromatogram()

        gcms.plot_smoothed_chromatogram()

        gcms.plot_baseline_subtraction()

        gcms.plot_detected_baseline()

        gcms.chromatogram_settings.use_deconvolution = False
        gcms.process_chromatogram()
        

def test_run_gcms_pipeline():

    out_put_file_name = 'test_gcms'
    calibrate_and_search(out_put_file_name)

if __name__ == "__main__":
    #test_sql_database()
    test_run_gcms_pipeline()