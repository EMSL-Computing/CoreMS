
__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append(".")

import pytest
from pathlib import Path

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch

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

    gcms.process_chromatogram()

    return gcms

def get_reference_dict():

    file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"
    
        
    ref_file_path = Path.cwd() / "tests/tests_data/gcms/" / "FAMES_REF.MSL"

    gcms = get_gcms(file_path)

    lowResSearch = LowResMassSpectralMatch(gcms, ref_file_path, calibration=True)

    lowResSearch.run()

    dict_ri_rt = {}

    list_of_compound_obj = {}

    for gcms_peak in gcms:

        # has a compound matched
        if gcms_peak:
            
            compound_obj = gcms_peak.highest_ss_compound
            
            if not compound_obj.ri in dict_ri_rt.keys():
                
                dict_ri_rt[compound_obj.ri] = [(gcms_peak.mass_spectrum.rt, compound_obj)]

            else:
                
                dict_ri_rt[compound_obj.ri].append((gcms_peak.mass_spectrum.rt, compound_obj))
            
            print(compound_obj.name, gcms_peak.mass_spectrum.rt, compound_obj.spectral_similarity_score)
    
    ris = [i for i in  dict_ri_rt.keys()]
    rts = [max(i, key = lambda c: c[1].spectral_similarity_score)[0] for i in dict_ri_rt.values()]
    
    rt_ri_pairs = list(zip(rts, ris)) 
    
    print(rt_ri_pairs)

    return rt_ri_pairs

def run(args):
    
    file_path, ref_file_path, ref_dict = args
    gcms = get_gcms(file_path)
            
    for gcms_peak in gcms:
        
        gcms_peak.calc_ri(ref_dict)
        
    lowResSearch = LowResMassSpectralMatch(gcms, ref_file_path)
    
    lowResSearch.run()

    return gcms

def calibrate_and_search(out_put_file_name):
    
    import csv
    
    ref_dict = get_reference_dict()
    
    if ref_dict:
        
        file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"
        ref_file_path = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"
        
        gcms = run((file_path, ref_file_path, ref_dict))
        
        gcms.to_csv(out_put_file_name)
        gcms.to_excel(out_put_file_name, highest_score=False)
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

def test_run_gcms_pipeline():

    out_put_file_name = 'STD_Mix1'
    calibrate_and_search(out_put_file_name)

if __name__ == "__main__":
    test_run_gcms_pipeline()