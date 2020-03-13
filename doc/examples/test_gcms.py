import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("C:\\Users\\eber373\\Desenvolvimento\\Projects-Python\\CoreMS")

from pathlib import Path

from multiprocessing import Pool
from matplotlib import pyplot
from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch

from numpy import array, polyfit, poly1d
from corems.encapsulation.settings.processingSetting import CompoundSearchSettings

def sql_database(file_location):
    
    print(file_location)

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

    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_path = file_dialog.getOpenFileName(None, "FAMES REF FILE", filter="*.cdf")[0]
    app.exit()
    
    ref_file_path = Path.cwd() / "tests/tests_data/gcms/" / "FAMES_REF.MSL"

    gcms = get_gcms(file_path)

    lowResSearch = LowResMassSpectralMatch(gcms, ref_file_path, calibration=True)

    lowResSearch.run()

    dict_rt_ri = {}

    for gcms_peak in gcms:

        # has a compound matched
        if gcms_peak:
            
            compound_obj = gcms_peak.highest_score_compound
            
            #print(compound_obj.name, gcms_peak.mass_spectrum.rt, compound_obj.similarity_score)
            dict_rt_ri[gcms_peak.mass_spectrum.rt] = compound_obj.ri
    

    rts = list(dict_rt_ri.keys() )
    ris = list(dict_rt_ri.values() )
    
    # retention time calibration curve RT vs RI
    poli = poly1d(polyfit(ris,rts, 1))
    
    # ensures all datapoint are represented up to C40 
    RI =  array(list(range(0,4000+CompoundSearchSettings.ri_spacing, CompoundSearchSettings.ri_spacing)))
    RT = (poli(RI))
    
    dict_RT_RI = dict(zip(RT, RI))

    #This will retains the experimental retention index
    dict_RT_RI.update(dict_rt_ri)
    
    return dict(sorted(dict_RT_RI.items()))

def run(args):
    
    file_path, ref_file_path, ref_dict = args
    gcms = get_gcms(file_path)
            
    for gcms_peak in gcms:
        
        gcms_peak.calc_ri(ref_dict)
        
    lowResSearch = LowResMassSpectralMatch(gcms, ref_file_path)
    
    lowResSearch.run()

    return gcms

def calibrate_and_search(out_put_file_name, cores):
    
    import csv

    with open(out_put_file_name, mode='w', newline='') as results_file:
        
        results_writer = csv.writer(results_file, delimiter=',')

        results_writer.writerow(["sample_name", "retention_time", "retention_time_ref",  "abundance", "area", "retention_index", "retention_index_ref", "cosine_correlation", "compound_name" ])
        
        ref_dict = get_reference_dict()
        
        #app = QApplication(sys.argv)
        file_dialog = QFileDialog()
        file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
        
        file_locations = file_dialog.getOpenFileNames(None, "Standard Files", filter="*.cdf")
        
        ref_file_path = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"

        p = Pool(cores)

        args = [(file_path, ref_file_path, ref_dict) for file_path in file_locations[0]]
        
        gcmss = p.map(run, args)
        
        for gcms in gcmss:
            
            #gcms.plot_processed_chromatogram()

            #gcms.plot_gc_peaks()

            #gcms.plot_chromatogram()

            #gcms.plot_smoothed_chromatogram()

            #gcms.plot_baseline_subtraction()

            #gcms.plot_detected_baseline()

            #pyplot.show()

            for gcms_peak in gcms:
                
                if gcms_peak:
                    
                    for compound_obj in gcms_peak:
                        
                        
                        results_writer.writerow([gcms.sample_name, gcms_peak.rt, compound_obj.rt, gcms_peak.tic, gcms_peak.area, gcms_peak.ri, compound_obj.ri, compound_obj.similarity_score, compound_obj.name])    
                
                else:
                    
                        results_writer.writerow([gcms.sample_name, gcms_peak.rt, None, gcms_peak.tic, gcms_peak.area, gcms_peak.ri, None, None, None])    
            
if __name__ == '__main__':                           
    
    cores = 6
    out_put_file_name = 'Group1_Standards.csv'
    calibrate_and_search(out_put_file_name, cores)


