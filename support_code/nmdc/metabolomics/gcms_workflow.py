import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path
from multiprocessing import Pool
import cProfile


from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch
from corems.mass_spectra.calc.GC_RI_Calibration import get_rt_ri_pairs
from support_code.filefinder import get_dirname, get_filename

import glob
from nmdc.metadata.nmdc_registration import NMDC_Metadata

def start_sql_from_file():

    ref_lib_path = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"
    sql_obj = ReadNistMSI(ref_lib_path).get_sqlLite_obj()
    return sql_obj


def get_gcms(file_path):

    reader_gcms = ReadAndiNetCDF(file_path)

    reader_gcms.run()

    gcms = reader_gcms.get_gcms_obj()

    # gcms.process_chromatogram()

    return gcms

def get_reference_dict(calibration_file_path=False):

    from PySide2.QtWidgets import QFileDialog, QApplication
    from PySide2.QtCore import Qt

    if not calibration_file_path:
        app = QApplication(sys.argv)
        file_dialog = QFileDialog()
        file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
        file_path = file_dialog.getOpenFileName(None, "FAMES REF FILE", filter="*.cdf")[0]
        file_dialog.close()
        app.exit()
    else:
        file_path = calibration_file_path

    if not file_path:
        return None

    else:

        gcms_ref_obj = get_gcms(file_path)
        # sql_obj = start_sql_from_file()
        rt_ri_pairs = get_rt_ri_pairs(gcms_ref_obj)  # sql_obj=sql_obj)
        # !!!!!! READ !!!!! use the previous two lines if db/pnnl_lowres_gcms_compounds.sqlite does not exist
        # and comment the next line
        # rt_ri_pairs = get_rt_ri_pairs(gcms_ref_obj)

        return rt_ri_pairs, file_path

def run(args):

    file_path, ref_dict, cal_file_path = args

    gcms = get_gcms(file_path)

    gcms.process_chromatogram()

    gcms.calibrate_ri(ref_dict, cal_file_path)

    # sql_obj = start_sql_from_file()

    lowResSearch = LowResMassSpectralMatch(gcms)  # sql_obj=sql_obj)
    # !!!!!! READ !!!!! use the previous two lines if db/pnnl_lowres_gcms_compounds.sqlite does not exist
    # and comment the next line
    # lowResSearch = LowResMassSpectralMatch(gcms)
    lowResSearch.run()

    return gcms

def calibrate_and_search(out_put_file_name, jobs, dms_file_path):

    from PySide2.QtWidgets import QFileDialog
    from PySide2.QtCore import Qt

    ref_dict, cal_file_path = get_reference_dict()
    if ref_dict:

        file_dialog = QFileDialog()
        file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)

        if file_dialog:

            file_locations = file_dialog.getOpenFileNames(None, "Standard Compounds Files", filter="*.cdf")
            file_dialog.close()

            # run in multiprocessing mode
            pool = Pool(jobs)
            args = [(file_path, ref_dict, cal_file_path) for file_path in file_locations[0]]
            gcmss = pool.map(run, args)
            pool.close()
            pool.join()

            for file_index, gcms in enumerate(gcmss):

                file_path = Path(file_locations[0][file_index])
                # print(out_put_file_name)

                gcms.to_csv(file_path, write_metadata=False, id_label="emsl:")
                nmdc = NMDC_Metadata(file_path, cal_file_path, file_path.with_suffix(".cvs"), dms_file_path)
                nmdc.create_nmdc_gcms_metadata(gcms)

                # gcms.to_excel(out_put_file_name)
                # gcms.to_pandas(out_put_file_name)
                # gcms.to_hdf()

                # df = gcms.get_dataframe()
                # json_data = gcms.to_json()

                # print(json_data)

                # gcms.plot_processed_chromatogram()

                # gcms.plot_gc_peaks()

                # gcms.plot_chromatogram()

                # gcms.plot_smoothed_chromatogram()

                # gcms.plot_baseline_subtraction()

                # gcms.plot_detected_baseline()

                # matplotlib.pyplot.show()

def worker(args):

    cProfile.runctx('run(args)', globals(), locals(), 'gc-ms.prof')



if __name__ == '__main__':
    # import matplotlib
    # matplotlib.use('TkAgg')
    # %%
    cores = 4
    out_put_file_group_name = 'sql_test'
    dms_file_path="db/GC-MS Metabolomics Experiments to Process Final.xlsx"
    calibrate_and_search(out_put_file_group_name, cores, dms_file_path)
   