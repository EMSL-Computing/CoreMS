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
    sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30), (17.111, 18.111))

def stand_alone():

    file_path = get_filename()

    reader_gcms = ReadAndiNetCDF(file_path)

    reader_gcms.run()

    gcms = reader_gcms.get_gcms_obj()

    gcms.process_chromatogram()

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

def auto_calibrate_and_search(file_locations, output_file_name, jobs, calibration_file_path):

    ref_dict, cal_file_path = get_reference_dict(calibration_file_path=calibration_file_path)

    if ref_dict:

        # run in multiprocessing mode
        pool = Pool(jobs)
        args = [(file_path, ref_dict, cal_file_path) for file_path in file_locations]
        gcmss = pool.map(run, args)
        pool.close()
        pool.join()
        for gcms in gcmss:

            gcms.to_hdf()
            gcms.to_csv(output_file_name)
            # print(output_file_name)


def calibrate_and_search(out_put_file_name, jobs):

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

                gcms.to_csv(out_put_file_name, write_metadata=True, id_label="emsl:")

                # gcms.to_excel(out_put_file_name)
                # gcms.to_pandas(out_put_file_name)
                gcms.to_hdf()

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

def auto_process(jobs):

    import os
    rootdir = get_dirname()

    out_put_file_names = list(os.walk(rootdir))[0][1]

    # print(out_put_file_names[0])
    for out_put_file_name in out_put_file_names:
        # print(out_put_file_name)

        file_locations = glob.glob(str((rootdir / out_put_file_name)) + "/*.cdf")
        calibration_file_path = ''
        for file_path in file_locations:
            if "FAME" in file_path:
                calibration_file_path = file_path
        if calibration_file_path:
            auto_calibrate_and_search(file_locations, out_put_file_name, jobs, calibration_file_path)
        else:
            print("Could not find a calibration experimental file for {}".format(out_put_file_name))


if __name__ == '__main__':
    # import matplotlib
    # matplotlib.use('TkAgg')
    # %%
    cores = 4
    out_put_file_group_name = 'sql_test'
    calibrate_and_search(out_put_file_group_name, cores)
    # start_sql_from_file()
    # auto_process(cores)
    # stand_alone()
