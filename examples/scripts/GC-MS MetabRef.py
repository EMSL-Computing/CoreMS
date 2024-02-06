import warnings
warnings.filterwarnings("ignore")

import os
from pathlib import Path
from multiprocessing import Pool
import sys

import corems
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch
from corems.mass_spectra.calc.GC_RI_Calibration import get_rt_ri_pairs

import metabrefapi


def start_gcms_metabref_sql(normalize=False, url='sqlite://'):
    # Get MetabRef GCMS library
    metabref_lib = metabrefapi.get_metabref_gcms_library()

    # Convert to CoreMS format
    corems_lib = metabrefapi.metabref_to_corems(metabref_lib, normalize=normalize)

    # Convert to SQLite
    corems_lib = metabrefapi.corems_to_sqlite(corems_lib, url=url)

    return corems_lib


def start_fames_metabref_sql(normalize=False, url='sqlite://'):
    # Get MetabRef GCMS library
    metabref_lib = metabrefapi.get_metabref_fames_library()

    # Convert to CoreMS format
    corems_lib = metabrefapi.metabref_to_corems(metabref_lib, normalize=normalize)

    # Convert to SQLite
    corems_lib = metabrefapi.corems_to_sqlite(corems_lib, url=url)

    return corems_lib


def get_gcms(filepath):
    # Initialize reader
    reader_gcms = ReadAndiNetCDF(filepath)

    # Run reader
    reader_gcms.run()

    # Get GCMS object
    gcms = reader_gcms.get_gcms_obj()

    # # Process chromatogram
    # gcms.process_chromatogram()

    return gcms

def stand_alone():
    # Determine filename
    filepath = get_filename()

    # Initialize reader
    reader_gcms = ReadAndiNetCDF(filepath)

    # Run
    reader_gcms.run()

    # Get GCMS object
    gcms = reader_gcms.get_gcms_obj()

    # Process chromatogram
    gcms.process_chromatogram()


def get_reference_dict(calibration_filepath=False):
    from PySide2.QtWidgets import QFileDialog, QApplication
    from PySide2.QtCore import Qt

    # Open dialog to select calibration file
    if not calibration_filepath:
        app = QApplication(sys.argv)
        file_dialog = QFileDialog()
        file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
        filepath = file_dialog.getOpenFileName(None, "FAMES REF FILE", filter="*.cdf")[0]
        file_dialog.close()
        app.exit()

    # Calibration file specified directly
    else:
        filepath = calibration_filepath

    # No filepath
    if not filepath:
        raise ValueError("Must supply calibration file.")
        
    # Parse supplied calibration data
    gcms_ref_obj = get_gcms(filepath)

    # Build calibration SQLite database from MetabRef
    sql_obj = start_fames_metabref_sql()

    # Determine calibration pairs
    rt_ri_pairs = get_rt_ri_pairs(gcms_ref_obj, sql_obj=sql_obj)

    return rt_ri_pairs, filepath


def run(args):
    # Unpack arguments
    filepath, ref_dict, cal_filepath = args

    # Parse supplied file
    gcms = get_gcms(filepath)

    # Process chromatogram
    gcms.process_chromatogram()

    # Calibrate retention index
    gcms.calibrate_ri(ref_dict, cal_filepath)

    # Initialize GCMS refeerence database from MetabRef
    sql_obj = start_gcms_metabref_sql()

    # Initialize spectral match
    lowResSearch = LowResMassSpectralMatch(gcms, sql_obj=sql_obj)

    # Run spectral match
    lowResSearch.run()

    return gcms


def auto_calibrate_and_search(file_locations, output_filename, jobs, calibration_filepath):
    ref_dict, cal_filepath = get_reference_dict(calibration_filepath=calibration_filepath)

    if ref_dict:
        # run in multiprocessing mode
        pool = Pool(jobs)
        args = [(filepath, ref_dict, cal_filepath) for filepath in file_locations]
        gcmss = pool.map(run, args)
        pool.close()
        pool.join()
        for gcms in gcmss:
            gcms.to_hdf()
            gcms.to_csv(output_filename)
            # print(output_filename)


def calibrate_and_search(filename, jobs):
    from PySide2.QtWidgets import QFileDialog, QApplication
    from PySide2.QtCore import Qt
    
    import csv
    
    ref_dict, cal_filepath = get_reference_dict()
    
    if ref_dict:
        file_dialog = QFileDialog()
        file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    
        if file_dialog:
            file_locations = file_dialog.getOpenFileNames(None, "Standard Compounds Files", filter="*.cdf")
            file_dialog.close()
    
            # run in multiprocessing mode
            pool = Pool(jobs)
            args = [(filepath, ref_dict, cal_filepath) for filepath in file_locations[0]]
            gcmss = pool.map(run, args)
            pool.close()
            pool.join()
            for gcms in gcmss:
    
                gcms.to_csv(filename)
                gcms.to_hdf()


def worker(args):
    cProfile.runctx('run(args)', globals(), locals(), 'gc-ms.prof')


def auto_process(jobs):
    rootdir = get_dirname()

    output_filenames = list(os.walk(rootdir))[0][1]

    # print(output_filenames[0])
    for output_filename in output_filenames:
        # print(output_filename)

        file_locations = glob.glob(str((rootdir / output_filename)) + "/*.cdf")
        calibration_filepath = ''
        for filepath in file_locations:
            if "FAME" in filepath:
                calibration_filepath = filepath
        if calibration_filepath:
            auto_calibrate_and_search(file_locations, output_filename, jobs, calibration_filepath)
        else:
            print("Could not find a calibration experimental file for {}".format(output_filename))


if __name__ == '__main__':
    TOKEN_PATH = "metabref.token"
    metabrefapi.set_token(TOKEN_PATH)

    jobs = 6
    filename = 'json_test'

    calibrate_and_search(filename, jobs)
