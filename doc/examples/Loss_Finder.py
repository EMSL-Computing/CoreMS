import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd
import statistics as st
import cProfile

from matplotlib import pyplot as plt
from numpy import array, polyfit, poly1d
from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems import get_dirname, get_filename
from corems.mass_spectra.factory.GC_Class import GCMSBase
from corems.mass_spectra.calc.LF_Targeted import LossFinderTargeted
import glob

def run_targetedLF(file_path, ref_file):

        Loss_finder = LossFinderTargeted()

        reader_gcms = ReadAndiNetCDF(file_path)
	
        reader_gcms.run()
    
        gcms = reader_gcms.get_gcms_obj()

        gc_ms = gcms._ms

        Loss_finder.noise_cutoff = float(0.7)

        Loss_finder.tolerance = float(5)

        mz_dict, abund = Loss_finder.ms_info_get(gc_ms)

        range_ref = Loss_finder.loss_ref_get(ref_file, Loss_finder.tolerance)

        mz_filtered, abund_filtered = Loss_finder.threshold_filter(mz_dict, abund, Loss_finder.noise_cutoff)

        offset_hits = Loss_finder.findpeakoffset(range_ref, mz_filtered, abund_filtered)

        Loss_finder.LF_out(offset_hits, Loss_finder.mz_count)

        return offset_hits, Loss_finder.mz_count

if __name__ == '__main__':

    file_path = get_filename()

    ref_file = '/mnt/c/ubuntu_home/loss_finder/NeutralLossList.csv'

    output, mz_count = run_targetedLF(file_path, ref_file)

