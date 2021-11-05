##### CWD 2021-10-27
#### (1) import .csv with LC-ICPMS data; (2) for each metal, pick peaks; (3) return dictionary with RT as key, metal as value

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from tqdm import tqdm 

import os
import pandas as pd
import numpy as np

from pathlib import Path

import matplotlib.pyplot as plt
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

from corems.mass_spectra.calc import SignalProcessing as sp

#Import LC-ICPMS data and plot it:

import csv


def find_metal_peaks(metal, icp_data, sn, min_peak_datapoints):
	icpt="Time_" + metal
	rt = icp_data[icpt]
	icp_signal = icp_data[metal]
	max_icp_signal = np.max(icp_signal)
	max_height = max_icp_signal
	max_prominence = 100 
	peak_indices = sp.peak_picking_first_derivative(rt, icp_signal, max_height, max_prominence, max_icp_signal, min_peak_datapoints,
                                                                   signal_threshold=sn, correct_baseline=True)
	rt_aps = []
	for peak_index in peak_indices:
		rt_aps.append(rt[peak_index[1]])

	rt_dict = dict.fromkeys(rt_aps,metal)
	return rt_dict

icpfile = "tests/tests_data/cwd_211018_day7_8_c18_1uMcobalamin_10uL.csv"
icpdata = np.genfromtxt(icpfile, dtype=float, delimiter=',', names=True) 

metal = '59Co'

icp_signal = icpdata[metal]  
min_peak_datapoints = 20

rts = find_metal_peaks(metal,icpdata, 1,  min_peak_datapoints)
print(rts)

fig, host = plt.subplots()
host.plot(icpdata['Time_'+metal],icpdata[metal])
host.set_xlabel('Time (s)')
host.set_ylabel(metal +' intensity (counts)')
plt.show()