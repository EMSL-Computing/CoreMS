__author__ = "Joshua Sakai"
__date__ = "Jul 29, 2020"
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append(".")

import pytest
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import statistics as st
import csv

from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch
from corems.mass_spectra.calc.LF_Targeted import LossFinderTargeted


def run_targetedLF(file_path, ref_file):

    Loss_finder = LossFinderTargeted()

    reader_gcms = ReadAndiNetCDF(file_path)
    reader_gcms.run()

    gcms = reader_gcms.get_gcms_obj()
    gc_ms = gcms._ms

    Loss_finder.noise_cutoff = float(0.85)
    Loss_finder.tolerance = float(2)
    Loss_finder.ref_file = ref_file

    mz_dict, abund = Loss_finder.ms_info_get(gc_ms)

    range_ref = Loss_finder.loss_ref_get(ref_file, Loss_finder.tolerance)

    mz_filtered, abund_filtered = Loss_finder.threshold_filter(mz_dict, abund, Loss_finder.noise_cutoff)

    offset_hits = Loss_finder.findpeakoffset(range_ref, mz_filtered, abund_filtered)

    return offset_hits, Loss_finder.mz_count

def run_LF_pipeline():

    file_path = Path.cwd() / "tests/tests_data/gcms/" / "GCMS_FAMES_01_GCMS-01_20191023.cdf"

    ref_file = Path.cwd() / "tests/tests_data/" / "NeutralLossList.csv"

    offset_hits, mz_count = run_targetedLF(file_path, ref_file)

    data = pd.DataFrame()

    for chem, stats in offset_hits.items():
    
        data = data.append(pd.DataFrame(stats), ignore_index=True)

    data.columns = ['mz1','mz2','Intensity1','Intensity2','Name','scan_number']
    data.to_csv('test_LF.csv', encoding='utf-8', index=False, header = True, mode='w')
    chem = data['Name'].value_counts().keys().tolist()
    freq = data['Name'].value_counts().tolist()
    chem.reverse()
    freq.reverse()
    plt.barh(chem,freq)
    
   
    plt.title("Chemical Loss Frequency")
    plt.xlabel("Number of Loss found")
    
    plt.tight_layout()
    #ax = plt.subplot()
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    for n, i in enumerate(chem):
    #Create an axis text object
        plt.text(freq[n] + 12, #X location of text (with adjustment)
            n, #Y location
            s=f'{round(((freq[n])/mz_count)*100,3)}%', #Required label with formatting
            #s=f'{range_ref[n]}', #Required label with formatting
            va='center', #Vertical alignment
            color='#661D98', #Font colour and size
            fontsize=12)

    plt.savefig('test_LF.png')

if __name__ == '__main__':
    
    run_LF_pipeline()
