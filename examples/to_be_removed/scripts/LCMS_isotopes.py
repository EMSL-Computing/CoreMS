##### CWD 2021-09-29
####Replicate isotope pattern algorithm in python for a single file and peak pair.

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


#set file here
file_location = "tests/tests_data/icpms/rmb_161221_kansas_h2o_2"

#Set peak detection threshold method
MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'
MSParameters.mass_spectrum.noise_threshold_min_relative_abundance = 1

MSParameters.mass_spectrum.noise_threshold_method = 'log'
MSParameters.mass_spectrum.noise_threshold_min_s2n = 10

#Parser for thermo RAW files. 
parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location)

t_ion_chromatogram = parser.get_tic()

t_ion_subset=t_ion_chromatogram[(t_ion_chromatogram["Time"]>8) & (t_ion_chromatogram["Time"]<9)]

#print(t_ion_chromatogram['Time'])
#print(parser.start_scan)
#print(parser.end_scan)


#Function for obtaining extracted ion chromatogram 
def get_EIC(parser,mass,dmz,scanrange):
    
    EIC = np.zeros((len(scanrange),3))  ## scan, time, eic
    EIC[:,0] = scanrange
    for scan in tqdm(scanrange, desc = 'extracting ion chromatogram'):    

        scanStatistics = parser.iRawDataPlus.GetScanStatsForScanNumber(scan)
        EIC[np.where(EIC[:,0] == scan),1] = scanStatistics.StartTime

        current_dictionary = parser.get_data(scan,1,scan_type="Profile")
        current = np.zeros((len(current_dictionary['m/z']),2))
        current[:,0] = current_dictionary['m/z']
        current[:,1] = current_dictionary['Peak Height']
        subset = current[abs(current[:,0] - mass) < dmz]

        EIC[np.where(EIC[:,0] == scan),2] =sum(subset[:,1])
    return(EIC)

scanrange=range(parser.start_scan,parser.end_scan)
mass=677
dmz=1

EIC=get_EIC(parser,mass,dmz,scanrange)

fig, host = plt.subplots()
host.plot(EIC[:,1],EIC[:,2])
host.set_xlabel('Time (min)')
host.set_ylabel('Intensity (counts)')
plt.show()




#Get MS1 scans numbers only.
first_scan = parser.start_scan
final_scan =parser.end_scan
scanrange = range(first_scan, final_scan)

MSn = np.zeros((len(scanrange),3))  ## scan, time, n
MSn[:,0] = scanrange
for scan in tqdm(range(parser.start_scan, parser.end_scan), desc = 'extracting MS1 scans'):    

	scanStatistics = parser.iRawDataPlus.GetScanStatsForScanNumber(scan)
	MSn[np.where(MSn[:,0] == scan),1] = scanStatistics.StartTime
	MSn[np.where(MSn[:,0] == scan),2] = int(parser.get_scan_header(scan)['Master Scan Number:'])


MS1scans=MSn[np.where(MSn[:,2] == 0)]
#MS1scans={'scan':MS1scans[:,0], 'time':MS1scans[:,1], 'n':MS1scans[:,2]}
#MSn_dict={'scan':MSn[:,0], 'time':MSn[:,1], 'n':MSn[:,2]}
#scanrange=MS1scans['scan']
scanrange=MS1scans[:,0]

mass=677
dmz=1

EIC=get_EIC(parser,mass,dmz,scanrange)

fig, host = plt.subplots()
host.plot(EIC[:,1],EIC[:,2])
host.set_xlabel('Time (min)')
host.set_ylabel('Intensity (counts)')
plt.show()



#Import LC-ICPMS data and plot it:

import csv

icpfile = "tests/tests_data/icpms/161220_soils_hypercarb_3_kansas_qH2O.csv"

icpdata = np.genfromtxt(icpfile, dtype=float, delimiter=',', names=True) 

fig, host = plt.subplots()
host.plot(icpdata['Time_56Fe'],icpdata['56Fe'])
host.set_xlabel('Time (s)')
host.set_ylabel('56Fe intensity (counts)')
plt.show()

fig, host = plt.subplots()
host.plot(icpdata['Time_63Cu'],icpdata['63Cu'])
host.set_xlabel('Time (s)')
host.set_ylabel('63Cu intensity (counts)')
plt.show()

fig, host = plt.subplots()
host.plot(icpdata['Time_59Co'],icpdata['59Co'])
host.set_xlabel('Time (s)')
host.set_ylabel('59Co intensity (counts)')
plt.show()




#Pick out a time-slice of the ICPMS data that will be correlated w/ EIC data. 

timestart=535
timestop=575
icpi="59Co"
icpt="Time_" + icpi
icpslice = icpdata[np.where((icpdata[icpt] >= timestart) & (icpdata[icpt] <= timestop))]

fig, host = plt.subplots()
host.plot(icpslice[icpt],icpslice[icpi])
host.set_xlabel('Time (s)')
host.set_ylabel(icpi+' intensity (counts)')
plt.show()



#Correlate every ESIMS mz detected across the time range with the metal intensity. 
#This section obtains EIC's for every m/z over the time range. 

MS1scans=MSn[np.where(MSn[:,2] == 0)]
scans=MS1scans[np.where((MS1scans[:,1] >= timestart/60) & (MS1scans[:,1] <= timestop/60))][:,0].tolist()
parser.chromatogram_settings.scans = scans

AverageMS = parser.get_average_mass_spectrum()
AverageMS.plot_mz_domain_profile()

plt.show()

print(AverageMS.mz_exp.size)

scanrange=scans

EICdict = {}

for mz in AverageMS.mz_exp[1:20]:
    #EIC = pd.DataFrame(index=scanrange, columns=['Time', 'EIC'])
    EIC = np.zeros((len(scanrange),3))  ## scan, time, eic
    EIC[:,0] = scanrange   
    mass=mz
    dmz=0.002
    print('m/z: ',mz)
    EIC=get_EIC(parser,mass,dmz,scanrange)
    EICdict[mz]=EIC

# sums all the mass spectra
mass_spectrum = parser.get_average_mass_spectrum()


mass_spectrum.plot_mz_domain_profile()
plt.show()

mass_spectrum.plot_profile_and_noise_threshold()
plt.show()


mass_spectrum.molecular_search_settings.error_method = 'None'
mass_spectrum.molecular_search_settings.min_ppm_error = -5
mass_spectrum.molecular_search_settings.max_ppm_error = 5

mass_spectrum.molecular_search_settings.url_database = None
mass_spectrum.molecular_search_settings.min_dbe = 0
mass_spectrum.molecular_search_settings.max_dbe = 50

mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 100)
mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 200)
mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 30)
mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 0)

mass_spectrum.molecular_search_settings.isProtonated = True
mass_spectrum.molecular_search_settings.isRadical = False
mass_spectrum.molecular_search_settings.isAdduct = False

# mass_spectrum.filter_by_max_resolving_power(15, 2)
SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

mass_spectrum.percentile_assigned(report_error=True)
mass_spectrum.molecular_search_settings.score_method = "prob_score"
mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

# export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")

mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)

mass_spectrum_by_classes.plot_ms_assigned_unassigned()
plt.show()
mass_spectrum_by_classes.plot_mz_error()
plt.show()
mass_spectrum_by_classes.plot_ms_class("O2")
plt.show()




mass_spectrum.molecular_search_settings.error_method = 'None'
mass_spectrum.molecular_search_settings.min_ppm_error = -2
mass_spectrum.molecular_search_settings.max_ppm_error = 4

mass_spectrum.molecular_search_settings.url_database = None
mass_spectrum.molecular_search_settings.min_dbe = 0
mass_spectrum.molecular_search_settings.max_dbe = 50

mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 100)
mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 200)
mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 30)
mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 6)
mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 0)

mass_spectrum.molecular_search_settings.isProtonated = True
mass_spectrum.molecular_search_settings.isRadical = False
mass_spectrum.molecular_search_settings.isAdduct = False

# mass_spectrum.filter_by_max_resolving_power(15, 2)
SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

mass_spectrum.percentile_assigned(report_error=True)
mass_spectrum.molecular_search_settings.score_method = "prob_score"
mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

# export_calc_isotopologues(mass_spectrum, "15T_Neg_ESI_SRFA_Calc_Isotopologues")

mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)

mass_spectrum_by_classes.plot_ms_assigned_unassigned()
plt.show()
mass_spectrum_by_classes.plot_mz_error()
plt.show()
mass_spectrum_by_classes.plot_ms_class("O2")
plt.show()