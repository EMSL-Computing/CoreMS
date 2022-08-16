##### rMB 2022-04-11

####Replicate isotope pattern algorithm in python for batch processing

##### RMB 2022-04-11

####Replicate isotope pattern algorithm in python for a single file and pattern. 

###### Settings for pattern mining
timerange=(8,14) #in minutes
peakwidth=0.25 #in minutes
slope_filter=(0.5,2) # normalized slope (1=true)
correlation=0.6 #minimum r-squared correlation cut-off. 

mass_tolerance=0.001
ratio_tolerance=1.6

###### Set file folder and THERMO RAW file name here:
file_location = '/Users/boiteaur/Desktop/Logan Samples/'

file_name="rmb_20220627_pos_logan_7002_1_6.raw"

#####Set isotope pattern using atom.epattern. Just requires elements and max number of isotopes used
element='Zn'
nisotope_used=2
valence=2

save_file="_"+element

#### Alternative way to set isotope pattern using atom.ipattern
isotopes=['64Zn','66Zn','68Zn']
requirement=['Y','Y','N']

# Import the os module
import os
import pandas as pd
import numpy as np

#Isotope pattern matching algorithm for LC-FTMS Data
# Change the current working directory

os.chdir('/Users/boiteaur/Desktop/CoreMS_metallomics/CoreMS')

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt

from corems.mass_spectra.input import rawFileReader
#from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
#from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

import matplotlib.backends.backend_pdf
from support_code import isotope_pattern as ip
from support_code import AtomsDescription_standardized as atom

#Set peak detection threshold method
MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 2
MSParameters.ms_peak.peak_min_prominence_percent = 0.001

MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -2
MSParameters.molecular_search.max_ppm_error = 2

MSParameters.molecular_search.url_database = None
MSParameters.molecular_search.min_dbe = 0
MSParameters.molecular_search.max_dbe = 16

MSParameters.molecular_search.usedAtoms['C'] = (1,30)
MSParameters.molecular_search.usedAtoms['H'] = (4,60)
MSParameters.molecular_search.usedAtoms['O'] = (1,8)
MSParameters.molecular_search.usedAtoms['N'] = (0,6)
MSParameters.molecular_search.usedAtoms['S'] = (0,0)
MSParameters.molecular_search.usedAtoms[element] = (0,1)

MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = False

MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"

########################## Main work functions

#Create LCMS object
parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+file_name)

#Create isotope pattern
#pattern=atom.epattern(element,nisotope_used,atom.atoms)
pattern=atom.ipattern(isotopes,requirement,atom.atoms)

print(pd.DataFrame(pattern))

#Mine pattern.
results=parser.isotopehunter(pattern,timerange,mass_tolerance,ratio_tolerance,peakwidth,correlation,slope_filter)

#Assign molecular formula.
for result in results[2]:
    ip.metal_assignment(parser,result,pattern,file_name)
######################### Main reporting functions

pdf = matplotlib.backends.backend_pdf.PdfPages(file_location+file_name+save_file+'.pdf')

#Generate plots:
pdf.savefig(ip.isotopehunter_qc_plots(results,pattern,file_name,correlation),dpi=200,bbox_inches='tight')

#ip.isotopehunter_qc_plots(results,pattern,file_name,correlation)

for result in results[2]:
    fig, (ax1,ax2,ax3) = plt.subplots(3,1)
    fig.set_size_inches(8,18)
    ip.metal_chromatogram(ax1,parser,result,pattern,timerange,file_name)
    ip.apoplot(ax2,parser,result,pattern,timerange,file_name,2)
    ip.MS_pattern_plot(ax3,parser,result,pattern)
    pdf.savefig(fig,dpi=200,bbox_inches='tight')

pdf.close()

#Print data out. 
pd.DataFrame(results[2]).to_csv(file_location+file_name+save_file+'.csv')

print('end')