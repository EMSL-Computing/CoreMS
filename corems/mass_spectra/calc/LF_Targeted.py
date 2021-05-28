from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import statistics as st
import csv

from corems.mass_spectra.calc.GC_Calc import GC_Calculations
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecBase
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecCentroidLowRes
from corems.chroma_peak.factory.ChromaPeakClasses import GCPeak
from corems.mass_spectra.calc import SignalProcessing as sp

class LossFinderTargeted(GC_Calculations):

    def __init__(self, ref_file = 'Unknown', noise_cutoff = 'Unknown', tolerance = "Unknown"):

        self.tolerance = float()

        self.noise_cutoff = float()

        self.offset_hits = {}
        
        self.mz_filtered = {}

        self.abund_filtered = {}
        
        self.mz_count = float()
        
        #self._mz_exp = MassSpecBase._mspeaks
        #self._abundance = MassSpecBase._abundance

    def ms_info_get(self, mass_spectra_obj):

        mz_dict = {}
        abund = {}

        self.mz_count = sum([len(mz) for mz in mass_spectra_obj.values()])

        for scan_number, ms_obj in mass_spectra_obj.items():
            
            mz_dict.update({scan_number:ms_obj.mz_exp})
            abund.update({scan_number:ms_obj.abundance})
        
        return mz_dict, abund

    def loss_ref_get(self, file_location, tolerance):

        offset_ref = {}
        range_ref = {}

        with open(file_location) as ref:

            ref_reader = csv.reader(ref, delimiter=',')
            next(ref_reader)

            for name, mass_offset in ref_reader:
                offset_ref.setdefault(name,float(mass_offset))

        for key, val in offset_ref.items():

            range_ref.update({key:(val-tolerance, val+tolerance)})

        return range_ref
    
    def threshold_filter(self, mz_dict, Intensity, noise_cutoff):
        
        for scan_number, info in Intensity.items():

            cutoff = st.mean(Intensity[scan_number])*noise_cutoff

            noise = set([peak for peak in Intensity[scan_number] if peak <=cutoff])

            mz_dict[scan_number] = [mz for mz, peak in zip(mz_dict[scan_number], Intensity[scan_number]) if peak not in noise ]

            Intensity[scan_number] = [peak for peak in Intensity[scan_number] if peak >= cutoff]

        return mz_dict , Intensity
    
    def mz_pair_checker(self, chem, lower, upper, mz1, spectrum, Intensity, scan_number):

        for mz2 in spectrum:

            if mz1 > mz2 and lower <= abs(mz1-mz2) <= upper:

                if chem not in self.offset_hits.keys():
                    self.offset_hits.update( {chem:[((mz2, mz1, Intensity[spectrum.index(mz2)] , Intensity[spectrum.index(mz1)], chem, scan_number ))]} )

                else:
                    self.offset_hits[chem].append( ((mz2,mz1, Intensity[spectrum.index(mz2)] , Intensity[spectrum.index(mz1)], chem, scan_number )) )

        

    def findpeakoffset(self, range_ref, mz_filtered, abund_filtered):

        count = 0    

        for chem, ref_pair in range_ref.items():
        
            for scan_number in mz_filtered.keys():

                while count < len(mz_filtered[scan_number])-1:
                #while count < 2:    

                    self.mz_pair_checker(chem, ref_pair[0], ref_pair[1], mz_filtered[scan_number][count], mz_filtered[scan_number], abund_filtered[scan_number], scan_number)
            
                    count+=1

                #print(self.offset_hits)

                count=0

        return self.offset_hits

    def LF_out(self, LF_dict, mz_count):
        data = pd.DataFrame()

        for chem, stats in LF_dict.items():
        
            data = data.append(pd.DataFrame(stats), ignore_index=True)

        data.columns = ['mz1','mz2','Intensity1','Intensity2','Name','scan_number']

        data.to_csv('Loss_finder.csv', encoding='utf-8', index=False, header = True, mode='w')

        chem = data['Name'].value_counts().keys().tolist()
        freq = data['Name'].value_counts().tolist()

        chem.reverse()
        freq.reverse()

        plt.barh(chem,freq)
        
       
        plt.title("Chemical Loss Frequency")
        plt.xlabel("Number of Loss found")
        
        plt.tight_layout()

        ax = plt.subplot()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        for n, i in enumerate(chem):
        #Create an axis text object
            plt.text(freq[n] + 12, #X location of text (with adjustment)
                n, #Y location
                s=f'{round(((freq[n])/mz_count)*100,3)}%', #Required label with formatting
                #s=f'{range_ref[n]}', #Required label with formatting
                va='center', #Vertical alignment
                color='#661D98', #Font colour and size
                fontsize=12)

        # plt.show()

    def plot_offset(self):
        #MassSpecBase.plot_mz_domain_profile(MassSpecBase)
        #MassSpecBase.plot_profile_and_noise_threshold(MassSpecBase)

        #plt.show()

        out = MassSpecBase.get_mz_and_abundance_peaks_tuples(MassSpecBase)

        return out