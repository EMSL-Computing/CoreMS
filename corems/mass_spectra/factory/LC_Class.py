"""
Created on Oct 15, 2021
"""
__author__ = "Yuri E. Corilo"
__date__ = "Oct 15, 2021"


from collections.abc import Mapping
from dataclasses import field, dataclass
import logging
from pathlib import Path
from typing import List, Dict, List, Tuple
import sys
import site


from matplotlib import axes
import numpy as np
from corems.chroma_peak.factory.ChromaPeakClasses import DataDependentPeak
from corems.encapsulation.factory.parameters import LCMSParameters
from corems.encapsulation.factory.processingSetting import LiquidChromatographSetting, MolecularFormulaSearchSettings

from corems.mass_spectra.calc.LC_Calc import LC_Calculations
from corems.mass_spectra.factory.LC_Temp import TIC_Data, EIC_Data

import clr

sys.path.append(site.getsitepackages()[0] + '/ext_lib/dotnet')
sys.path.append('ext_lib/dotnet')

clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.Data')
clr.AddReference('ThermoFisher.CommonCore.MassPrecisionEstimator')

from ThermoFisher.CommonCore.Data import ToleranceUnits
from ThermoFisher.CommonCore.Data.Business import ChromatogramTraceSettings, TraceType, MassOptions
from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, Range
from ThermoFisher.CommonCore.Data.FilterEnums import MSOrderType
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings

class LCMSBase(Mapping, LC_Calculations):
    """
    classdocs
    """

    def __init__(self, file_location, analyzer='Unknown', instrument_label='Unknown', sample_name=None):
        
        '''
         # Parameters
		----------
		file_location: text,  pathlib.Path(), or s3path.S3Path 
            Path object from pathlib containing the file location
        '''
		
        if  isinstance(file_location, str):
			# if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():
        
            raise FileExistsError("File does not exist: " + file_location)
        
        self.file_location = file_location
        
        if sample_name: 
        
            self.sample_name = sample_name

        else: 
        
            self.sample_name = file_location.stem
        
        self.analyzer = analyzer
        self.instrument_label = instrument_label
        
        self._retention_time_list = []
        self._scans_number_list = []
        self._tic_list = []

        self._ms = {}
        """
        key is scan number; value is MassSpectrum Class
        """
    
    def __len__(self):
        
        return len(self._ms)
        
    def __getitem__(self, scan_number):
        
        return self._ms.get(scan_number)

    def __iter__(self):

         return iter(self._ms.values()) 

    def add_mass_spectrum(self, mass_spec):

        self._ms[mass_spec.scan_number] = mass_spec

    def set_tic_list_from_data(self):

        self.tic = [self._ms.get(i).tic for i in self.scans_number]
        
        # self.set_tic_list([self._ms.get(i).get_sumed_signal_to_noise() for i in self.get_scans_number()])

    def set_retention_time_from_data(self):

        retention_time_list = []

        for key_ms in sorted(self._ms.keys()):

            retention_time_list.append(self._ms.get(key_ms).retention_time)

        self.retention_time = retention_time_list 

        # self.set_retention_time_list(sorted(self._ms.keys()))

    def set_scans_number_from_data(self):
        
        self.scans_number = sorted(self._ms.keys())

    @property
    def scans_number(self):

        return self._scans_number_list
    
    @scans_number.setter
    def scans_number(self, l):

        self._scans_number_list = l

    @property
    def retention_time(self):

        return self._retention_time_list
    
    @property
    def tic(self):

        return self._tic_list

    @retention_time.setter
    def retention_time(self, l):
        # self._retention_time_list = linspace(0, 80, num=len(self._scans_number_list))
        self._retention_time_list = l

    @tic.setter
    def tic(self, l):

        self._tic_list = l    


class DataDependentLCMS(LC_Calculations):
    
    def __init__(self, file_location, target_mzs:List[float], parser, analyzer:str ='Unknown', 
                 instrument_label:str='Unknown', sample_name:str=None) -> None:
        
        self._parser = parser

        if  isinstance(file_location, str):
			# if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():
        
            raise FileExistsError("File does not exist: " + file_location)
        
        self.file_location = file_location
        
        if sample_name: 
        
            self.sample_name = sample_name

        else: 
        
            self.sample_name = file_location.stem

        # place holder for DataDependentPeak list
        self._lcmspeaks = []
        
        # place holders for Full Chromatogram data 
        self._retention_time_list = []
        self._scans_number_list = []
        self._tic_list = []
        

    def __len__(self):
        
        return len(self._lcmspeaks)
        
    def __getitem__(self, position) -> DataDependentPeak:
        
        return self._lcmspeaks[position]

    def add_peak(self, mass_spectrum, peak_scans, eic_data, possible_molform=None):
        
        lcms_peak  = DataDependentPeak(self, mass_spectrum, peak_scans, eic_data, possible_molform=possible_molform)
        self._lcmspeaks.append(lcms_peak)

    @property
    def ms1_molecular_search_settings(self) -> MolecularFormulaSearchSettings:  
        return self.parameters.ms1_molecular_search
    
    @property
    def chromatogram_settings(self) -> LiquidChromatographSetting:
        return self.parameters.lc_ms

    @property
    def parameters(self) -> LCMSParameters:  
        return self._parser.parameters

    
    @property
    def eic_scans_number(self):

        return [i.apex_scan for i in self._lcmspeaks]

    @property
    def scans_number(self):

        return self._scans_number_list

    @property
    def retention_time(self):

        return self._retention_time_list
    
    @property
    def tic(self):

        return self._tic_list

    @retention_time.setter
    def retention_time(self, l):
        # self._retention_time_list = linspace(0, 80, num=len(self._scans_number_list))
        self._retention_time_list = l

    

    @tic.setter
    def tic(self, l):

        self._tic_list = l        

    def get_tic(self, ms_type='MS !d', peak_detection=True, 
                smooth=True, plot=False, ax=None) -> Tuple[ TIC_Data, axes.Axes]:

        '''ms_type: str ('MS', MS2')
        start_scan: int default -1 will select the lowest available
        end_scan: int default -1 will select the highest available
        returns:
            chroma: dict
            {
            Scan: [int]
                original thermo scan numberMS
            Time: [floats]
                list of retention times
            TIC: [floats]
                total ion chromatogram
            Apexes: [int]    
                original thermo apex scan number after peak picking     
            }
        '''

        settings = ChromatogramTraceSettings(TraceType.TIC)
        settings.Filter = ms_type

        chroma_settings = IChromatogramSettings(settings)

        data = self._parser.iRawDataPlus.GetChromatogramData([chroma_settings],
                                                     self.parameters.lc_ms.start_scan, self.parameters.lc_ms.end_scan)

        trace = ChromatogramSignal.FromChromatogramData(data)
        
        data = TIC_Data(time= [], scans= [], tic= [], apexes= [])

        if trace[0].Length > 0:

            for i in range(trace[0].Length):
                # print(trace[0].HasBasePeakData,trace[0].EndTime )
                
                # print("  {} - {}, {}".format( i, trace[0].Times[i], trace[0].Intensities[i] ))
                data.time.append(trace[0].Times[i])
                data.tic.append(trace[0].Intensities[i])
                data.scans.append(trace[0].Scans[i])

                #print(trace[0].Scans[i])
            if smooth:
                
                data.tic= self.smooth_tic(data.tic)

            else:
                
                data.tic= np.array(data.tic)

            if peak_detection:
                
                centroid_peak_indexes = [i for i in self.centroid_detector(data.time, data.tic)]
                data.apexes = centroid_peak_indexes

            if plot:
                if not ax:
                    import matplotlib.pyplot as plt
                    ax = plt.gca()
                    # fig, ax = plt.subplots(figsize=(6, 3))

                ax.plot(data.time, data.tic, label=' TIC')    
                ax.set_xlabel('Time (min)')
                ax.set_ylabel('a.u.')
                if peak_detection:
                    for peak_indexes in data.apexes:
                        
                        apex_index = peak_indexes[1]
                        ax.plot(data.time[apex_index], data.tic[apex_index], marker='x', linewidth=0)
                    
                
                # plt.show()
                return data, ax
            
            return data, None

        else:
            return None, None

    def get_eics(self, tic_data: dict, ms_type='MS !d', 
                 peak_detection=True, smooth=True, plot=False, 
                 ax=None, legend=False) -> Tuple[Dict[float, EIC_Data], axes.Axes]:

        '''ms_type: str ('MS', MS2')
        start_scan: int default -1 will select the lowest available
        end_scan: int default -1 will select the highest available

        returns:

            chroma: dict{target_mz: EIC_Data(
                                        Scans: [int]
                                            original thermo scan numbers
                                        Time: [floats]
                                            list of retention times
                                        TIC: [floats]
                                            total ion chromatogram
                                        Apexes: [int]    
                                            original thermo apex scan number after peak picking 
                                        )
                                        
        '''
        target_mzs = self._parser.selected_mzs
        
        options = MassOptions()
        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = self.chromatogram_settings.eic_tolerance_ppm

        all_chroma_settings = []

        for target_mz in target_mzs:

            settings = ChromatogramTraceSettings(TraceType.MassRange)
            settings.Filter = ms_type
            settings.MassRanges = [Range(target_mz, target_mz)]

            chroma_settings = IChromatogramSettings(settings)

            all_chroma_settings.append(chroma_settings)

        # chroma_settings2 = IChromatogramSettings(settings)
        # print(chroma_settings.FragmentMass)
        # print(chroma_settings.FragmentMass)
        # print(chroma_settings)
        # print(chroma_settings)

        data = self._parser.iRawDataPlus.GetChromatogramData(all_chroma_settings,
                                                     self.chromatogram_settings.start_scan, self.chromatogram_settings.end_scan, options)

        traces = ChromatogramSignal.FromChromatogramData(data)

        chroma = {}    
       
        if plot:
            from matplotlib.transforms import Bbox
            import matplotlib.pyplot as plt
            if not ax:
                # ax = plt.gca()
                # ax.clear()
                fig, ax = plt.subplots()
                
            else:
                fig = plt.gcf()    

                
            # plt.show()
        
        for i, trace in enumerate(traces):
            if trace.Length > 0:
                rt, eic, scans  = self._parser.get_rt_time_from_trace(trace)
                if smooth:
                    eic= self.smooth_tic(eic)
    
                chroma[target_mzs[i]] = EIC_Data(scans=scans, time=rt, eic= eic)
                if plot:
                    ax.plot(rt, eic, label="{:.5f}".format(target_mzs[i]))
            
        if peak_detection:
            
            #max_eic = self.get_max_eic(chroma)
            max_signal = max(tic_data.tic)
            
            for eic_data in chroma.values():
               
                eic = eic_data.eic
                time = eic_data.time

                if len(eic) != len(tic_data.tic):
                    logging.warn("The software assumes same lenth of TIC and EIC, this does not seems to be the case and the results mass spectrum selected by the scan number might not be correct")
                
                centroid_eics = self.eic_centroid_detector(time, eic, max_signal)
                eic_data.apexes = [i for i in centroid_eics]
                
                if plot:
                    for peak_indexes in eic_data.apexes:
                        
                        apex_index = peak_indexes[1]
                        ax.plot(time[apex_index], eic[apex_index], marker='x', linewidth=0)

        if plot:
            
            ax.set_xlabel('Time (min)')
            ax.set_ylabel('a.u.')
            ax.set_title(ms_type + ' EIC')
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.axes.spines['top'].set_visible(False)
            ax.axes.spines['right'].set_visible(False)

            if legend:
                legend = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 0, 0.07, 1))
                fig.subplots_adjust(right=0.76)
                #ax.set_prop_cycle(color=plt.cm.gist_rainbow(np.linspace(0, 1, len(traces))))

                d = {"down": 30, "up": -30}
                def func(evt):
                    if legend.contains(evt):
                        bbox = legend.get_bbox_to_anchor()
                        bbox = Bbox.from_bounds(bbox.x0, bbox.y0 + d[evt.button], bbox.width, bbox.height)
                        tr = legend.axes.transAxes.inverted()
                        legend.set_bbox_to_anchor(bbox.transformed(tr))
                        fig.canvas.draw_idle()

                fig.canvas.mpl_connect("scroll_event", func)    
            return chroma, ax
        else:
            return chroma, None     
            rt = []
            tic = []
            scans = []
            for i in range(traces[0].Length):
                # print(trace[0].HasBasePeakData,trace[0].EndTime )

                # print("  {} - {}, {}".format( i, trace[0].Times[i], trace[0].Intensities[i] ))
                rt.append(traces[0].Times[i])
                tic.append(traces[0].Intensities[i])
                scans.append(traces[0].Scans[i])

            return traces
            # plot_chroma(rt, tic)
            # plt.show()
    
    def process_ms1(self, dict_tarrget_mzs):

        target_mzs = self._parser.selected_mzs
        
        tic_data, ax_tic = self.get_tic(ms_type='MS !d', peak_detection=True, 
                                      smooth=True, plot=False)

        eics_data, axes = self.get_eics(tic_data)
        
        results_list = []

        scan_number_mass_spectrum = {}

        for mz, eic_data in eics_data.items():
            
            #all possible m/z from the same mix, should be one per m/z as per current lib
            possible_mf = dict_tarrget_mzs.get(mz)
            
            if eic_data.apexes:                    
                
                #dict_res = {}

                #names = [mf_obj.name for mf_obj in possible_mf]
                #molecular_formulae = [mf_obj.string for mf_obj in possible_mf]
                #rts = [eic_data.time[apex[1]] for apex in eic_data.apexes]
                #scans = [eic_data.scans[apex[1]] for apex in eic_data.apexes]
                #peak_height = [eic_data.eic[apex[1]] for apex in eic_data.apexes]

                #print("m/z =  {}, formulas = {}, names = {}, peaks indexes = {}, retention times = {}, abundance = {}".format(mz,
                #                                                                        molecular_formulae,
                #                                                                        names,
                #                                                                        scans,
                #                                                                        rts,
                #                                                                        peak_height) )
                #dict_res["Mix Name"] = current_mix
                #dict_res["Dataset"] = self.sample_name
                #dict_res["Compound Name"] = names[0]
                #dict_res["Neutral Formula"] = molecular_formulae[0]
                #dict_res["Target m/z (de)protonated"] = round(mz,6)
                #dict_res["Retention Times"] = rts
                #dict_res["Scans"] = scans
                #dict_res["Peak Height"] = peak_height
                
                #results_list.append(dict_res)

                for peak_indexes in eic_data.apexes:
            
                    apex_index = peak_indexes[1]
                    
                    original_scan = eic_data.scans[apex_index]
                    
                    original_indexes = [eic_data.scans[i] for i in peak_indexes] 
                    
                    eics_original_scans = self.eic_scans_number

                    if original_scan in self.eic_scans_number:

                        index = eics_original_scans.index(original_scan)
                        
                        scan_number_mass_spectrum[original_scan][1].extend(possible_mf)

                        self._lcmspeaks[index].add_molecular_formula(possible_mf)

                    else:

                        mass_spec = self._parser.get_average_mass_spectrum(auto_process=False)
                        
                        self._parser.chromatogram_settings.scans = [original_scan]
                        
                        mass_spec = self._parser.get_average_mass_spectrum()
                        
                        #mass_spec.min_ppm_error = - 5
                        #mass_spec.max_ppm_error = 5

                        scan_number_mass_spectrum[original_scan] = [mass_spec, [i for i in possible_mf]]
                        
                        #broadcast ms1 molecular formula search
                        mass_spec.molecular_search_settings = self.parameters.ms1_molecular_search
                        
                        self.add_peak(mass_spec, original_indexes, eic_data, possible_molform=possible_mf)
                        # create data dependent chromapeak here
                        # add mass spectrum and eic data 
                        # lcms will return chromapeaks and check on if if is already existing to add new molecular formulas
                        # lcms needs to be dict based to be able to perform scan checks
                        
                        #mass_spec.plot_mz_domain_profile()
                        #plt.show()
        
           