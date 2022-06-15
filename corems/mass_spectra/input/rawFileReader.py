__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"


from dataclasses import dataclass, field
from warnings import warn

from matplotlib import axes

from corems.mass_spectra.calc.LC_Calc import LC_Calculations
import numpy as np
import sys
import site
from pathlib import Path

import clr
import pandas as pd
from s3path import S3Path
from tqdm import tqdm

from typing import Dict, List, Tuple
from corems.encapsulation.constant import Labels
from corems.mass_spectra.factory.LC_Class import DataDependentLCMS
from corems.mass_spectra.factory.LC_Temp import EIC_Data, TIC_Data
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.mass_spectra.calc.MZSearch import MZSearch
from corems.encapsulation.factory.parameters import LCMSParameters, default_parameters


# do not change the order from the imports statements and reference ThermoFisher below
sys.path.append(site.getsitepackages()[0] + '/ext_lib')
sys.path.append('ext_lib')

clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.Data')
clr.AddReference('ThermoFisher.CommonCore.MassPrecisionEstimator')

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data import ToleranceUnits, Extensions
from ThermoFisher.CommonCore.Data.Business import ChromatogramTraceSettings, TraceType, MassOptions
from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, Range
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings
from ThermoFisher.CommonCore.Data.FilterEnums import MSOrderType
from System.Collections.Generic import List

class ThermoBaseClass():

    def __init__(self, file_location):
        ''' file_location: srt pathlib.Path or s3path.S3Path
                Thermo Raw file path
        '''
        # Thread.__init__(self)
        if isinstance(file_location, str):
            file_path = Path(file_location)

        elif isinstance(file_location, S3Path):

            temp_dir = Path('tmp/')
            temp_dir.mkdir(exist_ok=True)

            file_path = temp_dir / file_location.name
            with open(file_path, 'wb') as fh:
                fh.write(file_location.read_bytes())

        else:
            file_path = file_location

        self.iRawDataPlus = RawFileReaderAdapter.FileFactory(str(file_path))

        self.res = self.iRawDataPlus.SelectInstrument(0, 1)
        
        self.file_path = file_location

        self._init_settings()

    def _init_settings(self):

        self._parameters = LCMSParameters()

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, instance_LCMSParameters):
        self._parameters = instance_LCMSParameters

    @property
    def chromatogram_settings(self): return self.parameters.lc_ms

    @chromatogram_settings.setter
    def chromatogram_settings(self, instance_LiquidChromatographSetting):

        self.parameters.lc_ms =  instance_LiquidChromatographSetting

    @property
    def start_scan(self):
        
        if self.chromatogram_settings.start_scan == -1:
            return self.iRawDataPlus.RunHeaderEx.FirstSpectrum
        else:
           return self.chromatogram_settings.start_scan
      

    @property
    def end_scan(self):
          
        if self.chromatogram_settings.end_scan == -1:
            return self.iRawDataPlus.RunHeaderEx.LastSpectrum   
        else:        
            return self.chromatogram_settings.end_scan
        

    def remove_temp_file(self):
        '''if the path is from S3Path data cannot be serialized to io.ByteStream and
           a temporary copy is stored at the temp dir
           use this function only at the end of your execution scrip
           some LCMS class methods depend on this file
        '''

        self.file_path.unlink()

    def get_polarity_mode(self, scan_number: int):

        polarity_symbol = self.get_filter_for_scan_num(scan_number)[1]

        if polarity_symbol == '+':

            return 1
            # return 'POSITIVE_ION_MODE'

        elif polarity_symbol == '-':

            return -1

        else:

            raise Exception('Polarity Mode Unknown, please set it manually')

    def get_filter_for_scan_num(self, scan_number: int):
        '''
        Returns the closest matching run time that corresponds to scan_number for the current
        controller. This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']
        '''
        scan_label = self.iRawDataPlus.GetScanEventStringForScanNumber(
            scan_number)

        return str(scan_label).split()

    def check_full_scan(self, scan_number: int):
        # scan_filter.ScanMode 0 = FULL
        scan_filter = self.iRawDataPlus.GetFilterForScanNumber(scan_number)

        return scan_filter.ScanMode == MSOrderType.Ms

    def get_all_filters(self):
        '''
        Get all scan filters.
        This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']
        '''
        
        scanrange = range(self.start_scan, self.end_scan + 1)
        scanfiltersdic = {}
        scanfilterslist = []
        for scan_number in scanrange:
            scan_label = self.iRawDataPlus.GetScanEventStringForScanNumber(scan_number)
            scanfiltersdic[scan_number] = scan_label
            scanfilterslist.append(scan_label)
        scanfilterset = list(set(scanfilterslist))
        return scanfiltersdic, scanfilterset

    def get_scan_header(self, scan: int):
        '''
        Get full dictionary of scan header meta data, i.e. AGC status, ion injection time, etc.
        '''
        header = self.iRawDataPlus.GetTrailerExtraInformation(scan)
        header_dic = {}
        for i in np.arange(header.Length):
            header_dic.update({header.Labels[i]: header.Values[i]})
        return header_dic

    @staticmethod
    def get_rt_time_from_trace(trace) -> List[List[float]]:
        '''trace: ThermoFisher.CommonCore.Data.Business.ChromatogramSignal'''
        return list(trace.Times), list(trace.Intensities), list(trace.Scans)
        

    def get_eics(self, target_mzs: list, tic_data: dict, ms_type='MS !d', 
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

        data = self.iRawDataPlus.GetChromatogramData(all_chroma_settings,
                                                     self.start_scan, self.end_scan, options)

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
                rt, eic, scans  = self.get_rt_time_from_trace(trace)
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
                    warn("The software assumes same lenth of TIC and EIC, this does not seems to be the case and the results mass spectrum selected by the scan number might not be correct")
                
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

        data = self.iRawDataPlus.GetChromatogramData([chroma_settings],
                                                     self.start_scan, self.end_scan)

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

    def get_average_mass_spectrum_by_scanlist(self, scans_list: List[int], auto_process: bool = True,
                                              ppm_tolerance: float = 5.0) -> MassSpecProfile:

        '''
        Averages selected scans mass spectra using Thermo's AverageScans method
        scans_list: list[int]
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object
        Returns:
            MassSpecProfile
        '''

        """
        Averages selected scans mass spectra using Thermo's AverageScans method
        scans_list: list[int]
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object
        Returns:
            MassSpecProfile
        """

        d_params = self.set_metadata(scans_list=scans_list)

        # assumes scans is full scan or reduced profile scan

        scans = List[int]()
        for scan in scans_list:
            scans.Add(scan)

        # Create the mass options object that will be used when averaging the scans
        options = MassOptions()
        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = ppm_tolerance

        # Get the scan filter for the first scan.  This scan filter will be used to located
        # scans within the given scan range of the same type

        averageScan = Extensions.AverageScans(self.iRawDataPlus, scans, options)

        len_data = averageScan.SegmentedScan.Positions.Length

        mz_list = list(averageScan.SegmentedScan.Positions)
        abund_list = list(averageScan.SegmentedScan.Intensities)

        data_dict = {Labels.mz: mz_list,
                     Labels.abundance: abund_list,
                     }

        mass_spec = MassSpecProfile(data_dict, d_params, auto_process=auto_process)

        return mass_spec

    def get_centroid_msms_data(self, scan):
        
        d_params = self.set_metadata(scans_list=[scan], label=Labels.thermo_centroid)

        centroidStream = self.iRawDataPlus.GetCentroidStream(scan, False)

        noise = list(centroidStream.Noises)

        baselines = list(centroidStream.Baselines)

        rp = list(centroidStream.Resolutions)

        magnitude = list(centroidStream.Intensities)

        mz = list(centroidStream.Masses)

        # charge = scans_labels[5]
        array_noise_std = (np.array(noise) - np.array(baselines)) / 3
        l_signal_to_noise = np.array(magnitude) / array_noise_std

        d_params['baselise_noise'] = np.average(array_noise_std)

        d_params['baselise_noise_std'] = np.std(array_noise_std)

        data_dict = {
            Labels.mz: mz,
            Labels.abundance: magnitude,
            Labels.rp: rp,
            Labels.s2n: list(l_signal_to_noise)
        }

        
        mass_spec = MassSpecCentroid(data_dict, d_params, auto_process=False)
        mass_spec.settings.threshold_method = 'relative_abundance'
        mass_spec.settings.relative_abundance_threshold = 1
        mass_spec.process_mass_spec()
        return mass_spec


    def get_average_mass_spectrum_in_scan_range(self, auto_process: bool = True, ppm_tolerance: float = 5.0,
                                                ms_type: int = 0) -> MassSpecProfile:

        '''
        Averages mass spectra over a scan range using Thermo's AverageScansInScanRange method
        start_scan: int
        end_scan: int
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object
        ms_type: MSOrderType.MS
            Type of mass spectrum scan, default for full scan acquisition
         Returns:
            MassSpecProfile
        '''

        d_params = self.set_metadata(firstScanNumber=self.start_scan,
                                     lastScanNumber=self.end_scan)

        # Create the mass options object that will be used when averaging the scans
        options = MassOptions()
        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = ppm_tolerance

        # Get the scan filter for the first scan.  This scan filter will be used to located
        # scans within the given scan range of the same type
        scanFilter = self.iRawDataPlus.GetFilterForScanNumber(self.start_scan)

        # force it to only look for the MSType
        scanFilter.MSOrder = ms_type

        averageScan = Extensions.AverageScansInScanRange(self.iRawDataPlus, 
                                                         self.start_scan, self.end_scan,
                                                         scanFilter, options)
        
        if averageScan:
            
            mz_list = list(averageScan.SegmentedScan.Positions)
            abund_list = list(averageScan.SegmentedScan.Intensities)

            data_dict = {Labels.mz: mz_list,
                         Labels.abundance: abund_list,
                         }

            
            mass_spec = MassSpecProfile(data_dict, d_params, auto_process=auto_process)

            return mass_spec
        else:
            raise Exception('no data found for the MSOrderType = {}'.format(ms_type))

   
    def set_metadata(self, firstScanNumber=0, lastScanNumber=0, scans_list=False, label=Labels.thermo_profile):
        '''
        Collect metadata to be ingested in the mass spectrum object

        scans_list: list[int] or false
        lastScanNumber: int
        firstScanNumber: int
        '''

        d_params = default_parameters(self.file_path)

        # assumes scans is full scan or reduced profile scan

        d_params['label'] = label

        if scans_list:
            d_params['scan_number'] = scans_list

            d_params['polarity'] = self.get_polarity_mode(scans_list[0])

        else:

            d_params['scan_number'] = '{}-{}'.format(firstScanNumber, lastScanNumber)

            d_params['polarity'] = self.get_polarity_mode(firstScanNumber)

        d_params['analyzer'] = self.iRawDataPlus.GetInstrumentData().Model

        d_params['instrument_label'] = self.iRawDataPlus.GetInstrumentData().Name

        return d_params    

class ImportDataDependentThermoMSFileReader(ThermoBaseClass, LC_Calculations):

    '''  Collection of methdos to import LC data dependent acquisition from Thermo's raw file
         Intended do create the LCMS object --> ChromaPeaks --> MSobj FullScan --> Dependent MS/MS Obj
    '''

    def __init__(self, file_location: str, selected_mzs: List[float] = None):
        '''
        target_mzs: list[float] monoisotopic target m/z  or None
            Details: None will defalt to depends scans selected m/
        file_location: str, Path, or S3Path
        
        '''
        super().__init__(file_location)
        
        eic_tolerance_ppm = self.chromatogram_settings.eic_tolerance_ppm
        enforce_target_ms2 = self.chromatogram_settings.enforce_target_ms2
        average_target_mz = self.chromatogram_settings.average_target_mz
        
        print('TOLERANCE = {} ppm'.format(eic_tolerance_ppm))
        self._selected_mzs = self._init_target_mz(selected_mzs, enforce_target_ms2, 
                                                 eic_tolerance_ppm, average_target_mz)

        self.lcms = DataDependentLCMS(file_location, self._selected_mzs, self)

    @property
    def selected_mzs(self) -> List[float]:
        return list(self._selected_mzs)

    def get_lcms_obj(self):
        
        return self.lcms

    def get_precursors_list(self, precision_decimals=5):
        '''returns a set of unique precursors m/z
        precision_decimals: int
            change this parameters does not seem to affect the number of dependent scans selected
            needs more investigation
        '''

        precursors_mzs = set()
        
        for scan in range(self.start_scan, self.end_scan):

            scan_filter = self.iRawDataPlus.GetFilterForScanNumber(scan)

            MSOrder = scan_filter.MSOrder
            
            if MSOrder == MSOrderType.Ms:

                scanDependents = self.iRawDataPlus.GetScanDependents(scan, precision_decimals)

                for scan_dependent_detail in scanDependents.ScanDependentDetailArray:

                    for precursor_mz in scan_dependent_detail.PrecursorMassArray:
                        
                        precursors_mzs.add(precursor_mz)
            
        return precursors_mzs

    def _init_target_mz(self, selected_mzs: List[float], enforce_target_ms2: bool, 
                        tolerance_ppm: float, average_target_mz: bool):

        precursors_mzs = self.get_precursors_list()
        
        if selected_mzs is None:
            # no selected m/z list provided, default to use the precursos m/z
            if average_target_mz:
                searchmz = MZSearch(precursors_mzs, precursors_mzs, tolerance_ppm, average_target_mz=average_target_mz)
                return searchmz.averaged_target_mz
            else:
                return precursors_mzs
            #searchmz.start()
            #searchmz.join()


        elif selected_mzs and enforce_target_ms2 is False:
            # selected m/z list provided, and not enforcing being selected as precursor
            return selected_mzs

        elif selected_mzs and enforce_target_ms2:
            # search the selected m/z list in the precursors m/z with a ms/ms experiment
            print("YEAHHHHH")
            searchmz = MZSearch(precursors_mzs, selected_mzs, tolerance_ppm, average_target_mz=average_target_mz)
            searchmz.start()
            searchmz.join()
            return sorted(searchmz.results.keys())

class ImportMassSpectraThermoMSFileReader(ThermoBaseClass, LC_Calculations):

    '''  Collection of methdos to import Summed/Averaged mass spectrum from Thermo's raw file
         Currently only for profile mode data
         Returns MassSpecProfile object
    '''

    def get_icr_transient_times(self):
        '''
        Return a list for transient time targets for all scans, or selected scans range
        Resolving Power and Transient time targets based on 7T FT-ICR MS system
        '''
        
        res_trans_time = {'50': 0.384,
                          '100000': 0.768,
                          '200000': 1.536,
                          '400000': 3.072,
                          '750000': 6.144,
                          '1000000': 12.288
                          }

        firstScanNumber = self.start_scan

        lastScanNumber = self.end_scan

        transient_time_list = []

        for scan in range(firstScanNumber, lastScanNumber):

            scan_header = self.get_scan_header(scan)

            rp_target = scan_header['FT Resolution:']

            transient_time = res_trans_time.get(rp_target)

            transient_time_list.append(transient_time)

            # print(transient_time, rp_target)

        return transient_time_list

    def get_summed_mass_spectrum(self, auto_process=True, pd_method=True, pd_merge_n=100) -> MassSpecProfile:

        '''
        Manually sum mass spectrum over a scan range
        start_scan: int
        end_scan: int
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object 
        pd_method: bool
            If true uses pandas to align and sum data
            Else: Assumes data is aligned and sum each data point across all mass spectra
        Returns:
            MassSpecProfile
        '''
        
        d_params = default_parameters(self.file_path)

        # assumes scans is full scan or reduced profile scan

        d_params['label'] = Labels.thermo_profile

        if type(self.start_scan) is list:
            d_params['polarity'] = self.get_polarity_mode(self.start_scan[0])

            scanrange = self.start_scan
        else:
            d_params['polarity'] = self.get_polarity_mode(self.start_scan)

            scanrange = range(self.start_scan, self.end_scan + 1)

        if pd_method:

            def sort_sum_df(df):
                '''
                Nested function to sort dataframe and sum rows with exact matching indexes (m/z)
                '''
                df = df.sort_index()
                df = df.groupby(level=0).sum()
                return df

            # initialise empty Pandas series
            big_df = pd.Series(index=[], dtype='float64')

            for scan_number in tqdm(scanrange):
                scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(scan_number)
                segmentedScan = self.iRawDataPlus.GetSegmentedScanFromScanNumber(scan_number, scanStatistics)

                tmp_df = pd.Series(index=list(segmentedScan.Positions),
                                   dtype='float64', data=list(segmentedScan.Intensities))
                big_df = big_df.append(tmp_df)

                # this allows you to merge/sum the values earlier, however it slows down a lot
                # limited benefit unless running into memory issues
                # for complex data it is necessary to stop the iterations getting too slow
                if scan_number % pd_merge_n == 0:
                    big_df = sort_sum_df(big_df)

            big_df = sort_sum_df(big_df)
            data_dict = {Labels.mz: list(big_df.index.values),
                         Labels.abundance: list(big_df.values),
                         }
        else:
            all_mz = dict()

            for scan_number in tqdm(scanrange):

                scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(scan_number)

                segmentedScan = self.iRawDataPlus.GetSegmentedScanFromScanNumber(scan_number, scanStatistics)

                len_data = segmentedScan.Positions.Length

                for i in range(len_data):

                    mz = segmentedScan.Positions[i]
                    abundance = segmentedScan.Intensities[i]

                    if mz in all_mz:
                        all_mz[mz] = all_mz[mz] + abundance
                    else:
                        all_mz[mz] = abundance

            mz_all = []
            abun_all = []

            for mz in sorted(all_mz):
                
                mz_all.append(mz)
                abun_all.append(all_mz[mz])

            data_dict = {Labels.mz: mz_all,
                         Labels.abundance: abun_all,
                         }

        print('Summed. Now Processing.')

        mass_spec = MassSpecProfile(data_dict, d_params, auto_process=auto_process)

        return mass_spec

    def get_data(self, scan: int, d_parameter: dict, scan_type: str):

        if scan_type == 'Centroid':

            centroidStream = self.iRawDataPlus.GetCentroidStream(scan, False)

            noise = list(centroidStream.Noises)

            baselines = list(centroidStream.Baselines)

            rp = list(centroidStream.Resolutions)

            magnitude = list(centroidStream.Intensities)

            mz = list(centroidStream.Masses)

            # charge = scans_labels[5]
            array_noise_std = (np.array(noise) - np.array(baselines)) / 3
            l_signal_to_noise = np.array(magnitude) / array_noise_std

            d_parameter['baselise_noise'] = np.average(array_noise_std)

            d_parameter['baselise_noise_std'] = np.std(array_noise_std)

            data_dict = {
                Labels.mz: mz,
                Labels.abundance: magnitude,
                Labels.rp: rp,
                Labels.s2n: l_signal_to_noise,
            }

        else:

            scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(scan)

            profileStream = self.iRawDataPlus.GetSegmentedScanFromScanNumber(
                scan, scanStatistics)

            magnitude = list(profileStream.Intensities)

            mz = list(profileStream.Positions)

            data_dict = {
                Labels.mz: mz,
                Labels.abundance: magnitude,
            }

        return data_dict

   

    
    def get_best_scans_idx(self, stdevs=2, method='mean', plot=False):
        '''
        Method to determine the best scan indexes for selective co-addition
        Based on calculating the mean (default) of the TIC values
        and setting an upper limit above/below that within X standard deviations.
        Mean or median makes limited difference, it seems.
        Empirically, 1-2 stdevs enough to filter out the worst datapoints.
        Optionally, plot the TIC with horizontal lines for the standard dev cutoffs.
        '''
        tic = pd.Dataframe(self.get_tic(plot=plot))

        if method == 'median':
            tic_median = tic['TIC'].median()
        elif method == 'mean':
            tic_median = tic['TIC'].mean()
        else:
            print('Method ' + str(method) + ' undefined')

        tic_std = tic['TIC'].std()

        upperlimit = tic_median - (stdevs * tic_std)
        lowerlimit = tic_median + (stdevs * tic_std)

        tic_filtered = tic[(tic['TIC'] > upperlimit) & (tic['TIC'] < lowerlimit)]
        scans = list(tic_filtered.Scans.values)

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(tic['Time'], tic['TIC'])
            ax.axhline(y=upperlimit, c='r')
            ax.axhline(y=lowerlimit, c='r')
            return fig, scans
        else:
            return scans

    def isotopehunter(self,pattern,timerange,mass_tolerance,ratio_tolerance,peakwidth,correlation,slope_filter):
        #Function matches required ('Y') peaks in pattern to spectra in self
        #Requires 'timerange','mass_tolerance','ratio_tolerance' for pattern matching
        #Requires 'pattern','peakwidth','correlation','slope_filter' for QC filtering

        #Define pattern boundaries
        umass=pattern.mdiff[pattern.requirement=='Y']+mass_tolerance
        lmass=pattern.mdiff[pattern.requirement=='Y']-mass_tolerance
        uratio=pattern.ratio[pattern.requirement=='Y']*ratio_tolerance
        lratio=pattern.ratio[pattern.requirement=='Y']/ratio_tolerance
        nisotope=len(umass)

        interval=peakwidth

        times=np.arange(timerange[0],timerange[1],interval).tolist()

        #Retrieve TIC for MS1 scans only within timerange
        tic=self.get_tic(ms_type='MS')[0]
        tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
        #scans=tic_df[tic_df.time.between(timerange[0],timerange[1])].scan.tolist()

        #Create  empty results dictionaries. Will be filed with {Scan: {1st peak:{mz,intense}, 2nd peak:{mz,intense}, npeak:{mz,intense}}}
        results=[]
        clean_results=[]
        final_results=[]

        #Determine 2 most abundant isotopologues for correlation analysis.  
        req_isotopes=pattern[pattern.requirement=='Y'].isotope
        isotope1=req_isotopes[0]
        isotope2=req_isotopes[1]

        #for s in scans:
        for timestart in times:
            scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
            s=scans[0]
            print('Scan:'+str(s)+' Time (min): '+ str(tic_df.time[tic_df.scan==s].round(2).max()))
            ms=self.get_average_mass_spectrum_by_scanlist(scans)
            if(ms.mspeaks):

                spectrum=pd.DataFrame({'mz':ms.mz_exp, 'intense':ms.abundance})
                print(len(spectrum))

                for j in spectrum.index:
                    k=1
                    hitlist=[]
                    result={}

                    while(k<nisotope):
                        hits=spectrum[(spectrum.mz > spectrum.mz[j]+lmass[k]) & (spectrum.mz < spectrum.mz[j]+umass[k]) & (spectrum.intense > spectrum.intense[j]*lratio[k]) & (spectrum.intense < spectrum.intense[j]*uratio[k])].index.tolist()
                        if hits:
                            hitlist.append(hits)
                            k=k+1
                        else:
                            k=nisotope+1
                    if k==(nisotope):
                        result['scan']=scans
                        result['time']=tic_df[tic_df.scan==s].time.iloc[0]

                        result[pattern.isotope[0]]={'mz':spectrum.mz[j],'intense':spectrum.intense[j]}

                        for i, iso in enumerate(req_isotopes[1:]):
                            result[iso]={'mz':spectrum.mz[hitlist[i][0]],'intense':spectrum.intense[hitlist[i][0]]}
                        mass1=result[isotope1]['mz']
                        mass2=result[isotope2]['mz']

                        EIC=self.get_eics(target_mzs=[mass1,mass2],tic_data={},peak_detection=False,smooth=False)
                        df=pd.DataFrame({'mz1':EIC[0][mass1].eic,'mz2':EIC[0][mass2].eic,'time':EIC[0][mass1].time})
                        df_sub_a=df[df['time'].between(timestart,timestart+interval)]
                        peakmax_i=df_sub_a.mz1.max()
                        peakmax_t=df_sub_a.time[df_sub_a.mz1==peakmax_i].max()
                        
                        df_sub=df[df['time'].between(peakmax_t-interval,peakmax_t+interval)]

                        #Calculate correlation and slope between two main isotopologues.        
                        corr=df_sub.corr(method='pearson').iat[0,1]**2
                        slope=np.polyfit(df_sub.mz1,df_sub.mz2,1)[0]/(pattern.sort_values(by='ratio',ascending=False).ratio[1]/pattern.sort_values(by='ratio',ascending=False).ratio[0])

                        result['corr']=corr 
                        result['file']=self.file_path
                        result['slope']=slope
                        result['mass']=round(mass1,3)
                        result['abundance']=result[isotope1]['intense']
                        result['time_peak']=peakmax_t
                        result['abundance_peak']=peakmax_i
                        result['qc']='match'
                        result['dmz']=result[isotope2]['mz']-result[isotope1]['mz']

                        if corr>correlation:
                            if ((slope > slope_filter[0]) & (slope < slope_filter[1])):
                                clean_results.append(result)
                                result['qc']='qc'

                        results.append(result)

        clean_results_df=pd.DataFrame(clean_results)

        for result in clean_results:
            masses=clean_results_df[(abs(clean_results_df.mass-result['mass']) < mass_tolerance)& (abs(clean_results_df.time - result['time']) < peakwidth*2)]
            max_value=max(masses['abundance_peak'])

            if (result['abundance_peak']==max_value):
                final_results.append(result)
                result['qc']='max'
                
        print("Results:")
        print(len(results))
        print("Final Results:")
        print(len(final_results))

        return(results,clean_results,final_results)   