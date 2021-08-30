__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"

from corems.mass_spectra.calc.LC_Calc import LC_Calculations
import numpy as np
import sys
import site
from pathlib import Path

import clr
import pandas as pd
from s3path import S3Path
from tqdm import tqdm

from typing import List
from corems.encapsulation.constant import Labels
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

        self._start_scan = self.iRawDataPlus.RunHeaderEx.FirstSpectrum

        self._end_scan = self.iRawDataPlus.RunHeaderEx.LastSpectrum

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
        return self._start_scan

    @property
    def end_scan(self):
        return self._end_scan

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
        scanrange = range(self._start_scan, self._end_scan + 1)
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
        

    def get_eics(self, target_mzs: list, ppm_tolerance=5,
                 start_scan=-1, end_scan=-1, ms_type='MS', plot=False, ax=None):

        '''ms_type: str ('MS', MS2')
        start_scan: int default -1 will select the lowest available
        end_scan: int default -1 will select the highest available

        returns:

            chroma: dict{target_mz: dict{
                                        Scan: [int]
                                            original thermo scan number
                                        Time: [floats]
                                            list of retention times
                                        TIC: [floats]
                                            total ion chromatogram
                                        }
        '''

        options = MassOptions()
        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = ppm_tolerance

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
                                                     start_scan, end_scan, options)

        traces = ChromatogramSignal.FromChromatogramData(data)

        chroma = {}    
       
        if plot:
            from matplotlib.transforms import Bbox
            import matplotlib.pyplot as plt
            if not ax:
                # ax = plt.gca()
                # ax.clear()
                fig, ax = plt.subplots()
                ax.set_prop_cycle(color=plt.cm.gist_rainbow(np.linspace(0, 1, len(traces))))
            else:
                fig = plt.gcf()    

            ax.set_xlabel('Time (min)')
            ax.set_ylabel('a.u.')
            ax.set_title(ms_type + ' EIC')
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.axes.spines['top'].set_visible(False)
            ax.axes.spines['right'].set_visible(False)

            legend = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 0, 0.07, 1))
            fig.subplots_adjust(right=0.76)
            d = {"down": 30, "up": -30}

            def func(evt):
                if legend.contains(evt):
                    bbox = legend.get_bbox_to_anchor()
                    bbox = Bbox.from_bounds(bbox.x0, bbox.y0 + d[evt.button], bbox.width, bbox.height)
                    tr = legend.axes.transAxes.inverted()
                    legend.set_bbox_to_anchor(bbox.transformed(tr))
                    fig.canvas.draw_idle()

            fig.canvas.mpl_connect("scroll_event", func)        
            # plt.show()
        
        for i, trace in enumerate(traces):
            if trace.Length > 0:
                rt, ic, scans  = self.get_rt_time_from_trace(trace)
                chroma[target_mzs[i]] = {'Scans' : scans, 'Time': rt, 'ECI': ic}
                if plot:
                    ax.plot(rt, ic, label="{:.5f}".format(target_mzs[i]))    
                    
        if plot:
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

    def get_tic(self, start_scan=-1, end_scan=-1, ms_type='MS', smooth=False, plot=False, ax=None) -> dict:

        '''ms_type: str ('MS', MS2')
        start_scan: int default -1 will select the lowest available
        end_scan: int default -1 will select the highest available
        returns:
            chroma: dict
            {
            Scan: [int]
                original thermo scan number
            Time: [floats]
                list of retention times
            TIC: [floats]
                total ion chromatogram
            }
        '''

        settings = ChromatogramTraceSettings(TraceType.TIC)
        settings.Filter = ms_type

        chroma_settings = IChromatogramSettings(settings)

        data = self.iRawDataPlus.GetChromatogramData([chroma_settings],
                                                     start_scan, end_scan)

        trace = ChromatogramSignal.FromChromatogramData(data)
        rt = []
        tic = []
        data = {'Time': [], 'Scan': [], 'TIC': []}

        if trace[0].Length > 0:

            for i in range(trace[0].Length):
                # print(trace[0].HasBasePeakData,trace[0].EndTime )

                # print("  {} - {}, {}".format( i, trace[0].Times[i], trace[0].Intensities[i] ))
                data['Time'].append(trace[0].Times[i])
                data['TIC'].append(trace[0].Intensities[i])
                data['Scan'].append(trace[0].Scans[i])

            chroma = pd.DataFrame(data)
            if smooth:
                
                chroma['TIC']= self.smooth_tic(data['TIC'])

            if plot:
                if not ax:
                    import matplotlib.pyplot as plt
                    ax = plt.gca()
                    # fig, ax = plt.subplots(figsize=(6, 3))

                ax.plot(chroma['Time'], chroma['TIC'], label=ms_type + ' TIC')    
                ax.set_xlabel('Time (min)')
                ax.set_ylabel('a.u.')
                plt.legend()
                # plt.show()
                return chroma, ax
            
            return chroma, None

        else:
            return None, None

class ImportDataDependentThermoMSFileReader(ThermoBaseClass, LC_Calculations):

    '''  Collection of methdos to import LC data dependent acquisition from Thermo's raw file
         Intended do create the LCMS object --> ChromaPeaks --> MSobj FullScan --> Dependent MS/MS Obj
    '''

    def __init__(self, file_location: str, start_scan: int = -1, end_scan: int = -1,
                 selected_mzs: List[float] = None, enforce_target_ms2: bool = True, 
                 eic_tolerance_ppm: float = 5.0):
        '''
        target_mzs: list[float] monoisotopic target m/z  or None
            Details: None will defalt to depends scans selected m/
        file_location: str, Path, or S3Path
        enforce_target_ms2: bool
            only perform EIC for target_mz if the m/z was selected as precursor for ms2
        start_scan: int
            default -1 will select the lowest available
        end_scan: int
            default -1 will select the highest available
        '''
        super().__init__(file_location)

        self._selected_mzs = self._init_target_mz(selected_mzs, enforce_target_ms2, eic_tolerance_ppm)

    @property
    def selected_mzs(self) -> List[float]:
        return list(self._selected_mzs)

    def get_precursors_list(self, precision_decimals=6):
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

    def _init_target_mz(self, selected_mzs: List[float], enforce_target_ms2: List[float], tolerance_ppm: float):

        precursors_mzs = self.get_precursors_list()

        if selected_mzs is None:
            # no selected m/z list provided, default to use the precursos m/z

            return precursors_mzs

        elif selected_mzs and enforce_target_ms2 is False:
            # selected m/z list provided, and not enforcing being selected as precursor
            return selected_mzs

        elif selected_mzs and enforce_target_ms2:
            # search the selected m/z list in the precursors m/z with a ms/ms experiment

            searchmz = MZSearch(precursors_mzs, selected_mzs, tolerance_ppm)
            searchmz.start()
            searchmz.join()
            return searchmz.results.keys()

class ImportMassSpectraThermoMSFileReader(ThermoBaseClass):

    '''  Collection of methdos to import Summed/Averaged mass spectrum from Thermo's raw file
         Currently only for profile mode data
         Returns MassSpecProfile object
    '''

    def get_icr_transient_times(self, first_scan: int = None, last_scan: int = None):
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

        firstScanNumber = self._start_scan if first_scan is None else first_scan

        lastScanNumber = self._end_scan if last_scan is None else last_scan

        transient_time_list = []

        for scan in range(firstScanNumber, lastScanNumber):

            scan_header = self.get_scan_header(scan)

            rp_target = scan_header['FT Resolution:']

            transient_time = res_trans_time.get(rp_target)

            transient_time_list.append(transient_time)

            # print(transient_time, rp_target)

        return transient_time_list

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

    def set_metadata(self, firstScanNumber=0, lastScanNumber=0, scans_list=False):
        '''
        Collect metadata to be ingested in the mass spectrum object

        scans_list: list[int] or false
        lastScanNumber: int
        firstScanNumber: int
        '''

        d_params = default_parameters(self.file_path)

        # assumes scans is full scan or reduced profile scan

        d_params['label'] = Labels.thermo_profile

        if scans_list:
            d_params['scan_number'] = scans_list

            d_params['polarity'] = self.get_polarity_mode(scans_list[0])

        else:

            d_params['scan_number'] = '{}-{}'.format(firstScanNumber, lastScanNumber)

            d_params['polarity'] = self.get_polarity_mode(firstScanNumber)

        d_params['analyzer'] = self.iRawDataPlus.GetInstrumentData().Model

        d_params['instrument_label'] = self.iRawDataPlus.GetInstrumentData().Name

        return d_params

    def get_average_mass_spectrum_by_scanlist(self, scans_list: List[int], auto_process: bool = True,
                                              ppm_tolerance: float = 5.0):

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

    def get_average_mass_spectrum_in_scan_range(self, first_scan: int = None, last_scan: int = None,
                                                auto_process: bool = True, ppm_tolerance: float = 5.0,
                                                ms_type: str = 0):

        '''
        Averages mass spectra over a scan range using Thermo's AverageScansInScanRange method
        first_scan: int
        last_scan: int
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object
        ms_type: MSOrderType.MS
            Type of mass spectrum scan, default for full scan acquisition
         Returns:
            MassSpecProfile
        '''

        firstScanNumber = self._start_scan if first_scan is None else first_scan

        lastScanNumber = self._end_scan if last_scan is None else last_scan

        d_params = self.set_metadata(firstScanNumber=firstScanNumber, lastScanNumber=lastScanNumber)

        # Create the mass options object that will be used when averaging the scans
        options = MassOptions()

        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = ppm_tolerance

        # Get the scan filter for the first scan.  This scan filter will be used to located
        # scans within the given scan range of the same type
        scanFilter = self.iRawDataPlus.GetFilterForScanNumber(firstScanNumber)

        # force it to only look for the MSType
        scanFilter.MSOrder = ms_type

        averageScan = Extensions.AverageScansInScanRange(self.iRawDataPlus, firstScanNumber, lastScanNumber,
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

    def get_summed_mass_spectrum(self, start_scan: int, end_scan: int = None,
                                 auto_process=True, pd_method=True, pd_merge_n=100):

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

        if type(start_scan) is list:
            d_params['polarity'] = self.get_polarity_mode(start_scan[0])

            scanrange = start_scan
        else:
            d_params['polarity'] = self.get_polarity_mode(start_scan)

            if end_scan is None:
                end_scan = self._end_scan

            scanrange = range(start_scan, end_scan + 1)

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
