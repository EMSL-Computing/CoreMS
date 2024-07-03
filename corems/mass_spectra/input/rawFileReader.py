__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"


from warnings import warn
import warnings

from matplotlib import axes
from corems.encapsulation.factory.processingSetting import LiquidChromatographSetting

from corems.mass_spectra.calc.LC_Calc import LC_Calculations
import numpy as np
import sys
import site
from pathlib import Path
import datetime

import clr
import pandas as pd
from s3path import S3Path


from typing import Any, Dict, List, Optional, Tuple
from corems.encapsulation.constant import Labels
from corems.mass_spectra.factory.LC_Class import DataDependentLCMS
from corems.mass_spectra.factory.LC_Temp import EIC_Data, TIC_Data
from corems.mass_spectrum.factory.MassSpectrumClasses import (
    MassSpecProfile,
    MassSpecCentroid,
)
from corems.mass_spectra.calc.MZSearch import MZSearch
from corems.encapsulation.factory.parameters import LCMSParameters, default_parameters


# do not change the order from the imports statements and reference ThermoFisher below
sys.path.append(site.getsitepackages()[0] + "/ext_lib/dotnet/")
sys.path.append("ext_lib/dotnet/")

clr.AddReference("ThermoFisher.CommonCore.RawFileReader")
clr.AddReference("ThermoFisher.CommonCore.Data")
clr.AddReference("ThermoFisher.CommonCore.MassPrecisionEstimator")

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data import ToleranceUnits, Extensions
from ThermoFisher.CommonCore.Data.Business import (
    ChromatogramTraceSettings,
    TraceType,
    MassOptions,
)
from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, Range
from ThermoFisher.CommonCore.Data.Business import Device
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings
from ThermoFisher.CommonCore.Data.Business import MassOptions, FileHeaderReaderFactory
from ThermoFisher.CommonCore.Data.FilterEnums import MSOrderType
from System.Collections.Generic import List


class ThermoBaseClass:
    """
    Class for reading Thermo Raw files and extracting information from them.

    Parameters:
    -----------
    file_location : str or pathlib.Path or s3path.S3Path
        Thermo Raw file path or S3 path.

    Attributes:
    -----------
    file_path : str or pathlib.Path or s3path.S3Path
        The file path of the Thermo Raw file.
    parameters : LCMSParameters
        The LCMS parameters for the Thermo Raw file.
    chromatogram_settings : LiquidChromatographSetting
        The chromatogram settings for the Thermo Raw file.
    scans : list or tuple
        The selected scans for the Thermo Raw file.
    start_scan : int
        The starting scan number for the Thermo Raw file.
    end_scan : int
        The ending scan number for the Thermo Raw file.

    Methods:
    --------
    * set_msordertype(scanFilter, mstype: str = 'ms1') -> scanFilter
        Convert the user-passed MS Type string to a Thermo MSOrderType object.
    * get_creation_time() -> datetime.datetime
        Extract the creation date stamp from the .RAW file and return it as a formatted datetime object.
    * remove_temp_file()
        Remove the temporary file if the path is from S3Path.
    * get_polarity_mode(scan_number: int) -> int
        Get the polarity mode for the given scan number.
    * get_filter_for_scan_num(scan_number: int) -> List[str]
        Get the filter for the given scan number.
    * check_full_scan(scan_number: int) -> bool
        Check if the given scan number is a full scan.
    * get_all_filters() -> Tuple[Dict[int, str], List[str]]
        Get all scan filters for the Thermo Raw file.
    * get_scan_header(scan: int) -> Dict[str, Any]
        Get the full dictionary of scan header metadata for the given scan number.
    * get_rt_time_from_trace(trace) -> Tuple[List[float], List[float], List[int]]
        Get the retention time, intensity, and scan number from the given trace.
    * get_eics(target_mzs: List[float], tic_data: Dict[str, Any], ms_type: str = 'MS !d',
             peak_detection: bool = True, smooth: bool = True, plot: bool = False,
             ax: Optional[matplotlib.axes.Axes] = None, legend: bool = False) -> Tuple[Dict[float, EIC_Data], matplotlib.axes.Axes]
        Get the extracted ion chromatograms (EICs) for the target m/z values.

    """

    def __init__(self, file_location):
        """file_location: srt pathlib.Path or s3path.S3Path
        Thermo Raw file path
        """
        # Thread.__init__(self)
        if isinstance(file_location, str):
            file_path = Path(file_location)

        elif isinstance(file_location, S3Path):
            temp_dir = Path("tmp/")
            temp_dir.mkdir(exist_ok=True)

            file_path = temp_dir / file_location.name
            with open(file_path, "wb") as fh:
                fh.write(file_location.read_bytes())

        else:
            file_path = file_location

        self.iRawDataPlus = RawFileReaderAdapter.FileFactory(str(file_path))
        print((self.iRawDataPlus).IsOpen)

        if not self.iRawDataPlus.IsOpen:
            raise FileNotFoundError(
                "Unable to access the RAW file using the RawFileReader class!"
            )

        # Check for any errors in the RAW file
        if self.iRawDataPlus.IsError:
            raise IOError(
                "Error opening ({}) - {}".format(self.iRawDataPlus.FileError, file_path)
            )

        self.res = self.iRawDataPlus.SelectInstrument(Device.MS, 1)

        self.file_path = file_location
        self.iFileHeader = FileHeaderReaderFactory.ReadFile(str(file_path))

        # removing tmp file

        self._init_settings()

    def _init_settings(self):
        """
        Initialize the LCMSParameters object.
        """
        self._parameters = LCMSParameters()

    @property
    def parameters(self) -> LCMSParameters:
        """
        Get or set the LCMSParameters object.
        """
        return self._parameters

    @parameters.setter
    def parameters(self, instance_LCMSParameters: LCMSParameters):
        self._parameters = instance_LCMSParameters

    @property
    def chromatogram_settings(self) -> LiquidChromatographSetting:
        """
        Get or set the LiquidChromatographSetting object.
        """
        return self.parameters.lc_ms

    @chromatogram_settings.setter
    def chromatogram_settings(
        self, instance_LiquidChromatographSetting: LiquidChromatographSetting
    ):
        self.parameters.lc_ms = instance_LiquidChromatographSetting

    @property
    def scans(self) -> list | tuple:
        """scans : list or tuple
        If list uses Thermo AverageScansInScanRange for selected scans, ortherwise uses Thermo AverageScans for a scan range
        """
        return self.chromatogram_settings.scans

    @property
    def start_scan(self) -> int:
        """
        Get the starting scan number for the Thermo Raw file.
        """
        if self.scans[0] == -1:
            return self.iRawDataPlus.RunHeaderEx.FirstSpectrum
        else:
            return self.scans[0]

    @property
    def end_scan(self) -> int:
        """
        Get the ending scan number for the Thermo Raw file.
        """
        if self.scans[-1] == -1:
            return self.iRawDataPlus.RunHeaderEx.LastSpectrum
        else:
            return self.scans[-1]

    def set_msordertype(self, scanFilter, mstype: str = "ms1"):
        """
        Function to convert user passed string MS Type to Thermo MSOrderType object
        Limited to MS1 through MS10.

        Parameters:
        -----------
        scanFilter : Thermo.ScanFilter
            The scan filter object.
        mstype : str, optional
            The MS Type string, by default 'ms1'

        """
        mstype = mstype.upper()
        # Check that a valid mstype is passed
        if (int(mstype.split("MS")[1]) > 10) or (int(mstype.split("MS")[1]) < 1):
            warn("MS Type not valid, must be between MS1 and MS10")

        msordertypedict = {
            "MS1": MSOrderType.Ms,
            "MS2": MSOrderType.Ms2,
            "MS3": MSOrderType.Ms3,
            "MS4": MSOrderType.Ms4,
            "MS5": MSOrderType.Ms5,
            "MS6": MSOrderType.Ms6,
            "MS7": MSOrderType.Ms7,
            "MS8": MSOrderType.Ms8,
            "MS9": MSOrderType.Ms9,
            "MS10": MSOrderType.Ms10,
        }
        scanFilter.MSOrder = msordertypedict[mstype]
        return scanFilter

    def get_creation_time(self) -> datetime.datetime:
        """
        Extract the creation date stamp from the .RAW file
        Return formatted creation date stamp.

        """
        credate = self.iRawDataPlus.CreationDate.get_Ticks()
        credate = datetime.datetime(1, 1, 1) + datetime.timedelta(
            microseconds=credate / 10
        )
        return credate

    def remove_temp_file(self) -> None:
        """if the path is from S3Path data cannot be serialized to io.ByteStream and
        a temporary copy is stored at the temp dir
        use this function only at the end of your execution scrip
        some LCMS class methods depend on this file
        """

        self.file_path.unlink()

    def get_polarity_mode(self, scan_number: int) -> int:
        """
        Get the polarity mode for the given scan number.

        Parameters:
        -----------
        scan_number : int
            The scan number.

        Raises:
        -------
        Exception
            If the polarity mode is unknown.

        """
        polarity_symbol = self.get_filter_for_scan_num(scan_number)[1]

        if polarity_symbol == "+":
            return 1
            # return 'POSITIVE_ION_MODE'

        elif polarity_symbol == "-":
            return -1

        else:
            raise Exception("Polarity Mode Unknown, please set it manually")

    def get_filter_for_scan_num(self, scan_number: int) -> List[str]:
        """
        Returns the closest matching run time that corresponds to scan_number for the current
        controller. This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']

        Parameters:
        -----------
        scan_number : int
            The scan number.

        """
        scan_label = self.iRawDataPlus.GetScanEventStringForScanNumber(scan_number)

        return str(scan_label).split()

    def check_full_scan(self, scan_number: int) -> bool:
        # scan_filter.ScanMode 0 = FULL
        scan_filter = self.iRawDataPlus.GetFilterForScanNumber(scan_number)

        return scan_filter.ScanMode == MSOrderType.Ms

    def get_all_filters(self) -> Tuple[Dict[int, str], List[str]]:
        """
        Get all scan filters.
        This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']

        """

        scanrange = range(self.start_scan, self.end_scan + 1)
        scanfiltersdic = {}
        scanfilterslist = []
        for scan_number in scanrange:
            scan_label = self.iRawDataPlus.GetScanEventStringForScanNumber(scan_number)
            scanfiltersdic[scan_number] = scan_label
            scanfilterslist.append(scan_label)
        scanfilterset = list(set(scanfilterslist))
        return scanfiltersdic, scanfilterset

    def get_scan_header(self, scan: int) -> Dict[str, Any]:
        """
        Get full dictionary of scan header meta data, i.e. AGC status, ion injection time, etc.

        Parameters:
        -----------
        scan : int
            The scan number.

        """
        header = self.iRawDataPlus.GetTrailerExtraInformation(scan)

        header_dic = {}
        for i in range(header.Length):
            header_dic.update({header.Labels[i]: header.Values[i]})
        return header_dic

    @staticmethod
    def get_rt_time_from_trace(trace) -> Tuple[List[float], List[float], List[int]]:
        """trace: ThermoFisher.CommonCore.Data.Business.ChromatogramSignal"""
        return list(trace.Times), list(trace.Intensities), list(trace.Scans)

    def get_eics(
        self,
        target_mzs: List[float],
        tic_data: Dict[str, Any],
        ms_type="MS !d",
        peak_detection=True,
        smooth=True,
        plot=False,
        ax: Optional[axes.Axes] = None,
        legend=False,
    ) -> Tuple[Dict[float, EIC_Data], axes.Axes]:
        """ms_type: str ('MS', MS2')
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

        """

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

        data = self.iRawDataPlus.GetChromatogramData(
            all_chroma_settings, self.start_scan, self.end_scan, options
        )

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
                rt, eic, scans = self.get_rt_time_from_trace(trace)
                if smooth:
                    eic = self.smooth_tic(eic)

                chroma[target_mzs[i]] = EIC_Data(scans=scans, time=rt, eic=eic)
                if plot:
                    ax.plot(rt, eic, label="{:.5f}".format(target_mzs[i]))

        if peak_detection:
            # max_eic = self.get_max_eic(chroma)
            max_signal = max(tic_data.tic)

            for eic_data in chroma.values():
                eic = eic_data.eic
                time = eic_data.time

                if len(eic) != len(tic_data.tic):
                    warn(
                        "The software assumes same lenth of TIC and EIC, this does not seems to be the case and the results mass spectrum selected by the scan number might not be correct"
                    )

                if eic.max() > 0:
                    centroid_eics = self.eic_centroid_detector(time, eic, max_signal)
                    eic_data.apexes = [i for i in centroid_eics]

                    if plot:
                        for peak_indexes in eic_data.apexes:
                            apex_index = peak_indexes[1]
                            ax.plot(
                                time[apex_index],
                                eic[apex_index],
                                marker="x",
                                linewidth=0,
                            )

        if plot:
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("a.u.")
            ax.set_title(ms_type + " EIC")
            ax.tick_params(axis="both", which="major", labelsize=12)
            ax.axes.spines["top"].set_visible(False)
            ax.axes.spines["right"].set_visible(False)

            if legend:
                legend = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 0, 0.07, 1))
                fig.subplots_adjust(right=0.76)
                # ax.set_prop_cycle(color=plt.cm.gist_rainbow(np.linspace(0, 1, len(traces))))

                d = {"down": 30, "up": -30}

                def func(evt):
                    if legend.contains(evt):
                        bbox = legend.get_bbox_to_anchor()
                        bbox = Bbox.from_bounds(
                            bbox.x0, bbox.y0 + d[evt.button], bbox.width, bbox.height
                        )
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

    def get_tic(
        self, ms_type="MS !d", peak_detection=True, smooth=True, plot=False, ax=None,trace_type='TIC',
    ) -> Tuple[TIC_Data, axes.Axes]:
        """ms_type: str ('MS !d', 'MS2', None)
            if you use None you get all scans.
        peak_detection: bool
        smooth: bool
        plot: bool
        ax: matplotlib axis object
        trace_type: str ('TIC','BPC')

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
        """
        if trace_type == 'TIC':
            settings = ChromatogramTraceSettings(TraceType.TIC)
        elif trace_type == 'BPC':
            settings = ChromatogramTraceSettings(TraceType.BasePeak)
        else:
            print(f'{trace_type} undefined')
        if ms_type == "all":
            settings.Filter = None
        else:
            settings.Filter = ms_type

        chroma_settings = IChromatogramSettings(settings)

        data = self.iRawDataPlus.GetChromatogramData(
            [chroma_settings], self.start_scan, self.end_scan
        )

        trace = ChromatogramSignal.FromChromatogramData(data)

        data = TIC_Data(time=[], scans=[], tic=[], bpc=[], apexes=[])

        if trace[0].Length > 0:
            for i in range(trace[0].Length):
                # print(trace[0].HasBasePeakData,trace[0].EndTime )

                # print("  {} - {}, {}".format( i, trace[0].Times[i], trace[0].Intensities[i] ))
                data.time.append(trace[0].Times[i])
                data.tic.append(trace[0].Intensities[i])
                data.scans.append(trace[0].Scans[i])

                # print(trace[0].Scans[i])
            if smooth:
                data.tic = self.smooth_tic(data.tic)

            else:
                data.tic = np.array(data.tic)

            if peak_detection:
                centroid_peak_indexes = [
                    i for i in self.centroid_detector(data.time, data.tic)
                ]

                data.apexes = centroid_peak_indexes

            if plot:
                if not ax:
                    import matplotlib.pyplot as plt

                    ax = plt.gca()
                    # fig, ax = plt.subplots(figsize=(6, 3))

                ax.plot(data.time, data.tic, label=trace_type)
                ax.set_xlabel("Time (min)")
                ax.set_ylabel("a.u.")
                if peak_detection:
                    for peak_indexes in data.apexes:
                        apex_index = peak_indexes[1]
                        ax.plot(
                            data.time[apex_index],
                            data.tic[apex_index],
                            marker="x",
                            linewidth=0,
                        )

                # plt.show()
                if trace_type == 'BPC':
                    data.bpc = data.tic
                    data.tic = []
                return data, ax
            if trace_type == 'BPC':
                data.bpc = data.tic
                data.tic = []
            return data, None

        else:
            return None, None

    def get_average_mass_spectrum(
        self,
        spectrum_mode: str = "profile",
        auto_process: bool = True,
        ppm_tolerance: float = 5.0,
        ms_type: str = "MS1",
    ) -> MassSpecProfile | MassSpecCentroid:
        """
        Averages mass spectra over a scan range using Thermo's AverageScansInScanRange method
        or a scan list using Thermo's AverageScans method
        spectrum_mode: str
            centroid or profile mass spectrum
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object
        ms_type: str
            String of form 'ms1' or 'ms2' or 'MS3' etc. Valid up to MS10.
            Internal function converts to Thermo MSOrderType class.

        """

        def get_profile_mass_spec(averageScan, d_params: dict, auto_process: bool):
            mz_list = list(averageScan.SegmentedScan.Positions)
            abund_list = list(averageScan.SegmentedScan.Intensities)

            data_dict = {
                Labels.mz: mz_list,
                Labels.abundance: abund_list,
            }

            return MassSpecProfile(data_dict, d_params, auto_process=auto_process)

        def get_centroid_mass_spec(averageScan, d_params: dict):
            noise = list(averageScan.centroidScan.Noises)

            baselines = list(averageScan.centroidScan.Baselines)

            rp = list(averageScan.centroidScan.Resolutions)

            magnitude = list(averageScan.centroidScan.Intensities)

            mz = list(averageScan.centroidScan.Masses)

            array_noise_std = (np.array(noise) - np.array(baselines)) / 3
            l_signal_to_noise = np.array(magnitude) / array_noise_std

            d_params["baseline_noise"] = np.average(array_noise_std)

            d_params["baseline_noise_std"] = np.std(array_noise_std)

            data_dict = {
                Labels.mz: mz,
                Labels.abundance: magnitude,
                Labels.rp: rp,
                Labels.s2n: list(l_signal_to_noise),
            }

            mass_spec = MassSpecCentroid(data_dict, d_params, auto_process=False)

            return mass_spec

        d_params = self.set_metadata(
            firstScanNumber=self.start_scan, lastScanNumber=self.end_scan
        )

        # Create the mass options object that will be used when averaging the scans
        options = MassOptions()
        options.ToleranceUnits = ToleranceUnits.ppm
        options.Tolerance = ppm_tolerance

        # Get the scan filter for the first scan.  This scan filter will be used to located
        # scans within the given scan range of the same type
        scanFilter = self.iRawDataPlus.GetFilterForScanNumber(self.start_scan)

        # force it to only look for the MSType
        scanFilter = self.set_msordertype(scanFilter, ms_type)

        if isinstance(self.scans, tuple):
            averageScan = Extensions.AverageScansInScanRange(
                self.iRawDataPlus, self.start_scan, self.end_scan, scanFilter, options
            )

            if averageScan:
                if spectrum_mode == "profile":
                    mass_spec = get_profile_mass_spec(
                        averageScan, d_params, auto_process
                    )

                    return mass_spec

                elif spectrum_mode == "centroid":
                    if averageScan.HasCentroidStream:
                        mass_spec = get_centroid_mass_spec(averageScan, d_params)

                        return mass_spec

                    else:
                        raise ValueError(
                            "No Centroind data available for the selected scans"
                        )
                else:
                    raise ValueError("spectrum_mode must be 'profile' or centroid")
            else:
                raise ValueError("No data found for the selected scans")

        elif isinstance(self.scans, list):
            d_params = self.set_metadata(scans_list=self.scans)

            scans = List[int]()
            for scan in self.scans:
                scans.Add(scan)

            averageScan = Extensions.AverageScans(self.iRawDataPlus, scans, options)

            if averageScan:
                if spectrum_mode == "profile":
                    mass_spec = get_profile_mass_spec(
                        averageScan, d_params, auto_process
                    )

                    return mass_spec

                elif spectrum_mode == "centroid":
                    if averageScan.HasCentroidStream:
                        mass_spec = get_centroid_mass_spec(averageScan, d_params)

                        return mass_spec

                    else:
                        raise ValueError(
                            "No Centroind data available for the selected scans"
                        )

                else:
                    raise ValueError("spectrum_mode must be 'profile' or centroid")

            else:
                raise ValueError("No data found for the selected scans")

        else:
            raise ValueError("scans must be a list intergers or a tuple if integers")

    def set_metadata(
        self,
        firstScanNumber=0,
        lastScanNumber=0,
        scans_list=False,
        label=Labels.thermo_profile,
    ):
        """
        Collect metadata to be ingested in the mass spectrum object

        scans_list: list[int] or false
        lastScanNumber: int
        firstScanNumber: int
        """

        d_params = default_parameters(self.file_path)

        # assumes scans is full scan or reduced profile scan

        d_params["label"] = label

        if scans_list:
            d_params["scan_number"] = scans_list

            d_params["polarity"] = self.get_polarity_mode(scans_list[0])

        else:
            d_params["scan_number"] = "{}-{}".format(firstScanNumber, lastScanNumber)

            d_params["polarity"] = self.get_polarity_mode(firstScanNumber)

        d_params["analyzer"] = self.iRawDataPlus.GetInstrumentData().Model

        d_params["acquisition_time"] = self.get_creation_time()

        d_params["instrument_label"] = self.iRawDataPlus.GetInstrumentData().Name

        return d_params

    def get_centroid_msms_data(self, scan):
        """
        .. deprecated:: 2.0
            This function will be removed in CoreMS 2.0. Please use `get_average_mass_spectrum()` instead for similar functionality.
        """

        warnings.warn(
            "The `get_centroid_msms_data()` is deprecated as of CoreMS 2.0 and will be removed in a future version. "
            "Please use `get_average_mass_spectrum()` instead.",
            DeprecationWarning,
        )

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

        d_params["baseline_noise"] = np.average(array_noise_std)

        d_params["baseline_noise_std"] = np.std(array_noise_std)

        data_dict = {
            Labels.mz: mz,
            Labels.abundance: magnitude,
            Labels.rp: rp,
            Labels.s2n: list(l_signal_to_noise),
        }

        mass_spec = MassSpecCentroid(data_dict, d_params, auto_process=False)
        mass_spec.settings.noise_threshold_method = "relative_abundance"
        mass_spec.settings.noise_threshold_min_relative_abundance = 1
        mass_spec.process_mass_spec()
        return mass_spec

    def get_average_mass_spectrum_by_scanlist(
        self,
        scans_list: List[int],
        auto_process: bool = True,
        ppm_tolerance: float = 5.0,
    ) -> MassSpecProfile:
        """
        Averages selected scans mass spectra using Thermo's AverageScans method
        scans_list: list[int]
        auto_process: bool
            If true performs peak picking, and noise threshold calculation after creation of mass spectrum object
        Returns:
            MassSpecProfile

         .. deprecated:: 2.0
        This function will be removed in CoreMS 2.0. Please use `get_average_mass_spectrum()` instead for similar functionality.
        """

        warnings.warn(
            "The `get_average_mass_spectrum_by_scanlist()` is deprecated as of CoreMS 2.0 and will be removed in a future version. "
            "Please use `get_average_mass_spectrum()` instead.",
            DeprecationWarning,
        )

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

        data_dict = {
            Labels.mz: mz_list,
            Labels.abundance: abund_list,
        }

        mass_spec = MassSpecProfile(data_dict, d_params, auto_process=auto_process)

        return mass_spec


class ImportDataDependentThermoMSFileReader(ThermoBaseClass, LC_Calculations):

    """Collection of methdos to import LC data dependent acquisition from Thermo's raw file
    Intended do create the LCMS object --> ChromaPeaks --> MSobj FullScan --> Dependent MS/MS Obj
    """

    def __init__(self, file_location: str, selected_mzs: List[float] = None):
        """
        target_mzs: list[float] monoisotopic target m/z  or None
            Details: None will defalt to depends scans selected m/
        file_location: str, Path, or S3Path

        """
        super().__init__(file_location)

        eic_tolerance_ppm = self.chromatogram_settings.eic_tolerance_ppm
        enforce_target_ms2 = self.chromatogram_settings.enforce_target_ms2
        average_target_mz = self.chromatogram_settings.average_target_mz

        print("TOLERANCE = {} ppm".format(eic_tolerance_ppm))
        self._selected_mzs = self._init_target_mz(
            selected_mzs, enforce_target_ms2, eic_tolerance_ppm, average_target_mz
        )

        self.lcms = DataDependentLCMS(file_location, self._selected_mzs, self)

    @property
    def selected_mzs(self) -> List[float]:
        return list(self._selected_mzs)

    def get_lcms_obj(self):
        return self.lcms

    def get_precursors_list(self, precision_decimals=5):
        """returns a set of unique precursors m/z
        precision_decimals: int
            change this parameters does not seem to affect the number of dependent scans selected
            needs more investigation
        """

        precursors_mzs = set()

        for scan in range(self.start_scan, self.end_scan):
            scan_filter = self.iRawDataPlus.GetFilterForScanNumber(scan)

            MSOrder = scan_filter.MSOrder

            if MSOrder == MSOrderType.Ms:
                scanDependents = self.iRawDataPlus.GetScanDependents(
                    scan, precision_decimals
                )

                for scan_dependent_detail in scanDependents.ScanDependentDetailArray:
                    for precursor_mz in scan_dependent_detail.PrecursorMassArray:
                        precursors_mzs.add(precursor_mz)

        return precursors_mzs

    def _init_target_mz(
        self,
        selected_mzs: List[float],
        enforce_target_ms2: bool,
        tolerance_ppm: float,
        average_target_mz: bool,
    ):
        precursors_mzs = self.get_precursors_list()

        if selected_mzs is None:
            # no selected m/z list provided, default to use the precursos m/z
            if average_target_mz:
                searchmz = MZSearch(
                    precursors_mzs,
                    precursors_mzs,
                    tolerance_ppm,
                    average_target_mz=average_target_mz,
                )
                return searchmz.averaged_target_mz
            else:
                return precursors_mzs
            # searchmz.start()
            # searchmz.join()

        elif selected_mzs and enforce_target_ms2 is False:
            # selected m/z list provided, and not enforcing being selected as precursor
            return selected_mzs

        elif selected_mzs and enforce_target_ms2:
            # search the selected m/z list in the precursors m/z with a ms/ms experiment
            # print("YEAHHHHH")
            searchmz = MZSearch(
                precursors_mzs,
                selected_mzs,
                tolerance_ppm,
                average_target_mz=average_target_mz,
            )
            searchmz.start()
            searchmz.join()
            return sorted(searchmz.results.keys())


class ImportMassSpectraThermoMSFileReader(ThermoBaseClass, LC_Calculations):

    """Collection of methdos to import Summed/Averaged mass spectrum from Thermo's raw file
    Currently only for profile mode data
    Returns MassSpecProfile object
    """

    def get_icr_transient_times(self):
        """
        Return a list for transient time targets for all scans, or selected scans range
        Resolving Power and Transient time targets based on 7T FT-ICR MS system
        """

        res_trans_time = {
            "50": 0.384,
            "100000": 0.768,
            "200000": 1.536,
            "400000": 3.072,
            "750000": 6.144,
            "1000000": 12.288,
        }

        firstScanNumber = self.start_scan

        lastScanNumber = self.end_scan

        transient_time_list = []

        for scan in range(firstScanNumber, lastScanNumber):
            scan_header = self.get_scan_header(scan)

            rp_target = scan_header["FT Resolution:"]

            transient_time = res_trans_time.get(rp_target)

            transient_time_list.append(transient_time)

            # print(transient_time, rp_target)

        return transient_time_list

    import numpy as np

    class RawFileReader:
        def get_data(self, scan: int, d_parameter: dict, scan_type: str):
            """
            Retrieve mass spectrometry data for a given scan.

            Parameters
            ----------
            scan : int
                The scan number.
            d_parameter : dict
                A dictionary to store additional parameters.
            scan_type : str
                The type of scan ("Centroid" or "Profile").

            """
            if scan_type == "Centroid":
                centroidStream = self.iRawDataPlus.GetCentroidStream(scan, False)

                noise = list(centroidStream.Noises)
                baselines = list(centroidStream.Baselines)
                rp = list(centroidStream.Resolutions)
                magnitude = list(centroidStream.Intensities)
                mz = list(centroidStream.Masses)

                array_noise_std = (np.array(noise) - np.array(baselines)) / 3
                l_signal_to_noise = np.array(magnitude) / array_noise_std

                d_parameter["baseline_noise"] = np.average(array_noise_std)
                d_parameter["baseline_noise_std"] = np.std(array_noise_std)

                data_dict = {
                    Labels.mz: mz,
                    Labels.abundance: magnitude,
                    Labels.rp: rp,
                    Labels.s2n: l_signal_to_noise,
                }

            else:
                scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(scan)
                profileStream = self.iRawDataPlus.GetSegmentedScanFromScanNumber(scan, scanStatistics)
                magnitude = list(profileStream.Intensities)

            

                mz = list(profileStream.Positions)

                data_dict = {
                    Labels.mz: mz,
                    Labels.abundance: magnitude,
                }

            return data_dict

    def get_best_scans_idx(self, stdevs=2, method="mean", plot=False):
            """
            Method to determine the best scan indexes for selective co-addition.

            Parameters
            ----------
            stdevs : int, optional
                The number of standard deviations to use as the cutoff for filtering out datapoints. Default is 2.
            method : str, optional
                The method to calculate the mean or median of the TIC values. Default is "mean".
            plot : bool, optional
                Whether to plot the TIC with horizontal lines for the standard deviation cutoffs. Default is False.

            Notes
            -----
            This method calculates the mean (default) or median of the TIC values and sets an upper and lower limit
            based on a specified number of standard deviations. The scans with TIC values outside of this range are
            considered the best scans for selective co-addition.

            Empirically, using 1-2 standard deviations is enough to filter out the worst datapoints.

            If `plot` is True, a matplotlib figure is returned along with the list of scan indexes.

            Examples
            --------
            >>> reader = RawFileReader()
            >>> scans = reader.get_best_scans_idx(stdevs=2, method="mean", plot=True)
            """
            tic = pd.DataFrame(self.get_tic(plot=plot))

            if method == "median":
                tic_median = tic["TIC"].median()
            elif method == "mean":
                tic_median = tic["TIC"].mean()
            else:
                print("Method " + str(method) + " undefined")

            tic_std = tic["TIC"].std()

            upperlimit = tic_median - (stdevs * tic_std)
            lowerlimit = tic_median + (stdevs * tic_std)

            tic_filtered = tic[(tic["TIC"] > upperlimit) & (tic["TIC"] < lowerlimit)]
            scans = list(tic_filtered.Scans.values)

            if plot:
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots(figsize=(8, 4))
                ax.plot(tic["Time"], tic["TIC"])
                ax.axhline(y=upperlimit, c="r")
                ax.axhline(y=lowerlimit, c="r")
                return fig, scans
            else:
                return scans
