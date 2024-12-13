__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"


from warnings import warn
import warnings
from collections import defaultdict

from matplotlib import axes
from corems.encapsulation.factory.processingSetting import LiquidChromatographSetting

import numpy as np
import sys
import site
from pathlib import Path
import datetime
import importlib.util
import os

import clr
import pandas as pd
from s3path import S3Path


from typing import Any, Dict, List, Optional, Tuple
from corems.encapsulation.constant import Labels
from corems.mass_spectra.factory.lc_class import MassSpectraBase, LCMSBase
from corems.mass_spectra.factory.chromat_data import EIC_Data, TIC_Data
from corems.mass_spectrum.factory.MassSpectrumClasses import (
    MassSpecProfile,
    MassSpecCentroid,
)
from corems.encapsulation.factory.parameters import LCMSParameters, default_parameters
from corems.mass_spectra.input.parserbase import SpectraParserInterface

# Add the path of the Thermo .NET libraries to the system path
spec = importlib.util.find_spec("corems")
sys.path.append(str(Path(os.path.dirname(spec.origin)).parent) + "/ext_lib/dotnet/")

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
    """Class for parsing Thermo Raw files and extracting information from them.

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
        self,
        ms_type="MS !d",
        peak_detection=True,
        smooth=True,
        plot=False,
        ax=None,
        trace_type="TIC",
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
        if trace_type == "TIC":
            settings = ChromatogramTraceSettings(TraceType.TIC)
        elif trace_type == "BPC":
            settings = ChromatogramTraceSettings(TraceType.BasePeak)
        else:
            raise ValueError(f"{trace_type} undefined")
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
                if trace_type == "BPC":
                    data.bpc = data.tic
                    data.tic = []
                return data, ax
            if trace_type == "BPC":
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


class ImportMassSpectraThermoMSFileReader(ThermoBaseClass, SpectraParserInterface):
    """A class for parsing Thermo RAW mass spectrometry data files and instatiating MassSpectraBase or LCMSBase objects

    Parameters
    ----------
    file_location : str or Path
        The path to the RAW file to be parsed.
    analyzer : str, optional
        The type of mass analyzer used in the instrument. Default is "Unknown".
    instrument_label : str, optional
        The name of the instrument used to acquire the data. Default is "Unknown".
    sample_name : str, optional
        The name of the sample being analyzed. If not provided, the stem of the file_location path will be used.

    Attributes
    ----------
    file_location : Path
        The path to the RAW file being parsed.
    analyzer : str
        The type of mass analyzer used in the instrument.
    instrument_label : str
        The name of the instrument used to acquire the data.
    sample_name : str
        The name of the sample being analyzed.

    Methods
    -------
    * run(spectra=True).
        Parses the RAW file and returns a dictionary of mass spectra dataframes and a scan metadata dataframe.
    * get_mass_spectrum_from_scan(scan_number, polarity, auto_process=True)
        Parses the RAW file and returns a MassSpecBase object from a single scan.
    * get_mass_spectra_obj().
        Parses the RAW file and instantiates a MassSpectraBase object.
    * get_lcms_obj().
        Parses the RAW file and instantiates an LCMSBase object.
    * get_icr_transient_times().
        Return a list for transient time targets for all scans, or selected scans range

    Inherits from ThermoBaseClass and SpectraParserInterface
    """

    def __init__(
        self,
        file_location,
        analyzer="Unknown",
        instrument_label="Unknown",
        sample_name=None,
    ):
        super().__init__(file_location)
        if isinstance(file_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)
        if not file_location.exists():
            raise FileExistsError("File does not exist: " + str(file_location))

        self.file_location = file_location
        self.analyzer = analyzer
        self.instrument_label = instrument_label

        if sample_name:
            self.sample_name = sample_name
        else:
            self.sample_name = file_location.stem

    def load(self):
        pass

    def get_scan_df(self):
        # This automatically brings in all the data
        self.chromatogram_settings.scans = (-1, -1)

        # Get scan df info; starting with bulk ms1 and ms2 scans
        ms1_tic_data, _ = self.get_tic(ms_type="MS", peak_detection=False, smooth=False)
        ms1_scan_dict = {
            "scan": ms1_tic_data.scans,
            "scan_time": ms1_tic_data.time,
            "tic": ms1_tic_data.tic,
        }
        ms1_tic_df = pd.DataFrame.from_dict(ms1_scan_dict)
        ms1_tic_df["ms_level"] = "ms1"

        ms2_tic_data, _ = self.get_tic(
            ms_type="MS2", peak_detection=False, smooth=False
        )
        ms2_scan_dict = {
            "scan": ms2_tic_data.scans,
            "scan_time": ms2_tic_data.time,
            "tic": ms2_tic_data.tic,
        }
        ms2_tic_df = pd.DataFrame.from_dict(ms2_scan_dict)
        ms2_tic_df["ms_level"] = "ms2"

        scan_df = (
            pd.concat([ms1_tic_df, ms2_tic_df], axis=0).sort_values(by="scan").reindex()
        )

        # get scan text
        scan_filter_df = pd.DataFrame.from_dict(
            self.get_all_filters()[0], orient="index"
        )
        scan_filter_df.reset_index(inplace=True)
        scan_filter_df.rename(columns={"index": "scan", 0: "scan_text"}, inplace=True)

        scan_df = scan_df.merge(scan_filter_df, on="scan", how="left")
        scan_df["scan_window_lower"] = scan_df.scan_text.str.extract(
            r"\[(\d+\.\d+)-\d+\.\d+\]"
        )
        scan_df["scan_window_upper"] = scan_df.scan_text.str.extract(
            r"\[\d+\.\d+-(\d+\.\d+)\]"
        )
        scan_df["polarity"] = np.where(
            scan_df.scan_text.str.contains(" - "), "negative", "positive"
        )
        scan_df["precursor_mz"] = scan_df.scan_text.str.extract(r"(\d+\.\d+)@")
        scan_df["precursor_mz"] = scan_df["precursor_mz"].astype(float)
        scan_df["ms_level"] = scan_df["ms_level"].str.replace("ms", "").astype(int)

        # Assign each scan as centroid or profile
        scan_df["ms_format"] = None
        for i in scan_df.scan.to_list():
            if self.iRawDataPlus.IsCentroidScanFromScanNumber(i):
                scan_df.loc[scan_df.scan == i, "ms_format"] = "centroid"
            else:
                scan_df.loc[scan_df.scan == i, "ms_format"] = "profile"

        return scan_df

    def get_ms_raw(self, spectra, scan_df):
        if spectra == "all":
            scan_df_forspec = scan_df
        elif spectra == "ms1":
            scan_df_forspec = scan_df[scan_df.ms_level == 1]
        elif spectra == "ms2":
            scan_df_forspec = scan_df[scan_df.ms_level == 2]
        else:
            raise ValueError("spectra must be 'none', 'all', 'ms1', or 'ms2'")

        # Result container
        res = {}

        # Row count container
        counter = {}

        # Column name container
        cols = {}

        # set at float32
        dtype = np.float32

        # First pass: get nrows
        N = defaultdict(lambda: 0)
        for i in scan_df_forspec.scan.to_list():
            level = scan_df_forspec.loc[scan_df_forspec.scan == i, "ms_level"].values[0]
            scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(i)
            profileStream = self.iRawDataPlus.GetSegmentedScanFromScanNumber(
                i, scanStatistics
            )
            abun = list(profileStream.Intensities)
            abun = np.array(abun)[np.where(np.array(abun) > 0)[0]]

            N[level] += len(abun)

        # Second pass: parse
        for i in scan_df_forspec.scan.to_list():
            scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(i)
            profileStream = self.iRawDataPlus.GetSegmentedScanFromScanNumber(
                i, scanStatistics
            )
            abun = list(profileStream.Intensities)
            mz = list(profileStream.Positions)

            # Get index of abun that are > 0
            inx = np.where(np.array(abun) > 0)[0]
            mz = np.array(mz)[inx]
            mz = np.float32(mz)
            abun = np.array(abun)[inx]
            abun = np.float32(abun)

            level = scan_df_forspec.loc[scan_df_forspec.scan == i, "ms_level"].values[0]

            # Number of rows
            n = len(mz)

            # No measurements
            if n == 0:
                continue

            # Dimension check
            if len(mz) != len(abun):
                warnings.warn("m/z and intensity array dimension mismatch")
                continue

            # Scan/frame info
            id_dict = i

            # Columns
            cols[level] = ["scan", "mz", "intensity"]
            m = len(cols[level])

            # Subarray init
            arr = np.empty((n, m), dtype=dtype)
            inx = 0

            # Populate scan/frame info
            arr[:, inx] = i
            inx += 1

            # Populate m/z
            arr[:, inx] = mz
            inx += 1

            # Populate intensity
            arr[:, inx] = abun
            inx += 1

            # Initialize output container
            if level not in res:
                res[level] = np.empty((N[level], m), dtype=dtype)
                counter[level] = 0

            # Insert subarray
            res[level][counter[level] : counter[level] + n, :] = arr
            counter[level] += n

        # Construct ms1 and ms2 mz dataframes
        for level in res.keys():
            res[level] = pd.DataFrame(res[level])
            res[level].columns = cols[level]
        # rename keys in res to add 'ms' prefix
        res = {f"ms{key}": value for key, value in res.items()}

        return res

    def run(self, spectra="all", scan_df=None):
        """
        Extracts mass spectra data from a raw file.

        Parameters
        ----------
        spectra : str, optional
            Which mass spectra data to include in the output. Default is all.  Other options: none, ms1, ms2.
        scan_df : pandas.DataFrame, optional
            Scan dataframe.  If not provided, the scan dataframe is created from the mzML file.

        Returns
        -------
        tuple
            A tuple containing two elements:
            - A dictionary containing mass spectra data, separated by MS level.
            - A pandas DataFrame containing scan information, including scan number, scan time, TIC, MS level,
                scan text, scan window lower and upper bounds, polarity, and precursor m/z (if applicable).
        """
        # Prepare scan_df
        if scan_df is None:
            scan_df = self.get_scan_df()

        # Prepare mass spectra data
        if spectra != "none":
            res = self.get_ms_raw(spectra=spectra, scan_df=scan_df)
        else:
            res = None

        return res, scan_df

    def get_mass_spectrum_from_scan(
        self, scan_number, spectrum_mode, auto_process=True
    ):
        """Instatiate a MassSpecBase object from a single scan number from the binary file, currently only supports profile mode.

        Parameters
        ----------
        scan_number : int
            The scan number to extract the mass spectrum from.
        polarity : int
            The polarity of the scan.  1 for positive mode, -1 for negative mode.
        spectrum_mode : str
            The type of mass spectrum to extract.  Must be 'profile' or 'centroid'.
        auto_process : bool, optional
            If True, perform peak picking and noise threshold calculation after creating the mass spectrum object. Default is True.

        Returns
        -------
        MassSpecProfile | MassSpecCentroid
            The MassSpecProfile or MassSpecCentroid object containing the parsed mass spectrum.
        """

        if spectrum_mode == "profile":
            scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(scan_number)
            profileStream = self.iRawDataPlus.GetSegmentedScanFromScanNumber(
                scan_number, scanStatistics
            )
            abun = list(profileStream.Intensities)
            mz = list(profileStream.Positions)
            data_dict = {
                Labels.mz: mz,
                Labels.abundance: abun,
            }
            d_params = self.set_metadata(
                firstScanNumber=scan_number,
                lastScanNumber=scan_number,
                scans_list=False,
                label=Labels.thermo_profile,
            )
            mass_spectrum_obj = MassSpecProfile(
                data_dict, d_params, auto_process=auto_process
            )

        elif spectrum_mode == "centroid":
            centroid_scan = self.iRawDataPlus.GetCentroidStream(scan_number, False)
            if centroid_scan.Masses is not None:
                mz = list(centroid_scan.Masses)
                abun = list(centroid_scan.Intensities)
                rp = list(centroid_scan.Resolutions)
                magnitude = list(centroid_scan.Intensities)
                noise = list(centroid_scan.Noises)
                baselines = list(centroid_scan.Baselines)
                array_noise_std = (np.array(noise) - np.array(baselines)) / 3
                l_signal_to_noise = np.array(magnitude) / array_noise_std
                data_dict = {
                    Labels.mz: mz,
                    Labels.abundance: abun,
                    Labels.rp: rp,
                    Labels.s2n: list(l_signal_to_noise),
                }
            else:  # For CID MS2, the centroid data are stored in the profile data location, they do not have any associated rp or baseline data, but they should be treated as centroid data
                scanStatistics = self.iRawDataPlus.GetScanStatsForScanNumber(
                    scan_number
                )
                profileStream = self.iRawDataPlus.GetSegmentedScanFromScanNumber(
                    scan_number, scanStatistics
                )
                abun = list(profileStream.Intensities)
                mz = list(profileStream.Positions)
                data_dict = {
                    Labels.mz: mz,
                    Labels.abundance: abun,
                    Labels.rp: [np.nan] * len(mz),
                    Labels.s2n: [np.nan] * len(mz),
                }
            d_params = self.set_metadata(
                firstScanNumber=scan_number,
                lastScanNumber=scan_number,
                scans_list=False,
                label=Labels.thermo_centroid,
            )
            mass_spectrum_obj = MassSpecCentroid(
                data_dict, d_params, auto_process=auto_process
            )

        return mass_spectrum_obj

    def get_mass_spectra_obj(self):
        """Instatiate a MassSpectraBase object from the binary data file file.

        Returns
        -------
        MassSpectraBase
            The MassSpectra object containing the parsed mass spectra.  The object is instatiated with the mzML file, analyzer, instrument, sample name, and scan dataframe.
        """
        _, scan_df = self.run(spectra="none")
        mass_spectra_obj = MassSpectraBase(
            self.file_location,
            self.analyzer,
            self.instrument_label,
            self.sample_name,
            self,
        )
        scan_df = scan_df.set_index("scan", drop=False)
        mass_spectra_obj.scan_df = scan_df

        return mass_spectra_obj

    def get_lcms_obj(self, spectra="all"):
        """Instatiates a LCMSBase object from the mzML file.

        Parameters
        ----------
        verbose : bool, optional
            If True, print progress messages. Default is True.
        spectra : str, optional
            Which mass spectra data to include in the output. Default is "all".  Other options: "none", "ms1", "ms2".

        Returns
        -------
        LCMSBase
            LCMS object containing mass spectra data. The object is instatiated with the file location, analyzer, instrument, sample name, scan info, mz dataframe (as specifified), polarity, as well as the attributes holding the scans, retention times, and tics.
        """
        _, scan_df = self.run(spectra="none")  # first run it to just get scan info
        res, scan_df = self.run(
            scan_df=scan_df, spectra=spectra
        )  # second run to parse data
        lcms_obj = LCMSBase(
            self.file_location,
            self.analyzer,
            self.instrument_label,
            self.sample_name,
            self,
        )
        if spectra != "none":
            for key in res:
                key_int = int(key.replace("ms", ""))
                res[key] = res[key][res[key].intensity > 0]
                res[key] = (
                    res[key].sort_values(by=["scan", "mz"]).reset_index(drop=True)
                )
                lcms_obj._ms_unprocessed[key_int] = res[key]
        lcms_obj.scan_df = scan_df.set_index("scan", drop=False)
        # Check if polarity is mixed
        if len(set(scan_df.polarity)) > 1:
            raise ValueError("Mixed polarities detected in scan data")
        lcms_obj.polarity = scan_df.polarity[0]
        lcms_obj._scans_number_list = list(scan_df.scan)
        lcms_obj._retention_time_list = list(scan_df.scan_time)
        lcms_obj._tic_list = list(scan_df.tic)

        return lcms_obj

    def get_icr_transient_times(self):
        """Return a list for transient time targets for all scans, or selected scans range

        Notes
        --------
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
