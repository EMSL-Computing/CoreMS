__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 20201"

import clr
from dataclasses import dataclass
import sys
import site

from matplotlib import pyplot as plt

@dataclass
class ScanDependentDetail:
    scanIndex: int
    filterString: int
    precursorMass: float
    isolationWidth: float
    activation: str
    collision_energy: float
    tic: float
    rt: float
    mz: list
    abundance: list

controllerType = {-1: 'No device',
                  0: 'MS',
                  1: 'Analog',
                  2: 'A/D card',
                  3: 'PDA',
                  4: 'UV',
                  'No device': -1,
                  'MS': 0,
                  'Analog': 1,
                  'A/D card': 2,
                  'PDA': 3,
                  'UV': 4}

massAnalyzerType = {'ITMS': 0,
                    'TQMS': 1,
                    'SQMS': 2,
                    'TOFMS': 3,
                    'FTMS': 4,
                    'Sector': 5,
                    0: 'ITMS',
                    1: 'TQMS',
                    2: 'SQMS',
                    3: 'TOFMS',
                    4: 'FTMS',
                    5: 'Sector'}

activationType = {'CID': 0,
                  'MPD': 1,
                  'ECD': 2,
                  'PQD': 3,
                  'ETD': 4,
                  'HCD': 5,
                  'Any activation type': 6,
                  'SA': 7,
                  'PTR': 8,
                  'NETD': 9,
                  'NPTR': 10,
                  0: 'CID',
                  1: 'MPD',
                  2: 'ECD',
                  3: 'PQD',
                  4: 'ETD',
                  5: 'HCD',
                  6: 'Any activation type',
                  7: 'SA',
                  8: 'PTR',
                  9: 'NETD',
                  10: 'NPTR'}

detectorType = {'CID': 0,
                'PQD': 1,
                'ETD': 2,
                'HCD': 3,
                0: 'CID',
                1: 'PQD',
                2: 'ETD',
                3: 'HCD'}

scanType = {'ScanTypeFull': 0,
            'ScanTypeSIM': 1,
            'ScanTypeZoom': 2,
            'ScanTypeSRM': 3,
            0: 'ScanTypeFull',
            1: 'ScanTypeSIM',
            2: 'ScanTypeZoom',
            3: 'ScanTypeSRM'}

scanType = {'ScanTypeFull': 0,
            'ScanTypeSIM': 1,
            'ScanTypeZoom': 2,
            'ScanTypeSRM': 3,
            0: 'ScanTypeFull',
            1: 'ScanTypeSIM',
            2: 'ScanTypeZoom',
            3: 'ScanTypeSRM'}

controllerType = {-1: 'No device',
                      0: 'MS',
                      1: 'Analog',
                      2: 'A/D card',
                      3: 'PDA',
                      4: 'UV',
                      'No device': -1,
                      'MS': 0,
                      'Analog': 1,
                      'A/D card': 2,
                      'PDA': 3,
                      'UV': 4}

sys.path.append(site.getsitepackages()[0] + "/ext_lib")
sys.path.append("ext_lib")

clr.AddReference("ThermoFisher.CommonCore.RawFileReader")
clr.AddReference("ThermoFisher.CommonCore.Data")
clr.AddReference("ThermoFisher.CommonCore.MassPrecisionEstimator")

clr.AddReference("ThermoFisher.CommonCore.RawFileReader")

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data.FilterEnums import MSOrderType
from ThermoFisher.CommonCore.Data.Business import ChromatogramTraceSettings, TraceType, MassOptions
from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, Range
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings
from ThermoFisher.CommonCore.Data import ToleranceUnits, Extensions

dirpath = "C:\\Users\\eber373\\Desktop\\Data\\LCMS\\RAW Files\\HILIC"
filepath = "\\NEG\\LCMS_5191_CapDev_HILIC_Mix1_NEG_30Apr2021.raw"
iRawDataPlus = RawFileReaderAdapter.FileFactory(dirpath + filepath)

print("FileName", iRawDataPlus.FileName)

print("InstrumentCount", iRawDataPlus.InstrumentCount)

device = iRawDataPlus.GetInstrumentType(0)

print("device", controllerType.get(device))

print()

controler = (iRawDataPlus.SelectInstrument(controllerType.get('MS'), 1))

fileHeader = iRawDataPlus.FileHeader

print("CreationDate", fileHeader.CreationDate)

print("WhoCreatedId", fileHeader.WhoCreatedId)

print("Revision", fileHeader.Revision)

print()


instrumentData = iRawDataPlus.GetInstrumentData()

print("Model", instrumentData.Model)

print("Name", instrumentData.Name)

print("SerialNumber", instrumentData.SerialNumber)

print("SoftwareVersion", instrumentData.SoftwareVersion)

print("Flags", instrumentData.Flags)

print("AxisLabelX", instrumentData.AxisLabelX)

print("AxisLabelY", instrumentData.AxisLabelY)

print()

iRunHeader = iRawDataPlus.RunHeaderEx

print("SpectraCount", iRunHeader.SpectraCount)

print("FirstSpectrum", iRunHeader.FirstSpectrum)

print("LastSpectrum", iRunHeader.LastSpectrum)

print("StartTime", iRunHeader.StartTime)

print("EndTime", iRunHeader.EndTime)

print("LowMass", iRunHeader.LowMass)

print("HighMass", iRunHeader.HighMass)

print("MassResolution", iRunHeader.MassResolution)

print(iRunHeader.FirstSpectrum, iRunHeader.LastSpectrum)

def plot_chroma(x, y, ax=None, c='g'):

    if ax is None:
        ax = plt.gca()

    ax.plot(x, y, color=c)

    ax.set_xlabel("$\t{RT}$", fontsize=12)
    ax.set_ylabel('TIC', fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.axes.spines['top'].set_visible(False)
    ax.axes.spines['right'].set_visible(False)

    # ax.get_yaxis().set_visible(False)
    # ax.spines['left'].set_visible(False)
    return ax

def get_centroid(scan, iRawDataPlus):

    centroidStream = iRawDataPlus.GetCentroidStream(scan, False)

    noises = list(centroidStream.Noises)
    baselines = centroidStream.Baselines

    resolutions = list(centroidStream.Resolutions)

    coefficientsCount = centroidStream.CoefficientsCount
    coefficients = centroidStream.Coefficients

    abundance = list(centroidStream.Intensities)
    mz = list(centroidStream.Masses)

    return mz, abundance

def get_profile(scan, iRawDataPlus, scanStatistics):

    segmentedStream = iRawDataPlus.GetSegmentedScanFromScanNumber(scan, scanStatistics)

    segmentCount = (segmentedStream.SegmentCount)

    ms_length = segmentedStream.SegmentLengths[0]

    abundance = list(segmentedStream.Intensities)

    mz = (list(segmentedStream.Positions))

    flags = segmentedStream.Flags

    massRanges = segmentedStream.MassRanges

    # print([i.High for i in massRanges])
    return mz, abundance

def plot_ms(x, y, ax=None, c='g', is_centroid=True):

    if ax is None:
        ax = plt.gca()

    if is_centroid:
        markerline_a, stemlines_a, baseline_a = ax.stem(x, y, linefmt='-', markerfmt=" ")
        plt.setp(markerline_a, 'color', c, 'linewidth', 2)
        plt.setp(stemlines_a, 'color', c, 'linewidth', 2)
        plt.setp(baseline_a, 'color', c, 'linewidth', 2)
    else:
        ax.plot(x, y)

    ax.set_xlabel("$\t{m/z}$", fontsize=12)
    ax.set_ylabel('Abundance', fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)

    ax.axes.spines['top'].set_visible(False)
    ax.axes.spines['right'].set_visible(False)

    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)

    return ax

def get_ei_chromatogram(IRawDataPlus, target_mz, ppm_tolerance=1000, start_scan=-1, end_scan=-1, ms_type='MS'):

    '''ms_type: str ('MS', MS2')
    start_scan: int default -1 will select the lowest available
    end_scan: int default -1 will select the highest available
    '''

    options = MassOptions()
    options.ToleranceUnits = ToleranceUnits.ppm
    options.Tolerance = ppm_tolerance

    settings = ChromatogramTraceSettings(TraceType.MassRange)
    settings.Filter = ms_type
    settings.MassRanges = [Range(target_mz, target_mz)]

    chroma_settings = IChromatogramSettings(settings)

    settings2 = ChromatogramTraceSettings(TraceType.MassRange)
    settings2.Filter = ms_type
    settings2.MassRanges = [Range(target_mz, target_mz)]
    settings2.MassRanges = [Range(404.5, 404.5)]
    chroma_settings2 = IChromatogramSettings(settings2)
    # chroma_settings2 = IChromatogramSettings(settings)
    # print(chroma_settings.FragmentMass)
    # print(chroma_settings.FragmentMass)
    # print(chroma_settings)
    # print(chroma_settings)

    data = IRawDataPlus.GetChromatogramData([chroma_settings, chroma_settings2], start_scan, end_scan, options)

    trace = ChromatogramSignal.FromChromatogramData(data)

    if trace[1].Length > 0:

        print("Base Peak chromatogram ({} points)".format(trace[0].Length))
        rt = []
        tic = []
        scan = []
        for i in range(trace[1].Length):
            # print(trace[0].HasBasePeakData,trace[0].EndTime )

            # print("  {} - {}, {}".format( i, trace[0].Times[i], trace[0].Intensities[i] ))
            rt.append(trace[1].Times[i])
            tic.append(trace[1].Intensities[i])
            scan.append(trace[1].Scans[i])

        plot_chroma(rt, tic)
        plt.show()

    if trace[0].Length > 0:

        print("Base Peak chromatogram ({} points)".format(trace[0].Length))

        rt = []
        tic = []
        for i in range(trace[0].Length):
            # print(trace[0].HasBasePeakData,trace[0].EndTime )

            # print("  {} - {}, {}".format( i, trace[0].Times[i], trace[0].Intensities[i] ))
            rt.append(trace[0].Times[i])
            tic.append(trace[0].Intensities[i])
        plot_chroma(rt, tic)
        plt.show()

def parse_scan_dependent(iRawDataPlus, iScanDependentDetailArray) -> [ScanDependentDetail]:

    data = []

    for scan_dependent_detail in iScanDependentDetailArray:

        dscan = scan_dependent_detail.ScanIndex

        scan_filter = iRawDataPlus.GetFilterForScanNumber(scan_dependent_detail.ScanIndex)

        scanStatistics = iRawDataPlus.GetScanStatsForScanNumber(dscan)

        tic = scanStatistics.TIC

        rt = iRawDataPlus.RetentionTimeFromScanNumber(dscan)

        reaction = scan_filter.GetReaction(0)

        activation = activationType.get(reaction.ActivationType)

        mz, abun = get_centroid(dscan, iRawDataPlus)

        scanDependentDetail = ScanDependentDetail(scan_dependent_detail.ScanIndex,
                                                  scan_dependent_detail.FilterString,
                                                  scan_dependent_detail.PrecursorMassArray[0],
                                                  scan_dependent_detail.IsolationWidthArray[0],
                                                  activation,
                                                  reaction.CollisionEnergy,
                                                  tic,
                                                  rt,
                                                  mz,
                                                  abun
                                                  )

        data.append(scanDependentDetail)
        # print(scanDependentDetail)

    return data

def get_data():
    ms1_rt = []
    ms2_rt = []
    ms1_tic = []
    ms2_tic = []

    for scan in range(iRunHeader.FirstSpectrum, iRunHeader.LastSpectrum):

        scan_filter = iRawDataPlus.GetFilterForScanNumber(scan)

        ionizationMode = scan_filter.IonizationMode

        MSOrder = scan_filter.MSOrder

        massAnalyzer = massAnalyzerType.get(scan_filter.MassAnalyzer)

        scanMode = scanType.get(scan_filter.ScanMode)

        # print("scanMode", scan_filter.ScanMode)

        rt = iRawDataPlus.RetentionTimeFromScanNumber(scan)

        scanStatistics = iRawDataPlus.GetScanStatsForScanNumber(scan)

        isCentroid = scanStatistics.IsCentroidScan

        tic = scanStatistics.TIC

        scanEvent = iRawDataPlus.GetScanEventForScanNumber(scan)
        time = iRawDataPlus.RetentionTimeFromScanNumber(scan)

        if MSOrder == MSOrderType.Ms2:
            pass
            # ms2_rt.append(rt)
            # ms2_tic.append(tic)

        if MSOrder == MSOrderType.Ms:

            ms1_rt.append(rt)
            ms1_tic.append(tic)
            scanDependents = iRawDataPlus.GetScanDependents(scan, 5)

            print("Scan number {} @ time {} - Instrument type={}, Number dependent scans={}".format(
                  scan, time, scanDependents.RawFileInstrumentType, scanDependents.ScanDependentDetailArray.Length))

            dependent_list = parse_scan_dependent(iRawDataPlus, scanDependents.ScanDependentDetailArray)

            total_tic = 0

            for dependent in dependent_list:

                # print(dependent.precursorMass)
                total_tic = total_tic + dependent.tic
                # print(dependent)
                # trailerData = iRawDataPlus.GetTrailerExtraInformation(dscan)

                # for i, label in enumerate(trailerData.Labels):
                #   print(label, trailerData.Values[i])

            ms2_rt.append(rt)
            ms2_tic.append(total_tic)

            mz, abun = get_profile(scan, iRawDataPlus, scanStatistics)

            # mz, abun = get_centroid(scan, iRawDataPlus)
            # plot_ms(mz, abun, is_centroid=False)
            # plt.show()
            # print(dependent_list)
            # parse_scan_dependent()

        # print(scanEvent)

        # scan_label = (iRawDataPlus.GetScanEventStringForScanNumber(scan))

        # print(scan_label)

        # get centroid data
        # get profile data

    ax = plot_chroma(ms1_rt, ms1_tic)
    plot_chroma(ms2_rt, ms2_tic, ax=ax, c='red')
plt.show()

# get_data()
get_ei_chromatogram(iRawDataPlus, 335.35)
