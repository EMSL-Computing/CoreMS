__author__ = "Yuri E. Corilo"
__date__ = "Nov 26, 2019"


import clr

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

import sys

sys.path.append("./lib")

clr.AddReference("ThermoFisher.CommonCore.RawFileReader")

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter

iRawDataPlus = RawFileReaderAdapter.FileFactory("C:/Users/eber373/Desktop/data/WK_ps_lignin_190301112616.raw")

print("FileName", iRawDataPlus.FileName)

print("InstrumentCount", iRawDataPlus.InstrumentCount)

device = iRawDataPlus.GetInstrumentType(0)

print("device", controllerType.get(device))

print()

controller = (iRawDataPlus.SelectInstrument(controllerType.get('MS'), 1))

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

for scan in range( iRunHeader.FirstSpectrum, iRunHeader.LastSpectrum):

        scan_filter = iRawDataPlus.GetFilterForScanNumber(scan)
        
        ionizationMode = scan_filter.IonizationMode 
        
        mSOrder = scan_filter.MSOrder
        
        massAnalyzer = massAnalyzerType.get(scan_filter.MassAnalyzer)
        
        activation = activationType.get(massAnalyzerType.get(scan_filter.GetActivation))
        
        scanMode = scanType.get(scan_filter.ScanMode)

        rt = iRawDataPlus.RetentionTimeFromScanNumber(scan)
        
        scanStatistics = iRawDataPlus.GetScanStatsForScanNumber(scan)
        
        isCentroid = scanStatistics.IsCentroidScan
        
        tic = scanStatistics.TIC

        scanEvent = iRawDataPlus.GetScanEventForScanNumber(scan)
        
        scan_label = (iRawDataPlus.GetScanEventStringForScanNumber(scan))

        #get centroid data
        centroidStream = iRawDataPlus.GetCentroidStream(scan, False)
        
        noises= list(centroidStream.Noises)
        baselines= centroidStream.Baselines
        
        resolutions = list(centroidStream.Resolutions)
        
        coefficientsCount = centroidStream.CoefficientsCount 
        coefficients= centroidStream.Coefficients

        abundance = list(centroidStream.Intensities)
        mz = list(centroidStream.Masses)

        #get profile data
        segmentedStream = iRawDataPlus.GetSegmentedScanFromScanNumber(scan, scanStatistics)
        
        segmentCount = (segmentedStream.SegmentCount)

        ms_length = segmentedStream.SegmentLengths[0]

        abundance = list(segmentedStream.Intensities)
        
        mz = (list(segmentedStream.Positions))
        
        flags = segmentedStream.Flags
        
        massRanges =  segmentedStream.MassRanges
      
       