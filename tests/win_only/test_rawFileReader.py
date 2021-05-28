import sys
import clr

sys.path.append("ext_lib")

clr.AddReference("ThermoFisher.CommonCore.RawFileReader")
clr.AddReference("ThermoFisher.CommonCore.Data")
clr.AddReference("ThermoFisher.CommonCore.MassPrecisionEstimator")

from matplotlib import pyplot

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data import ToleranceUnits, Extensions

from ThermoFisher.CommonCore.Data.Business import Scan
from ThermoFisher.CommonCore.Data.Business import Device
from ThermoFisher.CommonCore.Data.Business import MassOptions
from ThermoFisher.CommonCore.MassPrecisionEstimator import PrecisionEstimate
from System.Collections.Generic import List
from System import String 

def getRawFile(filename):

     rawFile = RawFileReaderAdapter.FileFactory(filename)

     print("The RAW file has data from {0} instruments".format(rawFile.InstrumentCount))

     rawFile.SelectInstrument(Device.MS, 1)

     return rawFile

def readAllSpectra( rawFile, firstScanNumber:int, lastScanNumber:int, outputData:bool):

    for scanNumber in range(firstScanNumber,lastScanNumber+1):
            
        scanFilter = rawFile.GetFilterForScanNumber(firstScanNumber)
        
        print( scanFilter.ToString() )
        
        # Get the scan from the RAW file.  This method uses the Scan.FromFile method which returns a
        # Scan object that contains both the segmented and centroid (label) data from an FTMS scan
        # or just the segmented data in non-FTMS scans.  The GetSpectrum method demonstrates an
        # alternative method for reading scans.
        
        scan = Scan.FromFile(rawFile, scanNumber)
        
        # If that scan contains FTMS data then Centroid stream will be populated so check to see if it is present.
        labelSize = 0

        if scan.HasCentroidStream:
        
            labelSize = scan.CentroidScan.Length
        

        # for non-FTMS data, the preferred data will be populated
        dataSize = scan.PreferredMasses.Length

        if outputData:
        
            print("Spectrum {0} - {1}: normal {2}, label {3} points".format(scanNumber, scanFilter.ToString(), dataSize, labelSize) )

def AnalyzeAllScans(rawFile, firstScanNumber:int, lastScanNumber:int):

            # Test the preferred (normal) data and centroid (high resolution/label) data
            failedCentroid = 0
            failedPreferred = 0

            for scanNumber in range(firstScanNumber,lastScanNumber+1):
                # Get each scan from the RAW file
                scan = Scan.FromFile(rawFile, scanNumber)

                # Check to see if the RAW file contains label (high-res) data and if it is present
                # then look for any data that is out of order
                if scan.HasCentroidStream:
                
                    if scan.CentroidScan.Length > 0:
                    
                        currentMass = scan.CentroidScan.Masses[0]

                        for index in range(1, scan.CentroidScan.Length):
                        
                            if scan.CentroidScan.Masses[index] > currentMass:
                            
                                currentMass = scan.CentroidScan.Masses[index]
                            
                            else:
                            
                                if failedCentroid == 0:
                                
                                    print("First failure: Failed in scan data at: Scan: " + scanNumber + " Mass: "
                                        + currentMass.ToString("F4"))
                                

                                failedCentroid += 1
                   

                # Check the normal (non-label) data in the RAW file for any out-of-order data
                if scan.PreferredMasses.Length > 0:
                
                    currentMass = scan.PreferredMasses[0]
                    
                    # print(scan.PreferredMasses.Length, scan.CentroidScan.Length, scan.SegmentedScan.Positions.Length)
                    
                    for index in range(1, scan.PreferredMasses.Length):
                    
                        if scan.PreferredMasses[index] > currentMass:
                            
                            currentMass = scan.PreferredMasses[index]
                            
                        else:
                        
                            if (failedPreferred == 0):
                            
                                print("First failure: Failed in scan data at: Scan: " + str(scanNumber) + " Mass: "
                                    + currentMass.ToString("F2"))
                            
                            failedPreferred += 1
            
            # Display a message indicating if any of the scans had data that was "out of order"
            if failedPreferred == 0 and failedCentroid == 0:
            
                print("Analysis completed: No out of order data found")
            
            else:
            
                print("Analysis completed: Preferred data failed: " + str(failedPreferred) + " Centroid data failed: " + str(failedCentroid) )
        

def CalculateMassPrecision(rawFile, scanNumber:int):
        
            # Get the scan from the RAW file
            scan = Scan.FromFile(rawFile, scanNumber)

            # Get the scan event and from the scan event get the analyzer type for this scan
            scanEvent = rawFile.GetScanEventForScanNumber(scanNumber)

            scanFilter = rawFile.GetFilterForScanNumber(scanNumber)

            print(scanFilter.MassAnalyzer)
            print(scanEvent)

            # Get the trailer extra data to get the ion time for this file
            logEntry = rawFile.GetTrailerExtraInformation(scanNumber)

            print(logEntry.Labels)    

            trailerHeadings = List[String]()
            trailerValues = List[String]()
            for i in range(logEntry.Length):
                
                trailerHeadings.Add(String(logEntry.Labels[i]))
                trailerValues.Add(String(logEntry.Values[i]))
            
            # create the mass precision estimate object
            precisionEstimate = PrecisionEstimate()

            # Get the ion time from the trailer extra data values
            ionTime = precisionEstimate.GetIonTime(scanFilter.MassAnalyzer, scan, trailerHeadings, trailerValues)

            # Calculate the mass precision for the scan
            listResults = precisionEstimate.GetMassPrecisionEstimate(scan, scanFilter.MassAnalyzer, ionTime, rawFile.RunHeader.MassResolution)
            
            # Output the mass precision results
            if len(listResults) > 0:

                print("Mass Precision Results:")

                for result in listResults:

                    print("Mass {}, mmu = {}, ppm = {}".format(result.Mass, result.MassAccuracyInMmu, result.MassAccuracyInPpm) )
        

def GetAverageSpectrum(rawFile, firstScanNumber:int, lastScanNumber:int, outputData:bool):
        
            # Create the mass options object that will be used when averaging the scans
            options = MassOptions()

            options.ToleranceUnits = ToleranceUnits.ppm
            options.Tolerance = 5.0

            # Get the scan filter for the first scan.  This scan filter will be used to located
            # scans within the given scan range of the same type
            scanFilter = rawFile.GetFilterForScanNumber(firstScanNumber)
            
            print(scanFilter.ScanMode)
            
            # Get the average mass spectrum for the provided scan range. In addition to getting the
            # average scan using a scan range, the library also provides a similar method that takes
            # a time range.
            averageScan = Extensions.AverageScansInScanRange(rawFile, firstScanNumber, lastScanNumber, scanFilter, options)

            #average= ScanAveragerFactory.GetScanAverager(rawFile)

            #averageScan = rawFile.AverageScansInScanRange(firstScanNumber, lastScanNumber, scanFilter, options)
            if averageScan.HasCentroidStream:
            
                print("Average spectrum ({0} points)".format( averageScan.CentroidScan.Length) )

                # Print the spectral data (mass, intensity values)
                #if outputData:
                
                #    for i in range(averageScan.CentroidScan.Length):
                    
                #        print("  {}\t{}".format(averageScan.CentroidScan.Masses[i], averageScan.CentroidScan.Intensities[i]))
            
            # This example uses a different method to get the same average spectrum that was calculated in the
            # previous portion of this method.  Instead of passing the start and end scan, a list of scans will
            # be passed to the GetAveragedMassSpectrum function.
            scans = List[int]()
            for scan in (1, 6, 7, 9, 11, 12, 14):
                scans.Add(scan)

            averageScan = Extensions.AverageScans(rawFile,scans, options)
            
            len_data = averageScan.SegmentedScan.Positions.Length
            
            mz_list = list(averageScan.SegmentedScan.Positions)
            abund_list = list(averageScan.SegmentedScan.Intensities)
            
            #for i in range(len_data):
                    
            # mz_list.append(averageScan.SegmentedScan.Positions[i])
                    
            # abund_list.append(averageScan.SegmentedScan.Intensities[i])

            pyplot.plot(mz_list, abund_list)    
            #pyplot.show()
            
            centroid_mz_list = []
            abundance_mz_list = []
            if averageScan.HasCentroidStream:
            
                print("Average spectrum ({0} points)".format(averageScan.CentroidScan.Length))  

                # Print the spectral data (mass, intensity values)
                if outputData:
                
                    for i in range(averageScan.CentroidScan.Length):
                        centroid_mz_list.append(averageScan.CentroidScan.Masses[i])
                        averageScan.CentroidScan.Resolutions
                        abundance_mz_list.append(averageScan.CentroidScan.Intensities[i])
                        #print("  {}\t{}".format(averageScan.CentroidScan.Masses[i], averageScan.CentroidScan.Intensities[i]) )

            pyplot.plot(centroid_mz_list, abundance_mz_list, linewidth=0, marker='o' )    
            pyplot.show()
                    
            print()
        

def get_metadata(rawFile):
    
    firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum
    lastScanNumber = rawFile.RunHeaderEx.LastSpectrum

    startTime = rawFile.RunHeaderEx.StartTime
    endTime = rawFile.RunHeaderEx.EndTime

    print()
    print("General File Information:")
    print("   RAW file: " + str( rawFile.FileName) )
    print("   RAW file version: " + str( rawFile.FileHeader.Revision))
    print("   Creation date: " + str( rawFile.FileHeader.CreationDate))
    print("   Operator: " + str( rawFile.FileHeader.WhoCreatedId))
    print("   Number of instruments: " + str( rawFile.InstrumentCount))
    print("   Description: " + str( rawFile.FileHeader.FileDescription))
    print("   Instrument model: " + str( rawFile.GetInstrumentData().Model))
    print("   Instrument name: " + str( rawFile.GetInstrumentData().Name))
    print("   Serial number: " + str( rawFile.GetInstrumentData().SerialNumber))
    print("   Software version: " + str( rawFile.GetInstrumentData().SoftwareVersion))
    print("   Firmware version: " + str( rawFile.GetInstrumentData().HardwareVersion))
    print("   Units: " + str( rawFile.GetInstrumentData().Units))
    print("   Mass resolution: {} ".format(rawFile.RunHeaderEx.MassResolution) )
    print("   Number of scans: {}".format( rawFile.RunHeaderEx.SpectraCount) )
    print("   Scan range: {} - {}".format( firstScanNumber, lastScanNumber) )
    print("   Time range: {} - {}".format( startTime, endTime) )
    print("   Mass range: {} - {}".format( rawFile.RunHeaderEx.LowMass, rawFile.RunHeaderEx.HighMass) )
    print()

    #Get information related to the sample that was processed
    print("Sample Information:")
    print("   Sample name: " + str( rawFile.SampleInformation.SampleName))
    print("   Sample id: " + str( rawFile.SampleInformation.SampleId))
    print("   Sample type: " + str( rawFile.SampleInformation.SampleType))
    print("   Sample comment: " + str( rawFile.SampleInformation.Comment))
    print("   Sample vial: " + str( rawFile.SampleInformation.Vial))
    print("   Sample volume: " + str( rawFile.SampleInformation.SampleVolume))
    print("   Sample injection volume: " + str( rawFile.SampleInformation.InjectionVolume))
    print("   Sample row number: " + str( rawFile.SampleInformation.RowNumber))
    print("   Sample dilution factor: " + str( rawFile.SampleInformation.DilutionFactor))

if __name__ == "__main__":
    
    filename = "tests/tests_data/SRFA_NEG_ESI_ORB.raw"
    
    rawFile = getRawFile(filename)
    
    get_metadata(rawFile)

    firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum
    lastScanNumber = rawFile.RunHeaderEx.LastSpectrum

    # readAllSpectra(rawFile, firstScanNumber, lastScanNumber, True)

    #AnalyzeAllScans(rawFile, firstScanNumber, lastScanNumber)

    # CalculateMassPrecision(rawFile, firstScanNumber)
    
    GetAverageSpectrum(rawFile, firstScanNumber, lastScanNumber, True)