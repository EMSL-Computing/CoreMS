import numpy
import multiprocessing
from threading import Thread
import sys
import site
from pathlib import Path
from io import BytesIO

import clr
from threading import Thread
import multiprocessing
import numpy
import pandas as pd
from s3path import S3Path
from tqdm import tqdm


from corems.encapsulation.constant import Labels
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.mass_spectra.factory.LC_Class import LCMSBase
from corems.encapsulation.factory.parameters import default_parameters


# do not change the order from the imports statements and reference below 
sys.path.append(site.getsitepackages()[0]+ "/ext_lib")
# sys.path.append("ext_lib")

clr.AddReference("ThermoFisher.CommonCore.RawFileReader")
clr.AddReference("ThermoFisher.CommonCore.Data")
clr.AddReference("ThermoFisher.CommonCore.MassPrecisionEstimator")

from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data import ToleranceUnits, Extensions
from ThermoFisher.CommonCore.Data.Business import MassOptions
from ThermoFisher.CommonCore.Data.FilterEnums import MSOrderType 
from System.Collections.Generic import List



class ImportMassSpectraThermoMSFileReader():

    """     Read FULL mode spectra only from raw file data and store it return a LC-MS class
    *  Default behavior is to load all scans numbers

    *  set start_scan_number  and final_scan_number to change it before calling start(), or run()
    """

    def __init__(self, file_location):

        # Thread.__init__(self)
        if isinstance(file_location, str):
            file_path = Path(file_location)

        if isinstance(file_location, S3Path):
            
            temp_dir = Path('tmp/')
            temp_dir.mkdir(exist_ok=True)

            file_path = temp_dir / file_location.name 
            with open(file_path,'wb') as fh:
                fh.write(file_location.read_bytes())
        
        self.iRawDataPlus = RawFileReaderAdapter.FileFactory(str(file_path))
        
        #removing tmp file
        
        if isinstance(file_location, S3Path):
            file_path.unlink()

        self.res = self.iRawDataPlus.SelectInstrument(0, 1)

        self._initial_scan_number = self.iRawDataPlus.RunHeaderEx.FirstSpectrum

        self._final_scan_number = self.iRawDataPlus.RunHeaderEx.LastSpectrum

        self.file_location = file_location

    @property
    def initial_scan_number(self):
        return self._initial_scan_number

    @property
    def final_scan_number(self):
        return self._final_scan_number

    def get_filter_for_scan_num(self, scan_number):
        """Returns the closest matching run time that corresponds to scan_number for the current
        controller. This function is only supported for MS device controllers.
        e.g.  ['FTMS', '-', 'p', 'NSI', 'Full', 'ms', '[200.00-1000.00]']
        """
        scan_label = self.iRawDataPlus.GetScanEventStringForScanNumber(
            scan_number)

        return str(scan_label).split()

    def check_full_scan(self, scan_number):
        # scan_filter.ScanMode 0 = FULL
        scan_filter = self.iRawDataPlus.GetFilterForScanNumber(scan_number)

        return scan_filter.ScanMode == 0

    def get_polarity_mode(self, scan_number):

        polarity_symbol = self.get_filter_for_scan_num(scan_number)[1]

        if polarity_symbol == "+":

            return 1
            # return "POSITIVE_ION_MODE"

        elif polarity_symbol == "-":

            return -1

        else:

            raise Exception("Polarity Mode Unknown, please set it manually")

    def get_scan_header(self, scan):
        '''
        Get full dictionary of scan header meta data, i.e. AGC status, ion injection time, etc.
        '''
        header = self.iRawDataPlus.GetTrailerExtraInformation(scan)
        header_dic = {}
        for i in numpy.arange(header.Length):
            header_dic.update({header.Labels[i]:header.Values[i]})
        return header_dic

    def get_data(self, scan, d_parameter, scan_type):

        if scan_type == "Centroid":

            centroidStream = self.iRawDataPlus.GetCentroidStream(scan, False)

            noise = list(centroidStream.Noises)

            baselines = list(centroidStream.Baselines)

            rp = list(centroidStream.Resolutions)

            magnitude = list(centroidStream.Intensities)

            mz = list(centroidStream.Masses)

            # charge = scans_labels[5]
            array_noise_std = (numpy.array(noise) - numpy.array(baselines)) / 3
            l_signal_to_noise = numpy.array(magnitude) / array_noise_std

            d_parameter["baselise_noise"] = numpy.average(array_noise_std)

            d_parameter["baselise_noise_std"] = numpy.std(array_noise_std)

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

        d_params = default_parameters(self.file_location)

        # assumes scans is full scan or reduced profile scan

        d_params["label"] = Labels.thermo_profile

        if scans_list:
            d_params['scan_number'] = scans_list

            d_params["polarity"] = self.get_polarity_mode(scans_list[0])

        else:

            d_params['scan_number'] = "{}-{}".format(firstScanNumber, lastScanNumber)

            d_params["polarity"] = self.get_polarity_mode(firstScanNumber)

        d_params['analyzer'] = self.iRawDataPlus.GetInstrumentData().Model

        d_params['instrument_label'] = self.iRawDataPlus.GetInstrumentData().Name

        return d_params

    def get_average_mass_spectrum_by_scanlist(self, scans_list, auto_process: bool = True, ppm_tolerance: float = 5.0):

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

    def get_average_mass_spectrum_in_scan_range(self, first_scan: int = None, last_scan: int = None, auto_process: bool = True, ppm_tolerance: float = 5.0, ms_type=MSOrderType.Ms):

        firstScanNumber = self._initial_scan_number if first_scan is None else first_scan

        lastScanNumber = self._final_scan_number if last_scan is None else last_scan

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

        averageScan = Extensions.AverageScansInScanRange(self.iRawDataPlus, firstScanNumber, lastScanNumber, scanFilter, options)

        if averageScan:
            mz_list = list(averageScan.SegmentedScan.Positions)
            abund_list = list(averageScan.SegmentedScan.Intensities)        

            data_dict = {
                        Labels.mz: mz_list,
                        Labels.abundance: abund_list,
                    }

            mass_spec = MassSpecProfile(data_dict, d_params, auto_process=auto_process)

            return mass_spec
        else:
            raise Exception('no data found for the MSOrderType = {}'.format(ms_type))

    def get_summed_mass_spectrum(self, initial_scan_number, final_scan_number=None,
                                 auto_process=True, pd_method=True, pd_merge_n=100):

        d_params = default_parameters(self.file_location)

        # assumes scans is full scan or reduced profile scan

        d_params["label"] = Labels.thermo_profile

        if type(initial_scan_number) is list:
            d_params["polarity"] = self.get_polarity_mode(initial_scan_number[0])

            scanrange = initial_scan_number
        else:
            d_params["polarity"] = self.get_polarity_mode(initial_scan_number)

            if final_scan_number is None:
                final_scan_number = self._final_scan_number

            scanrange = range(initial_scan_number, final_scan_number + 1)

        if pd_method:

            def sort_sum_df(df):
                """
                Nested function to sort dataframe and sum rows with exact matching indexes (m/z)
                """
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
            data_dict = {
                        Labels.mz: list(big_df.index.values),
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

            data_dict = {
                        Labels.mz: mz_all,
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
        tic = self.get_tic()

        if method == 'median':
            tic_median = tic['TIC'].median()
        elif method == 'mean':
            tic_median = tic['TIC'].mean()
        else:
            print("Method " + str(method) + " undefined")

        tic_std = tic['TIC'].std()

        upperlimit = tic_median - (stdevs * tic_std)
        lowerlimit = tic_median + (stdevs * tic_std)

        tic_filtered = tic[(tic['TIC'] > upperlimit) & (tic['TIC'] < lowerlimit)]
        scans = list(tic_filtered.index.values)

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(tic['Time'], tic['TIC'])
            ax.axhline(y=upperlimit, c='r')
            ax.axhline(y=lowerlimit, c='r')
            return fig, scans
        else:
            return scans
