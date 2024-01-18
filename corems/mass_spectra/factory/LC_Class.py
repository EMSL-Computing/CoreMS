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
import numpy as np
import pandas as pd

from matplotlib import axes
import numpy as np
from corems.chroma_peak.factory.ChromaPeakClasses import DataDependentPeak
from corems.encapsulation.factory.parameters import LCMSParameters
from corems.encapsulation.factory.processingSetting import LiquidChromatographSetting, MolecularFormulaSearchSettings

from corems.mass_spectra.calc.LC_Calc import LC_Calculations
from corems.mass_spectra.factory.LC_Temp import TIC_Data, EIC_Data
from corems.mass_spectrum.input.numpyArray import ms_from_array_profile, ms_from_array_centroid

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

class MassSpectraBase():
    """Base class for mass spectra objects.  

    Parameters
    -----------
    file_location : str or Path
        The location of the file containing the mass spectra data.
    analyzer : str, optional
        The type of analyzer used to generate the mass spectra data. Defaults to 'Unknown'.
    instrument_label : str, optional
        The type of instrument used to generate the mass spectra data. Defaults to 'Unknown'.
    sample_name : str, optional
        The name of the sample; defaults to the file name if not provided to the parser. Defaults to None.
    spectra_parser : object, optional
        The spectra parser object used to create the mass spectra object. Defaults to None.

    Attributes
    -----------
    spectra_parser_class : class
        The class of the spectra parser used to create the mass spectra object.
    file_location : str or Path
        The location of the file containing the mass spectra data.
    sample_name : str
        The name of the sample; defaults to the file name if not provided to the parser.
    analyzer : str
        The type of analyzer used to generate the mass spectra data. Derived from the spectra parser.
    instrument_label : str
        The type of instrument used to generate the mass spectra data. Derived from the spectra parser.
    _scan_info : dict
        A dictionary containing the scan data with columns for scan number, scan time, ms level, precursor m/z, scan text, and scan window (lower and upper). Associated with the property scan_df, which returns a pandas DataFrame or can set this attribute from a pandas DataFrame.
    _ms : dict
        A dictionary containing mass spectra for the dataset, keys of dictionary are scan numbers. Initialized as an empty dictionary.
    _ms_unprocessed: dictionary of pandas.DataFrames or None
        A dictionary of unprocssed mass spectra data, as an (optional) intermediate data product for peak picking.  Key is ms_level, and value is dataframe with columns for scan number, m/z, and intensity. Default is None.
    
    Methods
    --------
    * add_mass_spectrum(MassSpecBase).
        Adds a MassSpecBase object to the dataset within the _ms dictionary.
    * add_mass_spectra(scan_list, spectrum_mode: str = 'profile', use_parser = True, auto_process=True).
        Add mass spectra to _ms slot, from a list of scans
    """
    def __init__(
            self, 
            file_location, 
            analyzer='Unknown', 
            instrument_label='Unknown', 
            sample_name=None,
            spectra_parser = None
            ):

        if isinstance(file_location, str):
            file_location = Path(file_location)
        else:
            file_location = file_location
        if not file_location.exists():
            raise FileExistsError("File does not exist: " + str(file_location))
        
        self.file_location = file_location
        self.sample_name = sample_name
        self.analyzer = analyzer
        self.instrument_label = instrument_label

        # Add the spectra parser class to the object if it is not None
        if spectra_parser is not None:
            self.spectra_parser_class = spectra_parser.__class__
            self.spectra_parser = spectra_parser
            # Check that spectra_pasrser.sample_name is same as sample_name etc, raise warning if not
            if self.sample_name is not None and self.sample_name != self.spectra_parser.sample_name:
                logging.warning("sample_name provided to MassSpectraBase object does not match sample_name provided to spectra parser object")
            if self.analyzer != self.spectra_parser.analyzer:
                logging.warning("analyzer provided to MassSpectraBase object does not match analyzer provided to spectra parser object")
            if self.instrument_label != self.spectra_parser.instrument_label:
                logging.warning("instrument provided to MassSpectraBase object does not match instrument provided to spectra parser object")
            if self.file_location != self.spectra_parser.file_location:
                logging.warning("file_location provided to MassSpectraBase object does not match file_location provided to spectra parser object")  

        # Instantiate empty dictionaries for scan information and mass spectra
        self._scan_info = {}
        self._ms = {}
        self._ms_unprocessed = None
    

    def add_mass_spectra(self, scan_list, spectrum_mode: str = 'profile', ms_level = 1, use_parser = True, auto_process=True):
        """Add mass spectra to _ms dictionary, from a list of scans or single scan

        Notes
        -----
        The mass spectra will inherit the mass_spectrum, ms_peak, and molecular_search parameters from the LCMSBase object.

        
        Parameters
        -----------
        scan_list : list of ints
            List of scans to use to populate _ms slot
        spectrum_mode : str, optional
            The spectrum mode to use for the mass spectra.  Only profile mode is currently supported. Defaults to 'profile'.
        ms_level : int, optional
            The MS level to use for the mass spectra.  This is used to pass the molecular_search parameters from the LCMS object to the individual MassSpectrum objects. Defaults to 1.
        using_parser : bool
            Whether to use the mass spectra parser to get the mass spectra.  Defaults to True.
        auto_process : bool
            Whether to auto-process the mass spectra.  Defaults to True.

        Raises
        ------
        TypeError
            If scan_list is not a list of ints
        ValueError
            If polarity is not 'positive' or 'negative' or if it is a mix
            If ms_level is not 1 or 2
        """
        def add_mass_spectrum(self, mass_spec):
            """Adds a mass spectrum to the dataset.

            Parameters
            -----------
            mass_spec : MassSpectrum
                The corems MassSpectrum object to be added to the dataset.

            Notes
            -----
            This is a helper function for the add_mass_spectra() method, and is not intended to be called directly.
            """
            # check if mass_spec has a scan_number attribute
            if not hasattr(mass_spec, 'scan_number'):
                raise ValueError("Mass spectrum must have a scan_number attribute to be added to the dataset correctly")
            self._ms[mass_spec.scan_number] = mass_spec
        
        # check if scan_list is a list or a single int; if single int, convert to list
        if isinstance(scan_list, int):
            scan_list = [scan_list]
        if not isinstance(scan_list, list):
            raise TypeError("scan_list must be a list of integers")
        for scan in scan_list:
            if not isinstance(scan, int):
                raise TypeError("scan_list must be a list of integers")
            
        # set polarity to -1 if negative mode, 1 if positive mode (for mass spectrum creation)
        if self.polarity == set(['negative']):
            polarity = -1
        elif self.polarity == set(['positive']):
            polarity = 1
        else:
            raise ValueError("Polarity not set for dataset, must be a set containing either 'positive' or 'negative'")
        
        # is not using_parser, check that ms1 and ms2 are not None
        if use_parser == False:
            if ms_level not in self._ms_unprocessed.keys():
                raise ValueError("ms_level {} not found in _ms_unprocessed dictionary".format(ms_level))
        
        scan_list = list(set(scan_list))
        scan_list.sort()
        for scan in scan_list: 
            ms = None
            # Instantiate the mass spectrum object using the parser or the unprocessed data
            if use_parser == False:
                if self._ms_unprocessed[ms_level] is not None and scan in self._ms_unprocessed[ms_level].scan.tolist():
                    my_ms_df = self._ms_unprocessed[ms_level][self._ms_unprocessed[ms_level].scan == scan]
                    ms = ms_from_array_profile(my_ms_df.mz, my_ms_df.intensity, self.file_location, polarity=polarity, auto_process=False)
            if use_parser == True:
                ms = self.spectra_parser.get_mass_spectrum_from_scan(scan_number = scan, spectrum_mode = spectrum_mode, polarity = polarity, auto_process=auto_process)
            
            # Set the mass spectrum parameters, auto-process if auto_process is True, and add to the dataset
            if ms is not None:
                ms.parameters.mass_spectrum = self.parameters.mass_spectrum
                ms.parameters.ms_peak = self.parameters.ms_peak
                if ms_level == 1:
                    ms.parameters.molecular_search = self.parameters.ms1_molecular_search
                elif ms_level == 2:
                    ms.parameters.molecular_search = self.parameters.ms2_molecular_search
                else:
                    raise ValueError("ms_level must be 1 or 2")
                ms.scan_number = scan
                if len(ms.mz_exp_profile) > 3:
                    if auto_process:
                        ms.process_mass_spec()
                self.add_mass_spectrum(ms)
    
    @property
    def scan_df(self):
        """
        pandas.DataFrame : A pandas DataFrame containing the scan info data with columns for scan number, scan time, ms level, precursor m/z, scan text, and scan window (lower and upper).
        """
        scan_df = pd.DataFrame.from_dict(self._scan_info)
        return scan_df
    
    @scan_df.setter
    def scan_df(self, df):
        """
        Sets the scan data for the dataset.

        Parameters
        -----------
        df : pandas.DataFrame
            A pandas DataFrame containing the scan data with columns for scan number, scan time, ms level, precursor m/z, scan text, and scan window (lower and upper).
        """
        self._scan_info = df.to_dict()    

    def __getitem__(self, scan_number):
        
        return self._ms.get(scan_number)


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
        
           