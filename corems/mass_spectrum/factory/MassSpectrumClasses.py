from pathlib import Path
from copy import deepcopy


#from matplotlib import rcParamsDefault, rcParams
from numpy import array, power, float64, where

from corems.mass_spectrum.calc.MassSpectrumCalc import MassSpecCalc
from corems.mass_spectrum.calc.KendrickGroup import KendrickGrouping
from corems.encapsulation.constant import Labels
from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak as MSPeak
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.input.parameter_from_json import load_and_set_parameters_ms

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

def overrides(interface_class):
    def overrider(method):
        assert method.__name__ in dir(interface_class)
        return method
    return overrider

class MassSpecBase(MassSpecCalc, KendrickGrouping):
    '''
    - A iterative mass spectrum base class, stores the profile data and instrument settings
    - Iteration over a list of MSPeaks classes stored at the _mspeaks attributes
    - _mspeaks is populated under the hood by calling process_mass_spec method
    - iteration is null if _mspeaks is empty

    Parameters
    ----------
    mz_exp : list(float)
        list containing the imported experiemental m/z 
        (default is to store profile mode data, but it depends on the input type)
    abundance: list(float)
        list containing the imported abundance 
        (default is to store profile mode data, but it depends on the input type)
    d_params : dict{'str': float,int and str}
        The keyword arguments are used for ...

    Attributes
    ----------
    _mz_exp : list(float)
        This is where we store mz_exp,
    _abundance : list(float)     
        This is where we store _abundance,
    _mspeaks : list(MSPeak)
        store MSpeaks objects identified by a peak picking algorithm

    Relevant Methods
    ----------
    process_mass_spec()
        find or set the noise threshold base on the setting encapsulated at settings.input.ProcessingSetting.MassSpectrumSetting
        - run the peak peaking algorithm and use the method addMSPeaks() to populate _mspeaks attribute

    see also: MassSpecCentroid(), MassSpecfromFreq(), MassSpecProfile()
    '''
    def __init__(self, mz_exp, abundance, d_params, **kwargs):

        self._abundance = array(abundance, dtype=float64)
        self._mz_exp = array(mz_exp, dtype=float64)

        # objects created after process_mass_spec() function
        self._mspeaks = list()
        self._dict_nominal_masses_indexes = dict()
        self._baselise_noise = 0.001
        self._baselise_noise_std = 0.001
        self._dynamic_range = None
        # set to None: initialization occurs inside subclass MassSpecfromFreq
        self._transient_settings = None
        self._frequency_domain = None
        self._mz_cal_profile = None
        self.is_calibrated = False

        self._set_parameters_objects(d_params)
        self._init_settings()

        self.is_centroid = False
        self.has_frequency = False

        self.calibration_order = None
        self.calibration_points = None
        self.calibration_RMS = None

    def _init_settings(self):

        self._parameters = MSParameters()

    def __len__(self):

        return len(self.mspeaks)

    def __getitem__(self, position):

        return self.mspeaks[position]

    def set_indexes(self, list_indexes):
        ''' set the mass spectrum to interate over only the selected MSpeaks indexes'''
        self.mspeaks = [self._mspeaks[i] for i in list_indexes]

        for i, mspeak in  enumerate(self.mspeaks): mspeak.index = i

        self._set_nominal_masses_start_final_indexes()

    def reset_indexes(self):
        ''' reset the mass spectrum to interate over all MSpeaks obj'''
        self.mspeaks = self._mspeaks

        for i, mspeak in  enumerate(self.mspeaks): mspeak.index = i

        self._set_nominal_masses_start_final_indexes()

    def add_mspeak(self, ion_charge, mz_exp,
                    abundance,
                    resolving_power,
                    signal_to_noise,
                    massspec_indexes,
                    exp_freq=None,
                    ms_parent=None
                ):

        mspeak = MSPeak(
                ion_charge,
                mz_exp,
                abundance,
                resolving_power,
                signal_to_noise,
                massspec_indexes,
                len(self._mspeaks),
                exp_freq=exp_freq,
                ms_parent=ms_parent,
        )

        self._mspeaks.append(mspeak)

    def _set_parameters_objects(self, d_params):

        self._calibration_terms = (
            d_params.get("Aterm"),
            d_params.get("Bterm"),
            d_params.get("Cterm"),
        )

        self.label = d_params.get(Labels.label)

        self.analyzer = d_params.get('analyzer')

        self.instrument_label = d_params.get('instrument_label')

        self.polarity = int(d_params.get("polarity"))

        self.scan_number = d_params.get("scan_number")

        self.rt = d_params.get("rt")

        self.mobility_rt = d_params.get("mobility_rt")

        self.mobility_scan = d_params.get("mobility_scan")

        self._filename = d_params.get("filename_path")

        self._dir_location = d_params.get("dir_location")

        self._baselise_noise = d_params.get("baselise_noise")

        self._baselise_noise_std = d_params.get("baselise_noise_std")

        if d_params.get('sample_name') != 'Unknown':

            self.sample_name = d_params.get('sample_name')
            if not self.sample_name:
                self.sample_name = self.filename.stem
        else:

            self.sample_name = self.filename.stem

    def reset_cal_therms(self, Aterm, Bterm, C, fas= 0):

        self._calibration_terms = (Aterm, Bterm, C)

        self._mz_exp = self._f_to_mz()
        self._abundance = self._abundance
        self.find_peaks()
        self.reset_indexes()
        # self.reset_indexes()

    def clear_molecular_formulas(self):

        self.check_mspeaks()
        return array([mspeak.clear_molecular_formulas() for mspeak in self.mspeaks])

    def process_mass_spec(self, keep_profile=True, auto_noise=True, noise_bayes_est=False):
        
        # if runned mannually make sure to rerun filter_by_noise_threshold     
        # calculates noise threshold 
        # do peak picking( create mspeak_objs) 
        # reset mspeak_obj the indexes
         
        self.cal_noise_threshold(auto=auto_noise, bayes=noise_bayes_est)

        self.find_peaks()
        self.reset_indexes()

        if self.mspeaks:
            self._dynamic_range = self.max_abundance / self.min_abundance
        else:
            self._dynamic_range = 0
        if not keep_profile:

            self._abundance *= 0
            self._mz_exp *= 0
            

    def cal_noise_threshold(self, auto=True, bayes=False):

        if self.label == Labels.simulated_profile:

            self._baselise_noise, self._baselise_noise_std = 0.1, 1

        else:

            self._baselise_noise, self._baselise_noise_std = self.run_noise_threshold_calc(auto, bayes=bayes)
        
    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, instance_MSParameters):
        self._parameters = instance_MSParameters

    def set_parameter_from_json(self, parameters_path):
        load_and_set_parameters_ms(self, parameters_path=parameters_path)    

    @property
    def mspeaks_settings(self): return self.parameters.ms_peak

    @mspeaks_settings.setter
    def mspeaks_settings(self, instance_MassSpecPeakSetting):

            self.parameters.ms_peak = instance_MassSpecPeakSetting

    @property
    def settings(self): return self.parameters.mass_spectrum

    @settings.setter
    def settings(self, instance_MassSpectrumSetting):

        self.parameters.mass_spectrum =  instance_MassSpectrumSetting

    @property
    def molecular_search_settings(self):  return self.parameters.molecular_search

    @molecular_search_settings.setter
    def molecular_search_settings(self, instance_MolecularFormulaSearchSettings):

        self.parameters.molecular_search =  instance_MolecularFormulaSearchSettings

    @property
    def mz_cal_profile(self):

        return self._mz_cal_profile

    @mz_cal_profile.setter
    def mz_cal_profile(self, mz_cal_list):

        if len(mz_cal_list) == len(self.mz_exp_profile):
            self._mz_cal_profile = mz_cal_list
        else:
            raise Exception( "calibrated array (%i) is not of the same size of the data (%i)" % (len(mz_cal_list),  len(self.mz_exp_profile)))    

    @property
    def mz_cal(self):
        return array([mspeak.mz_cal for mspeak in self.mspeaks])

    @mz_cal.setter
    def mz_cal(self, mz_cal_list):

            if  len(mz_cal_list) == len(self.mspeaks):
                self.is_calibrated = True
                for index, mz_cal in enumerate(mz_cal_list):
                    self.mspeaks[index].mz_cal = mz_cal
            else: 
                raise Exception( "calibrated array (%i) is not of the same size of the data (%i)" % (len(mz_cal_list),  len(self._mspeaks)))    

    @property
    def mz_exp(self):

        self.check_mspeaks()

        if self.is_calibrated:

            return array([mspeak.mz_cal for mspeak in self.mspeaks])

        else:

            return array([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def freq_exp_profile(self):
        return self._frequency_domain
    
    @freq_exp_profile.setter
    def freq_exp_profile(self, new_data): self._frequency_domain = array(new_data)

    @property
    def mz_exp_profile(self): return self._mz_exp

    @mz_exp_profile.setter
    def mz_exp_profile(self, new_data ): self._mz_exp = array(new_data)

    @property
    def abundance_profile(self): return self._abundance

    @abundance_profile.setter
    def abundance_profile(self, new_data): self._abundance = array(new_data)

    @property
    def abundance(self):
        self.check_mspeaks()
        return array([mspeak.abundance for mspeak in self.mspeaks])

    def freq_exp(self):
        self.check_mspeaks()
        return array([mspeak.freq_exp for mspeak in self.mspeaks])

    @property
    def resolving_power(self):
        self.check_mspeaks()
        return array([mspeak.resolving_power for mspeak in self.mspeaks])

    @property
    def signal_to_noise(self):
        self.check_mspeaks()
        return array([mspeak.signal_to_noise for mspeak in self.mspeaks])

    @property
    def nominal_mz(self):

        if self._dict_nominal_masses_indexes:

            return sorted(list(self._dict_nominal_masses_indexes.keys()))

        else:

            raise ValueError("Nominal indexes not yet set")    

    def get_mz_and_abundance_peaks_tuples(self):

        self.check_mspeaks()
        return [(mspeak.mz_exp, mspeak.abundance) for mspeak in self.mspeaks]

    @property
    def kmd(self):
        self.check_mspeaks()
        return array([mspeak.kmd for mspeak in self.mspeaks])

    @property
    def kendrick_mass(self):
        self.check_mspeaks()
        return array([mspeak.kendrick_mass for mspeak in self.mspeaks])

    @property
    def max_mz_exp(self):
        return max([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def min_mz_exp(self):
        return min([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def max_abundance(self):
        
        return max([mspeak.abundance for mspeak in self.mspeaks])

    @property
    def max_signal_to_noise(self):
        return max([mspeak.signal_to_noise for mspeak in self.mspeaks])

    @property
    def most_abundant_mspeak(self):
        return max(self.mspeaks, key=lambda m: m.abundance)

    @property
    def min_abundance(self):
        return min([mspeak.abundance for mspeak in self.mspeaks])

    # takes too much cpu time 
    @property
    def dynamic_range(self):
        return self._dynamic_range

    @property
    def baselise_noise(self):
        if self._baselise_noise:
            return self._baselise_noise
        else:     
            return None

    @property
    def baselise_noise_std(self):
        if self._baselise_noise_std:
            return self._baselise_noise_std
        else:     
            return None

    @property
    def Aterm(self):
        return self._calibration_terms[0]

    @property
    def Bterm(self):
        return self._calibration_terms[1]

    @property
    def Cterm(self):
        return self._calibration_terms[2]


    @property
    def filename(self):
        return Path(self._filename)

    @property
    def dir_location(self):
        return self._dir_location

    def sort_by_mz(self):
        return sorted(self, key=lambda m: m.mz_exp)

    def sort_by_abundance(self, reverse=False):
        return sorted(self, key=lambda m: m.abundance, reverse=reverse)

    @property
    def tic(self):

        return sum(self.abundance_profile)

    def check_mspeaks_warning(self):
        import warnings
        if self.mspeaks:
            pass
        else:
            warnings.warn(
                "mspeaks list is empty, continuing without filtering data"
            )

    def check_mspeaks(self):
        if self.mspeaks:
            pass
        else:
            raise Exception(
                "mspeaks list is empty, please run process_mass_spec() first"
            )

    def remove_assignment_by_index(self, indexes):
        for i in indexes: self.mspeaks[i].clear_molecular_formulas()

    def filter_by_index(self, list_indexes):

        self.mspeaks = [self.mspeaks[i] for i in range(len(self.mspeaks)) if i not in list_indexes]

        for i, mspeak in  enumerate(self.mspeaks): mspeak.index = i

        self._set_nominal_masses_start_final_indexes()

    def filter_by_mz(self, min_mz, max_mz):

        self.check_mspeaks_warning()
        indexes = [index for index, mspeak in enumerate(self.mspeaks) if min_mz <= mspeak.mz_exp <= max_mz]
        self.filter_by_index(indexes)

    def filter_by_s2n(self, min_s2n, max_s2n=False):

        self.check_mspeaks_warning()
        if not max_s2n:
            max_s2n = self.max_signal_to_noise

        self.check_mspeaks_warning()
        indexes = [index for index, mspeak in enumerate(self.mspeaks) if min_s2n <= mspeak.signal_to_noise <= max_s2n ]
        self.filter_by_index(indexes)

    def filter_by_abundance(self, min_abund, max_abund=False):

        self.check_mspeaks_warning()
        if not max_abund:
            max_abund = self.max_abundance
        indexes = [index for index, mspeak in enumerate(self.mspeaks) if min_abund <= mspeak.abundance <= max_abund]
        self.filter_by_index(indexes)

    def filter_by_max_resolving_power(self, B, T):

        rpe = lambda m, z: (1.274e7 * z * B * T)/(m*z)

        self.check_mspeaks_warning()

        indexes_to_remove = [index for index, mspeak in enumerate(self.mspeaks) if  mspeak.resolving_power >= rpe(mspeak.mz_exp,mspeak.ion_charge)]
        self.filter_by_index(indexes_to_remove)

    def filter_by_min_resolving_power(self, B, T):

        rpe = lambda m, z: (1.274e7 * z * B * T)/(m*z)

        self.check_mspeaks_warning()

        indexes_to_remove = [index for index, mspeak in enumerate(self.mspeaks) if  mspeak.resolving_power <= rpe(mspeak.mz_exp,mspeak.ion_charge)]
        self.filter_by_index(indexes_to_remove)

    def filter_by_noise_threshold(self):
        
        threshold = self.get_noise_threshold()[1][0]
        
        self.check_mspeaks_warning()
        
        indexes_to_remove = [index for index, mspeak in enumerate(self.mspeaks) if  mspeak.abundance <= threshold]
        self.filter_by_index(indexes_to_remove)

    
    def find_peaks(self):
        """needs to clear previous results from peak_picking"""
        self._mspeaks = list()
        """then do peak picking"""

        self.do_peak_picking()
        # print("A total of %i peaks were found" % len(self._mspeaks))

    def change_kendrick_base_all_mspeaks(self, kendrick_dict_base):
        """kendrick_dict_base = {"C": 1, "H": 2} or {{"C": 1, "H": 1, "O":1} etc """

        self.parameters.ms_peak.kendrick_base = kendrick_dict_base

        for mspeak in self.mspeaks:

            mspeak.change_kendrick_base(kendrick_dict_base)

    def get_nominal_mz_first_last_indexes(self, nominal_mass):

        if self._dict_nominal_masses_indexes:

            if nominal_mass in self._dict_nominal_masses_indexes.keys():

                return (self._dict_nominal_masses_indexes.get(nominal_mass)[0], self._dict_nominal_masses_indexes.get(nominal_mass)[1]+1)

            else:
                # import warnings
                # uncomment warn to distribution
                # warnings.warn("Nominal mass not found in _dict_nominal_masses_indexes, returning (0, 0) for nominal mass %i"%nominal_mass)
                return (0, 0)
        else:
            raise Exception("run process_mass_spec() function before trying to access the data")

    def get_masses_count_by_nominal_mass(self):

        dict_nominal_masses_count = {}

        all_nominal_masses = list(set([i.nominal_mz_exp for i in self.mspeaks]))

        for nominal_mass in all_nominal_masses:
            if nominal_mass not in dict_nominal_masses_count:
                dict_nominal_masses_count[nominal_mass] = len(list(self.get_nominal_mass_indexes(nominal_mass)))

        return dict_nominal_masses_count

    def datapoints_count_by_nominal_mz(self, mz_overlay=0.1):

        dict_nominal_masses_count ={}

        all_nominal_masses = list(set([i.nominal_mz_exp for i in self.mspeaks]))

        for nominal_mass in all_nominal_masses:

            if nominal_mass not in dict_nominal_masses_count:

                min_mz = nominal_mass - mz_overlay

                max_mz = nominal_mass + 1 + mz_overlay

                indexes = indexes = where((self.mz_exp_profile > min_mz) & (self.mz_exp_profile < max_mz)) 

                dict_nominal_masses_count[nominal_mass] = indexes[0].size

        return dict_nominal_masses_count

    def get_nominal_mass_indexes(self, nominal_mass, overlay=0.1):

        min_mz_to_look = nominal_mass - overlay
        max_mz_to_look = nominal_mass + 1 + overlay

        return (i for i in range(len(self.mspeaks)) if min_mz_to_look <= self.mspeaks[i].mz_exp <= max_mz_to_look)

        # indexes = (i for i in range(len(self.mspeaks)) if min_mz_to_look <= self.mspeaks[i].mz_exp <= max_mz_to_look)
        # return indexes

    def _set_nominal_masses_start_final_indexes(self):

        '''return ms peaks objs indexes(start and end) on the mass spectrum for all nominal masses'''
        dict_nominal_masses_indexes ={}

        all_nominal_masses = set(i.nominal_mz_exp for i in self.mspeaks)

        for nominal_mass in all_nominal_masses:

            indexes = self.get_nominal_mass_indexes(nominal_mass)

            defaultvalue = None
            first = last = next(indexes, defaultvalue)
            for last in indexes:
                pass

            dict_nominal_masses_indexes[nominal_mass] = (first, last)

        self._dict_nominal_masses_indexes = dict_nominal_masses_indexes

    def plot_centroid(self, ax=None, c='g'):

        import matplotlib.pyplot as plt
        if self._mspeaks:

            if ax is None:
                ax = plt.gca()

            markerline_a, stemlines_a, baseline_a = ax.stem(self.mz_exp, self.abundance, linefmt='-', markerfmt=" ", use_line_collection=True)

            plt.setp(markerline_a, 'color', c, 'linewidth', 2)
            plt.setp(stemlines_a, 'color', c, 'linewidth', 2)
            plt.setp(baseline_a, 'color', c, 'linewidth', 2)

            ax.set_xlabel("$\t{m/z}$", fontsize=12)
            ax.set_ylabel('Abundance', fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

            ax.axes.spines['top'].set_visible(False)
            ax.axes.spines['right'].set_visible(False)

            ax.get_yaxis().set_visible(False)
            ax.spines['left'].set_visible(False)

        else:

            raise Exception("No centroid data found, please run process_mass_spec")

        return ax


    def plot_profile_and_noise_threshold(self, ax=None): 

        import matplotlib.pyplot as plt
        if self.baselise_noise_std and self.baselise_noise_std:

            # x = (self.mz_exp_profile.min(), self.mz_exp_profile.max())
            baseline = (self.baselise_noise, self.baselise_noise)

            # std = self.parameters.mass_spectrum.noise_threshold_std
            # threshold = self.baselise_noise_std + (std * self.baselise_noise_std)
            x, y = self.get_noise_threshold()    
            
            if ax is None:
                ax = plt.gca()
            
            ax.plot(self.mz_exp_profile, self.abundance_profile, color="green")
            ax.plot(x, (baseline, baseline), color="yellow")
            ax.plot(x, y, color="red")

            ax.set_xlabel("$\t{m/z}$", fontsize=12)
            ax.set_ylabel('Abundance', fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

            ax.axes.spines['top'].set_visible(False)
            ax.axes.spines['right'].set_visible(False)

            ax.get_yaxis().set_visible(False)
            ax.spines['left'].set_visible(False)

        else:

            raise Exception("Calculate noise threshold first")

        return ax

    def plot_mz_domain_profile(self, color='green', ax=None): #pragma: no cover

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()
        ax.plot(self.mz_exp_profile, self.abundance_profile, color=color)
        ax.set(xlabel='m/z', ylabel='abundance')

        return ax

    def to_excel(self, out_file_path):

        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_excel()

    def to_hdf(self, out_file_path):
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_hdf()

    def to_csv(self, out_file_path):
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_csv()

    def to_pandas(self, out_file_path):
        # pickle dataframe (pkl extension)
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_pandas()

    def to_dataframe(self,):
        # returns pandas dataframe

        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.get_pandas_df()

    def to_json(self):
        # returns pandas dataframe

        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.to_json()

    def parameters_json(self):

        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.parameters_to_json()

class MassSpecProfile(MassSpecBase):
    '''
    - A iterative mass spectrum class when the entry point is on profile format
    - Stores the profile data and instrument settings
    - Iteration over a list of MSPeaks classes stored at the _mspeaks attributes
    - _mspeaks is populated under the hood by calling process_mass_spec method
    - iteration is null if _mspeaks is empty

    Parameters
    ----------
    dataframe : pandas Dataframe(Series(floats))
        contains columns [m/z, Abundance, Resolving Power, S/N] 
    d_params : dict{'str': float, int or str}

    Attributes
    ----------
    _mz_exp : list(float)
        This is where we store mz_exp,
    _abundance : list(float)     
        This is where we store _abundance,
    _mspeaks : list(MSPeak)
        store MSpeaks objects identified by a peak picking algorithm

    Relevant Methods
    ----------
    process_mass_spec()
        find or set the noise threshold base on the setting encapsulated at settings.input.ProcessingSetting.MassSpectrumSetting
        - run the peak peaking algorithm and use the method addMSPeaks() to populate _mspeaks attribute

    see also: MassSpecBase(), MassSpecfromFreq(), MassSpecProfile()
    '''

    def __init__(self, data_dict, d_params, auto_process=True, auto_noise=True, noise_bayes_est=False):
        """
        method docs
        """
        # print(data_dict.keys())
        super().__init__(data_dict.get(Labels.mz), data_dict.get(Labels.abundance), d_params)
       
        if auto_process:
            self.process_mass_spec(auto_noise=auto_noise, noise_bayes_est=noise_bayes_est)

class MassSpecfromFreq(MassSpecBase):
    '''
    - A iterative mass spectrum class when data entry is on frequency(Hz) domain
    - Transform to m/z based on the settings stored at d_params
    - Stores the profile data and instrument settings
    - Iteration over a list of MSPeaks classes stored at the _mspeaks attributes
    - _mspeaks is populated under the hood by calling process_mass_spec method
    - iteration is null if _mspeaks is empty

    Parameters
    ----------
    frequency_domain : list(float)
        all datapoints in frequency domain in Hz
    magnitude :  frequency_domain : list(float)
        all datapoints in for magnitude of each frequency datapoint

    Attributes
    ----------
    _mz_exp : list(float)
        This is where we store mz_exp,
    _frequency_domain : list(float)
        This is where we store _frequency_domain,
    _abundance : list(float)     
        This is where we store _abundance,
    _mspeaks : list(MSPeak)
        store MSpeaks objects identified by a peak picking algorithm     
    label : str
        store label (Bruker, Midas Transient, see Labels class ). It across distinct processing points
    
    Relevant Methods
    ----------
    _set_mz_domain()
        calculates the m_z based on the setting of d_params

    process_mass_spec()
        find or set the noise threshold base on the setting encapsulated at settings.input.ProcessingSetting.MassSpectrumSetting
        - run the peak peaking algorithm and use the method addMSPeaks() to populate _mspeaks attribute
    
    see also: MassSpecBase(), MassSpecfromFreq(), MassSpecProfile()
    '''

    def __init__(self, frequency_domain, magnitude, d_params, 
                auto_process=True, keep_profile=True, auto_noise=False, noise_bayes_est=False):
        """
        method docs
        """
        super().__init__(None, magnitude, d_params)

        self._frequency_domain = frequency_domain
        self.has_frequency = True
        self._set_mz_domain()

        
        """ use this call to automatically process data as the object is created, Setting need to be changed before initiating the class to be in effect"""
        
        if auto_process:
            self.process_mass_spec(keep_profile=keep_profile, auto_noise=auto_noise, noise_bayes_est=noise_bayes_est)


    def _set_mz_domain(self):

        if self.label == Labels.bruker_frequency:

            self._mz_exp = self._f_to_mz_bruker()

        else:

            self._mz_exp = self._f_to_mz()

    @property
    def transient_settings(self): return self.parameters.transient

    @transient_settings.setter
    def transient_settings(self, instance_TransientSetting):
     
        self.parameters.transient = instance_TransientSetting  


class MassSpecCentroid(MassSpecBase):

    '''
    - A iterative mass spectrum class when data entry is centroid mode
    - Stores the centroid data and instrument settings
    - Simulate profile data based on Gaussian or Lorentzian peak shape
    - Iteration over a list of MSPeaks classes stored at the _mspeaks attributes
    - _mspeaks is populated under the hood by calling process_mass_spec method
    - iteration is null if _mspeaks is empty

    Parameters
    ----------
    data_dict : dict {string: numpy array float64 )
        contains keys [m/z, Abundance, Resolving Power, S/N] 
    d_params : dict{'str': float, int or str}
        
    Attributes
    ----------
    _mz_exp : list(float)
        This is where we store mz_exp,
    _abundance : list(float)     
        This is where we store _abundance,
    _mspeaks : list(MSPeak)
        store MSpeaks objects identified by a peak picking algorithm  

    Attributes
    ----------
    _mz_exp : list(float)
        This is where we store mz_exp,
    _frequency_domain : list(float)
        This is where we store _frequency_domain,
    _abundance : list(float)     
        This is where we store _abundance,
    _mspeaks : list(MSPeak)
        store MSpeaks objects identified by a peak picking algorithm     
    label : str
        store label (Bruker, Midas Transient, see Labels class)
    
    Relevant Methods
    ----------
    _set_mz_domain()
        calculates the m_z based on the setting of d_params

    process_mass_spec()
        - overrides the base class function
        - Populates _mspeaks list with MSpeaks class using the centroid date

    see also: MassSpecBase(), MassSpecfromFreq(), MassSpecProfile()
    '''

    def __init__(self, data_dict, d_params):

        """needs to simulate peak shape and pass as mz_exp and magnitude."""

        super().__init__([], [], d_params)

        self._set_parameters_objects(d_params)
        
        if self.label == Labels.thermo_centroid:
            self._baselise_noise = d_params.get("baselise_noise")
            self._baselise_noise_std = d_params.get("baselise_noise_std")

        self.is_centroid = True

        self.process_mass_spec(data_dict)

    def __simulate_profile__data__(self, exp_mz_centroid, magnitude_centroid):
        '''needs theoretical resolving power calculation and define peak shape
        this is a quick fix to trick a line plot be able to plot as sticks'''

        x, y = [], []
        for i in range(len(exp_mz_centroid)):
            x.append(exp_mz_centroid[i] - 0.0000001)
            x.append(exp_mz_centroid[i])
            x.append(exp_mz_centroid[i] + 0.0000001)
            y.append(0)
            y.append(magnitude_centroid[i])
            y.append(0)
        return x, y

    @property
    def mz_exp_profile(self):

        mz_list = []
        for mz in self.mz_exp:
            mz_list.append(mz - 0.0000001)
            mz_list.append(mz)
            mz_list.append(mz + 0.0000001)
        return mz_list
    
    @mz_exp_profile.setter
    def mz_exp_profile(self, _mz_exp ): self._mz_exp = _mz_exp

    @property
    def abundance_profile(self):
        ab_list = []
        for ab in self.abundance:
            ab_list.append(0)
            ab_list.append(ab)
            ab_list.append(0)
        return ab_list

    @abundance_profile.setter
    def abundance_profile(self, abundance ): self._abundance = abundance

    @property
    def tic(self):

        return sum(self.abundance)

    def process_mass_spec(self, data_dict, auto_noise=True, noise_bayes_est=False):
        import tqdm
        # overwrite process_mass_spec 
        # mspeak objs are usually added inside the PeaKPicking class 
        # for profile and freq based data
        

        s2n = True
        ion_charge = self.polarity
        #l_exp_mz_centroid = data_dict.get(Labels.mz)
        #l_intes_centr = data_dict.get(Labels.abundance)
        #l_peak_resolving_power = data_dict.get(Labels.rp)
        l_s2n = data_dict.get(Labels.s2n)
        
        if not l_s2n: s2n = False
        
        print("Loading mass spectrum object")
        
        for index, mz in enumerate(data_dict.get(Labels.mz)):
            
            # centroid peak does not have start and end peak index pos
            massspec_indexes = (index, index, index)
            
            if s2n:
                
                self.add_mspeak(
                    ion_charge,
                    mz,
                    data_dict.get(Labels.abundance)[index],
                    data_dict.get(Labels.rp)[index],
                    l_s2n[index],
                    massspec_indexes,
                    ms_parent=self
                )

            else:

                self.add_mspeak(
                    ion_charge,
                    mz,
                    data_dict.get(Labels.abundance)[index],
                    data_dict.get(Labels.rp)[index],
                    -999,
                    massspec_indexes,
                    ms_parent=self
                )

        self.mspeaks = self._mspeaks
        self._dynamic_range = self.max_abundance / self.min_abundance
        self._set_nominal_masses_start_final_indexes()
        
        if self.label != Labels.thermo_centroid:
            self._baselise_noise, self._baselise_noise_std = self.run_noise_threshold_calc(auto=auto_noise, bayes=noise_bayes_est)
            
class MassSpecCentroidLowRes(MassSpecCentroid,):
    
    '''Does not store MSPeak Objs, will iterate over mz, abundance pairs instead'''
    
    def __init__(self, data_dict, d_params):
    
        self._set_parameters_objects(d_params)
        self._mz_exp = array(data_dict.get(Labels.mz))
        self._abundance = array(data_dict.get(Labels.abundance))
        self._processed_tic = None
    
    def __len__(self):
        
        return len(self.mz_exp)
        
    def __getitem__(self, position):
        
        return (self.mz_exp[position], self.abundance[position])

    @property
    def mz_exp(self):
        
        return self._mz_exp 

    @property
    def abundance(self):

        return self._abundance

    @property
    def processed_tic(self):
    
        return sum(self._processed_tic)
    
    @property
    def tic(self):
    
        if self._processed_tic:
            return self._processed_tic
        else:
            return sum(self.abundance)
    
    @property
    def mz_abun_tuples(self):

        r = lambda x: ( int(round(x[0],0), int(round(x[1],0))) )

        return [r(i) for i in self]
    
    @property
    def mz_abun_dict(self):
        
        r = lambda x: int(round(x,0))
            
        return { r(i[0]):r(i[1]) for i in self}
    
    '''
    def to_excel(self, out_file_path):
    
    
    exportMS= HighResMassSpecExport(out_file_path, self)
    exportMS.to_excel()


    def to_hdf(self, out_file_path):
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS= HighResMassSpecExport(out_file_path, self)
        exportMS.to_hdf()

    def to_csv(self, out_file_path):
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS= HighResMassSpecExport(out_file_path, self)
        exportMS.to_csv()
        
    def to_pandas(self, out_file_path):
        #pickle dataframe (pkl extension)
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS= HighResMassSpecExport(out_file_path, self)
        exportMS.to_pandas()

    def to_dataframe(self,):
        #returns pandas dataframe
        
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS= HighResMassSpecExport(self.filename, self)
        return exportMS.get_pandas_df()
       ''' 


