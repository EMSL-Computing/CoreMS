from pathlib import Path
from copy import deepcopy


#from matplotlib import rcParamsDefault, rcParams
from numpy import array, power, float64, where, histogram, trapz

from pandas import DataFrame
from lmfit.models import GaussianModel

from corems.mass_spectrum.calc.MassSpectrumCalc import MassSpecCalc
from corems.mass_spectrum.calc.KendrickGroup import KendrickGrouping
from corems.encapsulation.constant import Labels
from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak as MSPeak
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.input.parameter_from_json import load_and_set_parameters_ms, load_and_set_toml_parameters_ms
from corems.mass_spectrum.calc.MeanResolvingPowerFilter import MeanResolvingPowerFilter

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

def overrides(interface_class):
    """Checks if the method overrides a method from an interface class."""
    def overrider(method):
        assert method.__name__ in dir(interface_class)
        return method
    return overrider

class MassSpecBase(MassSpecCalc, KendrickGrouping):
    """A mass spectrum base class, stores the profile data and instrument settings.

    Iteration over a list of MSPeaks classes stored at the _mspeaks attributes.
    _mspeaks is populated under the hood by calling process_mass_spec method.
    Iteration is null if _mspeaks is empty.

    Parameters
    ----------
    mz_exp : array_like
        The m/z values of the mass spectrum.
    abundance : array_like
        The abundance values of the mass spectrum.
    d_params : dict
        A dictionary of parameters for the mass spectrum.
    **kwargs
        Additional keyword arguments.

    Attributes
    ----------

    mspeaks : list
        A list of mass peaks.
    is_calibrated : bool
        Whether the mass spectrum is calibrated.
    is_centroid : bool
        Whether the mass spectrum is centroided.
    has_frequency : bool
        Whether the mass spectrum has a frequency domain.
    calibration_order : None or int
        The order of the mass spectrum's calibration.
    calibration_points : None or ndarray
        The calibration points of the mass spectrum.
    calibration_RMS : None or float
        The root mean square of the mass spectrum's calibration.
    calibration_segment : None or CalibrationSegment
        The calibration segment of the mass spectrum.
    _abundance : ndarray
        The abundance values of the mass spectrum.
    _mz_exp : ndarray
        The m/z values of the mass spectrum.
    _mspeaks : list
        A list of mass peaks.
    _dict_nominal_masses_indexes : dict
        A dictionary of nominal masses and their indexes.
    _baseline_noise : float
        The baseline noise of the mass spectrum.
    _baseline_noise_std : float
        The standard deviation of the baseline noise of the mass spectrum.
    _dynamic_range : float or None
        The dynamic range of the mass spectrum.
    _transient_settings : None or TransientSettings
        The transient settings of the mass spectrum.
    _frequency_domain : None or FrequencyDomain
        The frequency domain of the mass spectrum.
    _mz_cal_profile : None or MzCalibrationProfile
        The m/z calibration profile of the mass spectrum.

    Methods
    -------
    * process_mass_spec(). Main function to process the mass spectrum, 
    including calculating the noise threshold, peak picking, and resetting the MSpeak indexes.

    See also: MassSpecCentroid(), MassSpecfromFreq(), MassSpecProfile()
    """
    def __init__(self, mz_exp, abundance, d_params, **kwargs):
        
        self._abundance = array(abundance, dtype=float64)
        self._mz_exp = array(mz_exp, dtype=float64)
                    
        # objects created after process_mass_spec() function
        self._mspeaks = list()
        self.mspeaks = list()
        self._dict_nominal_masses_indexes = dict()
        self._baseline_noise = 0.001
        self._baseline_noise_std = 0.001
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
        self.calibration_segment = None

    def _init_settings(self):
        """Initializes the settings for the mass spectrum."""
        self._parameters = MSParameters()

    def __len__(self):

        return len(self.mspeaks)

    def __getitem__(self, position) -> MSPeak:

        return self.mspeaks[position]

    def set_indexes(self, list_indexes):
        """Set the mass spectrum to iterate over only the selected MSpeaks indexes.

        Parameters
        ----------
        list_indexes : list of int
            A list of integers representing the indexes of the MSpeaks to iterate over.

        """
        self.mspeaks = [self._mspeaks[i] for i in list_indexes]

        for i, mspeak in  enumerate(self.mspeaks): mspeak.index = i

        self._set_nominal_masses_start_final_indexes()

    def reset_indexes(self):
        """Reset the mass spectrum to iterate over all MSpeaks objects.

        This method resets the mass spectrum to its original state, allowing iteration over all MSpeaks objects.
        It also sets the index of each MSpeak object to its corresponding position in the mass spectrum.

        """
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
        """Add a new MSPeak object to the MassSpectrum object.

        Parameters
        ----------
        ion_charge : int
            The ion charge of the MSPeak.
        mz_exp : float
            The experimental m/z value of the MSPeak.
        abundance : float
            The abundance of the MSPeak.
        resolving_power : float
            The resolving power of the MSPeak.
        signal_to_noise : float
            The signal-to-noise ratio of the MSPeak.
        massspec_indexes : list
            A list of indexes of the MSPeak in the MassSpectrum object.
        exp_freq : float, optional
            The experimental frequency of the MSPeak. Defaults to None.
        ms_parent : MSParent, optional
            The MSParent object associated with the MSPeak. Defaults to None.
        """
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
        """Set the parameters of the MassSpectrum object.

        Parameters
        ----------
        d_params : dict
            A dictionary containing the parameters to set.

        Notes
        -----
        This method sets the following parameters of the MassSpectrum object:
        - _calibration_terms
        - label
        - analyzer
        - acquisition_time
        - instrument_label
        - polarity
        - scan_number
        - retention_time
        - mobility_rt
        - mobility_scan
        - _filename
        - _dir_location
        - _baseline_noise
        - _baseline_noise_std
        - sample_name
        """
        self._calibration_terms = (
            d_params.get("Aterm"),
            d_params.get("Bterm"),
            d_params.get("Cterm"),
        )

        self.label = d_params.get(Labels.label)

        self.analyzer = d_params.get('analyzer')

        self.acquisition_time = d_params.get('acquisition_time')

        self.instrument_label = d_params.get('instrument_label')

        self.polarity = int(d_params.get("polarity"))

        self.scan_number = d_params.get("scan_number")

        self.retention_time = d_params.get("rt")

        self.mobility_rt = d_params.get("mobility_rt")

        self.mobility_scan = d_params.get("mobility_scan")

        self._filename = d_params.get("filename_path")

        self._dir_location = d_params.get("dir_location")

        self._baseline_noise = d_params.get("baseline_noise")

        self._baseline_noise_std = d_params.get("baseline_noise_std")

        if d_params.get('sample_name') != 'Unknown':

            self.sample_name = d_params.get('sample_name')
            if not self.sample_name:
                self.sample_name = self.filename.stem
        else:

            self.sample_name = self.filename.stem

    def reset_cal_therms(self, Aterm, Bterm, C, fas=0):
        """Reset calibration terms and recalculate the mass-to-charge ratio and abundance.

        Parameters
        ----------
        Aterm : float
            The A-term calibration coefficient.
        Bterm : float
            The B-term calibration coefficient.
        C : float
            The C-term calibration coefficient.
        fas : float, optional
            The frequency amplitude scaling factor. Default is 0.
        """
        self._calibration_terms = (Aterm, Bterm, C)

        self._mz_exp = self._f_to_mz()
        self._abundance = self._abundance
        self.find_peaks()
        self.reset_indexes()

    def clear_molecular_formulas(self):
        """Clear the molecular formulas for all mspeaks in the MassSpectrum.

        Returns
        -------
        numpy.ndarray
            An array of the cleared molecular formulas for each mspeak in the MassSpectrum.
        """
        self.check_mspeaks()
        return array([mspeak.clear_molecular_formulas() for mspeak in self.mspeaks])

    def process_mass_spec(self, keep_profile=True):
        """Process the mass spectrum.

        Parameters
        ----------
        keep_profile : bool, optional
            Whether to keep the profile data after processing. Defaults to True.

        Notes
        -----
        This method does the following:
        - calculates the noise threshold
        - does peak picking (creates mspeak_objs)
        - resets the mspeak_obj indexes
        """
        
        # if runned mannually make sure to rerun filter_by_noise_threshold     
        # calculates noise threshold 
        # do peak picking( create mspeak_objs) 
        # reset mspeak_obj the indexes
         
        self.cal_noise_threshold()

        self.find_peaks()
        self.reset_indexes()

        if self.mspeaks:
            self._dynamic_range = self.max_abundance / self.min_abundance
        else:
            self._dynamic_range = 0
        if not keep_profile:

            self._abundance *= 0
            self._mz_exp *= 0
            

    def cal_noise_threshold(self):
        """Calculate the noise threshold of the mass spectrum.

        """

        if self.label == Labels.simulated_profile:

            self._baseline_noise, self._baseline_noise_std = 0.1, 1

        if self.settings.noise_threshold_method == 'log':

            self._baseline_noise, self._baseline_noise_std = self.run_log_noise_threshold_calc()

        else:
            self._baseline_noise, self._baseline_noise_std = self.run_noise_threshold_calc()

    @property
    def parameters(self):
        """Return the parameters of the mass spectrum."""
        return self._parameters

    @parameters.setter
    def parameters(self, instance_MSParameters):
        self._parameters = instance_MSParameters

    def set_parameter_from_json(self, parameters_path):
        """Set the parameters of the mass spectrum from a JSON file.
        
        Parameters
        ----------
        parameters_path : str
            The path to the JSON file containing the parameters.
        """
        load_and_set_parameters_ms(self, parameters_path=parameters_path)    

    def set_parameter_from_toml(self, parameters_path):
        load_and_set_toml_parameters_ms(self, parameters_path=parameters_path)    

    @property
    def mspeaks_settings(self): 
        """Return the MS peak settings of the mass spectrum."""
        return self.parameters.ms_peak

    @mspeaks_settings.setter
    def mspeaks_settings(self, instance_MassSpecPeakSetting):

            self.parameters.ms_peak = instance_MassSpecPeakSetting

    @property
    def settings(self): 
        """Return the settings of the mass spectrum."""
        return self.parameters.mass_spectrum

    @settings.setter
    def settings(self, instance_MassSpectrumSetting):

        self.parameters.mass_spectrum =  instance_MassSpectrumSetting

    @property
    def molecular_search_settings(self):  
        """Return the molecular search settings of the mass spectrum."""
        return self.parameters.molecular_search

    @molecular_search_settings.setter
    def molecular_search_settings(self, instance_MolecularFormulaSearchSettings):

        self.parameters.molecular_search =  instance_MolecularFormulaSearchSettings

    @property
    def mz_cal_profile(self):
        """Return the calibrated m/z profile of the mass spectrum."""
        return self._mz_cal_profile

    @mz_cal_profile.setter
    def mz_cal_profile(self, mz_cal_list):
        
        if len(mz_cal_list) == len(self._mz_exp):
            self._mz_cal_profile = mz_cal_list
        else:
            raise Exception( "calibrated array (%i) is not of the same size of the data (%i)" % (len(mz_cal_list),  len(self.mz_exp_profile)))    

    @property
    def mz_cal(self):
        """Return the calibrated m/z values of the mass spectrum."""
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
        """Return the experimental m/z values of the mass spectrum."""
        self.check_mspeaks()

        if self.is_calibrated:

            return array([mspeak.mz_cal for mspeak in self.mspeaks])

        else:

            return array([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def freq_exp_profile(self):
        """Return the experimental frequency profile of the mass spectrum."""
        return self._frequency_domain
    
    @freq_exp_profile.setter
    def freq_exp_profile(self, new_data): self._frequency_domain = array(new_data)

    @property
    def mz_exp_profile(self): 
        """Return the experimental m/z profile of the mass spectrum."""
        if self.is_calibrated: 
            return self.mz_cal_profile
        else:
            return self._mz_exp

    @mz_exp_profile.setter
    def mz_exp_profile(self, new_data ): self._mz_exp = array(new_data)

    @property
    def abundance_profile(self): 
        """Return the abundance profile of the mass spectrum."""
        return self._abundance

    @abundance_profile.setter
    def abundance_profile(self, new_data): self._abundance = array(new_data)

    @property
    def abundance(self):
        """Return the abundance values of the mass spectrum."""
        self.check_mspeaks()
        return array([mspeak.abundance for mspeak in self.mspeaks])

    def freq_exp(self):
        """Return the experimental frequency values of the mass spectrum."""
        self.check_mspeaks()
        return array([mspeak.freq_exp for mspeak in self.mspeaks])

    @property
    def resolving_power(self):
        """Return the resolving power values of the mass spectrum."""
        self.check_mspeaks()
        return array([mspeak.resolving_power for mspeak in self.mspeaks])

    @property
    def signal_to_noise(self):
        self.check_mspeaks()
        return array([mspeak.signal_to_noise for mspeak in self.mspeaks])

    @property
    def nominal_mz(self):
        """Return the nominal m/z values of the mass spectrum."""
        if self._dict_nominal_masses_indexes:
            return sorted(list(self._dict_nominal_masses_indexes.keys()))
        else:
            raise ValueError("Nominal indexes not yet set")    

    def get_mz_and_abundance_peaks_tuples(self):
        """Return a list of tuples containing the m/z and abundance values of the mass spectrum."""
        self.check_mspeaks()
        return [(mspeak.mz_exp, mspeak.abundance) for mspeak in self.mspeaks]

    @property
    def kmd(self):
        """Return the Kendrick mass defect values of the mass spectrum."""
        self.check_mspeaks()
        return array([mspeak.kmd for mspeak in self.mspeaks])

    @property
    def kendrick_mass(self):
        """Return the Kendrick mass values of the mass spectrum."""
        self.check_mspeaks()
        return array([mspeak.kendrick_mass for mspeak in self.mspeaks])

    @property
    def max_mz_exp(self):
        """Return the maximum experimental m/z value of the mass spectrum."""
        return max([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def min_mz_exp(self):
        """Return the minimum experimental m/z value of the mass spectrum."""
        return min([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def max_abundance(self):
        """Return the maximum abundance value of the mass spectrum."""        
        return max([mspeak.abundance for mspeak in self.mspeaks])

    @property
    def max_signal_to_noise(self):
        """Return the maximum signal-to-noise ratio of the mass spectrum."""
        return max([mspeak.signal_to_noise for mspeak in self.mspeaks])

    @property
    def most_abundant_mspeak(self):
        """Return the most abundant MSpeak object of the mass spectrum."""
        return max(self.mspeaks, key=lambda m: m.abundance)

    @property
    def min_abundance(self):
        """Return the minimum abundance value of the mass spectrum."""
        return min([mspeak.abundance for mspeak in self.mspeaks])

    # takes too much cpu time 
    @property
    def dynamic_range(self):
        """Return the dynamic range of the mass spectrum."""
        return self._dynamic_range

    @property
    def baseline_noise(self):
        """Return the baseline noise of the mass spectrum."""
        if self._baseline_noise:
            return self._baseline_noise
        else:     
            return None

    @property
    def baseline_noise_std(self):
        """Return the standard deviation of the baseline noise of the mass spectrum."""
        if self._baseline_noise_std == 0:
            return self._baseline_noise_std
        if self._baseline_noise_std:
            return self._baseline_noise_std
        else:     
            return None

    @property
    def Aterm(self):
        """Return the A-term calibration coefficient of the mass spectrum."""
        return self._calibration_terms[0]

    @property
    def Bterm(self):
        """Return the B-term calibration coefficient of the mass spectrum."""
        return self._calibration_terms[1]

    @property
    def Cterm(self):
        """Return the C-term calibration coefficient of the mass spectrum."""
        return self._calibration_terms[2]

    @property
    def filename(self):
        """Return the filename of the mass spectrum."""
        return Path(self._filename)

    @property
    def dir_location(self):
        """Return the directory location of the mass spectrum."""
        return self._dir_location

    def sort_by_mz(self):
        """Sort the mass spectrum by m/z values."""
        return sorted(self, key=lambda m: m.mz_exp)

    def sort_by_abundance(self, reverse=False):
        """Sort the mass spectrum by abundance values."""
        return sorted(self, key=lambda m: m.abundance, reverse=reverse)

    @property
    def tic(self):
        """Return the total ion current of the mass spectrum."""
        return trapz(self.abundance_profile, self.mz_exp_profile)

    def check_mspeaks_warning(self):
        """Check if the mass spectrum has MSpeaks objects.
        
        Raises
        ------
        Warning
            If the mass spectrum has no MSpeaks objects.
        """
        import warnings
        if self.mspeaks:
            pass
        else:
            warnings.warn(
                "mspeaks list is empty, continuing without filtering data"
            )

    def check_mspeaks(self):
        """Check if the mass spectrum has MSpeaks objects.

        Raises
        ------
        Exception
            If the mass spectrum has no MSpeaks objects.
        """
        if self.mspeaks:
            pass
        else:
            raise Exception(
                "mspeaks list is empty, please run process_mass_spec() first"
            )

    def remove_assignment_by_index(self, indexes):
        """Remove the molecular formula assignment of the MSpeaks objects at the specified indexes.

        Parameters
        ----------
        indexes : list of int
            A list of indexes of the MSpeaks objects to remove the molecular formula assignment from.
        """
        for i in indexes: self.mspeaks[i].clear_molecular_formulas()

    def filter_by_index(self, list_indexes):
        """Filter the mass spectrum by the specified indexes.

        Parameters
        ----------
        list_indexes : list of int
            A list of indexes of the MSpeaks objects to keep.

        """

        self.mspeaks = [self.mspeaks[i] for i in range(len(self.mspeaks)) if i not in list_indexes]

        for i, mspeak in  enumerate(self.mspeaks): mspeak.index = i

        self._set_nominal_masses_start_final_indexes()

    def filter_by_mz(self, min_mz, max_mz):
        """Filter the mass spectrum by the specified m/z range.

        Parameters
        ----------
        min_mz : float
            The minimum m/z value to keep.
        max_mz : float
            The maximum m/z value to keep.

        """      
        self.check_mspeaks_warning()
        indexes = [index for index, mspeak in enumerate(self.mspeaks) if not min_mz <= mspeak.mz_exp <= max_mz]
        self.filter_by_index(indexes)

    def filter_by_s2n(self, min_s2n, max_s2n=False):
        """Filter the mass spectrum by the specified signal-to-noise ratio range.

        Parameters
        ----------
        min_s2n : float
            The minimum signal-to-noise ratio to keep.
        max_s2n : float, optional
            The maximum signal-to-noise ratio to keep. Defaults to False (no maximum).

        """
        self.check_mspeaks_warning()
        if max_s2n:
            indexes = [index for index, mspeak in enumerate(self.mspeaks) if not min_s2n <= mspeak.signal_to_noise <= max_s2n ]
        else:
            indexes = [index for index, mspeak in enumerate(self.mspeaks) if mspeak.signal_to_noise <= min_s2n ]
        self.filter_by_index(indexes)

    def filter_by_abundance(self, min_abund, max_abund=False):
        """Filter the mass spectrum by the specified abundance range.

        Parameters
        ----------
        min_abund : float
            The minimum abundance to keep.
        max_abund : float, optional
            The maximum abundance to keep. Defaults to False (no maximum).

        """
        self.check_mspeaks_warning()
        if max_abund:
            indexes = [index for index, mspeak in enumerate(self.mspeaks) if not min_abund <= mspeak.abundance <= max_abund]
        else:
            indexes = [index for index, mspeak in enumerate(self.mspeaks) if mspeak.abundance <= min_abund]
        self.filter_by_index(indexes)

    def filter_by_max_resolving_power(self, B, T):
        """Filter the mass spectrum by the specified maximum resolving power.
        
        Parameters
        ----------
        B : float
        T : float
        
        """

        rpe = lambda m, z: (1.274e7 * z * B * T)/(m*z)

        self.check_mspeaks_warning()

        indexes_to_remove = [index for index, mspeak in enumerate(self.mspeaks) if  mspeak.resolving_power >= rpe(mspeak.mz_exp,mspeak.ion_charge)]
        self.filter_by_index(indexes_to_remove)

    def filter_by_mean_resolving_power(self, ndeviations=3,plot=False,guess_pars=False):
        """Filter the mass spectrum by the specified mean resolving power.

        Parameters
        ----------
        ndeviations : float, optional
            The number of standard deviations to use for filtering. Defaults to 3.
        plot : bool, optional
            Whether to plot the resolving power distribution. Defaults to False.
        guess_pars : bool, optional
            Whether to guess the parameters for the Gaussian model. Defaults to False.

        """
        self.check_mspeaks_warning()
        indexes_to_remove = MeanResolvingPowerFilter(self,ndeviations,plot,guess_pars).main()
        self.filter_by_index(indexes_to_remove)


    def filter_by_min_resolving_power(self, B, T):
        """Filter the mass spectrum by the specified minimum resolving power.

        Parameters
        ----------
        B : float
        T : float

        """
        rpe = lambda m, z: (1.274e7 * z * B * T)/(m*z)

        self.check_mspeaks_warning()

        indexes_to_remove = [index for index, mspeak in enumerate(self.mspeaks) if  mspeak.resolving_power <= rpe(mspeak.mz_exp,mspeak.ion_charge)]
        self.filter_by_index(indexes_to_remove)

    def filter_by_noise_threshold(self):
        """Filter the mass spectrum by the noise threshold."""
        
        threshold = self.get_noise_threshold()[1][0]
        
        self.check_mspeaks_warning()
        
        indexes_to_remove = [index for index, mspeak in enumerate(self.mspeaks) if  mspeak.abundance <= threshold]
        self.filter_by_index(indexes_to_remove)

    
    def find_peaks(self):
        """Find the peaks of the mass spectrum."""
        #needs to clear previous results from peak_picking
        self._mspeaks = list()

        #then do peak picking
        self.do_peak_picking()
        # print("A total of %i peaks were found" % len(self._mspeaks))

    def change_kendrick_base_all_mspeaks(self, kendrick_dict_base):
        """Change the Kendrick base of all MSpeaks objects.

        Parameters
        ----------
        kendrick_dict_base : dict
            A dictionary of the Kendrick base to change to.

        Notes
        -----
        Example of kendrick_dict_base parameter: kendrick_dict_base = {"C": 1, "H": 2} or {"C": 1, "H": 1, "O":1} etc
        """
        self.parameters.ms_peak.kendrick_base = kendrick_dict_base

        for mspeak in self.mspeaks:

            mspeak.change_kendrick_base(kendrick_dict_base)

    def get_nominal_mz_first_last_indexes(self, nominal_mass):
        """Return the first and last indexes of the MSpeaks objects with the specified nominal mass.

        Parameters
        ----------
        nominal_mass : int
            The nominal mass to get the indexes for.

        Returns
        -------
        tuple
            A tuple containing the first and last indexes of the MSpeaks objects with the specified nominal mass.
        """
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
        """Return a dictionary of the nominal masses and their counts."""

        dict_nominal_masses_count = {}

        all_nominal_masses = list(set([i.nominal_mz_exp for i in self.mspeaks]))

        for nominal_mass in all_nominal_masses:
            if nominal_mass not in dict_nominal_masses_count:
                dict_nominal_masses_count[nominal_mass] = len(list(self.get_nominal_mass_indexes(nominal_mass)))

        return dict_nominal_masses_count

    def datapoints_count_by_nominal_mz(self, mz_overlay=0.1):
        """Return a dictionary of the nominal masses and their counts.

        Parameters
        ----------
        mz_overlay : float, optional
            The m/z overlay to use for counting. Defaults to 0.1.

        Returns
        -------
        dict
            A dictionary of the nominal masses and their counts.
        """
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
        """Return the indexes of the MSpeaks objects with the specified nominal mass.

        Parameters
        ----------
        nominal_mass : int
            The nominal mass to get the indexes for.
        overlay : float, optional
            The m/z overlay to use for counting. Defaults to 0.1.

        Returns
        -------
        generator
            A generator of the indexes of the MSpeaks objects with the specified nominal mass.
        """       
        min_mz_to_look = nominal_mass - overlay
        max_mz_to_look = nominal_mass + 1 + overlay

        return (i for i in range(len(self.mspeaks)) if min_mz_to_look <= self.mspeaks[i].mz_exp <= max_mz_to_look)

        # indexes = (i for i in range(len(self.mspeaks)) if min_mz_to_look <= self.mspeaks[i].mz_exp <= max_mz_to_look)
        # return indexes

    def _set_nominal_masses_start_final_indexes(self):
        """Set the start and final indexes of the MSpeaks objects for all nominal masses."""
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
        """Plot the centroid data of the mass spectrum.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The matplotlib axes to plot on. Defaults to None.
        c : str, optional
            The color to use for the plot. Defaults to 'g' (green).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Raises
        ------
        Exception
            If no centroid data is found.
        """

        import matplotlib.pyplot as plt
        if self._mspeaks:

            if ax is None:
                ax = plt.gca()

            markerline_a, stemlines_a, baseline_a = ax.stem(self.mz_exp, self.abundance, linefmt='-', markerfmt=" ")

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

    def plot_profile_and_noise_threshold(self, ax=None,legend=False): 
        """Plot the profile data and noise threshold of the mass spectrum.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The matplotlib axes to plot on. Defaults to None.
        legend : bool, optional
            Whether to show the legend. Defaults to False.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Raises
        ------
        Exception
            If no noise threshold is found.
        """
        import matplotlib.pyplot as plt
        if self.baseline_noise_std and self.baseline_noise_std:

            # x = (self.mz_exp_profile.min(), self.mz_exp_profile.max())
            baseline = (self.baseline_noise, self.baseline_noise)

            # std = self.parameters.mass_spectrum.noise_threshold_min_std
            # threshold = self.baseline_noise_std + (std * self.baseline_noise_std)
            x, y = self.get_noise_threshold()    
            
            if ax is None:
                ax = plt.gca()
            
            ax.plot(self.mz_exp_profile, self.abundance_profile, color="green",label="Spectrum")
            ax.plot(x, (baseline, baseline), color="yellow",label="Baseline Noise")
            ax.plot(x, y, color="red",label="Noise Threshold")

            ax.set_xlabel("$\t{m/z}$", fontsize=12)
            ax.set_ylabel('Abundance', fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

            ax.axes.spines['top'].set_visible(False)
            ax.axes.spines['right'].set_visible(False)

            ax.get_yaxis().set_visible(False)
            ax.spines['left'].set_visible(False)
            if legend:
                ax.legend()

        else:

            raise Exception("Calculate noise threshold first")

        return ax

    def plot_mz_domain_profile(self, color='green', ax=None): 
        """Plot the m/z domain profile of the mass spectrum.

        Parameters
        ----------
        color : str, optional
            The color to use for the plot. Defaults to 'green'.
        ax : matplotlib.axes.Axes, optional
            The matplotlib axes to plot on. Defaults to None.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """       

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()
        ax.plot(self.mz_exp_profile, self.abundance_profile, color=color)
        ax.set(xlabel='m/z', ylabel='abundance')

        return ax

    def to_excel(self, out_file_path, write_metadata=True):
        """Export the mass spectrum to an Excel file.

        Parameters
        ----------
        out_file_path : str
            The path to the Excel file to export to.
        write_metadata : bool, optional
            Whether to write the metadata to the Excel file. Defaults to True.

        Returns
        -------
        None
        """
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_excel(write_metadata=write_metadata)

    def to_hdf(self, out_file_path):
        """Export the mass spectrum to an HDF file.

        Parameters
        ----------
        out_file_path : str
            The path to the HDF file to export to.

        Returns
        -------
        None
        """
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_hdf()

    def to_csv(self, out_file_path, write_metadata=True):
        """Export the mass spectrum to a CSV file.
        
        Parameters
        ----------
        out_file_path : str
            The path to the CSV file to export to.
        write_metadata : bool, optional
            Whether to write the metadata to the CSV file. Defaults to True.
        
        """
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_csv(write_metadata=write_metadata)

    def to_pandas(self, out_file_path, write_metadata=True):
        """Export the mass spectrum to a Pandas dataframe with pkl extension.

        Parameters
        ----------
        out_file_path : str
            The path to the CSV file to export to.
        write_metadata : bool, optional
            Whether to write the metadata to the CSV file. Defaults to True.

        """
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(out_file_path, self)
        exportMS.to_pandas(write_metadata=write_metadata)

    def to_dataframe(self,):
        """Return the mass spectrum as a Pandas dataframe.
        
        Returns
        -------
        pandas.DataFrame
            The mass spectrum as a Pandas dataframe.
        """
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.get_pandas_df()

    def to_json(self):
        """Return the mass spectrum as a JSON file."""
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.to_json()

    def parameters_json(self):
        """Return the parameters of the mass spectrum as a JSON string."""
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.parameters_to_json()

    def parameters_toml(self):
        """Return the parameters of the mass spectrum as a TOML string."""
        from corems.mass_spectrum.output.export import HighResMassSpecExport
        exportMS = HighResMassSpecExport(self.filename, self)
        return exportMS.parameters_to_toml()

class MassSpecProfile(MassSpecBase):
    """A mass spectrum class when the entry point is on profile format
    
    Notes
    -----
    Stores the profile data and instrument settings. 
    Iteration over a list of MSPeaks classes stored at the _mspeaks attributes.
    _mspeaks is populated under the hood by calling process_mass_spec method.
    Iteration is null if _mspeaks is empty. Many more attributes and methods inherited from MassSpecBase().

    Parameters
    ----------
    data_dict : dict
        A dictionary containing the profile data.
    d_params : dict{'str': float, int or str}
        contains the instrument settings and processing settings
    auto_process : bool, optional
        Whether to automatically process the mass spectrum. Defaults to True.


    Attributes 
    ----------
    _abundance : ndarray
        The abundance values of the mass spectrum.
    _mz_exp : ndarray
        The m/z values of the mass spectrum.
    _mspeaks : list
        A list of mass peaks.

    Methods 
    ----------
    * process_mass_spec(). Process the mass spectrum.

    see also: MassSpecBase(), MassSpecfromFreq(), MassSpecCentroid()
    """

    def __init__(self, data_dict, d_params, auto_process=True):
        # print(data_dict.keys())
        super().__init__(data_dict.get(Labels.mz), data_dict.get(Labels.abundance), d_params)
       
        if auto_process:
            self.process_mass_spec()

class MassSpecfromFreq(MassSpecBase):
    """ A mass spectrum class when data entry is on frequency domain

    Notes
    -----
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
    d_params : dict{'str': float, int or str}
        contains the instrument settings and processing settings
    auto_process : bool, optional
        Whether to automatically process the mass spectrum. Defaults to True.
    keep_profile : bool, optional
        Whether to keep the profile data. Defaults to True.
  
    Attributes
    ----------
    has_frequency : bool
        Whether the mass spectrum has frequency data.
    _frequency_domain : list(float)
        Frequency domain in Hz
    label : str
        store label (Bruker, Midas Transient, see Labels class ). It across distinct processing points
    _abundance : ndarray
        The abundance values of the mass spectrum.
    _mz_exp : ndarray
        The m/z values of the mass spectrum.
    _mspeaks : list
        A list of mass peaks.
    See Also: all the attributes of MassSpecBase class
     
    Methods
    ----------
    * _set_mz_domain().
        calculates the m_z based on the setting of d_params
    * process_mass_spec().  Process the mass spectrum.
    
    see also: MassSpecBase(), MassSpecProfile(), MassSpecCentroid()
    """

    def __init__(self, frequency_domain, magnitude, d_params, 
                auto_process=True, keep_profile=True):

        super().__init__(None, magnitude, d_params)

        self._frequency_domain = frequency_domain
        self.has_frequency = True
        self._set_mz_domain()
        self._sort_mz_domain()
        
        self.magnetron_frequency = None
        self.magnetron_frequency_sigma = None

        #use this call to automatically process data as the object is created, Setting need to be changed before initiating the class to be in effect
        
        if auto_process:
            self.process_mass_spec(keep_profile=keep_profile)

    def _sort_mz_domain(self):
        """Sort the mass spectrum by m/z values."""

        if self._mz_exp[0] > self._mz_exp[-1]:
            self._mz_exp = self._mz_exp[::-1]
            self._abundance = self._abundance[::-1]
            self._frequency_domain = self._frequency_domain[::-1]

    def _set_mz_domain(self):
        """Set the m/z domain of the mass spectrum based on the settings of d_params."""
        if self.label == Labels.bruker_frequency:
            self._mz_exp = self._f_to_mz_bruker()

        else:

            self._mz_exp = self._f_to_mz()

    @property
    def transient_settings(self): 
        """Return the transient settings of the mass spectrum."""
        return self.parameters.transient

    @transient_settings.setter
    def transient_settings(self, instance_TransientSetting):
     
        self.parameters.transient = instance_TransientSetting  

    def calc_magnetron_freq(self, max_magnetron_freq=50,magnetron_freq_bins=300):
        """Calculates the magnetron frequency of the mass spectrum.

        Parameters
        ----------
        max_magnetron_freq : float, optional
            The maximum magnetron frequency. Defaults to 50.
        magnetron_freq_bins : int, optional
            The number of bins to use for the histogram. Defaults to 300.

        Returns
        -------
        None

        Notes
        -----
        Calculates the magnetron frequency by examining all the picked peaks and the distances between them in the frequency domain.
        A histogram of those values below the threshold 'max_magnetron_freq' with the 'magnetron_freq_bins' number of bins is calculated.
        A gaussian model is fit to this histogram - the center value of this (statistically probably) the magnetron frequency.
        This appears to work well or nOmega datasets, but may not work well for 1x datasets or those with very low magnetron peaks.
        """
        ms_df = DataFrame(self.freq_exp(),columns=['Freq'])
        ms_df['FreqDelta'] = ms_df['Freq'].diff()

        freq_hist = histogram(ms_df[ms_df['FreqDelta']<max_magnetron_freq]['FreqDelta'],bins=magnetron_freq_bins)
    
        mod = GaussianModel()
        pars = mod.guess(freq_hist[0], x=freq_hist[1][:-1])
        out = mod.fit(freq_hist[0], pars, x=freq_hist[1][:-1])
        self.magnetron_frequency = out.best_values['center']
        self.magnetron_frequency_sigma = out.best_values['sigma']
            

class MassSpecCentroid(MassSpecBase):

    """A mass spectrum class when the entry point is on centroid format

    Notes
    -----
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
        contains the instrument settings and processing settings
    auto_process : bool, optional
        Whether to automatically process the mass spectrum. Defaults to True.
        
    Attributes
    ----------
    label : str
        store label (Bruker, Midas Transient, see Labels class)
    _baseline_noise : float
        store baseline noise
    _baseline_noise_std : float
        store baseline noise std
    _abundance : ndarray
        The abundance values of the mass spectrum.
    _mz_exp : ndarray
        The m/z values of the mass spectrum.
    _mspeaks : list
        A list of mass peaks. 

    
    Methods
    ----------
    * process_mass_spec().
        Process the mass spectrum. Overriden from MassSpecBase. Populates the _mspeaks list with MSpeaks class using the centroid data.
    * __simulate_profile__data__().
        Simulate profile data based on Gaussian or Lorentzian peak shape. Needs theoretical resolving power calculation and define peak shape, intended for plotting and inspection purposes only.

    see also: MassSpecBase(), MassSpecfromFreq(), MassSpecProfile()
    """

    def __init__(self, data_dict, d_params, auto_process=True):

        super().__init__([], [], d_params)

        self._set_parameters_objects(d_params)
        
        if self.label == Labels.thermo_centroid:
            self._baseline_noise = d_params.get("baseline_noise")
            self._baseline_noise_std = d_params.get("baseline_noise_std")

        self.is_centroid = True
        self.data_dict = data_dict
        self._mz_exp = data_dict[Labels.mz]
        self._abundance = data_dict[Labels.abundance]

        if auto_process:
            self.process_mass_spec()
            

    def __simulate_profile__data__(self, exp_mz_centroid, magnitude_centroid):
        """Simulate profile data based on Gaussian or Lorentzian peak shape

        Notes
        -----
        Needs theoretical resolving power calculation and define peak shape.
        This is a quick fix to trick a line plot be able to plot as sticks for plotting and inspection purposes only.
        
        Parameters
        ----------
        exp_mz_centroid : list(float)
            list of m/z values
        magnitude_centroid : list(float)
            list of abundance values
            
            
        Returns
        -------
        x : list(float)
            list of m/z values
        y : list(float)
            list of abundance values
        """

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
        """Return the m/z profile of the mass spectrum."""
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
        """Return the abundance profile of the mass spectrum."""
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
        """Return the total ion current of the mass spectrum."""
        return sum(self.abundance)

    def process_mass_spec(self):
        """Process the mass spectrum.
       
        """
        import tqdm
        # overwrite process_mass_spec 
        # mspeak objs are usually added inside the PeaKPicking class 
        # for profile and freq based data
        
        data_dict = self.data_dict
        s2n = True
        ion_charge = self.polarity
        #l_exp_mz_centroid = data_dict.get(Labels.mz)
        #l_intes_centr = data_dict.get(Labels.abundance)
        #l_peak_resolving_power = data_dict.get(Labels.rp)
        l_s2n = list(data_dict.get(Labels.s2n))
        
        if not l_s2n: s2n = False
        
        print("Loading mass spectrum object")
        
        abun = array(data_dict.get(Labels.abundance)).astype(float)
        
        abundance_threshold, factor = self.get_threshold(abun)
        
        for index, mz in enumerate(data_dict.get(Labels.mz)):
            
            # centroid peak does not have start and end peak index pos
            massspec_indexes = (index, index, index)
            
            if s2n:
                
                if abun[index]/factor >= abundance_threshold:

                    self.add_mspeak(
                        ion_charge,
                        mz,
                        abun[index],
                        float(data_dict.get(Labels.rp)[index]),
                        float(l_s2n[index]),
                        massspec_indexes,
                        ms_parent=self
                    )

            else:

                if data_dict.get(Labels.abundance)[index]/factor >= abundance_threshold:

                    self.add_mspeak(
                        ion_charge,
                        mz,
                        abun[index],
                        float(data_dict.get(Labels.rp)[index]),
                        -999,
                        massspec_indexes,
                        ms_parent=self
                    )

        self.mspeaks = self._mspeaks
        self._dynamic_range = self.max_abundance / self.min_abundance
        self._set_nominal_masses_start_final_indexes()
        
        if self.label != Labels.thermo_centroid:
            
            if self.settings.noise_threshold_method == 'log':
                
                raise  Exception("log noise Not tested for centroid data")
                #self._baseline_noise, self._baseline_noise_std = self.run_log_noise_threshold_calc()
            
            else:
                self._baseline_noise, self._baseline_noise_std = self.run_noise_threshold_calc()
        
        del self.data_dict
    
class MassSpecCentroidLowRes(MassSpecCentroid):
    """A mass spectrum class when the entry point is on low resolution centroid format

    Notes
    -----
    Does not store MSPeak Objs, will iterate over mz, abundance pairs instead
    
    Parameters
    ----------
    data_dict : dict {string: numpy array float64 )
        contains keys [m/z, Abundance, Resolving Power, S/N]
    d_params : dict{'str': float, int or str}
        contains the instrument settings and processing settings

    Attributes
    ----------
    _processed_tic : float
        store processed total ion current
    _abundance : ndarray
        The abundance values of the mass spectrum.
    _mz_exp : ndarray
        The m/z values of the mass spectrum.
    """
    
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
        """Return the m/z values of the mass spectrum."""
        return self._mz_exp 

    @property
    def abundance(self):
        """Return the abundance values of the mass spectrum."""
        return self._abundance

    @property
    def processed_tic(self):
        """Return the processed total ion current of the mass spectrum."""
        return sum(self._processed_tic)
    
    @property
    def tic(self):
        """Return the total ion current of the mass spectrum."""
        if self._processed_tic:
            return self._processed_tic
        else:
            return sum(self.abundance)
    
    @property
    def mz_abun_tuples(self):
        """Return the m/z and abundance values of the mass spectrum as a list of tuples."""
        r = lambda x: ( int(round(x[0],0), int(round(x[1],0))) )

        return [r(i) for i in self]
    
    @property
    def mz_abun_dict(self):
        """Return the m/z and abundance values of the mass spectrum as a dictionary."""
        r = lambda x: int(round(x,0))
            
        return { r(i[0]):r(i[1]) for i in self}
