import time

from matplotlib import rcParamsDefault, rcParams
from numpy import array, flip
from enviroms.emsl.yec.encapsulation.constant.Constants import Labels
from enviroms.emsl.yec.encapsulation.settings.input.ProcessingSetting import MassSpectrumSetting
from enviroms.emsl.yec.mass_spectrum.calc.MassSpectrumCalc import MassSpecCalc
from enviroms.emsl.yec.mass_spectrum.factory.MSPeakClasses import MSPeak
import matplotlib.pyplot as plt

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

fig = plt.figure()

fig.patch.set_facecolor(None)

fig.patch.set_alpha(0)


def overrides(interface_class):
    def overrider(method):
        assert method.__name__ in dir(interface_class)
        return method
    return overrider


class MassSpecBase(MassSpecCalc):
    """
    Mass Spectrum class object with common features and functions
    """
    def __init__(self, mz_exp, abundance, d_params, **kwargs):

        self._abundance = array(abundance)
        self._mz_exp = array(mz_exp)
        
        #objects created after process_mass_spec() function
        self.mspeaks = []
        self._dict_nominal_masses_indexes  = {}

        self._set_parameters_objects(d_params)

        # frequency is set inside MassSpecfromFreq Class
        self._frequency_domain = None

        # for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value)

    def __len__(self):
        
        return len(self.mspeaks)
        
    def __getitem__(self, position):
        
        return self.mspeaks[position]

    def add_mspeak (self, ion_charge, mz_exp,
                    abundance,
                    resolving_power,
                    signal_to_noise,
                    massspec_index,
                    exp_freq=None,
                ):

        self.mspeaks.append(
            MSPeak(
                ion_charge,
                mz_exp,
                abundance,
                resolving_power,
                signal_to_noise,
                massspec_index,
                exp_freq=exp_freq,
            )
        )

    def _set_parameters_objects(self, d_params):

        self._calibration_terms = (
            d_params.get("Aterm"),
            d_params.get("Bterm"),
            d_params.get("Cterm"),
        )

        self.polarity = int(d_params.get("polarity"))

        self.scan_number = d_params.get("scan_number")

        self.rt = d_params.get("rt")

        self._filename = d_params.get("filename")

        self._dir_location = d_params.get("dir_location")

        self.location = 220

    def process_mass_spec(self):

        self.cal_noise_treshould()
        self.find_peaks()
        self._set_nominal_masses_start_final_indexes()
        

    def cal_noise_treshould(self, auto=True):

        self._baselise_noise, self._baselise_noise_std = self.run_noise_threshould_calc(
            auto
        )

    def scale_plot_size(self, factor=1.5):

        default_dpi = rcParamsDefault["figure.dpi"]
        rcParams["figure.dpi"] = default_dpi * factor

    @property
    def frequency_domain(self):
        return self._frequency_domain

    @property
    def mz_exp(self): return self._mz_exp
    
    @property
    def abundance(self): return self._abundance
    
    @property
    def mz_exp_centroide(self):
        self.check_mspeaks()
        return array([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def max_mz_exp(self):
        return max([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def min_mz_exp(self):
        return min([mspeak.mz_exp for mspeak in self.mspeaks])

    @property
    def abundance_centroid(self):
        self.check_mspeaks()
        return array([mspeak.abundance for mspeak in self.mspeaks])

    @property
    def max_abundance(self):
        return max([mspeak.abundance for mspeak in self.mspeaks])
    
    @property
    def most_abundant_mspeak(self):
        
        return max(self.mspeaks, key=lambda m: m.abundance)
    
    @property
    def min_abundance(self):
        return min([mspeak.abundance for mspeak in self.mspeaks])
    
    @property
    def baselise_noise(self):
        return self._baselise_noise

    @property
    def baselise_noise_std(self):
        return self._baselise_noise_std

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
        return self._filename

    @property
    def dir_location(self):
        return self._dir_location

    def check_mspeaks(self):
        if self.mspeaks:
            pass
        else:
            raise Exception(
                "mspeaks dictionary is empty, please run process_mass_spec() first"
            )

    def filter_by_s2n(self, s2n):

        self.check_mspeaks()
        return [mspeak for mspeak in self.mspeaks if mspeak.signal_to_noise > s2n]

    def get_mz_and_abundance_peaks_tuples(self):

        self.check_mspeaks()
        return [(mspeak.mz, mspeak.abundance) for mspeak in self.mspeaks]

    def find_peaks(self):
        """needs to clear previous results from peak_picking"""
        del self.mspeaks 
        self.mspeaks = list()
        """then do peak picking"""
        self.do_peak_picking()
        print("A total of %i peaks were found" % len(self.mspeaks))

    def change_kendrick_base_all_mspeaks(self, kendrick_dict_base):
        """kendrick_dict_base = {"C": 1, "H": 2} or {{"C": 1, "H": 1, "O":1} etc """
        for mspeak in self.mspeaks:
            mspeak.change_kendrick_base(kendrick_dict_base)

    def get_nominal_mz_frist_last_indexes(self, nominal_mass):
        
        if self._dict_nominal_masses_indexes:
            
            if nominal_mass in self._dict_nominal_masses_indexes.keys():
                
                    
                return (self._dict_nominal_masses_indexes.get(nominal_mass)[0], self._dict_nominal_masses_indexes.get(nominal_mass)[1]+1)
            
            else:
                #import warnings
                #uncomment warn to distribution
                #warnings.warn("Nominal mass not found in _dict_nominal_masses_indexes, returning (0, 0) for nominal mass %i"%nominal_mass)
                return (0,0)
        else:
            raise Exception("run process_mass_spec() function before trying to access the data")

    def get_nominal_mass_indexes(self, nominal_mass, overlay=0.1):
        min_mz_to_look = nominal_mass - overlay
        max_mz_to_look = nominal_mass+1+overlay
        indexes = [i for i in range(len(self.mspeaks)) if min_mz_to_look <= self.mspeaks[i].mz_exp <= max_mz_to_look]
        return indexes
    
    def get_masses_sum_for_nominal_mass(self):
        
        dict_nominal_masses_count ={}
        
        all_nominal_masses = list(set([i.nominal_mz_exp for i in self.mspeaks]))
        
        for nominal_mass in all_nominal_masses:
            
            dict_nominal_masses_count[nominal_mass] = len(self.get_nominal_mass_indexes(nominal_mass))

        return dict_nominal_masses_count

    def _set_nominal_masses_start_final_indexes(self):
        '''return ms peaks objs indexes(start and end) on the mass spectrum for all nominal masses'''
        dict_nominal_masses_indexes ={}
        
        all_nominal_masses = list(set([i.nominal_mz_exp for i in self.mspeaks]))
        
        for nominal_mass in all_nominal_masses:
            
            indexes = self.get_nominal_mass_indexes(nominal_mass)
                
            dict_nominal_masses_indexes[nominal_mass] = (indexes[0],indexes[-1]) 
          

        self._dict_nominal_masses_indexes = dict_nominal_masses_indexes

    def plot_mz_domain_profile_and_noise_threshold(self):

        if self.baselise_noise and self.baselise_noise:
            x = (self.mz_exp.min(), self.mz_exp.max())
            y = (self.baselise_noise, self.baselise_noise)

            stds = MassSpectrumSetting.noise_threshold_stds
            threshold = self.baselise_noise + (stds * self.baselise_noise_std)
            plt.plot(self.mz_exp, self.abundance, color="green")
            plt.plot(x, (threshold, threshold), color="yellow")
            plt.plot(x, y, color="red")
            plt.xlabel("m/z")
            plt.ylabel("abundance")
            plt.show()

        else:

            raise Exception("Calculate noise threshold first")

    def plot_mz_domain_profile(self):

        plt.plot(self.mz_exp, self.abundance, color="green")
        plt.xlabel("m/z")
        plt.ylabel("abundance")
        plt.show()


class MassSpecProfile(MassSpecBase):

    """
    class docs
    """

    def __init__(self, dataframe, d_params, auto_process=True):
        """
        method docs
        """
        self.label = d_params.get("label")

        mz_exp = dataframe["m/z"].values
        abundance = dataframe["Abundance"].values
        super().__init__(mz_exp, abundance, d_params)
        if auto_process:
            self.process_mass_spec()

    @overrides(MassSpecBase)
    def process_mass_spec(self):
        self.cal_noise_treshould()
        self.find_peaks()
        self._set_nominal_masses_start_final_indexes()


class MassSpecfromFreq(MassSpecBase):
    def __init__(self, frequency_domain, magnitude, d_params, auto_process=True):
        """
        method docs
        """
        super().__init__(None, flip(magnitude), d_params)

        self.label = d_params.get("label")
        self._frequency_domain = frequency_domain
        self._set_mz_domain()
        """ use this call to automatically process data as the object is created, Setting need to be changed before initiating the class to be in effect"""
        if auto_process:
            self.process_mass_spec()

    def _set_mz_domain(self):

        if self.label == Labels.bruker_frequency:

            self._mz_exp = flip(self._f_to_mz_bruker())

        else:

            self._mz_exp = flip(self._f_to_mz())


class MassSpecCentroid(MassSpecBase):

    """
    classdocs
    """

    def __init__(self, dataframe, d_params, auto_process=True):
        """

        """
        """needs to simulate peak shape and pass as mz_exp and magnitude."""
        exp_mz_centroid = dataframe["m/z"].values
        magnitude_centroid = dataframe["Abundance"].values
        # mz_exp, magnitude = self.__simulate_profile__data__(
        #    exp_mz_centroid, magnitude_centroid)

        # print( mz_exp)

        self.label = d_params.get("label")
        self.dataframe = dataframe
        super().__init__(exp_mz_centroid, magnitude_centroid, d_params)

        self._set_parameters_objects(d_params)
        if self.label == Labels.thermo_centroid:
            self._baselise_noise = d_params.get("baselise_noise")
            self._baselise_noise_std = d_params.get("baselise_noise_std")

        if auto_process:
            self.process_mass_spec(dataframe)
            del self.dataframe

    def __simulate_profile__data__(self, exp_mz_centroid, magnitude_centroid):
        """needs theoretical resolving power calculation and define peak shape
        this is a quick fix to be able to plot as lines
        peakshape = #Gaussian"""

        x, y = [], []
        for i in range(len(exp_mz_centroid)):
            x.append(exp_mz_centroid[i] - 0.0000001)
            x.append(exp_mz_centroid[i])
            x.append(exp_mz_centroid[i] + 0.0000001)
            y.append(0)
            y.append(magnitude_centroid[i])
            y.append(0)
        return x, y

    @overrides(MassSpecBase)
    def process_mass_spec(self, dataframe):

        # it is wasy too slow, it needs to be changed to other functional structure

        ion_charge = self.polarity
        l_exp_mz_centroid = dataframe["m/z"]
        l_intes_centr = dataframe["Abundance"]
        l_peak_resolving_power = dataframe["Resolving Power"]
        l_s_n = dataframe["S/N"]

        for index in range(dataframe["m/z"].size):
            self.add_mspeak(
                ion_charge,
                l_exp_mz_centroid[index],
                l_intes_centr[index],
                l_peak_resolving_power[index],
                l_s_n[index],
                index,
            )

        self._set_nominal_masses_start_final_indexes()    
