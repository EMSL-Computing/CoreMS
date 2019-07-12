import time

from matplotlib import rcParamsDefault, rcParams
from numpy import array, flip
from enviroms.emsl.yec.encapsulation.constant.Constants import Labels
from enviroms.emsl.yec.encapsulation.settings.ProcessingSetting import MassSpectrumSetting
from enviroms.emsl.yec.mass_spectrum.calc.MassSpectrumCalc import MassSpecCalc
from enviroms.emsl.yec.mass_spectrum.factory.MSPeakClasses import MSPeak
import matplotlib.pyplot as plt


__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


fig = plt.figure()

fig.patch.set_facecolor(None)

fig.patch.set_alpha(0)


class MassSpecBase(MassSpecCalc):
    '''
    Mass Spectrum class object with common features and functions
    '''

    def __init__(self, exp_mz, abundance, d_params, **kwargs):
        '''
        Constructor
        '''
        self._abundance = array(abundance)
        self._exp_mz = array(exp_mz)
        self.mspeaks = list()

        self._set_parameters_objects(d_params)

        # frequency is set inside MassSpecfromFreq Class
        self._frequency_domain = None

        # for (key, value) in kwargs.items():
        #    print(key, value)
        #    if hasattr(self, key):
        #        setattr(self, key, value)
        #        print(key, value)

    def _set_parameters_objects(self, d_params):

        self._calibration_terms = (d_params.get(
            "Aterm"), d_params.get("Bterm"), d_params.get("Cterm"))

        self.polarity = int(d_params.get("polarity"))

        self.scan_number = d_params.get("scan_number")

        self.rt = d_params.get("rt")

        self._filename = d_params.get("filename")

        self._dir_location = d_params.get("dir_location")

        self.location = 220

    def cal_noise_treshould(self, auto=True):

        self._baselise_noise, self._baselise_noise_std = self.run_noise_threshould_calc(
            auto)

    def scale_plot_size(self, factor=1.5):

        default_dpi = rcParamsDefault['figure.dpi']
        rcParams['figure.dpi'] = default_dpi*factor

    @property
    def frequency_domain(self): return self._frequency_domain

    @property
    def exp_mz(self): return self._exp_mz

    @property
    def exp_mz_centroide(self):
        self.check_mspeaks()
        return array([mspeak.exp_mz for mspeak in self.mspeaks])

    @property
    def abundance(self): return self._abundance

    @property
    def abundance_centroid(self):
        self.check_mspeaks()
        return array([mspeak.abundance for mspeak in self.mspeaks])

    @property
    def baselise_noise(self): return self._baselise_noise

    @property
    def baselise_noise_std(self): return self._baselise_noise_std

    @property
    def Aterm(self): return self._calibration_terms[0]

    @property
    def Bterm(self): return self._calibration_terms[1]

    @property
    def Cterm(self): return self._calibration_terms[2]

    @property
    def filename(self): return self._filename

    @property
    def dir_location(self): return self._dir_location

    def add_mspeak(self, ion_charge, exp_mz, abundance, resolving_power, signal_to_noise, massspec_index, exp_freq=None):

        self.mspeaks.append(MSPeak(ion_charge, exp_mz, abundance, resolving_power,
                                   signal_to_noise, massspec_index, exp_freq=exp_freq))

    def check_mspeaks(self):
        if self.mspeaks:
            pass
        else:
            raise Exception(
                "mspeaks dictionary is empty, please run find_peaks() first")

    def filter_by_s2n(self, s2n):

        self.check_mspeaks()
        return [mspeak for mspeak in self.mspeaks if mspeak.signal_to_noise > s2n]

    def get_peaks_as_list(self):

        self.check_mspeaks()
        return list(self.mspeaks)

    def get_mz_and_abundance_peaks_tuples(self):

        self.check_mspeaks()
        return [(mspeak.mz, mspeak.abundance) for mspeak in self.mspeaks]

    def find_peaks(self):
        '''needs to clear previous results from peak_picking'''
        self.mspeaks = list()
        '''then do peak picking'''
        self.do_peak_picking()
        print("A total of %i peaks were found" % len(self.mspeaks))

    def change_kendrick_base_all_mspeaks(self, kendrick_dict_base):
        '''kendrick_dict_base = {"C": 1, "H": 2} or {{"C": 1, "H": 1, "O":1} etc '''
        for mspeak in self.mspeaks:
            mspeak.change_kendrick_base(kendrick_dict_base)

    def plot_mz_domain_profile_and_noise_threshold(self):

        if self.baselise_noise and self.baselise_noise:
            x = (self.exp_mz.min(), self.exp_mz.max())
            y = (self.baselise_noise, self.baselise_noise)

            stds = MassSpectrumSetting.noise_threshold_stds
            threshold = (self.baselise_noise + (stds*self.baselise_noise_std))
            plt.plot(self.exp_mz, self.abundance, color='green')
            plt.plot(x, (threshold, threshold), color='yellow')
            plt.plot(x, y, color='red')
            plt.xlabel("m/z")
            plt.ylabel("abundance")
            plt.show()

        else:

            raise Exception("Calculate noise threshold first")

    def plot_mz_domain_profile(self):

        plt.plot(self.exp_mz, self.abundance, color='green')
        plt.xlabel("m/z")
        plt.ylabel("abundance")
        plt.show()


class MassSpecProfile(MassSpecBase):

    '''
    class docs
    '''

    def __init__(self, dataframe, d_params, **kwargs):
        '''
        Constructor
        '''
        self.label = d_params.get('label')

        exp_mz = dataframe['m/z'].values
        abundance = dataframe['Abundance'].values
        super().__init__(exp_mz, abundance, d_params)

        if self.label == Labels.thermo_profile:
            self.stn = dataframe["S/N"].values
            self.resolving_power = dataframe['Resolving Power'].values
            self._baselise_noise = d_params.get('baselise_noise')
            self._baselise_noise_std = d_params.get('baselise_noise_std')

    def process_mass_spec(self):

        if self.label != Labels.thermo_profile:

            self.cal_noise_treshould()

        self.find_peaks()


class MassSpecfromFreq(MassSpecBase):

    def __init__(self, frequency_domain, magnitude, d_params, **kwargs):
        '''
        Constructor
        '''
        super().__init__(None, flip(magnitude), d_params)

        self.label = d_params.get('label')
        self._frequency_domain = frequency_domain
        self._set_mz_domain()
        ''' use this call to automatically process data as the object is created, Setting need to be changed before initiating the class to be in effect'''
        # self.process_mass_spec()

    def _set_mz_domain(self):

        if self.label == Labels.bruker_frequency:

            self._exp_mz = flip(self._f_to_mz_bruker())

        else:

            self._exp_mz = flip(self._f_to_mz())

    def process_mass_spec(self):

        time0 = time.time()
        self.cal_noise_treshould()
        print(round(time.time() - time0, 2), "seconds to calc thresould")
        time1 = time.time()
        self.find_peaks()
        print(round(time.time() - time1, 2), "seconds to find the peaks")


class MassSpecCentroid(MassSpecBase):

    '''
    classdocs
    '''

    def __init__(self, dataframe, d_params, **kwargs):
        '''
        Constructor
        '''
        '''needs to simulate peak shape and pass as exp_mz and magnitude.'''
        exp_mz_centroid = dataframe['m/z'].values
        magnitude_centroid = dataframe['Abundance'].values
        # exp_mz, magnitude = self.__simulate_profile__data__(
        #    exp_mz_centroid, magnitude_centroid)

        # print( exp_mz)
        self.label = Labels.centroid
        super().__init__(exp_mz_centroid, magnitude_centroid, d_params)

        self._set_parameters_objects(d_params)
        self.__process__from__centroid(dataframe)

    def __simulate_profile__data__(self, exp_mz_centroid, magnitude_centroid):
        '''needs theoretical resolving power calculation and define peak shape
        this is a quick fix to be able to plot as lines
        peakshape = #Gaussian'''

        x, y = [], []
        for i in range(len(exp_mz_centroid)):
            x.append(exp_mz_centroid[i] - 0.0000001)
            x.append(exp_mz_centroid[i])
            x.append(exp_mz_centroid[i] + 0.0000001)
            y.append(0)
            y.append(magnitude_centroid[i])
            y.append(0)
        return x, y

    def __process__from__centroid(self, dataframe):

        # this need to change after mass spec deconvolution
        ion_charge = self.polarity
        for index in range(dataframe['m/z'].size):
            exp_mz_centroid = dataframe['m/z'][index]
            intes_centr = dataframe['Abundance'][index]
            s_n = dataframe['S/N'][index]
            peak_resolving_power = dataframe['Resolving Power'][index]
            self.add_mspeak(ion_charge, exp_mz_centroid,
                            intes_centr, peak_resolving_power, s_n, index)
