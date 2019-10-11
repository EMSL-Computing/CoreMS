from os.path import basename, dirname, normpath

from matplotlib import rcParamsDefault, rcParams
from numpy import linspace

from corems.encapsulation.settings.input.ProcessingSetting import TransientSetting
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq
from corems.transient.calc.TransientCalc import TransientCalculations
import matplotlib.pyplot as plt


__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"

'''
fig = plt.figure()

fig.patch.set_facecolor(None)

fig.patch.set_alpha(0)
'''

class Transient(TransientCalculations):
    def __init__(self, data, d_params):

        self._transient_data = data

        self.d_params = d_params

        self.frequency_domain = None

        self.magnitude = None

        self.__set__parameters__objects(d_params)

        self.__set__transient__time()

    def __set__parameters__objects(self, d_params):

        self._full_filename_path = d_params.get("filename_path")

        self.calibration_terms = (
            d_params.get("Aterm"),
            d_params.get("Bterm"),
            d_params.get("Cterm"),
        )

        self._exc_high_freq = d_params.get("exc_high_freq")

        self._exc_low_freq = d_params.get("exc_low_freq")

        self.bandwidth = d_params.get("bandwidth")

        self.number_data_points = d_params.get("number_data_points")

        self.polarity = int(d_params.get("polarity"))

        self.location = 220

    def scale_plot_size(self, factor=1.5):

        default_dpi = rcParamsDefault["figure.dpi"]
        rcParams["figure.dpi"] = default_dpi * factor

    def __set__transient__time(self):
        self.transient_time = self.cal_transient_time()

    def set_processing_parameter(
        self, apodization_method, number_of_truncations, number_of_zero_fills
    ):

        TransientSetting.apodization_method = apodization_method

        TransientSetting.number_of_truncations = number_of_truncations

        TransientSetting.number_of_zero_fills = number_of_zero_fills

    def get_frequency_domain(self, plot_result=True):

        if TransientSetting.number_of_truncations > 0:

            new_time_domain = self.truncation(self._transient_data)

        else:

            new_time_domain = self._transient_data

        if TransientSetting.apodization_method is not None:

            new_time_domain = self.apodization(new_time_domain)

        if plot_result:

            self._plot_transient(self._transient_data)

            self._plot_transient(new_time_domain)

        time_domain_y_zero_filled = self.zero_fill(new_time_domain)

        self.transient_time = self.transient_time * (
            TransientSetting.number_of_zero_fills + 1
        )

        if plot_result:

            self._plot_transient(time_domain_y_zero_filled)

        return self.perform_magniture_mode_ft(time_domain_y_zero_filled)
        # return frequency_domain, magnitude

    def get_mass_spectrum(self, auto_process=True, plot_result=True):

        #plt.figure(figsize=(13, 8))

        frequency_domain, magnitude = self.get_frequency_domain(plot_result=plot_result)

        if plot_result:

            self._plot_frequency_domain(frequency_domain, magnitude)

        self.d_params["filename"] = self.filename
        self.d_params["dir_location"] = self.dir_location
        
        return MassSpecfromFreq(
            frequency_domain, magnitude, self.d_params, auto_process=auto_process
        )

    @property
    def filename(self):

        # return dirname(self._full_filename_path)
        return basename(normpath(self._full_filename_path))

    @property
    def dir_location(self):

        return dirname(self._full_filename_path).strip(
            basename(normpath(self._full_filename_path))
        )

    @property
    def A_therm(self):

        return self.calibration_terms[0]

    @property
    def B_therm(self):

        return self.calibration_terms[1]

    @property
    def C_therm(self):

        return self.calibration_terms[2]

    def _plot_frequency_domain(self, frequency_domain, magnitude):

        self.location += 1
        plt.subplot(self.location)
        plt.plot(frequency_domain, magnitude, color="green")
        plt.xlabel("Hz")
        plt.ylabel("Magnitude")
        # reset grid location index to 0
        self.location = 220
        plt.show()

    def _plot_transient(self, transient_data):

        self.location += 1
        # print( self.location)
        time_axis = linspace(0, self.transient_time, num=len(transient_data))
        plt.subplot(self.location)
        plt.plot(time_axis, transient_data, color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        # plt.show()

    def plot_transient(self):

        # self.location +=1
        # print( self.location)
        time_axis = linspace(0, self.transient_time, num=len(self._transient_data))
        # plt.subplot(self.location)
        plt.plot(time_axis, self._transient_data, color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        plt.show()

    def plot_zerofilled_transient(self):

        new_time_domain = self.apodization(self._transient_data)
        time_domain_y_zero_filled = self.zero_fill(new_time_domain)
        self.transient_time = self.transient_time * (
            TransientSetting.number_of_zero_fills + 1
        )
        time_axis = linspace(0, self.transient_time, num=len(time_domain_y_zero_filled))
        # plt.subplot(self.location)
        plt.plot(time_axis, time_domain_y_zero_filled, color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        plt.show()

    def plot_apodized_transient(self):

        # self.location +=1
        # print( self.location)
        new_time_domain = self.apodization(self._transient_data)
        time_axis = linspace(0, self.transient_time, num=len(new_time_domain))
        # plt.subplot(self.location)
        plt.plot(time_axis, new_time_domain, color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        plt.show()

    def plot_frequency_domain(self):

        # self.location +=1
        # plt.subplot(self.location)
        frequency_domain, magnitude = self.get_frequency_domain(plot_result=False)
        plt.plot(frequency_domain / 1000, magnitude, color="green")
        plt.xlabel("KHz")
        plt.ylabel("Magnitude")
        plt.show()
