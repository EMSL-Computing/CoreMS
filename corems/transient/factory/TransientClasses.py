from os.path import basename, dirname, normpath

from matplotlib import rcParamsDefault, rcParams
from numpy import linspace

from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq
from corems.transient.calc.TransientCalc import TransientCalculations
import matplotlib.pyplot as plt
from copy import deepcopy
from corems.encapsulation.input.parameter_from_json import (
    load_and_set_parameters_class,
    load_and_set_toml_parameters_class,
)


__author__ = "Yuri E. Corilo"
__date__ = "Jun 19, 2019"


class Transient(TransientCalculations):
    """The Transient object contains the transient data and the parameters used to process it

    Parameters
    ----------
    data : numpy.ndarray
        Array with the transient data
    d_params : dict
        Dictionary with the parameters to be set

    Attributes
    ----------
    calibration_terms : tuple
        Tuple with the calibration terms (A, B, C)
    bandwidth : float
        The bandwidth of the transient (Hz)
    number_data_points : int
        The number of data points of the transient
    polarity : int
        The polarity of the transient
    transient_time : float
        The time domain length of the transient
    d_params : dict
        Dictionary with the parameters to be set
    frequency_domain : numpy.ndarray
        Array with the frequency domain
    magnitude : numpy.ndarray
        Array with the magnitude
    _full_filename_path : str
        The full path of the transient file
    _exc_high_freq : float
        The high frequency of the excitation (Hz)
    _exc_low_freq : float
        The low frequency of the excitation (Hz)
    _parameters : corems.transient.parameters.TransientParameters
        The transient parameters
    _transient_data : numpy.ndarray
        Array with the transient data


    Methods
    -------
    * get_frequency_domain(plot_result=True).
        Get the frequency domain and magnitude from the transient data
    * get_mass_spectrum(auto_process=True, plot_result=True, keep_profile=True).
        Get the mass spectrum from the transient data
    * set_processing_parameter(apodization_method, number_of_truncations, number_of_zero_fills).
        Set the processing parameters
    * scale_plot_size(factor=1.5).
        Scale the plot size by a factor
    * plot_transient(ax=None, c='k').
        Plot the transient data
    * plot_zerofilled_transient(ax=None, c='k').
        Plot the transient data with zero fill
    * plot_apodized_transient(ax=None, c='k').
        Plot the transient data with apodization
    * plot_frequency_domain(ax=None, c='k').
        Plot the frequency domain and magnitude
    * set_parameter_from_toml(parameters_path).
        Set the processing parameters from a toml file
    * set_parameter_from_json(parameters_path).
        Set the processing parameters from a json file



    """

    def __init__(self, data, d_params):
        self._transient_data = data

        self.d_params = d_params

        self.frequency_domain = None

        self.magnitude = None

        self.__set__parameters__objects(d_params)

        self.__set__transient__time()

    def __set__parameters__objects(self, d_params):
        """Set the parameters objects from the dictionary d_params

        Parameters
        ----------
        d_params : dict
            Dictionary with the parameters to be set

        """

        self._full_filename_path = d_params.get("filename_path")

        self.calibration_terms = (
            d_params.get("Aterm"),
            d_params.get("Bterm"),
            d_params.get("Cterm"),
        )

        self._exc_high_freq = d_params.get("exc_high_freq")

        self._exc_low_freq = d_params.get("exc_low_freq")

        self._qpd_enabled = d_params.get("qpd_enabled")  # Quadrupolar detection enabled

        self._mw_low = d_params.get("mw_low")  # low mass for detection

        self._mw_high = d_params.get("mw_high")  # high mass for detection

        self.bandwidth = d_params.get("bandwidth")

        self.number_data_points = d_params.get("number_data_points")

        self.polarity = int(d_params.get("polarity"))

        self.location = 220

        self._parameters = deepcopy(MSParameters.transient)

    def scale_plot_size(self, factor=1.5):
        """Scale the plot size by a factor

        Parameters
        ----------
        factor : float, optional
            The factor to scale the plot size, by default 1.5
        """

        default_dpi = rcParamsDefault["figure.dpi"]
        rcParams["figure.dpi"] = default_dpi * factor

    def __set__transient__time(self):
        """Set the transient time variable with the calculated length."""
        self.transient_time = self.cal_transient_time()

    def set_processing_parameter(
        self,
        apodization_method: str,
        number_of_truncations: int,
        number_of_zero_fills: int,
    ):
        """Set the processing parameters

        Parameters
        ----------
        apodization_method : str
            Apodization method to be used
        number_of_truncations : int
            Number of truncations to be used
        number_of_zero_fills : int
            Number of zero fills to be used
        """

        self.parameters.apodization_method = apodization_method

        self.parameters.number_of_truncations = number_of_truncations

        self.parameters.number_of_zero_fills = number_of_zero_fills

    @property
    def parameters(self):
        """The transient parameters"""
        return self._parameters

    @parameters.setter
    def parameters(self, instance_TransientParameters):
        self._parameters = instance_TransientParameters

    def set_parameter_from_toml(self, parameters_path):
        """Set the processing parameters from a toml file"""
        self._parameters = load_and_set_toml_parameters_class(
            "Transient", self._parameters, parameters_path=parameters_path
        )

    def set_parameter_from_json(self, parameters_path):
        """Set the processing parameters from a json file"""
        self._parameters = load_and_set_parameters_class(
            "Transient", self._parameters, parameters_path=parameters_path
        )

    def get_frequency_domain(self, plot_result=True):
        """Get the frequency domain and magnitude from the transient data

        Parameters
        ----------
        plot_result : bool, optional
            Plot the frequency domain and magnitude, by default True

        Returns
        -------
        frequency_domain : numpy.ndarray
            Array with the frequency domain
        magnitude : numpy.ndarray
            Array with the magnitude
        """

        if self.parameters.number_of_truncations > 0:
            new_time_domain = self.truncation(self._transient_data)

        else:
            new_time_domain = self._transient_data

        if self.parameters.apodization_method is not None:
            new_time_domain = self.apodization(new_time_domain)

        if plot_result:
            self._plot_transient(self._transient_data)

            self._plot_transient(new_time_domain)

        time_domain_y_zero_filled = self.zero_fill(new_time_domain)

        self.transient_time = self.transient_time * (
            self.parameters.number_of_zero_fills + 1
        )

        if plot_result:
            self._plot_transient(time_domain_y_zero_filled)

        return self.perform_magniture_mode_ft(time_domain_y_zero_filled)
        # return frequency_domain, magnitude

    def get_mass_spectrum(
        self,
        auto_process: bool = True,
        plot_result: bool = True,
        keep_profile: bool = True,
    ) -> MassSpecfromFreq:
        """Get the mass spectrum from the transient data

        Parameters
        ----------
        auto_process : bool, optional
            Process the transient data, by default True
        plot_result : bool, optional
            Plot the frequency domain and magnitude, by default True
        keep_profile : bool, optional
            Keep the profile data, by default True

        Returns
        -------
        MassSpecfromFreq
            Mass spectrum object
        """

        frequency_domain, magnitude = self.get_frequency_domain(plot_result=plot_result)

        if plot_result:
            self._plot_frequency_domain(frequency_domain, magnitude)

        self.d_params["filename"] = self.filename
        self.d_params["dir_location"] = self.dir_location

        return MassSpecfromFreq(
            frequency_domain,
            magnitude,
            self.d_params,
            auto_process=auto_process,
            keep_profile=keep_profile,
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

    def _plot_frequency_domain(self, frequency_domain, magnitude):  # pragma: no cover
        """Plot the frequency domain and magnitude

        Parameters
        ----------
        frequency_domain : numpy.ndarray
            Array with the frequency domain
        magnitude : numpy.ndarray
            Array with the magnitude
        """

        self.location += 1
        plt.subplot(self.location)
        plt.plot(frequency_domain, magnitude, color="green")
        plt.xlabel("Hz")
        plt.ylabel("Magnitude")
        # reset grid location index to 0
        self.location = 220
        # plt.show()

    def _plot_transient(self, transient_data):  # pragma: no cover
        """Plot the transient data

        Parameters
        ----------
        transient_data : numpy.ndarray
            Array with the transient data

        """

        self.location += 1
        # print( self.location)
        time_axis = linspace(0, self.transient_time, num=len(transient_data))
        plt.subplot(self.location)
        plt.plot(time_axis, transient_data, color="green")
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        # plt.show()

    def plot_transient(self, ax=None, c="k"):  # pragma: no cover
        """Plot the transient data

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes object, by default None
        c : str, optional
            Color, by default 'k'

        Returns
        -------
        matplotlib.axes
            Matplotlib axes object

        """

        # self.location +=1
        # print( self.location)
        if ax is None:
            ax = plt.gca()
        time_axis = linspace(0, self.transient_time, num=len(self._transient_data))
        # plt.subplot(self.location)
        ax.plot(time_axis, self._transient_data, color=c)
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        # plt.show()
        return ax

    def plot_zerofilled_transient(self, ax=None, c="k"):  # pragma: no cover
        """Plot the transient data with zero fill

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes object, by default None
        c : str, optional
            Color, by default 'k'

        Returns
        -------
        matplotlib.axes
            Matplotlib axes object

        """
        if ax is None:
            ax = plt.gca()
        new_time_domain = self.apodization(self._transient_data)
        time_domain_y_zero_filled = self.zero_fill(new_time_domain)
        self.transient_time = self.transient_time * (
            self.parameters.number_of_zero_fills + 1
        )
        time_axis = linspace(0, self.transient_time, num=len(time_domain_y_zero_filled))
        # plt.subplot(self.location)
        ax.plot(time_axis, time_domain_y_zero_filled, color=c)
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        # plt.show()
        return ax

    def plot_apodized_transient(self, ax=None, c="k"):  # pragma: no cover
        """Plot the transient data with apodization

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes object, by default None
        c : str, optional
            Color, by default 'k'

        Returns
        -------
        matplotlib.axes
            Matplotlib axes object

        """
        # self.location +=1
        # print( self.location)
        if ax is None:
            ax = plt.gca()
        new_time_domain = self.apodization(self._transient_data)
        time_axis = linspace(0, self.transient_time, num=len(new_time_domain))
        # plt.subplot(self.location)
        ax.plot(time_axis, new_time_domain, color=c)
        plt.xlabel("Time (s)")
        plt.ylabel("Magnitude")
        # plt.show()
        return ax

    def plot_frequency_domain(self, ax=None, c="k"):  # pragma: no cover
        """Plot the frequency domain and magnitude

        Parameters
        ----------
        ax : matplotlib.axes, optional
            Matplotlib axes object, by default None
        c : str, optional
            Color, by default 'k'

        Returns
        -------
        matplotlib.axes
            Matplotlib axes object

        """
        # self.location +=1
        # plt.subplot(self.location)
        if ax is None:
            ax = plt.gca()
        frequency_domain, magnitude = self.get_frequency_domain(plot_result=False)
        ax.plot(frequency_domain / 1000, magnitude, color=c)
        plt.xlabel("KHz")
        plt.ylabel("Magnitude")
        # plt.show()
        return ax
