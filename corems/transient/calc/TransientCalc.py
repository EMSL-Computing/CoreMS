import gc
import warnings

from numpy import (
    arange,
    blackman,
    ceil,
    fft,
    hamming,
    hanning,
    kaiser,
    linspace,
    log2,
    pi,
    power,
    sin,
    sqrt,
    where,
    zeros,
)

__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"


class TransientCalculations(object):
    """Transient Calculations

    Parameters
    ----------
    parameters : corems.transient.parameters.TransientParameters
        The transient parameters
    bandwidth : float
        The bandwidth of the transient (Hz)
    number_data_points : int
        The number of data points of the transient
    exc_low_freq : float
        The low frequency of the excitation (Hz)
    exc_high_freq : float
        The high frequency of the excitation (Hz)

    Attributes
    ----------
    parameters : corems.transient.parameters.TransientParameters
        The transient parameters
    bandwidth : float
        The bandwidth of the transient (Hz)
    number_data_points : int
        The number of data points of the transient
    exc_low_freq : float
        The low frequency of the excitation (Hz)
    exc_high_freq : float
        The high frequency of the excitation (Hz)

    Methods
    -------
    * cal_transient_time().
        Calculate the time domain length of the transient
    * zero_fill(transient).
        Zero fill the transient
    * truncation(transient).
        Truncate the transient
    * apodization(transient).
        Apodization of the transient
    * calculate_frequency_domain(number_data_points).
        Calculate the frequency domain (axis) of the transient
    * cut_freq_domain(freqdomain_X, freqdomain_Y).
        Cut the frequency domain of the transient
    * phase_and_absorption_mode_ft().
        [Not Functional] Produce a phased absorption mode FT spectrum
    * magnitude_mode_ft(transient).
        Perform magnitude mode FT of the transient
    * correct_dc_offset().
        [Not Yet Implemented] Correct the DC offset of the transient

    """

    def cal_transient_time(self):
        """Calculate the time domain length of the transient

        Returns
        -------
        float
            The time domain length of the transient (s)
        """
        return (1 / self.bandwidth) * ((self.number_data_points) / 2)

    def zero_fill(self, transient):
        """Zero fill the transient

        Parameters
        ----------
        transient : numpy.ndarray
            The transient data points

        Returns
        -------
        numpy.ndarray
            The transient data points zerofilled

        Notes
        -----
        The number of zero fills is defined by the transient parameter number_of_zero_fills.
        The function first calculate the next power of two of the transient length and zero fills to that length, to take advantage of FFT algorithm.
            If the parameter next_power_of_two is set to False, the function will zero fill to the length of the original transient times the number of zero fills

        """
        if self.parameters.next_power_of_two:
            exponent = int(
                ceil(log2(len(transient) * (self.parameters.number_of_zero_fills + 1)))
            )
            zeros_filled_transient = zeros(2**exponent)
        else:
            zeros_filled_transient = zeros(
                len(transient) * (self.parameters.number_of_zero_fills + 1)
            )

        zeros_filled_transient[0 : len(transient)] = transient

        del transient

        gc.collect()

        return zeros_filled_transient

    def truncation(self, transient):
        """Truncate the transient

        Parameters
        ----------
        transient : numpy.ndarray
            The transient data points

        Returns
        -------
        numpy.ndarray
            The truncated transient data points

        Notes
        -----
        The number of truncations is defined by the transient parameter number_of_truncations
        """

        data_count = len(transient)

        for _ in range(self.parameters.number_of_truncations):
            data_count = int(data_count / 2)

        time_domain_truncated = transient[0:data_count]

        del transient

        gc.collect()

        return time_domain_truncated

    def apodization(self, transient):
        """Apodization of the transient

        Parameters
        ----------
        transient : numpy.ndarray
            The transient data points

        Returns
        -------
        numpy.ndarray
            The apodized transient data points

        Notes
        -----
        The apodization method is defined by the transient parameter apodization_method.
        The following apodization methods are available:
            Hamming,
            Hanning,
            Blackman,
            Full-Sine,
            Half-Sine,
            Kaiser,
            Half-Kaiser,
            Rectangular/None

        For Kaiser and Half-Kaiser, an additional parameter 'beta' is required, set by the transient parameter kaiser_beta.

        """

        apodi_method = self.parameters.apodization_method
        beta = self.parameters.kaiser_beta

        length = len(transient)

        if apodi_method == "Hamming":
            H_function = hamming(length)
        elif apodi_method == "Hanning":
            H_function = hanning(length)
        elif apodi_method == "Blackman":
            H_function = blackman(length)
        elif apodi_method == "Full-Sine":
            H_function = sin(linspace(0, pi, num=length))
        elif apodi_method == "Half-Sine":
            H_function = sin(linspace((pi / 2), 0, num=length))
        elif apodi_method == "Kaiser":
            H_function = kaiser(length, beta)
        elif apodi_method == "Half-Kaiser":
            H_function = kaiser(length * 2, beta)[length:]
        elif apodi_method == 'Rectangular' or apodi_method is None:
            H_function = 1

        S_x = transient * H_function

        del transient
        gc.collect()

        return S_x

    def calculate_frequency_domain(self, number_data_points):
        """Calculate the frequency domain (axis) of the transient

        Parameters
        ----------
        number_data_points : int
            The number of data points of the transient

        Returns
        -------
        numpy.ndarray
            The frequency domain of the transient (Hz)


        """

        qntpoints = arange(0, (number_data_points))

        factor_distancy = (self.bandwidth) / (number_data_points)

        frequency_domain = qntpoints * factor_distancy

        del qntpoints
        del factor_distancy
        gc.collect()

        return frequency_domain

    def cut_freq_domain(self, freqdomain_X, freqdomain_Y):
        """Cut the frequency domain of the transient

        Parameters
        ----------
        freqdomain_X : numpy.ndarray
            The frequency domain of the transient (Hz)
        freqdomain_Y : numpy.ndarray
            The frequency domain of the transient (Hz)

        Returns
        -------
        numpy.ndarray
            The frequency domain of the transient (Hz)
        numpy.ndarray
            The frequency domain of the transient (Hz)


        """
        # If the mw_low and mw_high are set, the frequency domain is cut to the mw range
        # this accounts for the detection settings, not the excitation settings.
        # TODO: Implement this - right now the f to mz function is in the ms class, not the transient class, so it doesnt work.
        # if (self._mw_low != 0) & (self._mw_high != 0):
        #    high_freq = self._f_to_mz(self._mw_high)
        #    low_freq = self._f_to_mz(self._mw_low)
        #
        #    final =  where(freqdomain_X < high_freq)[-1][-1]
        #      start =  where(freqdomain_X > low_freq)[0][0]
        # else:
        if self._qpd_enabled == 1:
            low_freq = self._exc_low_freq * 2
            high_freq = self._exc_high_freq * 2
        else:
            low_freq = self._exc_low_freq
            high_freq = self._exc_high_freq

        if self._exc_low_freq > self._exc_high_freq:
            # TODO: This needs to be tested
            # I'm not sure that this is relevant anyway - the excitation pulse is ramped in frequency but the detection is simulatenous
            warnings.warn("This is not tested. Please check the results.")
            final = where(freqdomain_X > low_freq)[0][0]
            start = where(freqdomain_X > high_freq)[0][0]

        else:
            final = where(freqdomain_X < high_freq)[-1][-1]
            start = where(freqdomain_X > low_freq)[0][0]

        return freqdomain_X[start:final], freqdomain_Y[start:final]
        # del freqdomain_X, freqdomain_Y
        # gc.collect()

    def phase_and_absorption_mode_ft(self):
        """[Not Functional] Produce a phased absorption mode FT spectrum"""
        # anyone wants to play with this part please make yourself comfortable. I will:
        pass

    def perform_magniture_mode_ft(self, transient):
        """Perform magnitude mode FT of the transient

        Parameters
        ----------
        transient : numpy.ndarray
            The transient data points

        Returns
        -------
        numpy.ndarray
            The frequency domain of the transient (Hz)
        numpy.ndarray
            The magnitude of the transient (a.u.)


        """

        A = fft.rfft(transient)

        # A = fft.fft(transient)
        # A = A[0:int(len(A)/2)]

        factor = int(self.parameters.number_of_zero_fills - 1)
        if self.parameters.number_of_zero_fills:
            if self.parameters.number_of_zero_fills == 1:
                factor = 1 / 2

            else:
                factor = int(1 / self.parameters.number_of_zero_fills + 1)

            Max_index = int(len(A) / factor)

        else:
            Max_index = int(len(A))

        A = A[0:Max_index]

        datapoints = len(A)

        freqdomain_X = self.calculate_frequency_domain(datapoints)

        magnitude_Y = sqrt((power(A.real, 2)) + (power(A.imag, 2)))

        freqdomain_X_cut, magnitude_Y_cut = self.cut_freq_domain(
            freqdomain_X, magnitude_Y
        )

        del transient
        # del freqdomain_X
        # del magnitude_Y
        gc.collect()

        return freqdomain_X_cut, magnitude_Y_cut

    def correct_dc_offset(self):
        """[Not Yet Implemented] Correct the DC offset of the transient

        A simple baseline correction to compensate for a DC offset in the recorded transient.
        Not implemented.

        """
        pass
