import numpy as np

from scipy.signal import savgol_filter
from scipy.signal.windows import boxcar
from scipy import interpolate
from matplotlib import pyplot as plt
from numpy import abs
from numpy import array, polyfit, asarray


def peak_detector(tic, max_tic):  # TODO remove max_tic argument?
    """
    Find peaks by detecting minima in the first derivative of the data
    Used in LC/GC data processing

    Parameters
    ----------
    tic : array
        array of data points to find the peaks
    max_tic : float
        maximum value of the data points

    Returns
    -------
    tuple
        tuple of indexes of the start, apex and final points of the peak

    """
    dy = derivate(tic)

    indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]

    for index in indexes:
        start_index = find_minima(index, tic, right=False)
        final_index = find_minima(index, tic)

        yield (start_index, index, final_index)


def find_nearest_scan(data, nodes):
    """
    Find nearest data point in a list of nodes (derivated data)
    in LC/GC this is 'scan', in MS this is 'm/z' data point

    Parameters
    ----------
    data : float
        data point to find the nearest node
    nodes : array
        array of nodes to search for the nearest node

    Returns
    -------
    float
        nearest node to the data point
    """

    array_data = asarray(nodes)

    scan_index = (abs(array_data - data)).argmin()

    return nodes[scan_index]


def check_corrected_abundance(
    closest_left,
    closest_right,
    apex_index,
    signal,
    max_signal,
    signal_threshold,
    abun_norm,
):
    """
    Check the corrected abundance of the peak

    Parameters
    ----------
    closest_left : int
        index of the closest left node
    closest_right : int
        index of the closest right node
    apex_index : int
        index of the apex node
    signal : array
        array of data points to find the peaks
    max_signal : float
        maximum value of the data points
    signal_threshold : float
        threshold for the signal
    abun_norm : float
        abundance normalization factor

    Returns
    -------
    float
        corrected abundance of the peak


    """
    x = [closest_left, closest_right]
    y = [signal[closest_left], signal[closest_right]]

    pol = polyfit(x, y, 1)  # TODO replace with faster method in this file

    corrected_peak_height = signal[apex_index] - pol(apex_index)

    if (corrected_peak_height / max_signal) * abun_norm > signal_threshold:
        return corrected_peak_height
    else:
        return False


def peak_picking_first_derivative(
    domain,
    signal,
    max_height,
    max_prominence,
    max_signal,
    min_peak_datapoints,
    peak_derivative_threshold,
    signal_threshold=0.1,
    correct_baseline=True,
    plot_res=False,
    abun_norm=100,
    check_abundance=False,
    apex_indexes=[],
):
    """
    Find peaks by detecting minima in the first derivative of the data
    Used in LC/GC and MS data processing
    Optional baseline correction, then peak apex detection via 1st derivative.
    For each apex the peak datapoints surrounding the apex are determined.
    Some basic thresholding is applied (signal, number of datapoints, etc).

    Parameters
    ----------
    domain : array
        array of data points to find the peaks
    signal : array
        array of data points to find the peaks
    max_height : float
        maximum height of the peak
    max_prominence : float
        maximum prominence of the peak
    max_signal : float
        maximum signal of the peak
    min_peak_datapoints : int
        minimum number of data points in the peak
    peak_derivative_threshold : float
        threshold for the peak derivative
    signal_threshold : float
        threshold for the signal
    correct_baseline : bool
        flag to correct the baseline
    plot_res : bool
        flag to plot the results
    abun_norm : float
        abundance normalization factor
    check_abundance : bool
        flag to check the abundance


    Returns
    -------
    tuple
        tuple of indexes of the start, apex and final points of the peak


    """
    if correct_baseline:
        signal = signal - baseline_detector(signal, domain, max_height, max_prominence)

    domain = np.array(domain)
    signal = np.array(signal)

    dy = derivate(signal)
    if len(apex_indexes) == 0:
        # Find apexes
        apex_indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]
    else:
        apex_indexes = np.array(apex_indexes)

    if apex_indexes.size and apex_indexes is not None:
        apex_indexes = apex_indexes[
            signal[apex_indexes] / max_signal >= signal_threshold
        ]

    signal = signal / max(signal)
    start_peak = []
    end_peak = []

    pos_dy_threshold = peak_derivative_threshold  # max(dy) * peak_derivative_threshold
    neg_dy_threshold = -peak_derivative_threshold  # min(dy) * peak_derivative_threshold
    len_dy = len(dy)
    # take apex_index and move left to find start
    for index in apex_indexes:
        # catch for starting position

        if index == 0:
            index_start = index
        else:
            index_start = index - 1

        # catch for ending position
        if (index + 1) >= dy.shape[0]:
            index_end = index - 1
        else:
            index_end = index + 1

        # while dy[index_start-1] > 0 and index_start != 0:
        while dy[index_start - 1] > pos_dy_threshold and index_start > 0:
            index_start = index_start - 1
        start_peak.append(index_start)

        # while dy[index_end] < 0 and index_end != (len(dy) - 1):
        while dy[index_end] < neg_dy_threshold and index_end != (len_dy - 1):
            index_end = index_end + 1
        end_peak.append(index_end)

    start_peak = array(start_peak)
    end_peak = array(end_peak)

    for apex_index in apex_indexes:
        # index_gt_apex = np.where(end_peak >= apex_index)[0]
        # index_lt_apex = np.where(start_peak <= apex_index)[0]
        index_gt_apex = np.arange(np.searchsorted(end_peak, apex_index), len(end_peak))
        index_lt_apex = np.arange(
            0, np.searchsorted(start_peak, apex_index, side="right")
        )

        if not index_gt_apex.size == 0 and not index_lt_apex.size == 0:
            closest_right = find_nearest_scan(apex_index, end_peak[index_gt_apex])
            closest_left = find_nearest_scan(apex_index, start_peak[index_lt_apex])
            if check_abundance:
                corrected_peak_height = check_corrected_abundance(
                    closest_left,
                    closest_right,
                    apex_index,
                    signal,
                    max_signal,
                    signal_threshold,
                    abun_norm,
                )
            else:
                corrected_peak_height = signal[apex_index]

            if (closest_right - closest_left) >= min_peak_datapoints:
                if plot_res:
                    plt.plot(
                        domain[closest_left : closest_right + 1],
                        dy[closest_left : closest_right + 1],
                        c="red",
                    )
                    plt.plot(
                        domain[closest_left : closest_right + 1],
                        signal[closest_left : closest_right + 1],
                        c="black",
                    )
                    plt.title(str((corrected_peak_height / max_signal) * 100))
                    plt.show()

                yield (closest_left, apex_index, closest_right)


def find_minima(index, tic, right=True):
    """
    Find the index of the local minima in the given time-of-flight (TOF) intensity array.

    Parameters:
    -----------
    index: int
        The starting index to search for the minima.
    tic: list
        TIC data points
    right : bool, optional
        Determines the direction of the search. If True, search to the right of the index. If False, search to the left of the index. Default is True.

    Returns:
    --------
    int
        The index of the local minima in the TIC  array.
    """

    j = index
    # apex_abundance = tic[index]
    tic_len = len(tic)

    if right:
        minima = tic[j] >= tic[j + 1]
    else:
        minima = tic[j] >= tic[j - 1]

    while minima:
        if j == 1 or j == tic_len - 2:
            break

        if right:
            j += 1

            minima = tic[j] >= tic[j + 1]

        else:
            j -= 1
            minima = tic[j] >= tic[j - 1]

    if right:
        return j
    else:
        return j


def derivate(data_array):
    """
    Calculate derivative of the data points.
    Replaces nan with infinity

    Parameters
    ----------
    data_array : array
        array of data points

    Returns
    -------
    array
        array of the derivative of the data points
    """
    data_array = np.array(data_array)

    dy = data_array[1:] - data_array[:-1]

    # replaces nan for infinity
    indices_nan = np.where(np.isnan(data_array))[0]

    if indices_nan.size:
        data_array[indices_nan] = np.inf
        dy[np.where(np.isnan(dy))[0]] = np.inf

    return dy


def minima_detector(tic, max_tic, peak_height_max_percent, peak_max_prominence_percent):
    """
    Minima detector for the TIC data points.

    Parameters
    ----------
    tic : array
        array of data points to find the peaks
    max_tic : float
        maximum value of the data points
    peak_height_max_percent : float
        maximum height of the peak
    peak_max_prominence_percent : float
        maximum prominence of the peak

    Returns
    -------
    generator
        generator of the indexes of the minima in the TIC array

    """
    peak_height_diff = lambda hi, li: ((tic[hi] - tic[li]) / max_tic) * 100

    for start_index, index, final_index in peak_detector(tic, max_tic):
        # abundance max threshold
        if (tic[index] / max_tic) * 100 < peak_height_max_percent:
            # calculates prominence and filter
            if (
                peak_height_diff(index, start_index)
                and peak_height_diff(index, final_index) < peak_max_prominence_percent
            ):
                yield from (start_index, final_index)


def baseline_detector(
    tic, rt, peak_height_max_percent, peak_max_prominence_percent, do_interpolation=True
):
    """
    Baseline detector for the TIC data points.
    For LC/GC data processing

    Parameters
    ----------
    tic : array
        array of data points to find the peaks
    rt : array
        array of retention time data points
    peak_height_max_percent : float
        maximum height of the peak
    peak_max_prominence_percent : float
        maximum prominence of the peak
    do_interpolation : bool, optional
        flag to interpolate the data points. Default is True

    Returns
    -------
    array
        array of the baseline corrected data points

    """
    rt = np.array(rt)

    max_tic = max(tic)

    indexes = sorted(
        list(
            set(
                i
                for i in minima_detector(
                    tic, max_tic, peak_height_max_percent, peak_max_prominence_percent
                )
            )
        )
    )

    y = -tic

    x1 = rt[indexes]

    y1 = y[indexes]

    if len(x1) <= 5:
        return tic

    if not do_interpolation:
        y0 = np.zeros(tic.shape)
        y0[indexes] = y[indexes]

        return y0

    else:
        f1 = interpolate.interp1d(x1, y1, kind="quadratic", fill_value="extrapolate")

        ynew1 = f1(list(rt))

        # from matplotlib import pyplot as plt
        # if self.deconv_rt_list and  self.deconv_mz == 51:

        #   plt.plot(rt, tic-(-1* ynew1), color='green')

        # plt.plot(rt, -1* ynew1, c='black')

        # s = self.smooth(s, 10, 'blackman')

        # plt.plot(self.retention_time, -s)

        # plt.show()

        return -1 * ynew1


def peak_detector_generator(
    tic, stds, method, rt, max_height, min_height, max_prominence, min_datapoints
):
    """
    Peak detector generator for the TIC data points.

    Parameters
    ----------
    tic : array
        array of data points to find the peaks
    stds : float
        standard deviation
    method : str
        method to detect the peaks
        Available methods: 'manual_relative_abundance', 'auto_relative_abundance', 'second_derivative'
    rt : array
        array of retention time data points
    max_height : float
        maximum height of the peak
    min_height : float
        minimum height of the peak
    max_prominence : float
        maximum prominence of the peak
    min_datapoints : int
        minimum number of data points in the peak

    Returns
    -------
    generator
        generator of the indexes of the peaks in the TIC array

    """
    max_tic = max(tic)

    if method == "manual_relative_abundance":
        tic = tic - baseline_detector(tic, rt, max_height, max_prominence)

        norm_tic = (tic / max_tic) * 100

        remove_indexes = np.where(norm_tic < min_height)[0]

        # if self.deconv_rt_list and  self.deconv_mz == 51:
        #    plt.plot(self.deconv_rt_list, tic, label=self.deconv_mz)

    elif method == "auto_relative_abundance":
        tic = tic - baseline_detector(tic, rt, max_height, max_prominence)

        baseline = baseline_detector(tic, rt, max_height, max_prominence)

        peak_detect_threshold = np.nanmean(baseline) + (stds * np.std(baseline))

        remove_indexes = np.where(tic < peak_detect_threshold)[0]

    elif method == "second_derivative":
        remove_indexes = second_derivative_threshold(
            tic, stds, rt, max_height, max_prominence
        )

    else:
        NotImplemented(method)

    peak_height_diff = lambda hi, li: ((tic[hi] - tic[li]) / max_tic) * 100

    dy = derivate(tic)

    include_indexes = np.where((np.hstack((dy, 0)) < 0) & (np.hstack((0, dy)) > 0))[0]

    final_indexes = sorted(set(include_indexes) - set(remove_indexes))

    # from matplotlib import pyplot as plt

    # plt.plot(self.retention_time, tic, color='black')
    # plt.scatter(self.retention_time[remove_indexes], tic[remove_indexes], color='red')
    # plt.scatter(self.retention_time[include_indexes], tic[include_indexes], color='blue')
    # plt.scatter(self.retention_time[final_indexes], tic[final_indexes], color='blue')

    # plt.show()

    for index in final_indexes:
        start_index = find_minima(index, tic, right=False)
        final_index = find_minima(index, tic)

        if final_index - start_index > min_datapoints:
            # if min( peak_height_diff(index,start_index), peak_height_diff(index,final_index) )> self.chromatogram_settings.peak_min_prominence_percent :

            yield (start_index, index, final_index)


def smooth_signal(x, window_len, window, pol_order, implemented_smooth_method):
    """
    Smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x: array
        the input signal
    window_len: int
        the dimension of the smoothing window; should be an odd integer
    window: str
        the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    pol_order: int
        the order of the polynomial to fit the data
    implemented_smooth_method: list
        list of implemented smoothing methods

    Returns
    -------
    y: array
        the smoothed signal

    Notes:
    -----
    See also: numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.savgol_filter

    """
    x = np.array(x)

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input array needs to be bigger than window size")

    # if window_len < 3:
    #    return x

    if not window in implemented_smooth_method:
        raise ValueError(
            "Window method should be 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-1:-window_len:-1]]

    if window == "savgol":
        return savgol_filter(x, window_len, pol_order)

    elif window == "boxcar":  # moving average
        w = boxcar(window_len)

        y = np.convolve(w, s, mode="valid")

    elif window == "flat":  # moving average
        w = np.ones(window_len, "d")

        y = np.convolve(w / w.sum(), s, mode="valid")

    else:
        w = eval(window + "(window_len)")

        y = np.convolve(w / w.sum(), s, mode="valid")

    return y[int(window_len / 2 - 1) : int(-window_len / 2)]


def second_derivative_threshold(
    tic, stds, rt, peak_height_max_percent, peak_max_prominence_percent
):
    """
    Second derivative threshold for the TIC data points.
    For LC/GC data processing

    Parameters
    ----------
    tic : array
        array of data points to find the peaks
    stds : float
        standard deviation
    rt : array
        array of retention time data points
    peak_height_max_percent : float
        maximum height of the peak

    Returns
    -------
    array
        array of the indexes of the data points to remove

    """

    dy = derivate(tic)

    dydy = derivate(dy)
    dydy = np.hstack((dydy, 0))
    dydy = np.hstack((0, dydy))

    baseline = baseline_detector(
        dydy,
        rt,
        peak_height_max_percent,
        peak_max_prominence_percent,
        do_interpolation=False,
    )

    threshold_median = np.median(baseline) - (stds * np.std(baseline))

    remove_indexes = np.where(dydy > threshold_median)[0]

    return remove_indexes
