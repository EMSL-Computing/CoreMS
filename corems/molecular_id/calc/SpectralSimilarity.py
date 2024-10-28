__author__ = "Yuri E. Corilo"
__date__ = "Jun 09, 2021"

from numpy.fft import rfft
from scipy.stats import pearsonr, spearmanr, kendalltau
from numpy import (
    power,
    dot,
    absolute,
    sqrt,
)
from numpy import sum as np_sum
from numpy.linalg import norm
from pandas import DataFrame
import numpy as np

methods_name = {
    # "entropy_distance": "Entropy Distance",
    # "weighted_entropy_distance": "Dynamic weighted entropy Distance",
    "chebyshev_distance": "Chebyshev Distance",
    "squared_euclidean_distance": "Squared Euclidean Distance",
    "fidelity_similarity": "Fidelity Similarity",
    "matusita_distance": "Matusita Distance",
    "squared_chord_distance": "Squared-chord Distance",
    # "bhattacharya_1_distance": "Bhattacharya 1 Distance",
    # "bhattacharya_2_distance": "Bhattacharya 2 Distance",
    "harmonic_mean_similarity": "Harmonic mean Distance",
    "Pearson_chi_squared_distance": "Pearson Chi Squared Distance",
    "Neyman_chi_squared_distance": "Neyman Chi Squared Distance",
    "probabilistic_symmetric_chi_squared_distance": "Probabilistic symmetric X2 Distance",
    "topsoe_distance": "Topsoe Distance",
    "chernoff_distance": "Chernoff Distance",
    "ruzicka_distance": "Ruzicka Distance",
    "roberts_distance": "Roberts Distance",
    # "intersection_distance": "Intersection Distance",
    "motyka_distance": "Motyka Distance",
    "canberra_distance": "Canberra Distance",
    "canberra_metric": "Canberra Metric",
    "kulczynski_1_distance": "Kulczynski 1 Distance",
    # "baroni_urbani_buser_distance": "Baroni-Urbani-Buser Distance",
    # "penrose_size_distance": "Penrose size Distance",
    # "mean_character_distance": "Mean character Distance",
    "lorentzian_distance": "Lorentzian Distance",
    # "penrose_shape_distance": "Penrose shape Distance",
    "clark_distance": "Clark Distance",
    "hellinger_distance": "Hellinger Distance",
    "whittaker_index_of_association_distance": "Whittaker index of association Distance",
    # "similarity_index_distance": "Similarity Index Distance",
    # "improved_similarity_distance": "Improved Similarity",
    # "absolute_value_distance": "Absolute Value Distance",
    "spectral_contrast_angle_distance": "Spectral Contrast Angle",
    "wave_hedges_distance": "Wave Hedges Distance",
    "dice_similarity": "Dice Similarity",
    "inner_product_distance": "Inner Product Distance",
    "divergence_distance": "Divergence Distance",
    "jensen_difference_distance": "Jensen Differences Distance",
    "kumar_johnson_distance": "Kumar Johnson Distance",
    "avg_l_distance": "Avg (L1, L8) Distance",
    "vicis_wave_hadges_distance": "Vicis Wave Hadges Distance",
    "vicis_symmetric_chi_squared_1_distance": "Vicis-Symmetric X2 1 Distance",
    "vicis_symmetric_chi_squared_2_distance": "Vicis-Symmetric X2 2 Distance",
    "vicis_symmetric_chi_squared_3_distance": "Vicis-Symmetric X2 3 Distance",
    "max_symmetric_chi_squared_distance": "Max Symmetric Chi Squared Distance",
    "min_symmetric_chi_squared_distance": "Min Symmetric Chi Squared Distance",
    # "ms_for_id_v1": "MSforID Distance version 1",
    # "ms_for_id": "MSforID Distance",
    "additive_sym_chi_sq": "Additive Symmetric Chi Squared",
    "bhattacharya_distance": "Battacharya Distance",
    "generalized_ochiai_index": "Generalized Ochiai Index",
    "gower_distance": "Gower Distance",
    "impr_sqrt_cosine_sim": "Improved Square Root Cosine Similarity",
    "intersection_sim": "Intersection Similarity",
    "j_divergence": "J Divergence",
    "jensen_shannon_index": "Jensen Shannon Index",
    "k_divergence": "K Divergence",
    "VW6": "VW6",
    "VW5": "VW5",
    "VW4": "VW4",
    "VW3": "VW3",
    "VW2": "VW2",
    "VW1": "VW1",
    "taneja_divergence": "Taneja Divergence",
    "symmetric_chi_squared_distance": "Symmetric Chi Squared Distance",
    "squared_chi_squared_distance": "Squared Chi Squared Distance",
    "square_root_cosine_correlation": "Square Root Cosine Correlation",
    "sorensen_distance": "Sorensen Distance",
    "Minokowski_3": "Minokowski 3 Distance",
    "Minokowski_4": "Minokowski 4 Distance",
    "kumarjohnson_divergence": "Kumar Johnson Divergence",
    "kumarhassebrook_similarity": "Kumar Hassebrook Similarity",
    "kullbackleibler_divergence": "Kullback Leibler Divergence",
    "soergel_distance": "Soergel Distance",
}

methods_scale = {
    "entropy": [0, np.log(4)],
    "weighted_entropy": [0, np.log(4)],
    "absolute_value": [0, 2],
    "avg_l": [0, 1.5],
    "bhattacharya_1": [0, np.arccos(0) ** 2],
    "bhattacharya_2": [0, np.inf],
    "canberra": [0, np.inf],
    "clark": [0, np.inf],
    "divergence": [0, np.inf],
    "euclidean": [0, np.sqrt(2)],
    "hellinger": [0, np.inf],
    "improved_similarity": [0, np.inf],
    "lorentzian": [0, np.inf],
    "manhattan": [0, 2],
    "matusita": [0, np.sqrt(2)],
    "mean_character": [0, 2],
    "motyka": [-0.5, 0],
    "ms_for_id": [-np.inf, 0],
    "ms_for_id_v1": [0, np.inf],
    "pearson_correlation": [-1, 1],
    "penrose_shape": [0, np.sqrt(2)],
    "penrose_size": [0, np.inf],
    "probabilistic_symmetric_chi_squared": [0, 1],
    "similarity_index": [0, np.inf],
    "squared_chord": [0, 2],
    "squared_euclidean": [0, 2],
    "symmetric_chi_squared": [0, 0.5 * np.sqrt(2)],
    "topsoe": [0, np.sqrt(2)],
    "vicis_symmetric_chi_squared_3": [0, 2],
    "wave_hedges": [0, np.inf],
    "whittaker_index_of_association": [0, np.inf],
}


class SpectralSimilarity:
    """Class containing methods for calculating spectral similarity between two mass spectra.

    Parameters
    ----------
    ms_mz_abun_dict : dict
        Dictionary of mass to abundance values for the experimental mass spectrum.
    ref_obj : dict
        Dictionary of mass to abundance values for the reference mass spectrum.
    norm_func : function
        Function to normalize the abundance values.

    Attributes
    ----------
    normalize_func : function
        Function to normalize the abundance values.
    ms_mz_abun_dict : dict
        Dictionary of mass to abundance values for the experimental mass spectrum.
    ref_obj : dict
        Dictionary of mass to abundance values for the reference mass spectrum.
    exp_abun : list
        List of abundance values for the experimental mass spectrum.
    exp_mz : list
        List of mass values for the experimental mass spectrum.
    ref_mz : list
        List of mass values for the reference mass spectrum.
    ref_abun : list
        List of abundance values for the reference mass spectrum.
    ref_mz_abun_dict : dict
        Dictionary of mass to abundance values for the reference mass spectrum.
    df : DataFrame
        DataFrame containing the experimental and reference mass spectrum data.
    zero_filled_u_l : tuple
        Tuple containing the experimental and reference mass spectrum data after zero filling and normalization.
    common_mz_values : list
        List of common mass values between the experimental and reference mass spectra.
    n_x_y : int
        Number of common mass values between the experimental and reference mass spectra.

    Methods
    -------
    * nan_fill(df, fill_with=0).
        Fill missing mass values with a given value.
    * normalize(x, y, norm_func=sum).
        Normalize the abundance values.
    * weighted_cosine_correlation(a=0.5, b=1.3, nanfill=1e-10).
        Calculate the weighted cosine correlation between the experimental and reference mass spectra.
    * cosine_correlation().
        Calculate the cosine correlation between the experimental and reference mass spectra.
    * stein_scott().
        Calculate the Stein-Scott similarity between the experimental and reference mass spectra.
    * pearson_correlation().
        Calculate the Pearson correlation between the experimental and reference mass spectra.
    * spearman_correlation().
        Calculate the Spearman correlation between the experimental and reference mass spectra.


    """

    def __init__(self, ms_mz_abun_dict, ref_obj, norm_func=sum):
        self.normalize_func = norm_func
        self.ms_mz_abun_dict = ms_mz_abun_dict
        self.ref_obj = ref_obj

        self.exp_abun = list(self.ms_mz_abun_dict.values())
        self.exp_mz = list(self.ms_mz_abun_dict.keys())

        self.ref_mz = self.ref_obj.get("mz")
        self.ref_abun = self.ref_obj.get("abundance")

        self.ref_mz_abun_dict = dict(zip(self.ref_mz, self.ref_abun))

        # parse to dataframe, easier to zerofill and tranpose
        self.df = DataFrame([self.ms_mz_abun_dict, self.ref_mz_abun_dict])

        # fill missing mz with abundance 0
        x, y = self.nan_fill(self.df, fill_with=1e-10)

        self.zero_filled_u_l = self.normalize(x, y, norm_func=self.normalize_func)

        # filter out the mass values that have zero intensities in self.exp_abun
        exp_mz_filtered = set([k for k in self.exp_mz if self.ms_mz_abun_dict[k] != 0])

        # filter out the mass values that have zero intensities in self.ref_mz
        ref_mz_filtered = set([k for k in self.ref_mz if self.ref_mz_abun_dict[k] != 0])

        # find the intersection/common mass values of both ref and exp, and sort them
        self.common_mz_values = sorted(
            list(exp_mz_filtered.intersection(ref_mz_filtered))
        )

        # find the number of common mass values (after filtering 0s)
        self.n_x_y = len(self.common_mz_values)
        # print(self.n_x_y)

    def nan_fill(self, df, fill_with=0):
        """Fill missing mass values with a given value.

        Parameters
        ----------
        df : DataFrame
            DataFrame containing the experimental and reference mass spectrum data.
        fill_with : float
            Value to fill missing mass values with.

        Returns
        -------
        x : list
            List of abundance values for the experimental mass spectrum.
        y : list
            List of abundance values for the reference mass spectrum."""
        df.fillna(fill_with, inplace=True)

        return df.T[0].values, df.T[1].values

    def normalize(self, x, y, norm_func=sum):
        """Normalize the abundance values.

        Parameters
        ----------
        x : list
            List of abundance values for the experimental mass spectrum.
        y : list
            List of abundance values for the reference mass spectrum.
        norm_func : function
            Function to normalize the abundance values.
            Default is sum

        Returns
        -------
        u_l : tuple
            Tuple containing the experimental and reference mass spectrum data after zero filling and normalization.
        """
        if norm_func is not None:
            u_l = (x / norm_func(x), y / norm_func(y))
            return u_l
        else:
            return (x, y)

    def weighted_cosine_correlation(self, a=0.5, b=1.3, nanfill=1e-10):
        """Calculate the weighted cosine correlation between the experimental and reference mass spectra.

        Parameters
        ----------
        a : float
            Weighting factor for the abundance values.
            Default is 0.5
        b : float
            Weighting factor for the mass values.
            Default is 1.3
        nanfill : float
            Value to fill missing mass values with.
            Default is 1e-10

        Returns
        -------
        correlation : float
            Weighted cosine correlation between the experimental and reference mass spectra.
        """
        # create dict['mz'] = abundance, for experimental data
        # ms_mz_abun_dict = mass_spec.mz_abun_dict
        # weight exp data

        xc = power(self.exp_abun, a) * power(self.exp_mz, b)

        # track back to individual mz
        weighted_exp_dict = dict(zip(self.ms_mz_abun_dict.keys(), xc))

        # weight ref data
        yc = power(self.ref_obj.get("abundance"), a) * power(self.ref_obj.get("mz"), b)

        ref_mz_abun_dict = dict(zip(self.ref_obj.get("mz"), yc))

        # parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([weighted_exp_dict, ref_mz_abun_dict])

        # fill missing mz with weight {abun**a}{m/z**b} to 0
        x, y = self.nan_fill(df, fill_with=nanfill)

        # correlation = (1 - cosine(x, y))

        correlation = dot(x, y) / (norm(x) * norm(y))

        return correlation

    def cosine_correlation(self):
        """Calculate the cosine correlation between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            Cosine correlation between the experimental and reference mass spectra.

        """
        # calculate cosine correlation,
        x = self.zero_filled_u_l[0]
        y = self.zero_filled_u_l[1]

        # correlation = (1 - cosine(x, y))

        correlation = dot(x, y) / (norm(x) * norm(y))

        return correlation

    def stein_scott(self):
        """Calculate the Stein-Scott similarity between the experimental and reference mass spectra.

        Returns
        -------
        s_ss_x_y : float
            Stein-Scott similarity between the experimental and reference mass spectra.
        s_ss_x_y_nist : float
            Stein-Scott similarity between the experimental and reference mass spectra.
        """
        # TODO check this code
        if self.n_x_y == 0:
            return 0, 0

        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in self.exp_abun)

        s_r_x_y = 0

        a, b = 1, 0

        for i in range(1, self.n_x_y):
            current_value = self.common_mz_values[i]
            previous_value = self.common_mz_values[i - 1]

            y_i = self.ref_mz_abun_dict[current_value]
            y_i_minus1 = self.ref_mz_abun_dict[previous_value]

            lc_current = power(y_i, a) * power(current_value, b)
            lc_previous = power(y_i_minus1, a) * power(previous_value, b)

            x_i = self.ms_mz_abun_dict[current_value]
            x_i_minus1 = self.ms_mz_abun_dict[previous_value]

            uc_current = power(x_i, a) * power(current_value, b)
            uc_previous = power(x_i_minus1, a) * power(previous_value, b)

            T1 = lc_current / lc_previous

            T2 = uc_previous / uc_current

            temp_computation = T1 * T2

            n = 0
            if temp_computation <= 1:
                n = 1
            else:
                n = -1

            s_r_x_y = s_r_x_y + power(temp_computation, n)

        # finish the calculation of S_R(X,Y)

        s_r_x_y = s_r_x_y / self.n_x_y
        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation(a=0.5, b=3, nanfill=0)

        s_ss_x_y = ((n_x * s_wc_x_y) + (self.n_x_y * s_r_x_y)) / (n_x + self.n_x_y)

        s_wc_x_y_nist = self.weighted_cosine_correlation(a=0.5, b=1.3, nanfill=0)

        s_ss_x_y_nist = ((n_x * s_wc_x_y_nist) + (self.n_x_y * s_r_x_y)) / (
            n_x + self.n_x_y
        )
        # final step

        return s_ss_x_y, s_ss_x_y_nist

    def pearson_correlation(
        self,
    ):
        """Calculate the Pearson correlation between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            Pearson correlation between the experimental and reference mass spectra.
        """
        correlation = pearsonr(self.zero_filled_u_l[0], self.zero_filled_u_l[1])

        return correlation[0]

    def spearman_correlation(self):
        """Calculate the Spearman correlation between the experimental and reference mass spectra.

        Returns
        -------
        coorelation : float
            Spearman correlation between the experimental and reference mass spectra.
        """
        # calculate Spearman correlation
        # ## TODO - Check axis
        correlation = spearmanr(
            self.zero_filled_u_l[0], self.zero_filled_u_l[1], axis=0
        )

        return correlation[0]

    def kendall_tau(self):
        """Calculate the Kendall's tau correlation between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            Kendall's tau correlation between the experimental and reference mass spectra."""
        # create dict['mz'] = abundance, for experimental data
        # self.ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data

        # calculate Kendall's tau
        correlation = kendalltau(self.zero_filled_u_l[0], self.zero_filled_u_l[1])

        return correlation[0]

    def dft_correlation(self):
        """Calculate the DFT correlation between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            DFT correlation between the experimental and reference mass spectra.
        """
        if self.n_x_y == 0:
            return 0

        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in self.exp_abun)

        x, y = self.nan_fill(self.df, fill_with=0)

        x, y = self.normalize(x, y, norm_func=self.normalize_func)

        # get the Fourier transform of x and y
        x_dft = rfft(x).real
        y_dft = rfft(y).real

        s_dft_xy = dot(x_dft, y_dft) / (norm(x_dft) * norm(y_dft))

        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation(nanfill=0)

        # final step
        s_dft = (n_x * s_wc_x_y + self.n_x_y * s_dft_xy) / (n_x + self.n_x_y)

        return s_dft

    def dwt_correlation(self):
        """Calculate the DWT correlation between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            DWT correlation between the experimental and reference mass spectra.

        Notes
        -----
        This function requires the PyWavelets library to be installed.
            This is not a default requirement as this function is not widely used.
        """

        from pywt import dwt

        if self.n_x_y == 0:
            return 0

        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in self.exp_abun)

        # calculate cosine correlation,
        x, y = self.nan_fill(self.df, fill_with=0)

        x, y = self.normalize(x, y, norm_func=self.normalize_func)

        # Make x and y into an array
        x_a = list(x)
        y_a = list(y)

        # get the wavelet transform of x and y (Daubechies with a filter length of 4. Asymmetric. pywavelets function)
        # Will only use the detail dwt (dwtDd
        x_dwtD = dwt(x_a, "db2")[1]
        y_dwtD = dwt(y_a, "db2")[1]

        s_dwt_xy = dot(x_dwtD, y_dwtD) / (norm(x_dwtD) * norm(y_dwtD))

        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation(nanfill=0)

        # final step
        s_dwt = (n_x * s_wc_x_y + self.n_x_y * s_dwt_xy) / (n_x + self.n_x_y)

        return s_dwt

    def euclidean_distance(self):
        """Calculate the Euclidean distance between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            Euclidean distance between the experimental and reference mass spectra.
        """
        # correlation = euclidean_distance_manual(self.zero_filled_u_l[0], self.zero_filled_u_l[1])
        qlist = self.zero_filled_u_l[0]
        rlist = self.zero_filled_u_l[1]

        correlation = sqrt(np_sum(power(qlist - rlist, 2)))

        return correlation

    def manhattan_distance(self):
        """Calculate the Manhattan distance between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            Manhattan distance between the experimental and reference mass spectra.
        """
        qlist = self.zero_filled_u_l[0]
        rlist = self.zero_filled_u_l[1]

        return np_sum(absolute(qlist - rlist))

    def jaccard_distance(self):
        """Calculate the Jaccard distance between the experimental and reference mass spectra.

        Returns
        -------
        correlation : float
            Jaccard distance between the experimental and reference mass spectra.
        """

        def jaccard_similarity(list1, list2):
            intersection = len(list(set(list1).intersection(list2)))
            union = (len(list1) + len(list2)) - intersection
            return float(intersection) / union

        qlist = self.zero_filled_u_l[0]
        rlist = self.zero_filled_u_l[1]

        return np_sum(power(qlist - rlist, 2)) / (
            np_sum(power(qlist, 2)) + np_sum(power(rlist, 2)) - np_sum(qlist * rlist)
        )
        # correlation = jaccard_similarity(self.zero_filled_u_l[0], self.zero_filled_u_l[1])
        # @return correlation

    def extra_distances(self):
        """Function to calculate distances using additional metrics defined in math_distance.py

        Currently, calculates all distances.

        Returns
        -------
        dict_res : dict
            Dictionary containing the distances between the experimental and reference mass spectra.

        """
        from corems.molecular_id.calc import math_distance

        # qlist = self.zero_filled_u_l[2]
        # rlist = self.zero_filled_u_l[3]

        dict_res = {}

        for method in methods_name:
            # function_name = method + "_distance"
            function_name = method
            if hasattr(math_distance, function_name):
                f = getattr(math_distance, function_name)

                if function_name == "canberra_metric":
                    x, y = self.nan_fill(self.df, fill_with=0)

                    qlist, rlist = self.normalize(x, y, norm_func=self.normalize_func)
                    # print("qlist:")
                    # print(qlist)
                    # print("rlist:")
                    # print(rlist)

                else:
                    qlist = self.zero_filled_u_l[0]
                    rlist = self.zero_filled_u_l[1]

                dist = f(qlist, rlist)
                # if method == "Minokowski_3":
                #    print("qlist:")
                #    print(qlist)
                #    print("rlist")
                #    print(rlist)
                #    exit()
                # if dist == np.nan or dis == np.inf:
                # print(self.exp_abun)
                # print(self.exp_mz)
                # print(function_name)
                # print(len(self.exp_abun))
                # print(len(self.exp_mz))
                # print(self.zero_filled_u_l[1])
                dict_res[method] = dist

        return dict_res
