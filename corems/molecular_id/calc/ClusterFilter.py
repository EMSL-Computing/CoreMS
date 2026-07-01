import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

# import matplotlib.pyplot as plt


class ClusteringFilter:
    """Class for filtering and clustering mass spectra data using various algorithms.

    Attributes
    -------
    mass_spectrum : MassSpectrum
        Mass spectrum object.
    ms_peaks : list
        List of mass peaks.
    ms_peak_indexes : list
        List of peak indexes.
    min_samples : int
        Minimum number of samples in a cluster.
    eps : float
        The maximum distance between two samples for one to be considered as in the neighborhood of the other.
    bandwidth : float
        Bandwidth used in MeanShift algorithm.
    quantile : float
        Quantile used in estimate_bandwidth function.
    n_samples : int
        Number of samples used in estimate_bandwidth function.
    bin_seeding : bool
        If true, initial kernel locations are not locations of all points, but rather the location of the discretized version of points, where points are binned onto a grid whose coarseness corresponds to the bandwidth. Setting this option to True will speed up the algorithm because fewer seeds will be initialized.
    min_peaks_per_class : int
        Minimum number of peaks per class.

    Methods
    -------
    * get_mass_error_matrix_data(ms_peaks).
        Get the mass error matrix data from a list of mass peaks.
    * get_kendrick_matrix_data(mass_spectrum).
        Get the Kendrick matrix data from a mass spectrum.
    * filter_kendrick(mass_spectrum).
        Filter the mass spectrum data using the Kendrick algorithm.
    * filter_kendrick_by_index(ms_peak_indexes, mass_spectrum_obj).
        Filter the mass spectrum data using the Kendrick algorithm based on a list of peak indexes.
    * remove_assignment_by_mass_error(mass_spectrum).
        Remove assignments from the mass spectrum based on mass error.


    """

    def get_mass_error_matrix_data(self, ms_peaks):
        """Get the mass error matrix data from a list of mass peaks.

        Parameters
        ----------
        ms_peaks : list
            List of mass peaks.

        Returns
        -------
        matrix_data : ndarray
            Matrix data containing mass and error values.
        list_indexes_mass_spec : list
            List of indexes of mass peaks in the original mass spectrum.
        """
        mass_list = list()
        error_list = list()
        list_indexes_mass_spec = []

        for index, mspeak in enumerate(ms_peaks):
            if mspeak.is_assigned:
                # print(mspeak.mz_exp, len(mspeak))
                for mformula in mspeak:
                    mass_list.append(mspeak.mz_exp)
                    error_list.append(mformula.mz_error)
                    list_indexes_mass_spec.append(index)

        kendrick_dict = {"mass": mass_list, "error": error_list}
        df = pd.DataFrame(kendrick_dict)
        matrix_data = df.values.astype("float32", copy=False)
        return matrix_data, list_indexes_mass_spec

    def get_kendrick_matrix_data(self, mass_spectrum):
        """Get the Kendrick matrix data from a mass spectrum.

        Parameters
        ----------
        mass_spectrum : MassSpectrum
            Mass spectrum object.

        Returns
        -------
        matrix_data : ndarray
            Matrix data containing Kendrick mass and Kendrick mass defect values.
        """
        km = mass_spectrum.kendrick_mass
        kmd = mass_spectrum.kmd
        kendrick_dict = {"km": km, "kmd": kmd}
        df = pd.DataFrame(kendrick_dict)
        matrix_data = df.values.astype("float32", copy=False)
        return matrix_data

    def filter_kendrick(self, mass_spectrum):
        """Filter the mass spectrum data using the Kendrick algorithm.

        Parameters
        ----------
        mass_spectrum : MassSpectrum
            Mass spectrum object.

        """
        matrix_data = self.get_kendrick_matrix_data(mass_spectrum)

        stdscaler = StandardScaler().fit(matrix_data)

        matrix_data_scaled = stdscaler.transform(matrix_data)

        clusters = DBSCAN(eps=0.75, min_samples=50).fit_predict(matrix_data_scaled)

        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(clusters)) - (1 if -1 in clusters else 0)
        n_noise_ = list(clusters).count(-1)

        indexes = []
        for i in range(len(clusters)):
            if clusters[i] == -1:
                indexes.append(i)

        if mass_spectrum.parameters.mass_spectrum.verbose_processing:
            print("Estimated number of clusters: %d" % n_clusters_)
            print("Estimated number of noise points: %d" % n_noise_)
        mass_spectrum.filter_by_index(indexes)
        # from matplotlib import pyplot as plt
        # plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="jet")
        # plt.xlabel("km")
        # plt.ylabel("kmd")
        # plt.show()
        # plt.close()

    def filter_kendrick_by_index(self, ms_peak_indexes, mass_spectrum_obj):
        """Filter the mass spectrum data using the Kendrick algorithm based on a list of peak indexes.

        Parameters
        ----------
        ms_peak_indexes : list
            List of peak indexes.
        mass_spectrum_obj : MassSpectrum
            Mass spectrum object.

        Returns
        -------
        noise_idx : list
            List of indexes of noise points in the mass spectrum.
        """
        min_samples = mass_spectrum_obj.molecular_search_settings.min_peaks_per_class

        kendrick_dict = {"km": list(), "kmd": list()}

        if len(ms_peak_indexes) <= 1:
            return []

        for index, _ in ms_peak_indexes:
            kendrick_dict["km"].append(mass_spectrum_obj[index].kendrick_mass)
            kendrick_dict["kmd"].append(mass_spectrum_obj[index].kmd)

        # check min data points otherwise StandardScaler().fit(0 will fail

        df = pd.DataFrame(kendrick_dict)
        matrix_data = df.values.astype("float32", copy=False)

        stdscaler = StandardScaler().fit(matrix_data)
        matrix_data_scaled = stdscaler.transform(matrix_data)

        clusters = DBSCAN(eps=0.8, min_samples=min_samples).fit_predict(
            matrix_data_scaled
        )

        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(clusters)) - (1 if -1 in clusters else 0)
        n_noise_ = list(clusters).count(-1)

        if mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Estimated number of clusters: %d" % n_clusters_)
            print("Estimated number of noise points: %d" % n_noise_)

        noise_idx = []

        other_peaks_idx = []

        for i in range(len(clusters)):
            if clusters[i] == -1:
                noise_idx.append(ms_peak_indexes[i])

            else:
                other_peaks_idx.append(ms_peak_indexes[i])

        # mfs = [mass_spectrum_obj[index].best_molecular_formula_candidate.string for index in other_peaks_idx]

        # mfs_noise = [mass_spectrum_obj[index].best_molecular_formula_candidate.string for index in noise_idx]

        # print(mfs)
        # print(mfs_noise)

        # from matplotlib import pyplot as plt
        # plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="jet")
        # plt.xlabel("km")
        # plt.ylabel("kmd")
        # plt.show()
        # plt.close()

        return noise_idx

    def remove_assignment_by_mass_error(self, mass_spectrum):
        """Remove assignments from the mass spectrum based on mass error.

        Parameters
        ----------
        mass_spectrum : MassSpectrum
            Mass spectrum object.

        """
        # data need to be binned by mz unit or more to be able to use clustering

        matrix_data, list_indexes_mass_spec = self.get_mass_error_matrix_data(
            mass_spectrum
        )

        stdscaler = StandardScaler().fit(matrix_data)

        matrix_data_scaled = stdscaler.transform(matrix_data)

        # bandwidth = estimate_bandwidth(matrix_data_scaled, quantile=0.3, n_samples=int(len(ms_peaks)/3))

        # clusters = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit_predict(matrix_data_scaled)

        # eps and min_samp need to be optimized by precision and number of mspeaks
        clusters = DBSCAN(eps=0.15).fit_predict(matrix_data_scaled)

        indexes = []

        # from matplotlib import pyplot as plt
        # plt.scatter(matrix_data[:, 0], matrix_data[:, 1], c=clusters, cmap="plasma")
        # plt.xlabel("km")
        # plt.ylabel("kmd")
        # plt.show()
        # plt.close()

        for i in range(len(clusters)):
            if clusters[i] == -1:
                indexes.append(list_indexes_mass_spec[i])

        mass_spectrum.remove_assignment_by_index(indexes)
