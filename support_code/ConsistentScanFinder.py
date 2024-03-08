import pandas as pd

def get_best_scans_idx(thermo_parser, stdevs=2, method="mean", plot=False):
        """
        Method to determine the best scan indexes for selective co-addition.

        Parameters
        ----------
        thermo_parser : ImportMassSpectraThermoMSFileReader
            An instance of ImportMassSpectraThermoMSFileReader
        stdevs : int, optional
            The number of standard deviations to use as the cutoff for filtering out datapoints. Default is 2.
        method : str, optional
            The method to calculate the mean or median of the TIC values. Default is "mean".
        plot : bool, optional
            Whether to plot the TIC with horizontal lines for the standard deviation cutoffs. Default is False.

        Notes
        -----
        This method calculates the mean (default) or median of the TIC values and sets an upper and lower limit
        based on a specified number of standard deviations. The scans with TIC values outside of this range are
        considered the best scans for selective co-addition.

        Empirically, using 1-2 standard deviations is enough to filter out the worst datapoints.

        If `plot` is True, a matplotlib figure is returned along with the list of scan indexes.

        Examples
        --------
        >>> reader = ImportMassSpectraThermoMSFileReader()
        >>> scans = get_best_scans_idx(reader, stdevs=2, method="mean", plot=True)
        """
        tic = pd.DataFrame(thermo_parser.get_tic(plot=plot))

        if method == "median":
            tic_median = tic["TIC"].median()
        elif method == "mean":
            tic_median = tic["TIC"].mean()
        else:
            print("Method " + str(method) + " undefined")

        tic_std = tic["TIC"].std()

        upperlimit = tic_median - (stdevs * tic_std)
        lowerlimit = tic_median + (stdevs * tic_std)

        tic_filtered = tic[(tic["TIC"] > upperlimit) & (tic["TIC"] < lowerlimit)]
        scans = list(tic_filtered.Scans.values)

        if plot:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(tic["Time"], tic["TIC"])
            ax.axhline(y=upperlimit, c="r")
            ax.axhline(y=lowerlimit, c="r")
            return fig, scans
        else:
            return scans