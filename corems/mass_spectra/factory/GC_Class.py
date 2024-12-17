__author__ = "Yuri E. Corilo"
__date__ = "Feb 13, 2020"


from collections.abc import Mapping
from pathlib import Path
import json

from numpy import array


from corems.mass_spectra.calc.GC_Calc import GC_Calculations
from corems.mass_spectra.calc.GC_Deconvolution import MassDeconvolution
from corems.mass_spectra.calc import SignalProcessing as sp

from corems.chroma_peak.factory.chroma_peak_classes import GCPeak
from corems.mass_spectra.output.export import LowResGCMSExport
from corems.encapsulation.factory.parameters import GCMSParameters


class GCMSBase(GC_Calculations, MassDeconvolution):
    """Base class for GC-MS data processing.

    Parameters
    ----
    file_location : str, pathlib.Path, or s3path.S3Path
        Path object containing the file location.
    analyzer : str, optional
        Name of the analyzer. Defaults to 'Unknown'.
    instrument_label : str, optional
        Label of the instrument. Defaults to 'Unknown'.
    sample_name : str, optional
        Name of the sample. If not provided, it is derived from the file location.

    Attributes
    ------------
    file_location : pathlib.Path
        Path object containing the file location.
    sample_name : str
        Name of the sample.
    analyzer : str
        Name of the analyzer.
    instrument_label : str
        Label of the instrument.
    gcpeaks : list
        List of GCPeak objects.
    ri_pairs_ref : None
        Reference retention index pairs.
    cal_file_path : None
        Calibration file path.
    _parameters : GCMSParameters
        GC-MS parameters.
    _retention_time_list : list
        List of retention times.
    _scans_number_list : list
        List of scan numbers.
    _tic_list : list
        List of total ion chromatogram values.
    _ms : dict
        Dictionary containing all mass spectra.
    _processed_tic : list
        List of processed total ion chromatogram values.

    Methods
    -------
    * process_chromatogram(plot_res=False). Process the chromatogram.
    * plot_gc_peaks(ax=None, color='red'). Plot the GC peaks.
    """

    def __init__(
        self,
        file_location,
        analyzer="Unknown",
        instrument_label="Unknown",
        sample_name=None,
    ):
        if isinstance(file_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():
            raise FileExistsError("File does not exist: " + str(file_location))

        self.file_location = file_location

        if sample_name:
            self.sample_name = sample_name
        else:
            self.sample_name = file_location.stem

        self.analyzer = analyzer
        self.instrument_label = instrument_label
        self._init_settings()

        self._retention_time_list = []
        self._scans_number_list = []
        self._tic_list = []

        # all scans
        self._ms = {}

        # after peak detection
        self._processed_tic = []
        self.gcpeaks = []

        self.ri_pairs_ref = None
        self.cal_file_path = None

    def _init_settings(self):
        """Initialize the settings for GC_Class.

        This method initializes the settings for the GC_Class object using the GCMSParameters class.
        """
        self._parameters = GCMSParameters()

    def __len__(self):
        """Return the number of GC peaks in the GC_Class object."""
        return len(self.gcpeaks)

    def __getitem__(self, scan_number) -> GCPeak:
        """Return the GCPeak with the given scan number."""
        return self.gcpeaks[scan_number]

    # def __iter__(self):

    #     return iter(self.gcpeaks.values())

    def process_chromatogram(self, plot_res=False):
        """Process the chromatogram.

        This method processes the chromatogram.

        Parameters
        ----------
        plot_res : bool, optional
            If True, plot the results. Defaults to False.
        """

        # tic = self.tic - self.baseline_detector(self.tic)

        self._processed_tic = self.smooth_tic(self.tic)

        for index, tic in enumerate(self._processed_tic):
            self._ms[index]._processed_tic = tic

        # self.second_derivative_threshold(self._processed_tic)

        if self.chromatogram_settings.use_deconvolution:
            self.run_deconvolution(plot_res=False)

        else:
            peaks_index = self.centroid_detector(
                self._processed_tic, self.retention_time
            )

            for i in peaks_index:
                apex_index = i[1]

                gc_peak = GCPeak(self, self._ms[apex_index], i)

                gc_peak.calc_area(self._processed_tic, 1)

                self.gcpeaks.append(gc_peak)

                # self.gcpeaks[self.scans_number[apex_index]] = gc_peak

    def add_mass_spectrum(self, mass_spec):
        """Add a mass spectrum to the GC-MS object.

        This method adds a mass spectrum to the GC-MS object.

        Parameters
        ----------
        mass_spec : MassSpectrum
            Mass spectrum to be added.
        """

        self._ms[mass_spec.scan_number] = mass_spec

    def set_tic_list_from_data(self):
        """Set the total ion chromatogram list from the mass spectra data within the GC-MS data object."""

        self.tic = [self._ms.get(i).tic for i in self.scans_number]

        # self.set_tic_list([self._ms.get(i).get_sumed_signal_to_noise() for i in self.get_scans_number()])

    def set_retention_time_from_data(self):
        """Set the retention time list from the mass spectra data within the GC-MS data object."""

        retention_time_list = []

        for key_ms in sorted(self._ms.keys()):
            retention_time_list.append(self._ms.get(key_ms).retention_time)

        self.retention_time = retention_time_list

        # self.set_retention_time_list(sorted(self._ms.keys()))

    def set_scans_number_from_data(self):
        """Set the scan number list from the mass spectra data within the GC-MS data object."""

        self.scans_number = sorted(self._ms.keys())

    @property
    def parameters(self):
        """GCMS Parameters"""
        return self._parameters

    @parameters.setter
    def parameters(self, gcms_parameters_instance):
        self._parameters = gcms_parameters_instance

    # Note: maintaining `parameter` for backwards compatibility,
    # but proper usage would reference `parameters` to conform
    # to other classes.
    @property
    def parameter(self):
        """GCMS Parameters"""
        return self._parameters

    @parameter.setter
    def parameter(self, gcms_parameters_instance):
        self._parameters = gcms_parameters_instance

    @property
    def molecular_search_settings(self):
        """Molecular Search Settings"""
        return self.parameters.molecular_search

    @molecular_search_settings.setter
    def molecular_search_settings(self, settings_class_instance):
        self.parameters.molecular_search = settings_class_instance

    @property
    def chromatogram_settings(self):
        """Chromatogram Settings"""
        return self.parameters.gc_ms

    @chromatogram_settings.setter
    def chromatogram_settings(self, settings_class_instance):
        self.parameters.gc_ms = settings_class_instance

    @property
    def scans_number(self):
        """Scans Number"""
        return self._scans_number_list

    @property
    def retention_time(self):
        """Retention Time"""
        return self._retention_time_list

    @property
    def processed_tic(self):
        """Processed Total Ion Current"""
        return self._processed_tic

    @property
    def tic(self):
        """Total Ion Current"""
        return self._tic_list

    @property
    def max_tic(self):
        """Maximum Total Ion Current"""
        return max([gc_peak.tic for gc_peak in self])

    @property
    def min_tic(self):
        """Minimum Total Ion Current"""
        return min([gc_peak.tic for gc_peak in self])

    @property
    def dynamic_range(self):
        """Dynamic Range of the Total Ion Current"""
        return self.max_tic / self.min_tic

    @property
    def matched_peaks(self):
        """Matched Peaks"""
        return [gc_peak for gc_peak in self if gc_peak]

    @property
    def sorted_gcpeaks(self):
        """Sorted GC Peaks, by retention time"""
        return sorted(self, key=lambda g: g.retention_time)

    @property
    def unique_metabolites(self):
        """Unique Metabolites"""
        metabolites = set()
        for gc_peak in self:
            if gc_peak:
                for compound_obj in gc_peak:
                    metabolites.add(compound_obj.name)

        return metabolites

    @property
    def metabolites_data(self):
        """Metabolites Data"""
        metabolites = {}
        for gc_peak in self:
            if gc_peak:
                for compound_obj in gc_peak:
                    if compound_obj.name in metabolites.keys():
                        current_score = metabolites[compound_obj.name][
                            "highest_similarity_score"
                        ]
                        compound_score = compound_obj.spectral_similarity_score
                        metabolites[compound_obj.name]["highest_similarity_score"] = (
                            compound_score
                            if compound_score > current_score
                            else current_score
                        )

                    else:
                        if compound_obj.metadata:
                            metabolites[compound_obj.name] = {
                                "name": compound_obj.name,
                                "highest_similarity_score": compound_obj.spectral_similarity_score,
                                "casno": compound_obj.metadata.cas,
                                "kegg": compound_obj.metadata.kegg,
                                "inchi": compound_obj.metadata.inchi,
                                "inchi_key": compound_obj.metadata.inchikey,
                                "chebi": compound_obj.metadata.chebi,
                                "smiles": compound_obj.metadata.smiles,
                            }
                        else:
                            metabolites[compound_obj.name] = {
                                "name": compound_obj.name,
                                "highest_similarity_score": compound_obj.spectral_similarity_score,
                                "casno": "",
                                "kegg": "",
                                "inchi": "",
                                "inchikey": "",
                                "chebi": "",
                                "smiles": "",
                            }

        return list(metabolites.values())

    @property
    def no_matched_peaks(self):
        """Peaks with no Matched Metabolites"""
        return [peak for peak in self if not peak]

    @retention_time.setter
    def retention_time(self, alist):
        # self._retention_time_list = linspace(0, 80, num=len(self._scans_number_list))
        self._retention_time_list = alist

    @scans_number.setter
    def scans_number(self, alist):
        self._scans_number_list = alist

    @tic.setter
    def tic(self, alist):
        self._tic_list = array(alist)

    def plot_gc_peaks(self, ax=None, color="red"):  # pragma: no cover
        """Plot the GC peaks.

        This method plots the GC peaks.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot the GC peaks. Defaults to None.
        color : str, optional
            Color of the GC peaks. Defaults to 'red'.
        """

        import matplotlib.pyplot as plt

        fig = plt.gcf()
        if ax is None:
            ax = plt.gca()

        max_rts = [gc_peak.mass_spectrum.retention_time for gc_peak in self]
        max_tics = [gc_peak.mass_spectrum.tic for gc_peak in self]

        # min_rts = [self._ms[gc_peak.start_index].retention_time for gc_peak in self] + [self._ms[gc_peak.final_index].retention_time for gc_peak in self]
        # min_tics = [self._ms[gc_peak.start_index].tic for gc_peak in self] + [self._ms[gc_peak.final_index].tic for gc_peak in self]
        # sc = ax.scatter(min_rts, min_tics, color='yellow', linewidth=0, marker='v')

        sc = ax.scatter(max_rts, max_tics, color=color, marker="v")

        ax.set(xlabel="Retention Time (s)", ylabel="Total Ion Chromatogram")

        annot = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(20, 20),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w"),
            arrowprops=dict(arrowstyle="->"),
        )
        annot.set_visible(False)
        annot.get_bbox_patch().set_facecolor(("lightblue"))
        annot.get_bbox_patch().set_alpha(0.8)

        def update_annot(ind):
            pos = sc.get_offsets()[ind["ind"][0]]
            annot.xy = pos

            text = "RT: {}\nRT Ref: {}\nRI: {}\nRI Ref: {}\nSimilarity Score: {}\nName: {}".format(
                " ".join([str(round(self[n].retention_time, 2)) for n in ind["ind"]]),
                " ".join(
                    [
                        str(
                            round(self[n].highest_score_compound.retention_time, 2)
                            if self[n].highest_score_compound
                            else None
                        )
                        for n in ind["ind"]
                    ]
                ),
                " ".join(
                    [
                        str(round(self[n].ri, 2) if self[n].ri else None)
                        for n in ind["ind"]
                    ]
                ),
                " ".join(
                    [
                        str(
                            round(self[n].highest_score_compound.ri, 2)
                            if self[n].highest_score_compound
                            else None
                        )
                        for n in ind["ind"]
                    ]
                ),
                " ".join(
                    [
                        str(
                            round(self[n].highest_score_compound.similarity_score, 4)
                            if self[n].highest_score_compound
                            else None
                        )
                        for n in ind["ind"]
                    ]
                ),
                " ".join(
                    [
                        str(
                            self[n].highest_score_compound.name
                            if self[n].highest_score_compound
                            else None
                        )
                        for n in ind["ind"]
                    ]
                ),
            )
            annot.set_text(text)

        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = sc.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", hover)

        return ax

    def to_excel(
        self, out_file_path, write_mode="ab", write_metadata=True, id_label="corems:"
    ):
        """Export the GC-MS data to an Excel file.

        This method exports the GC-MS data to an Excel file.

        Parameters
        ----------
        out_file_path : str, pathlib.Path, or s3path.S3Path
            Path object containing the file location.
        write_mode : str, optional
            Write mode. Defaults to 'ab'.
        write_metadata : bool, optional
            If True, write the metadata. Defaults to True.
        id_label : str, optional
            Label of the ID. Defaults to 'corems:'.

        """

        if isinstance(out_file_path, str):
            out_file_path = Path(out_file_path)

        exportMS = LowResGCMSExport(out_file_path, self)
        exportMS.to_excel(
            id_label=id_label, write_mode=write_mode, write_metadata=write_metadata
        )

        return out_file_path.with_suffix(".xlsx")

    def to_csv(
        self,
        out_file_path,
        separate_output=False,
        write_metadata=True,
        id_label="corems:",
    ):
        """Export the GC-MS data to a CSV file.

        Parameters
        ----------
        out_file_path : str, pathlib.Path, or s3path.S3Path
            Path object containing the file location.
        separate_output : bool, optional
            If True, separate the output. Defaults to False.
        write_metadata : bool, optional
            If True, write the metadata. Defaults to True.

        """

        if isinstance(out_file_path, str):
            out_file_path = Path(out_file_path)

        exportMS = LowResGCMSExport(out_file_path, self)
        exportMS.to_csv(
            id_label=id_label,
            separate_output=separate_output,
            write_metadata=write_metadata,
        )

        return out_file_path.with_suffix(".csv")

    def to_pandas(self, out_file_path, write_metadata=True, id_label="corems:"):
        """Export the GC-MS data to a Pandas dataframe.

        Parameters
        ----------
        out_file_path : str, pathlib.Path, or s3path.S3Path
            Path object containing the file location.
        write_metadata : bool, optional
            If True, write the metadata. Defaults to True.
        id_label : str, optional
            Label of the ID. Defaults to 'corems:'.

        """

        if isinstance(out_file_path, str):
            out_file_path = Path(out_file_path)
        # pickle dataframe (pkl extension)
        exportMS = LowResGCMSExport(out_file_path, self)
        exportMS.to_pandas(id_label=id_label, write_metadata=write_metadata)

        return out_file_path.with_suffix(".pkl")

    def to_dataframe(self, id_label="corems:"):
        """Export the GC-MS data to a Pandas dataframe.

        Parameters
        ----------
        id_label : str, optional
            Label of the ID. Defaults to 'corems:'.

        """

        # returns pandas dataframe
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_pandas_df(id_label=id_label)

    def processing_stats(self):
        """Return the processing statistics."""

        # returns json string
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_data_stats(self)

    def parameters_json(self, id_label="corems:", output_path=" "):
        """Return the parameters in JSON format.

        Parameters
        ----------
        id_label : str, optional
            Label of the ID. Defaults to 'corems:'.
        output_path : str, optional
            Path object containing the file location. Defaults to " ".
        """

        # returns json string
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_parameters_json(self, id_label, output_path)

    def to_json(self, id_label="corems:"):
        """Export the GC-MS data to a JSON file.

        Parameters
        ----------
        id_label : str, optional
            Label of the ID. Defaults to 'corems:'.

        """

        # returns pandas dataframe
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_json(id_label=id_label)

    def to_hdf(self, id_label="corems:"):
        """Export the GC-MS data to a HDF file.

        Parameters
        ----------
        id_label : str, optional
            Label of the ID. Defaults to 'corems:'.

        """

        # returns pandas dataframe
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.to_hdf(id_label=id_label)

    def plot_chromatogram(self, ax=None, color="blue"):  # pragma: no cover
        """Plot the chromatogram.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot the chromatogram. Defaults to None.
        color : str, optional
            Color of the chromatogram. Defaults to 'blue'.

        """

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        ax.plot(self.retention_time, self.tic, color=color)
        ax.set(xlabel="Retention Time (s)", ylabel="Total Ion Chromatogram")

        return ax

    def plot_smoothed_chromatogram(self, ax=None, color="green"):  # pragma: no cover
        """Plot the smoothed chromatogram.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot the smoothed chromatogram. Defaults to None.
        color : str, optional
            Color of the smoothed chromatogram. Defaults to 'green'.

        """

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        ax.plot(self.retention_time, self.smooth_tic(self.tic), color=color)

        ax.set(xlabel="Retention Time (s)", ylabel="Total Ion Chromatogram")

        return ax

    def plot_detected_baseline(self, ax=None, color="blue"):  # pragma: no cover
        """Plot the detected baseline.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot the detected baseline. Defaults to None.
        color : str, optional
            Color of the detected baseline. Defaults to 'blue'.

        """

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        max_height = self.chromatogram_settings.peak_height_max_percent
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent

        baseline = sp.baseline_detector(
            self.tic, self.retention_time, max_height, max_prominence
        )
        ax.plot(self.retention_time, color=color)
        ax.set(xlabel="Retention Time (s)", ylabel="Total Ion Chromatogram")

        return ax

    def plot_baseline_subtraction(self, ax=None, color="black"):  # pragma: no cover
        """Plot the baseline subtraction.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot the baseline subtraction. Defaults to None.
        color : str, optional
            Color of the baseline subtraction. Defaults to 'black'.

        """

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        max_height = self.chromatogram_settings.peak_height_max_percent

        max_prominence = self.chromatogram_settings.peak_max_prominence_percent

        x = self.tic + sp.baseline_detector(
            self.tic, self.retention_time, max_height, max_prominence
        )

        ax.plot(self.retention_time, x, color=color)

        ax.set(xlabel="Retention Time (s)", ylabel="Total Ion Chromatogram")

        return ax

    def peaks_rt_tic(self, json_string=False):
        """Return the peaks, retention time, and total ion chromatogram.

        Parameters
        ----------
        json_string : bool, optional
            If True, return the peaks, retention time, and total ion chromatogram in JSON format. Defaults to False.

        """

        peaks_list = dict()

        all_candidates_data = {}

        all_peaks_data = {}

        for gcms_peak in self.sorted_gcpeaks:
            dict_data = {
                "rt": gcms_peak.rt_list,
                "tic": gcms_peak.tic_list,
                "mz": gcms_peak.mass_spectrum.mz_exp.tolist(),
                "abundance": gcms_peak.mass_spectrum.abundance.tolist(),
                "candidate_names": gcms_peak.compound_names,
            }

            peaks_list[gcms_peak.retention_time] = dict_data

            for compound in gcms_peak:
                if compound.name not in all_candidates_data.keys():
                    mz = array(compound.mz).tolist()
                    abundance = array(compound.abundance).tolist()
                    data = {"mz": mz, "abundance": abundance}
                    all_candidates_data[compound.name] = data

        all_peaks_data["peak_data"] = peaks_list
        all_peaks_data["ref_data"] = all_candidates_data

        if json_string:
            return json.dumps(all_peaks_data)

        else:
            return all_peaks_data

    def plot_processed_chromatogram(self, ax=None, color="black"):
        """Plot the processed chromatogram.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot the processed chromatogram. Defaults to None.
        color : str, optional
            Color of the processed chromatogram. Defaults to 'black'.

        """

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        ax.plot(self.retention_time, self.processed_tic, color=color)

        ax.set(xlabel="Retention Time (s)", ylabel="Total Ion Chromatogram")

        return ax
