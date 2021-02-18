__author__ = "Yuri E. Corilo"
__date__ = "Feb 13, 2020"


from collections.abc import Mapping
from pathlib import Path
import json

from numpy import array

from corems.mass_spectra.calc.GC_Calc import GC_Calculations
from corems.mass_spectra.calc.GC_Deconvolution import MassDeconvolution
from corems.mass_spectra.calc import SignalProcessing as sp

from corems.chroma_peak.factory.ChromaPeakClasses import GCPeak
from corems.mass_spectra.output.export import LowResGCMSExport
from corems.encapsulation.factory.parameters import GCMSParameters


class GCMSBase(GC_Calculations, MassDeconvolution):
    """
    classdocs
    """

    def __init__(self, file_location, analyzer='Unknown', instrument_label='Unknown', sample_name=None):

        """
            # Parameters
        ----------
        file_location: text,  pathlib.Path(), or s3path.S3Path
            Path object from pathlib containing the file location
        """
        if isinstance(file_location, str):
            # if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():

            raise FileExistsError("File does not exist: " + str(file_location))

        self.file_location = file_location

        if sample_name: self.sample_name = sample_name
        else: self.sample_name = file_location.stem

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

        self._parameters = GCMSParameters()

    def __len__(self):

        return len(self.gcpeaks)

    def __getitem__(self, scan_number):

        return self.gcpeaks[scan_number]

    # def __iter__(self):

    #     return iter(self.gcpeaks.values())

    def process_chromatogram(self, plot_res=False):

        # tic = self.tic - self.baseline_detector(self.tic)

        self._processed_tic = self.smooth_tic(self.tic)

        for index, tic in enumerate(self._processed_tic):

            self._ms[index]._processed_tic = tic

        # self.second_derivative_threshold(self._processed_tic)

        if self.chromatogram_settings.use_deconvolution:

            self.run_deconvolution(plot_res=False)

        else:

            peaks_index = self.centroid_detector(self._processed_tic, self.retention_time)

            for i in peaks_index:

                apex_index = i[1]

                gc_peak = GCPeak(self, self._ms[apex_index], i )

                gc_peak.calc_area(self._processed_tic, 1)

                self.gcpeaks.append(gc_peak)

                # self.gcpeaks[self.scans_number[apex_index]] = gc_peak

    def add_mass_spectrum(self, mass_spec):

        self._ms[mass_spec.scan_number] = mass_spec

    def set_tic_list_from_data(self):

        self.tic = [self._ms.get(i).tic for i in self.scans_number]

        # self.set_tic_list([self._ms.get(i).get_sumed_signal_to_noise() for i in self.get_scans_number()])

    def set_retention_time_from_data(self):

        retention_time_list = []

        for key_ms in sorted(self._ms.keys()):

            retention_time_list.append(self._ms.get(key_ms).rt)

        self.retention_time = retention_time_list

        # self.set_retention_time_list(sorted(self._ms.keys()))

    def set_scans_number_from_data(self):

        self.scans_number = sorted(self._ms.keys())

    @property
    def parameter(self):
        return self._parameters

    @parameter.setter
    def parameter(self, gcms_parameters_instance):
        self._parameters = gcms_parameters_instance

    @property
    def molecular_search_settings(self):
        return self.parameter.molecular_search

    @molecular_search_settings.setter
    def molecular_search_settings(self, settings_class_instance):

        self.parameter.molecular_search = settings_class_instance

    @property
    def chromatogram_settings(self):
        return self.parameter.gc_ms

    @chromatogram_settings.setter
    def chromatogram_settings(self, settings_class_instance):

        self.parameter.gc_ms = settings_class_instance

    @property
    def scans_number(self):

        return self._scans_number_list

    @property
    def retention_time(self):

        return self._retention_time_list

    @property
    def processed_tic(self):

        return self._processed_tic

    @property
    def tic(self):

        return self._tic_list

    @property
    def max_tic(self):

        return max([gc_peak.tic for gc_peak in self])

    @property
    def min_tic(self):

        return min([gc_peak.tic for gc_peak in self])

    @property
    def dynamic_range(self):

        return self.max_tic / self.min_tic

    @property
    def matched_peaks(self):
        return [gc_peak for gc_peak in self if gc_peak]

    @property
    def sorted_gcpeaks(self):
        return sorted(self, key=lambda g: g.rt)

    @property
    def unique_metabolites(self):

        metabolites = set()
        for gc_peak in self:
            if gc_peak:
                for compound_obj in gc_peak:
                    metabolites.add(compound_obj.name)

        return metabolites

    @property
    def no_matched_peaks(self):
        return [peak for peak in self if not peak]

    @retention_time.setter
    def retention_time(self, l):
        # self._retention_time_list = linspace(0, 80, num=len(self._scans_number_list))
        self._retention_time_list = l

    @scans_number.setter
    def scans_number(self, l):

        self._scans_number_list = l

    @tic.setter
    def tic(self,l):

        self._tic_list = array(l)    

    def plot_gc_peaks(self, ax=None, color="red"):  # pragma: no cover

        import matplotlib.pyplot as plt
        fig = plt.gcf()
        if ax is None:
            ax = plt.gca()

        max_rts = [gc_peak.mass_spectrum.rt for gc_peak in self]
        max_tics = [gc_peak.mass_spectrum.tic for gc_peak in self]

        # min_rts = [self._ms[gc_peak.start_index].rt for gc_peak in self] + [self._ms[gc_peak.final_index].rt for gc_peak in self]
        # min_tics = [self._ms[gc_peak.start_index].tic for gc_peak in self] + [self._ms[gc_peak.final_index].tic for gc_peak in self]
        # sc = ax.scatter(min_rts, min_tics, color='yellow', linewidth=0, marker='v')

        sc = ax.scatter(max_rts, max_tics, color=color, marker='v')

        ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')

        annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        annot.get_bbox_patch().set_facecolor(('lightblue'))
        annot.get_bbox_patch().set_alpha(0.8)

        def update_annot(ind):

            pos = sc.get_offsets()[ind["ind"][0]]
            annot.xy = pos

            text = "RT: {}\nRT Ref: {}\nRI: {}\nRI Ref: {}\nSimilarity Score: {}\nName: {}".format(" ".join([str(round(self[n].rt, 2)) for n in ind["ind"]]),
                           " ".join([str(round(self[n].highest_score_compound.rt, 2) if self[n].highest_score_compound else None) for n in ind["ind"]]),
                           " ".join([str(round(self[n].ri, 2) if self[n].ri else None) for n in ind["ind"]]),
                           " ".join([str(round(self[n].highest_score_compound.ri, 2) if self[n].highest_score_compound else None) for n in ind["ind"]]),                           
                           " ".join([str(round(self[n].highest_score_compound.similarity_score, 4) if self[n].highest_score_compound else None) for n in ind["ind"]]),
                           " ".join([str(self[n].highest_score_compound.name if self[n].highest_score_compound else None) for n in ind["ind"]])
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

    def to_excel(self, out_file_path, write_mode='ab', id_label="corems:"):

        exportMS = LowResGCMSExport(out_file_path, self)
        exportMS.to_excel(id_label=id_label)

    def to_csv(self, out_file_path, write_mode='ab', id_label="corems:"):

        exportMS = LowResGCMSExport(out_file_path, self)
        exportMS.to_csv(id_label=id_label)

    def to_pandas(self, out_file_path, id_label="corems:"):

        # pickle dataframe (pkl extension)
        exportMS = LowResGCMSExport(out_file_path, self)
        exportMS.to_pandas(id_label=id_label)

    def to_dataframe(self, id_label="corems:"):

        # returns pandas dataframe
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_pandas_df(id_label=id_label)

    def processing_stats(self):

        # returns json string
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_data_stats(self)

    def parameters_json(self, id_label="corems:", output_path=" "):

        # returns json string
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_parameters_json(self, id_label, output_path)

    def to_json(self, id_label="corems:"):

        # returns pandas dataframe
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.get_json(id_label=id_label)

    def to_hdf(self, id_label="corems:"):

        # returns pandas dataframe
        exportMS = LowResGCMSExport(self.sample_name, self)
        return exportMS.to_hdf(id_label=id_label)

    def plot_chromatogram(self, ax=None, color="blue"):  # pragma: no cover

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        ax.plot(self.retention_time, self.tic, color=color)
        ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')

        return ax

    def plot_smoothed_chromatogram(self, ax=None, color="green"):  # pragma: no cover

        import matplotlib.pyplot as plt

        if ax is None:

            ax = plt.gca()

        ax.plot(self.retention_time, self.smooth_tic(self.tic), color=color)

        ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')

        return ax

    def plot_detected_baseline(self, ax=None, color="blue"):  # pragma: no cover

        import matplotlib.pyplot as plt

        if ax is None:

            ax = plt.gca()

        max_height = self.chromatogram_settings.peak_height_max_percent
        max_prominence = self.chromatogram_settings.peak_max_prominence_percent

        baseline = sp.baseline_detector(self.tic, self.retention_time, max_height, max_prominence)
        ax.plot(self.retention_time, color=color)
        ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')

        return ax

    def plot_baseline_subtraction(self, ax=None, color="black"):  # pragma: no cover

        import matplotlib.pyplot as plt

        if ax is None:

            ax = plt.gca()

        max_height = self.chromatogram_settings.peak_height_max_percent

        max_prominence = self.chromatogram_settings.peak_max_prominence_percent

        x = self.tic + sp.baseline_detector(self.tic, self.retention_time, max_height, max_prominence)

        ax.plot(self.retention_time, x, color=color)

        ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')

        return ax

    def peaks_rt_tic(self, json_string=False):

        peaks_list = dict()

        all_candidates_data = {}

        all_peaks_data = {}

        for gcms_peak in self.sorted_gcpeaks:

            dict_data = {'rt': gcms_peak.rt_list,
                         'tic': gcms_peak.tic_list,
                         'mz': gcms_peak.mass_spectrum.mz_exp.tolist(),
                         'abundance': gcms_peak.mass_spectrum.abundance.tolist(),
                         'candidate_names': gcms_peak.compound_names,
                         }

            peaks_list[gcms_peak.rt] = dict_data

            for compound in gcms_peak:

                if compound.name not in all_candidates_data.keys():
                    mz = array(compound.mz).tolist()
                    abundance = array(compound.abundance).tolist()
                    data = {'mz': mz, "abundance": abundance}
                    all_candidates_data[compound.name] = data

        all_peaks_data["peak_data"] = peaks_list
        all_peaks_data["ref_data"] = all_candidates_data

        if json_string:

            return json.dumps(all_peaks_data)

        else:
            return all_peaks_data

    def plot_processed_chromatogram(self, ax=None, color="black"):

        import matplotlib.pyplot as plt

        if ax is None:

            ax = plt.gca()

        ax.plot(self.retention_time, self.processed_tic, color=color)

        ax.set(xlabel='Retention Time (s)', ylabel='Total Ion Chromatogram')

        return ax
