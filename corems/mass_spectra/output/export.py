__author__ = "Yuri E. Corilo"
__date__ = "Dec 14, 2010"


import csv
import json
import uuid
import warnings
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from numpy import concatenate
from openpyxl import load_workbook
from pandas import DataFrame, ExcelWriter, read_excel

from corems import __version__, corems_md5
from corems.encapsulation.output import parameter_to_dict
from corems.encapsulation.output.parameter_to_json import (
    dump_lcms_settings_json,
    dump_lcms_settings_toml,
)
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecfromFreq
from corems.mass_spectrum.output.export import HighResMassSpecExport
from corems.molecular_id.calc.SpectralSimilarity import methods_name


class LowResGCMSExport:
    """A class to export low resolution GC-MS data.

    This class provides methods to export low resolution GC-MS data to various formats such as Excel, CSV, HDF5, and Pandas DataFrame.

    Parameters:
    ----------
    out_file_path : str
        The output file path.
    gcms : object
        The low resolution GCMS object.

    Attributes:
    ----------
    output_file : Path
        The output file path as a Path object.
    gcms : object
        The low resolution GCMS object.

    Methods:
    -------
    * get_pandas_df(id_label="corems:"). Get the exported data as a Pandas DataFrame.
    * get_json(nan=False, id_label="corems:"). Get the exported data as a JSON string.
    * to_pandas(write_metadata=True, id_label="corems:"). Export the data to a Pandas DataFrame and save it as a pickle file.
    * to_excel(write_mode='a', write_metadata=True, id_label="corems:"),
        Export the data to an Excel file.
    * to_csv(separate_output=False, write_mode="w", write_metadata=True, id_label="corems:").
        Export the data to a CSV file.
    * to_hdf(id_label="corems:").
        Export the data to an HDF5 file.
    * get_data_stats(gcms).
        Get statistics about the GCMS data.

    """

    def __init__(self, out_file_path, gcms):
        self.output_file = Path(out_file_path)

        self.gcms = gcms

        self._init_columns()

    def _init_columns(self):
        """Initialize the column names for the exported data.

        Returns:
        -------
        list
            The list of column names.
        """

        columns = [
            "Sample name",
            "Peak Index",
            "Retention Time",
            "Retention Time Ref",
            "Peak Height",
            "Peak Area",
            "Retention index",
            "Retention index Ref",
            "Retention Index Score",
            "Similarity Score",
            "Spectral Similarity Score",
            "Compound Name",
            "Chebi ID",
            "Kegg Compound ID",
            "Inchi",
            "Inchi Key",
            "Smiles",
            "Molecular Formula",
            "IUPAC Name",
            "Traditional Name",
            "Common Name",
            "Derivatization",
        ]

        if self.gcms.molecular_search_settings.exploratory_mode:
            columns.extend(
                [
                    "Weighted Cosine Correlation",
                    "Cosine Correlation",
                    "Stein Scott Similarity",
                    "Pearson Correlation",
                    "Spearman Correlation",
                    "Kendall Tau Correlation",
                    "Euclidean Distance",
                    "Manhattan Distance",
                    "Jaccard Distance",
                    "DWT Correlation",
                    "DFT Correlation",
                ]
            )

            columns.extend(list(methods_name.values()))

        return columns

    def get_pandas_df(self, id_label="corems:"):
        """Get the exported data as a Pandas DataFrame.

        Parameters:
        ----------
        id_label : str, optional
            The ID label for the data. Default is "corems:".

        Returns:
        -------
        DataFrame
            The exported data as a Pandas DataFrame.
        """

        columns = self._init_columns()

        dict_data_list = self.get_list_dict_data(self.gcms)

        df = DataFrame(dict_data_list, columns=columns)

        df.name = self.gcms.sample_name

        return df

    def get_json(self, nan=False, id_label="corems:"):
        """Get the exported data as a JSON string.

        Parameters:
        ----------
        nan : bool, optional
            Whether to include NaN values in the JSON string. Default is False.
        id_label : str, optional
            The ID label for the data. Default is "corems:".

        """

        import json

        dict_data_list = self.get_list_dict_data(self.gcms)

        return json.dumps(
            dict_data_list, sort_keys=False, indent=4, separators=(",", ": ")
        )

    def to_pandas(self, write_metadata=True, id_label="corems:"):
        """Export the data to a Pandas DataFrame and save it as a pickle file.

        Parameters:
        ----------
        write_metadata : bool, optional
            Whether to write metadata to the output file.
        id_label : str, optional
            The ID label for the data.
        """

        columns = self._init_columns()

        dict_data_list = self.get_list_dict_data(self.gcms)

        df = DataFrame(dict_data_list, columns=columns)

        df.to_pickle(self.output_file.with_suffix(".pkl"))

        if write_metadata:
            self.write_settings(
                self.output_file.with_suffix(".pkl"), self.gcms, id_label="corems:"
            )

    def to_excel(self, write_mode="a", write_metadata=True, id_label="corems:"):
        """Export the data to an Excel file.

        Parameters:
        ----------
        write_mode : str, optional
            The write mode for the Excel file. Default is 'a' (append).
        write_metadata : bool, optional
            Whether to write metadata to the output file. Default is True.
        id_label : str, optional
            The ID label for the data. Default is "corems:".
        """

        out_put_path = self.output_file.with_suffix(".xlsx")

        columns = self._init_columns()

        dict_data_list = self.get_list_dict_data(self.gcms)

        df = DataFrame(dict_data_list, columns=columns)

        if write_mode == "a" and out_put_path.exists():
            writer = ExcelWriter(out_put_path, engine="openpyxl")
            # try to open an existing workbook
            writer.book = load_workbook(out_put_path)
            # copy existing sheets
            writer.sheets = dict((ws.title, ws) for ws in writer.book.worksheets)
            # read existing file
            reader = read_excel(out_put_path)
            # write out the new sheet
            df.to_excel(writer, index=False, header=False, startrow=len(reader) + 1)

            writer.close()
        else:
            df.to_excel(
                self.output_file.with_suffix(".xlsx"), index=False, engine="openpyxl"
            )

        if write_metadata:
            self.write_settings(out_put_path, self.gcms, id_label=id_label)

    def to_csv(
        self,
        separate_output=False,
        write_mode="w",
        write_metadata=True,
        id_label="corems:",
    ):
        """Export the data to a CSV file.

        Parameters:
        ----------
        separate_output : bool, optional
            Whether to separate the output into multiple files. Default is False.
        write_mode : str, optional
            The write mode for the CSV file. Default is 'w' (write).
        write_metadata : bool, optional
            Whether to write metadata to the output file. Default is True.
        id_label : str, optional
            The ID label for the data. Default is "corems:".
        """

        if separate_output:
            # set write mode to write
            # this mode will overwrite the file without warning
            write_mode = "w"
        else:
            # set write mode to append
            write_mode = "a"

        columns = self._init_columns()

        dict_data_list = self.get_list_dict_data(self.gcms)

        out_put_path = self.output_file.with_suffix(".csv")

        write_header = not out_put_path.exists()

        try:
            with open(out_put_path, write_mode, newline="") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                if write_header:
                    writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)

            if write_metadata:
                self.write_settings(out_put_path, self.gcms, id_label=id_label)

        except IOError as ioerror:
            print(ioerror)

    def to_hdf(self, id_label="corems:"):
        """Export the data to an HDF5 file.

        Parameters:
        ----------
        id_label : str, optional
            The ID label for the data. Default is "corems:".
        """

        # save sample at a time
        def add_compound(gc_peak, compound_obj):
            modifier = compound_obj.classify if compound_obj.classify else ""
            compound_group = compound_obj.name.replace("/", "") + " " + modifier

            if compound_group not in peak_group:
                compound_group = peak_group.create_group(compound_group)

                # compound_group.attrs["retention_time"] = compound_obj.retention_time
                compound_group.attrs["retention_index"] = compound_obj.ri
                compound_group.attrs["retention_index_score"] = compound_obj.ri_score
                compound_group.attrs["spectral_similarity_score"] = (
                    compound_obj.spectral_similarity_score
                )
                compound_group.attrs["similarity_score"] = compound_obj.similarity_score

                compond_mz = compound_group.create_dataset(
                    "mz", data=np.array(compound_obj.mz), dtype="f8"
                )
                compond_abundance = compound_group.create_dataset(
                    "abundance", data=np.array(compound_obj.abundance), dtype="f8"
                )

                if self.gcms.molecular_search_settings.exploratory_mode:
                    compound_group.attrs["Spectral Similarities"] = json.dumps(
                        compound_obj.spectral_similarity_scores,
                        sort_keys=False,
                        indent=4,
                        separators=(",", ":"),
                    )
            else:
                warnings.warn("Skipping duplicate reference compound.")

        import json
        from datetime import datetime, timezone

        import h5py
        import numpy as np

        output_path = self.output_file.with_suffix(".hdf5")

        with h5py.File(output_path, "w") as hdf_handle:
            timenow = str(datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z"))
            hdf_handle.attrs["time_stamp"] = timenow
            hdf_handle.attrs["data_structure"] = "gcms"
            hdf_handle.attrs["analyzer"] = self.gcms.analyzer
            hdf_handle.attrs["instrument_label"] = self.gcms.instrument_label

            hdf_handle.attrs["sample_id"] = "self.gcms.id"
            hdf_handle.attrs["sample_name"] = self.gcms.sample_name
            hdf_handle.attrs["input_data"] = str(self.gcms.file_location)
            hdf_handle.attrs["output_data"] = str(output_path)
            hdf_handle.attrs["output_data_id"] = id_label + uuid.uuid4().hex
            hdf_handle.attrs["corems_version"] = __version__

            hdf_handle.attrs["Stats"] = json.dumps(
                self.get_data_stats(self.gcms),
                sort_keys=False,
                indent=4,
                separators=(",", ": "),
            )
            hdf_handle.attrs["Calibration"] = json.dumps(
                self.get_calibration_stats(self.gcms, id_label),
                sort_keys=False,
                indent=4,
                separators=(",", ": "),
            )
            hdf_handle.attrs["Blank"] = json.dumps(
                self.get_blank_stats(self.gcms),
                sort_keys=False,
                indent=4,
                separators=(",", ": "),
            )

            corems_dict_setting = parameter_to_dict.get_dict_data_gcms(self.gcms)
            hdf_handle.attrs["CoreMSParameters"] = json.dumps(
                corems_dict_setting, sort_keys=False, indent=4, separators=(",", ": ")
            )

            scans_dataset = hdf_handle.create_dataset(
                "scans", data=np.array(self.gcms.scans_number), dtype="f8"
            )
            rt_dataset = hdf_handle.create_dataset(
                "rt", data=np.array(self.gcms.retention_time), dtype="f8"
            )
            tic_dataset = hdf_handle.create_dataset(
                "tic", data=np.array(self.gcms.tic), dtype="f8"
            )
            processed_tic_dataset = hdf_handle.create_dataset(
                "processed_tic", data=np.array(self.gcms.processed_tic), dtype="f8"
            )

            output_score_method = (
                self.gcms.molecular_search_settings.output_score_method
            )

            for gc_peak in self.gcms:
                # print(gc_peak.retention_time)
                # print(gc_peak.tic)

                # check if there is a compound candidate
                peak_group = hdf_handle.create_group(str(gc_peak.retention_time))
                peak_group.attrs["deconvolution"] = int(
                    self.gcms.chromatogram_settings.use_deconvolution
                )

                peak_group.attrs["start_scan"] = gc_peak.start_scan
                peak_group.attrs["apex_scan"] = gc_peak.apex_scan
                peak_group.attrs["final_scan"] = gc_peak.final_scan

                peak_group.attrs["retention_index"] = gc_peak.ri
                peak_group.attrs["retention_time"] = gc_peak.retention_time
                peak_group.attrs["area"] = gc_peak.area

                mz = peak_group.create_dataset(
                    "mz", data=np.array(gc_peak.mass_spectrum.mz_exp), dtype="f8"
                )
                abundance = peak_group.create_dataset(
                    "abundance",
                    data=np.array(gc_peak.mass_spectrum.abundance),
                    dtype="f8",
                )

                if gc_peak:
                    if output_score_method == "highest_sim_score":
                        compound_obj = gc_peak.highest_score_compound
                        add_compound(gc_peak, compound_obj)

                    elif output_score_method == "highest_ss":
                        compound_obj = gc_peak.highest_ss_compound
                        add_compound(gc_peak, compound_obj)

                    else:
                        for compound_obj in gc_peak:
                            add_compound(gc_peak, compound_obj)

    def get_data_stats(self, gcms):
        """Get statistics about the GCMS data.

        Parameters:
        ----------
        gcms : object
            The low resolution GCMS object.

        Returns:
        -------
        dict
            A dictionary containing the data statistics.
        """

        matched_peaks = gcms.matched_peaks
        no_matched_peaks = gcms.no_matched_peaks
        unique_metabolites = gcms.unique_metabolites

        peak_matchs_above_0p85 = 0
        unique_peak_match_above_0p85 = 0
        for match_peak in matched_peaks:
            gc_peak_above_85 = 0
            matches_above_85 = list(
                filter(lambda m: m.similarity_score >= 0.85, match_peak)
            )
            if matches_above_85:
                peak_matchs_above_0p85 += 1
            if len(matches_above_85) == 1:
                unique_peak_match_above_0p85 += 1

        data_stats = {}
        data_stats["average_signal_noise"] = "ni"
        data_stats["chromatogram_dynamic_range"] = gcms.dynamic_range
        data_stats["total_number_peaks"] = len(gcms)
        data_stats["total_peaks_matched"] = len(matched_peaks)
        data_stats["total_peaks_without_matches"] = len(no_matched_peaks)
        data_stats["total_matches_above_similarity_score_0.85"] = peak_matchs_above_0p85
        data_stats["single_matches_above_similarity_score_0.85"] = (
            unique_peak_match_above_0p85
        )
        data_stats["unique_metabolites"] = len(unique_metabolites)

        return data_stats

    def get_calibration_stats(self, gcms, id_label):
        """Get statistics about the GC-MS calibration.

        Parameters:
        ----------
        """
        calibration_parameters = {}

        calibration_parameters["calibration_rt_ri_pairs_ref"] = gcms.ri_pairs_ref
        calibration_parameters["data_url"] = str(gcms.cal_file_path)
        calibration_parameters["has_input"] = id_label + corems_md5(gcms.cal_file_path)
        calibration_parameters["data_name"] = str(gcms.cal_file_path.stem)
        calibration_parameters["calibration_method"] = ""

        return calibration_parameters

    def get_blank_stats(self, gcms):
        """Get statistics about the GC-MS blank."""
        blank_parameters = {}

        blank_parameters["data_name"] = "ni"
        blank_parameters["blank_id"] = "ni"
        blank_parameters["data_url"] = "ni"
        blank_parameters["has_input"] = "ni"
        blank_parameters["common_features_to_blank"] = "ni"

        return blank_parameters

    def get_instrument_metadata(self, gcms):
        """Get metadata about the GC-MS instrument."""
        instrument_metadata = {}

        instrument_metadata["analyzer"] = gcms.analyzer
        instrument_metadata["instrument_label"] = gcms.instrument_label
        instrument_metadata["instrument_id"] = uuid.uuid4().hex

        return instrument_metadata

    def get_data_metadata(self, gcms, id_label, output_path):
        """Get metadata about the GC-MS data.

        Parameters:
        ----------
        gcms : object
            The low resolution GCMS object.
        id_label : str
            The ID label for the data.
        output_path : str
            The output file path.

        Returns:
        -------
        dict
            A dictionary containing the data metadata.
        """
        if isinstance(output_path, str):
            output_path = Path(output_path)

        paramaters_path = output_path.with_suffix(".json")

        if paramaters_path.exists():
            with paramaters_path.open() as current_param:
                metadata = json.load(current_param)
                data_metadata = metadata.get("Data")
        else:
            data_metadata = {}
            data_metadata["data_name"] = []
            data_metadata["input_data_url"] = []
            data_metadata["has_input"] = []

        data_metadata["data_name"].append(gcms.sample_name)
        data_metadata["input_data_url"].append(str(gcms.file_location))
        data_metadata["has_input"].append(id_label + corems_md5(gcms.file_location))

        data_metadata["output_data_name"] = str(output_path.stem)
        data_metadata["output_data_url"] = str(output_path)
        data_metadata["has_output"] = id_label + corems_md5(output_path)

        return data_metadata

    def get_parameters_json(self, gcms, id_label, output_path):
        """Get the parameters as a JSON string.

        Parameters:
        ----------
        gcms : GCMS object
            The low resolution GCMS object.
        id_label : str
            The ID label for the data.
        output_path : str
            The output file path.

        Returns:
        -------
        str
            The parameters as a JSON string.
        """

        output_parameters_dict = {}
        output_parameters_dict["Data"] = self.get_data_metadata(
            gcms, id_label, output_path
        )
        output_parameters_dict["Stats"] = self.get_data_stats(gcms)
        output_parameters_dict["Calibration"] = self.get_calibration_stats(
            gcms, id_label
        )
        output_parameters_dict["Blank"] = self.get_blank_stats(gcms)
        output_parameters_dict["Instrument"] = self.get_instrument_metadata(gcms)
        corems_dict_setting = parameter_to_dict.get_dict_data_gcms(gcms)
        corems_dict_setting["corems_version"] = __version__
        output_parameters_dict["CoreMSParameters"] = corems_dict_setting
        output_parameters_dict["has_metabolite"] = gcms.metabolites_data
        output = json.dumps(
            output_parameters_dict, sort_keys=False, indent=4, separators=(",", ": ")
        )

        return output

    def write_settings(self, output_path, gcms, id_label="emsl:"):
        """Write the settings to a JSON file.

        Parameters:
        ----------
        output_path : str
            The output file path.
        gcms : GCMS object
            The low resolution GCMS object.
        id_label : str
            The ID label for the data. Default is "emsl:".

        """

        output = self.get_parameters_json(gcms, id_label, output_path)

        with open(
            output_path.with_suffix(".json"),
            "w",
            encoding="utf8",
        ) as outfile:
            outfile.write(output)

    def get_list_dict_data(self, gcms, include_no_match=True, no_match_inline=False):
        """Get the exported data as a list of dictionaries.

        Parameters:
        ----------
        gcms : object
            The low resolution GCMS object.
        include_no_match : bool, optional
            Whether to include no match data. Default is True.
        no_match_inline : bool, optional
            Whether to include no match data inline. Default is False.

        Returns:
        -------
        list
            The exported data as a list of dictionaries.
        """

        output_score_method = gcms.molecular_search_settings.output_score_method

        dict_data_list = []

        def add_match_dict_data():
            derivatization = "{}:{}:{}".format(
                compound_obj.classify,
                compound_obj.derivativenum,
                compound_obj.derivatization,
            )
            out_dict = {
                "Sample name": gcms.sample_name,
                "Peak Index": gcpeak_index,
                "Retention Time": gc_peak.retention_time,
                "Retention Time Ref": compound_obj.retention_time,
                "Peak Height": gc_peak.tic,
                "Peak Area": gc_peak.area,
                "Retention index": gc_peak.ri,
                "Retention index Ref": compound_obj.ri,
                "Retention Index Score": compound_obj.ri_score,
                "Spectral Similarity Score": compound_obj.spectral_similarity_score,
                "Similarity Score": compound_obj.similarity_score,
                "Compound Name": compound_obj.name,
                "Chebi ID": compound_obj.metadata.chebi,
                "Kegg Compound ID": compound_obj.metadata.kegg,
                "Inchi": compound_obj.metadata.inchi,
                "Inchi Key": compound_obj.metadata.inchikey,
                "Smiles": compound_obj.metadata.smiles,
                "Molecular Formula": compound_obj.formula,
                "IUPAC Name": compound_obj.metadata.iupac_name,
                "Traditional Name": compound_obj.metadata.traditional_name,
                "Common Name": compound_obj.metadata.common_name,
                "Derivatization": derivatization,
            }

            if self.gcms.molecular_search_settings.exploratory_mode:
                out_dict.update(
                    {
                        "Weighted Cosine Correlation": compound_obj.spectral_similarity_scores.get(
                            "weighted_cosine_correlation"
                        ),
                        "Cosine Correlation": compound_obj.spectral_similarity_scores.get(
                            "cosine_correlation"
                        ),
                        "Stein Scott Similarity": compound_obj.spectral_similarity_scores.get(
                            "stein_scott_similarity"
                        ),
                        "Pearson Correlation": compound_obj.spectral_similarity_scores.get(
                            "pearson_correlation"
                        ),
                        "Spearman Correlation": compound_obj.spectral_similarity_scores.get(
                            "spearman_correlation"
                        ),
                        "Kendall Tau Correlation": compound_obj.spectral_similarity_scores.get(
                            "kendall_tau_correlation"
                        ),
                        "DFT Correlation": compound_obj.spectral_similarity_scores.get(
                            "dft_correlation"
                        ),
                        "DWT Correlation": compound_obj.spectral_similarity_scores.get(
                            "dwt_correlation"
                        ),
                        "Euclidean Distance": compound_obj.spectral_similarity_scores.get(
                            "euclidean_distance"
                        ),
                        "Manhattan Distance": compound_obj.spectral_similarity_scores.get(
                            "manhattan_distance"
                        ),
                        "Jaccard Distance": compound_obj.spectral_similarity_scores.get(
                            "jaccard_distance"
                        ),
                    }
                )
                for method in methods_name:
                    out_dict[methods_name.get(method)] = (
                        compound_obj.spectral_similarity_scores.get(method)
                    )

            dict_data_list.append(out_dict)

        def add_no_match_dict_data():
            dict_data_list.append(
                {
                    "Sample name": gcms.sample_name,
                    "Peak Index": gcpeak_index,
                    "Retention Time": gc_peak.retention_time,
                    "Peak Height": gc_peak.tic,
                    "Peak Area": gc_peak.area,
                    "Retention index": gc_peak.ri,
                }
            )

        for gcpeak_index, gc_peak in enumerate(gcms.sorted_gcpeaks):
            # check if there is a compound candidate
            if gc_peak:
                if output_score_method == "highest_sim_score":
                    compound_obj = gc_peak.highest_score_compound
                    add_match_dict_data()

                elif output_score_method == "highest_ss":
                    compound_obj = gc_peak.highest_ss_compound
                    add_match_dict_data()

                else:
                    for compound_obj in gc_peak:
                        add_match_dict_data()  # add monoisotopic peak

            else:
                # include not_match
                if include_no_match and no_match_inline:
                    add_no_match_dict_data()

        if include_no_match and not no_match_inline:
            for gcpeak_index, gc_peak in enumerate(gcms.sorted_gcpeaks):
                if not gc_peak:
                    add_no_match_dict_data()

        return dict_data_list


class HighResMassSpectraExport(HighResMassSpecExport):
    """A class to export high resolution mass spectra data.

    This class provides methods to export high resolution mass spectra data to various formats such as Excel, CSV, HDF5, and Pandas DataFrame.

    Parameters
    ----------
    out_file_path : str | Path
        The output file path.
    mass_spectra : object
        The high resolution mass spectra object.
    output_type : str, optional
        The output type. Default is 'excel'.

    Attributes
    ----------
    output_file : Path
        The output file path without suffix
    dir_loc : Path
        The directory location for the output file, by default this will be the output_file + ".corems" and all output files will be written into this location
    mass_spectra : MassSpectraBase
        The high resolution mass spectra object.
    """

    def __init__(self, out_file_path, mass_spectra, output_type="excel"):
        super().__init__(
            out_file_path=out_file_path, mass_spectrum=None, output_type=output_type
        )

        self.dir_loc = Path(out_file_path + ".corems")
        self.dir_loc.mkdir(exist_ok=True)
        # Place the output file in the directory
        self.output_file = self.dir_loc / Path(out_file_path).name
        self._output_type = output_type  # 'excel', 'csv', 'pandas' or 'hdf5'
        self.mass_spectra = mass_spectra
        self._init_columns()

    def get_pandas_df(self):
        """Get the mass spectra as a list of Pandas DataFrames."""

        list_df = []

        for mass_spectrum in self.mass_spectra:
            columns = self.columns_label + self.get_all_used_atoms_in_order(
                mass_spectrum
            )

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number

            df.name = str(self.output_file) + "_" + str(scan_number)

            list_df.append(df)

        return list_df

    def to_pandas(self, write_metadata=True):
        """Export the data to a Pandas DataFrame and save it as a pickle file.

        Parameters:
        ----------
        write_metadata : bool, optional
            Whether to write metadata to the output file. Default is True.
        """

        for mass_spectrum in self.mass_spectra:
            columns = self.columns_label + self.get_all_used_atoms_in_order(
                mass_spectrum
            )

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number

            out_filename = Path(
                "%s_scan%s%s" % (self.output_file, str(scan_number), ".pkl")
            )

            df.to_pickle(self.dir_loc / out_filename)

            if write_metadata:
                self.write_settings(
                    self.dir_loc / out_filename.with_suffix(""), mass_spectrum
                )

    def to_excel(self, write_metadata=True):
        """Export the data to an Excel file.

        Parameters:
        ----------
        write_metadata : bool, optional
            Whether to write metadata to the output file. Default is True.
        """
        for mass_spectrum in self.mass_spectra:
            columns = self.columns_label + self.get_all_used_atoms_in_order(
                mass_spectrum
            )

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            df = DataFrame(dict_data_list, columns=columns)

            scan_number = mass_spectrum.scan_number

            out_filename = Path(
                "%s_scan%s%s" % (self.output_file, str(scan_number), ".xlsx")
            )

            df.to_excel(self.dir_loc / out_filename)

            if write_metadata:
                self.write_settings(
                    self.dir_loc / out_filename.with_suffix(""), mass_spectrum
                )

    def to_csv(self, write_metadata=True):
        """Export the data to a CSV file.

        Parameters:
        ----------
        write_metadata : bool, optional
            Whether to write metadata to the output file. Default is True.
        """
        import csv

        for mass_spectrum in self.mass_spectra:
            columns = self.columns_label + self.get_all_used_atoms_in_order(
                mass_spectrum
            )

            scan_number = mass_spectrum.scan_number

            dict_data_list = self.get_list_dict_data(mass_spectrum)

            out_filename = Path(
                "%s_scan%s%s" % (self.output_file, str(scan_number), ".csv")
            )

            with open(self.dir_loc / out_filename, "w", newline="") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=columns)
                writer.writeheader()
                for data in dict_data_list:
                    writer.writerow(data)

            if write_metadata:
                self.write_settings(
                    self.dir_loc / out_filename.with_suffix(""), mass_spectrum
                )

    def get_mass_spectra_attrs(self, mass_spectra):
        """Get the mass spectra attributes as a JSON string.

        Parameters:
        ----------
        mass_spectra : object
            The high resolution mass spectra object.

        Returns:
        -------
        str
            The mass spectra attributes as a JSON string.
        """
        dict_ms_attrs = {}
        dict_ms_attrs["analyzer"] = self.mass_spectra.analyzer
        dict_ms_attrs["instrument_label"] = self.mass_spectra.instrument_label
        dict_ms_attrs["sample_name"] = self.mass_spectra.sample_name

        return json.dumps(
            dict_ms_attrs, sort_keys=False, indent=4, separators=(",", ": ")
        )

    def to_hdf(self, overwrite=False):
        """Export the data to an HDF5 file."""
        if overwrite:
            if self.output_file.with_suffix(".hdf5").exists():
                self.output_file.with_suffix(".hdf5").unlink()

        with h5py.File(self.output_file.with_suffix(".hdf5"), "a") as hdf_handle:
            if not hdf_handle.attrs.get("date_utc"):
                # Set metadata
                timenow = str(
                    datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z")
                )
                hdf_handle.attrs["date_utc"] = timenow
                hdf_handle.attrs["filename"] = self.mass_spectra.file_location.name
                hdf_handle.attrs["data_structure"] = "mass_spectra"
                hdf_handle.attrs["analyzer"] = self.mass_spectra.analyzer
                hdf_handle.attrs["instrument_label"] = (
                    self.mass_spectra.instrument_label
                )
                hdf_handle.attrs["sample_name"] = self.mass_spectra.sample_name
                hdf_handle.attrs["polarity"] = self.mass_spectra.polarity

            if "mass_spectra" not in hdf_handle:
                mass_spectra_group = hdf_handle.create_group("mass_spectra")
            else:
                mass_spectra_group = hdf_handle.get("mass_spectra")

            for mass_spectrum in self.mass_spectra:
                list_results = self.list_dict_to_list(mass_spectrum, is_hdf5=True)
                dict_ms_attrs = self.get_mass_spec_attrs(mass_spectrum)
                setting_dicts = parameter_to_dict.get_dict_data_ms(mass_spectrum)
                columns_labels = json.dumps(
                    self.columns_label + self.get_all_used_atoms_in_order(mass_spectrum)
                )
                group_key = str(int(mass_spectrum.scan_number))
                if group_key not in mass_spectra_group.keys():
                    scan_group = mass_spectra_group.create_group(
                        str(int(mass_spectrum.scan_number))
                    )
                    if list(mass_spectrum.abundance_profile):
                        mz_abun_array = concatenate(
                            (
                                mass_spectrum.abundance_profile,
                                mass_spectrum.mz_exp_profile,
                            ),
                            axis=0,
                        )
                        raw_ms_dataset = scan_group.create_dataset(
                            "raw_ms", data=mz_abun_array, dtype="f8"
                        )
                    else:
                        # create empy dataset for missing raw data
                        raw_ms_dataset = scan_group.create_dataset("raw_ms", dtype="f8")
                    raw_ms_dataset.attrs["MassSpecAttrs"] = json.dumps(
                        dict_ms_attrs, sort_keys=False, indent=4, separators=(",", ": ")
                    )
                    if isinstance(mass_spectrum, MassSpecfromFreq):
                        raw_ms_dataset.attrs["TransientSetting"] = json.dumps(
                            setting_dicts.get("TransientSetting"),
                            sort_keys=False,
                            indent=4,
                            separators=(",", ": "),
                        )
                else:
                    scan_group = mass_spectra_group.get(group_key)
                    # if there is not processed data len = 0, otherwise len() will return next index
                index_processed_data = str(len(scan_group.keys()))
                timenow = str(
                    datetime.now(timezone.utc).strftime("%d/%m/%Y %H:%M:%S %Z")
                )
                processed_dset = scan_group.create_dataset(
                    index_processed_data, data=list_results
                )
                processed_dset.attrs["date_utc"] = timenow
                processed_dset.attrs["ColumnsLabels"] = columns_labels
                processed_dset.attrs["MoleculaSearchSetting"] = json.dumps(
                    setting_dicts.get("MoleculaSearch"),
                    sort_keys=False,
                    indent=4,
                    separators=(",", ": "),
                )
                processed_dset.attrs["MassSpecPeakSetting"] = json.dumps(
                    setting_dicts.get("MassSpecPeak"),
                    sort_keys=False,
                    indent=4,
                    separators=(",", ": "),
                )
                processed_dset.attrs["MassSpectrumSetting"] = json.dumps(
                    setting_dicts.get("MassSpectrum"),
                    sort_keys=False,
                    indent=4,
                    separators=(",", ": "),
                )


class LCMSExport(HighResMassSpectraExport):
    """A class to export high resolution LC-MS data.

    This class provides methods to export high resolution LC-MS data to HDF5.

    Parameters
    ----------
    out_file_path : str | Path
        The output file path, do not include the file extension.
    lcms_object : LCMSBase
        The high resolution lc-ms object.
    """

    def __init__(self, out_file_path, mass_spectra):
        super().__init__(out_file_path, mass_spectra, output_type="hdf5")

    def to_hdf(self, overwrite=False, save_parameters=True, parameter_format="toml"):
        """Export the data to an HDF5.

        Parameters
        ----------
        overwrite : bool, optional
            Whether to overwrite the output file. Default is False.
        save_parameters : bool, optional
            Whether to save the parameters as a separate json or toml file. Default is True.
        parameter_format : str, optional
            The format to save the parameters in. Default is 'toml'.

        Raises
        ------
        ValueError
            If parameter_format is not 'json' or 'toml'.
        """
        super().to_hdf(overwrite=overwrite)
        with h5py.File(self.output_file.with_suffix(".hdf5"), "a") as hdf_handle:
            # Add scan_info to hdf5 file
            if "scan_info" not in hdf_handle:
                scan_info_group = hdf_handle.create_group("scan_info")
                for k, v in self.mass_spectra._scan_info.items():
                    array = np.array(list(v.values()))
                    if array.dtype.str[0:2] == "<U":
                        array = array.astype("S")
                    scan_info_group.create_dataset(k, data=array)

            # Add ms_unprocessed to hdf5 file
            if self.mass_spectra._ms_unprocessed:
                if "ms_unprocessed" not in hdf_handle:
                    ms_unprocessed_group = hdf_handle.create_group("ms_unprocessed")
                else:
                    ms_unprocessed_group = hdf_handle.get("ms_unprocessed")
                for k, v in self.mass_spectra._ms_unprocessed.items():
                    array = np.array(v)
                    ms_unprocessed_group.create_dataset(str(k), data=array)

            # Add LCMS mass features to hdf5 file
            if len(self.mass_spectra.mass_features) > 0:
                if "mass_features" not in hdf_handle:
                    mass_features_group = hdf_handle.create_group("mass_features")
                else:
                    mass_features_group = hdf_handle.get("mass_features")

                # Create group for each mass feature, with key as the mass feature id
                for k, v in self.mass_spectra.mass_features.items():
                    mass_features_group.create_group(str(k))
                    # Loop through each of the mass feature attributes and add them as attributes (if single value) or datasets (if array)
                    for k2, v2 in v.__dict__.items():
                        if v2 is not None:
                            # Check if the attribute is an integer or float and set as an attribute in the mass feature group
                            if k2 not in [
                                "chromatogram_parent",
                                "ms2_mass_spectra",
                                "mass_spectrum",
                                "_eic_data",
                                "ms2_similarity_results",
                            ]:
                                if k2 == "ms2_scan_numbers":
                                    array = np.array(v2)
                                    mass_features_group[str(k)].create_dataset(
                                        str(k2), data=array
                                    )
                                elif (
                                    isinstance(v2, int)
                                    or isinstance(v2, float)
                                    or isinstance(v2, str)
                                    or isinstance(v2, np.integer)
                                ):
                                    mass_features_group[str(k)].attrs[str(k2)] = v2
                                else:
                                    raise TypeError(
                                        f"Attribute {k2} is not an integer, float, or string and cannot be added to the hdf5 file"
                                    )

            # Add EIC data to hdf5 file
            if len(self.mass_spectra.eics) > 0:
                if "eics" not in hdf_handle:
                    eic_group = hdf_handle.create_group("eics")
                else:
                    eic_group = hdf_handle.get("eics")

                # Create group for each eic
                for k, v in self.mass_spectra.eics.items():
                    eic_group.create_group(str(k))
                    eic_group[str(k)].attrs["mz"] = k
                    # Loop through each of the attributes and add them as datasets (if array)
                    for k2, v2 in v.__dict__.items():
                        if v2 is not None:
                            array = np.array(v2)
                            eic_group[str(k)].create_dataset(str(k2), data=array)

            # Add ms2_search results to hdf5 file
            if len(self.mass_spectra.spectral_search_results) > 0:
                if "spectral_search_results" not in hdf_handle:
                    spectral_search_results = hdf_handle.create_group(
                        "spectral_search_results"
                    )
                else:
                    spectral_search_results = hdf_handle.get("spectral_search_results")
                # Create group for each search result by ms2_scan / precursor_mz
                for k, v in self.mass_spectra.spectral_search_results.items():
                    spectral_search_results.create_group(str(k))
                    for k2, v2 in v.items():
                        spectral_search_results[str(k)].create_group(str(k2))
                        spectral_search_results[str(k)][str(k2)].attrs[
                            "precursor_mz"
                        ] = v2.precursor_mz
                        spectral_search_results[str(k)][str(k2)].attrs[
                            "ms2_spectrum_id"
                        ] = v2.ms2_spectrum_id
                        # Loop through each of the attributes and add them as datasets (if array)
                        for k3, v3 in v2.__dict__.items():
                            if v3 is not None and k3 not in [
                                "ms2_spectrum",
                                "precursor_mz",
                                "ms2_spectrum_id",
                            ]:
                                if k3 == "query_frags" or k3 == "lib_frags":
                                    v3 = [", ".join(x) for x in v3]
                                array = np.array(v3)
                                if array.dtype.str[0:2] == "<U":
                                    array = array.astype("S")
                                spectral_search_results[str(k)][str(k2)].create_dataset(
                                    str(k3), data=array
                                )

        # Save parameters as separate json
        if save_parameters:
            # Check if parameter_format is valid
            if parameter_format not in ["json", "toml"]:
                raise ValueError("parameter_format must be 'json' or 'toml'")

            if parameter_format == "json":
                dump_lcms_settings_json(
                    filename=self.output_file.with_suffix(".json"),
                    lcms_obj=self.mass_spectra,
                )
            elif parameter_format == "toml":
                dump_lcms_settings_toml(
                    filename=self.output_file.with_suffix(".toml"),
                    lcms_obj=self.mass_spectra,
                )
