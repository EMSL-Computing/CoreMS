from math import exp
from threading import Thread

from numpy import power

from corems.molecular_id.calc.SpectralSimilarity import SpectralSimilarity
from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite


class LowResMassSpectralMatch(Thread):
    """A class representing a low-resolution mass spectral match.

    Parameters
    -----------
    gcms_obj : object
        The GC-MS object.
    sql_obj : object, optional
        The SQL object for database operations. Default is None.
    calibration : bool, optional
        Flag indicating if the match is for calibration. Default is False.

    Attributes
    -----------
    gcms_obj : object
        The GC-MS object.
    sql_obj : object
        The SQL object for database operations.
    calibration : bool
        Flag indicating if the match is for calibration.

    Methods
    --------
    * metabolite_detector_score(gc_peak, ref_obj, spectral_simi).
        Calculates the spectral similarity scores and the similarity score for a given GC peak and reference object.
    * run().
        Runs the low-resolution mass spectral match.

    """

    def __init__(self, gcms_obj, sql_obj=None, calibration=False):
        Thread.__init__(self)

        self.gcms_obj = gcms_obj

        #  initiated at create_molecular_database()
        # self.dict_molecular_lookup_table = None
        self.calibration = calibration
        # reading local file for now,
        if not sql_obj:
            self.sql_obj = EI_LowRes_SQLite(
                url=self.gcms_obj.molecular_search_settings.url_database
            )
        else:
            self.sql_obj = sql_obj

    def metabolite_detector_score(self, gc_peak, ref_obj, spectral_simi):
        """
        Calculates the spectral similarity scores and the similarity score for a given GC peak and reference object.

        Parameters
        -----------
        gc_peak : object
            The GC peak object.
        ref_obj : object
            The reference object.
        spectral_simi : object
            The spectral similarity object.

        Returns
        --------
        tuple
            A tuple containing the spectral similarity scores, RI score, and similarity score.

        """
        spectral_similarity_scores = {}
        spectral_similarity_scores["cosine_correlation"] = (
            spectral_simi.cosine_correlation()
        )

        if self.gcms_obj.molecular_search_settings.exploratory_mode:
            spectral_similarity_scores["weighted_cosine_correlation"] = (
                spectral_simi.weighted_cosine_correlation()
            )
            ss, ss_nist = spectral_simi.stein_scott()
            spectral_similarity_scores["stein_scott_similarity"] = ss
            spectral_similarity_scores["stein_scott_similarity_nist"] = ss_nist

            spectral_similarity_scores["pearson_correlation"] = (
                spectral_simi.pearson_correlation()
            )
            spectral_similarity_scores["spearman_correlation"] = (
                spectral_simi.spearman_correlation()
            )
            spectral_similarity_scores["kendall_tau_correlation"] = (
                spectral_simi.kendall_tau()
            )
            spectral_similarity_scores["euclidean_distance"] = (
                spectral_simi.euclidean_distance()
            )
            spectral_similarity_scores["manhattan_distance"] = (
                spectral_simi.manhattan_distance()
            )
            spectral_similarity_scores["jaccard_distance"] = (
                spectral_simi.jaccard_distance()
            )
            spectral_similarity_scores["dft_correlation"] = (
                spectral_simi.dft_correlation()
            )
            spectral_similarity_scores["dwt_correlation"] = (
                spectral_simi.dwt_correlation()
            )
            spectral_similarity_scores.update(spectral_simi.extra_distances())
            # print(spectral_similarity_scores)
        # print(ref_obj.get('ri'), gc_peak.ri, self.gcms_obj.molecular_search_settings.ri_window)

        ri_score = exp(
            -1
            * (
                power((gc_peak.ri - ref_obj.get("ri")), 2)
                / (2 * power(self.gcms_obj.molecular_search_settings.ri_std, 2))
            )
        )

        similarity_score = (
            (spectral_similarity_scores.get("cosine_correlation") ** 2) * (ri_score)
        ) ** (1 / 3)

        return spectral_similarity_scores, ri_score, similarity_score

    # @timeit
    def run(self):
        """Runs the low-resolution mass spectral match."""
        # TODO select the best gcms peak
        import tqdm

        original_use_deconvolution = (
            self.gcms_obj.chromatogram_settings.use_deconvolution
        )

        if not self.gcms_obj:
            # Do not use deconvolution for the retention index calibration

            if self.calibration:
                self.gcms_obj.chromatogram_settings.use_deconvolution = False

            self.gcms_obj.process_chromatogram()

        self.gcms_obj.chromatogram_settings.use_deconvolution = (
            original_use_deconvolution
        )

        for gc_peak in tqdm.tqdm(self.gcms_obj):
            if not self.calibration:
                window = self.gcms_obj.molecular_search_settings.ri_search_range

                ri = gc_peak.ri

                min_mat_ri = (ri - window, ri + window)

                ref_objs = self.sql_obj.query_min_max_ri(min_mat_ri)

            else:
                compound_names = self.gcms_obj.molecular_search_settings.ri_calibration_compound_names

                window = self.gcms_obj.molecular_search_settings.rt_search_range

                rt = gc_peak.retention_time

                min_mat_rt = (rt - window, rt + window)

                ref_objs = self.sql_obj.query_names_and_rt(min_mat_rt, compound_names)

            for ref_obj in ref_objs:
                # uses spectral similarly and uses a threshold to only select peaks with high data correlation

                spectral_simi = SpectralSimilarity(
                    gc_peak.mass_spectrum.mz_abun_dict, ref_obj
                )

                if self.calibration:
                    spectral_similarity_scores = {}
                    spectral_similarity_scores["cosine_correlation"] = (
                        spectral_simi.cosine_correlation()
                    )

                    # print(w_correlation_value,correlation_value )
                    if (
                        spectral_similarity_scores["cosine_correlation"]
                        >= self.gcms_obj.molecular_search_settings.correlation_threshold
                    ):
                        gc_peak.add_compound(ref_obj, spectral_similarity_scores)

                # use score, usually a combination of Retention index and Spectral Similarity
                # Threshold is implemented by not necessarily used
                else:
                    # m/q developed methods will be implemented here
                    spectral_similarity_scores, ri_score, similarity_score = (
                        self.metabolite_detector_score(gc_peak, ref_obj, spectral_simi)
                    )

                    # TODO need to add similarity score option in the parameters encapsulation class

                    if (
                        similarity_score
                        >= self.gcms_obj.molecular_search_settings.score_threshold
                    ):
                        gc_peak.add_compound(
                            ref_obj,
                            spectral_similarity_scores,
                            ri_score,
                            similarity_score,
                        )

        self.sql_obj.session.close()
        self.sql_obj.engine.dispose()
