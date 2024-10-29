import re

import numpy as np

from corems.molecular_id.factory.spectrum_search_results import SpectrumSearchResults


class LCMSSpectralSearch:
    """
    Methods for searching LCMS spectra.

    This class is designed to be a mixin class for the :obj:`~corems.mass_spectra.factory.lc_class.LCMSBase` class.

    """

    @staticmethod
    def get_more_match_quals(
        query_mz_arr, lib_entry, mz_tol_da=0.1, include_fragment_types=False
    ):
        """
        Return additional match qualities between query and library entry.

        Parameters
        ----------
        query_mz_arr : np.array
            Array of query spectrum. Shape (N, 2), with m/z in the first column
            and abundance in the second.
        lib_entry : dict
            Library spectrum entry, with 'mz' key containing the spectrum in
            the format (mz, abundance),(mz, abundance), i.e. from MetabRef.
        mz_tol_da : float, optional
            Tolerance in Da for matching peaks (in MS2). Default is 0.1.
        include_fragment_types : bool, optional
            If True, include fragment type comparisons in output.
            Defaults to False.

        Returns
        -------
        tuple
            Tuple of (query_in_lib, query_in_lib_fract, lib_in_query, lib_in_query_fract, query_frags, lib_frags, lib_precursor_mz).

        Notes
        -----
        query_in_lib : int
            Number of peaks in query that are present in the library entry (within mz_tol_da).
        query_in_lib_fract : float
            Fraction of peaks in query that are present in the library entry (within mz_tol_da).
        lib_in_query : int
            Number of peaks in the library entry that are present in the query (within mz_tol_da).
        lib_in_query_fract : float
            Fraction of peaks in the library entry that are present in the query (within mz_tol_da).
        query_frags : list
            List of unique fragment types present in the query, generally 'MLF' or 'LSF' or both.
        lib_frags : list
            List of unique fragment types present in the library entry, generally 'MLF' or 'LSF' or both.

        Raises
        ------
        ValueError
            If library entry does not have 'fragment_types' key and include_fragment_types is True.

        """

        # Get the original mz values from the library entry
        lib_mzs = np.array(
            re.findall(r"\(([^,]+),([^)]+)\)", lib_entry["mz"]), dtype=float
        ).reshape(-1, 2)[:, 0]

        # Get count and fraction of peaks in query that are in lib entry
        query_in_lib = 0
        for peak in query_mz_arr:
            if np.any(np.isclose(lib_mzs, peak, atol=mz_tol_da)):
                query_in_lib += 1
        query_in_lib_fract = query_in_lib / len(query_mz_arr)

        # Get count and fraction of peaks in lib that are in query
        lib_in_query = 0
        for peak in lib_mzs:
            if np.any(np.isclose(query_mz_arr, peak, atol=mz_tol_da)):
                lib_in_query += 1
        lib_in_query_fract = lib_in_query / len(lib_mzs)

        if include_fragment_types:
            # Check that fragment types are present in the library entry
            if "fragment_types" not in lib_entry.keys():
                raise ValueError(
                    "Flash entropy library entry must have 'fragment_types' key to include fragment types in output."
                )

            # Get types of fragments in the lib entry
            lib_frags = lib_entry["fragment_types"].split(", ")
            # make list of the fragment types that are present in the query spectrum
            lib_in_query_ids = list(
                set(
                    [
                        ind
                        for ind, x in enumerate(lib_mzs)
                        if len(np.where(np.isclose(query_mz_arr, x, atol=mz_tol_da))[0])
                        > 0
                    ]
                )
            )
            query_frags = list(set([lib_frags[x] for x in lib_in_query_ids]))
            lib_frags = list(set(lib_frags))

        else:
            query_frags = None
            lib_frags = None

        return (
            query_in_lib,
            query_in_lib_fract,
            lib_in_query,
            lib_in_query_fract,
            query_frags,
            lib_frags,
        )

    def fe_search(
        self,
        scan_list,
        fe_lib,
        precursor_mz_list=[],
        use_mass_features=True,
        peak_sep_da=0.01,
        get_additional_metrics=True,
    ):
        """
        Search LCMS spectra using a FlashEntropy approach.

        Parameters
        ----------
        scan_list : list
            List of scan numbers to search.
        fe_lib : :obj:`~ms_entropy.FlashEntropySearch`
            FlashEntropy Search instance.
        precursor_mz_list : list, optional
            List of precursor m/z values to search, by default [], which implies
            matched with mass features; to enable this use_mass_features must be True.
        use_mass_features : bool, optional
            If True, use mass features to get precursor m/z values, by default True.
            If True, will add search results to mass features' ms2_similarity_results attribute.
        peak_sep_da : float, optional
            Minimum separation between m/z peaks spectra in Da. This needs match the
            approximate resolution of the search spectra and the FlashEntropySearch
            instance, by default 0.01.
        get_additional_metrics : bool, optional
            If True, get additional metrics from FlashEntropy search, by default True.

        Returns
        -------
        None, but adds results to self.spectral_search_results and associates these
        spectral_search_results with mass_features within the self.mass_features dictionary.

        """
        # Retrieve parameters from self
        # include_fragment_types should used for lipids queries only, not general metabolomics
        include_fragment_types = self.parameters.lc_ms.include_fragment_types
        min_match_score = self.parameters.lc_ms.ms2_min_fe_score

        # If precursor_mz_list is empty and use_mass_features is True, get precursor m/z values from mass features for each scan in scan_list
        if use_mass_features and len(precursor_mz_list) == 0:
            precursor_mz_list = []
            for scan in scan_list:
                mf_ids = [
                    key
                    for key, value in self.mass_features.items()
                    if scan in value.ms2_mass_spectra
                ]
                precursor_mz = [
                    value.mz
                    for key, value in self.mass_features.items()
                    if key in mf_ids
                ]
                precursor_mz_list.append(precursor_mz)

        # Check that precursor_mz_list same length as scan_list, if not, raise error
        if len(precursor_mz_list) != len(scan_list):
            raise ValueError("Length of precursor_mz_list is not equal to scan_list.")

        # Loop through each query spectrum / precursor match and save ids of db spectrum that are decent matches
        overall_results_dict = {}
        for i in np.arange(len(scan_list)):
            scan_oi = scan_list[i]
            if len(self._ms[scan_oi].mspeaks) > 0:
                precursor_mzs = precursor_mz_list[i]
                overall_results_dict[scan_oi] = {}
                for precursor_mz in precursor_mzs:
                    query_spectrum = fe_lib.clean_spectrum_for_search(
                        precursor_mz=precursor_mz,
                        peaks=np.vstack(
                            (self._ms[scan_oi].mz_exp, self._ms[scan_oi].abundance)
                        ).T,
                        precursor_ions_removal_da=None,
                        noise_threshold=self._ms[
                            scan_oi
                        ].parameters.mass_spectrum.noise_threshold_min_relative_abundance
                        / 100,
                        min_ms2_difference_in_da=peak_sep_da,
                    )
                    search_results = fe_lib.search(
                        precursor_mz=precursor_mz,
                        peaks=query_spectrum,
                        ms1_tolerance_in_da=self.parameters.mass_spectrum[
                            "ms1"
                        ].molecular_search.max_ppm_error
                        * 10**-6
                        * precursor_mz,
                        ms2_tolerance_in_da=peak_sep_da * 0.5,
                        method={"identity"},
                        precursor_ions_removal_da=None,
                        noise_threshold=self._ms[
                            scan_oi
                        ].parameters.mass_spectrum.noise_threshold_min_relative_abundance
                        / 100,
                        target="cpu",
                    )["identity_search"]
                    match_inds = np.where(search_results > min_match_score)[0]

                    # If any decent matches are found, add them to the results dictionary
                    if len(match_inds) > 0:
                        match_scores = search_results[match_inds]
                        ref_ms_ids = [fe_lib[x]["id"] for x in match_inds]
                        ref_mol_ids = [
                            fe_lib[x]["molecular_data_id"] for x in match_inds
                        ]
                        ref_precursor_mzs = [
                            fe_lib[x]["precursor_mz"] for x in match_inds
                        ]
                        ion_types = [fe_lib[x]["ion_type"] for x in match_inds]
                        overall_results_dict[scan_oi][precursor_mz] = {
                            "ref_mol_id": ref_mol_ids,
                            "ref_ms_id": ref_ms_ids,
                            "ref_precursor_mz": ref_precursor_mzs,
                            "precursor_mz_error_ppm": [
                                (precursor_mz - x) / precursor_mz * 10**6
                                for x in ref_precursor_mzs
                            ],
                            "entropy_similarity": match_scores,
                            "ref_ion_type": ion_types,
                        }
                        if get_additional_metrics:
                            more_match_quals = [
                                self.get_more_match_quals(
                                    self._ms[scan_oi].mz_exp,
                                    fe_lib[x],
                                    mz_tol_da=peak_sep_da,
                                    include_fragment_types=include_fragment_types,
                                )
                                for x in match_inds
                            ]
                            overall_results_dict[scan_oi][precursor_mz].update(
                                {
                                    "query_mz_in_ref_n": [
                                        x[0] for x in more_match_quals
                                    ],
                                    "query_mz_in_ref_fract": [
                                        x[1] for x in more_match_quals
                                    ],
                                    "ref_mz_in_query_n": [
                                        x[2] for x in more_match_quals
                                    ],
                                    "ref_mz_in_query_fract": [
                                        x[3] for x in more_match_quals
                                    ],
                                }
                            )
                            if include_fragment_types:
                                overall_results_dict[scan_oi][precursor_mz].update(
                                    {
                                        "query_frag_types": [
                                            x[4] for x in more_match_quals
                                        ],
                                        "ref_frag_types": [
                                            x[5] for x in more_match_quals
                                        ],
                                    }
                                )

        # Drop scans with no results from dictionary
        overall_results_dict = {k: v for k, v in overall_results_dict.items() if v}

        # Cast each entry as a MS2SearchResults object
        for scan_id in overall_results_dict.keys():
            for precursor_mz in overall_results_dict[scan_id].keys():
                ms2_spectrum = self._ms[scan_id]
                ms2_search_results = overall_results_dict[scan_id][precursor_mz]
                overall_results_dict[scan_id][precursor_mz] = SpectrumSearchResults(
                    ms2_spectrum, precursor_mz, ms2_search_results
                )

        # Add MS2SearchResults to the existing spectral search results dictionary
        self.spectral_search_results.update(overall_results_dict)

        # If there are mass features, associate the results with each mass feature
        if len(self.mass_features) > 0:
            for mass_feature_id, mass_feature in self.mass_features.items():
                scan_ids = mass_feature.ms2_scan_numbers
                for ms2_scan_id in scan_ids:
                    precursor_mz = mass_feature.mz
                    try:
                        self.spectral_search_results[ms2_scan_id][precursor_mz]
                    except KeyError:
                        pass
                    else:
                        self.mass_features[
                            mass_feature_id
                        ].ms2_similarity_results.append(
                            self.spectral_search_results[ms2_scan_id][precursor_mz]
                        )
