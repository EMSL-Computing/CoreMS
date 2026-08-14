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

        if "mz" in lib_entry.keys():
            # Get the original mz values from the library entry
            lib_mzs = np.array(
                re.findall(r"\(([^,]+),([^)]+)\)", lib_entry["mz"]), dtype=float
            ).reshape(-1, 2)[:, 0]
        elif "peaks" in lib_entry.keys() and lib_entry["peaks"] is not None:
            lib_mzs = lib_entry["peaks"][:, 0]

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

            # Get types of fragments in the lib entry and convert it to a list on ", "
            lib_frags = [x.strip() for x in lib_entry["fragment_types"].split(",")]
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
        peak_sep_da=None,
        get_additional_metrics=True,
        accumulate_results=False,
        ms_params_key: str = "ms2",
    ):
        """
        Search LC-MS MS2 spectra with FlashEntropy against a pre-built library.

        Annotation thresholds and default peak-separation for cleaning come from
        :class:`~corems.encapsulation.factory.processingSetting.SpectralSimilaritySearchSettings`
        on this object's parameter tree
        (``parameters.mass_spectrum[ms_params_key].spectral_similarity_search``).
        Build ``fe_lib`` with the **same** settings so library index tolerances
        match search cleaning (see
        :meth:`~corems.molecular_id.search.database_interfaces.MSPInterface.get_metabolomics_spectra_library`
        and
        :meth:`~corems.molecular_id.search.database_interfaces.LCLipidLibraryInterface.get_lipid_library`).

        Parameters
        ----------
        scan_list : list
            List of scan numbers to search (must already be loaded in ``self._ms``).
        fe_lib : :obj:`~ms_entropy.FlashEntropySearch`
            FlashEntropy search instance. Prefer constructing it with
            ``settings=settings_from_lcms(self)`` (or the same profile key as
            *ms_params_key*) so FE build knobs match this search.
        precursor_mz_list : list, optional
            List of precursor m/z values to search, by default [], which implies
            matched with mass features; to enable this use_mass_features must be True.
        use_mass_features : bool, optional
            If True, use mass features to get precursor m/z values, by default True.
            If True, will add search results to mass features' ms2_similarity_results attribute.
        peak_sep_da : float, optional
            Minimum separation between m/z peaks spectra in Da. This needs match the
            approximate resolution of the search spectra and the FlashEntropySearch
            instance. If None (default), uses
            ``parameters.mass_spectrum[ms_params_key].spectral_similarity_search.resolved_peak_sep_da``
            (typically ``2 * max_ms2_tolerance_in_da`` from that settings object).
        get_additional_metrics : bool, optional
            If True, get additional metrics from FlashEntropy search, by default True.
        accumulate_results : bool, optional
            If True, accumulate results with existing spectral_search_results instead of
            replacing them. This allows searching the same scans with multiple libraries
            without overwriting previous results, by default False.
        ms_params_key : str, optional
            Key in ``parameters.mass_spectrum`` whose ``spectral_similarity_search``
            settings are used (default ``"ms2"``). Use another key (e.g.
            ``"ms2_cid"``) when multiple MS2 parameter profiles are configured
            for different scan classes.

        Returns
        -------
        None
            Adds results to ``self.spectral_search_results`` and associates them
            with mass features in ``self.mass_features`` when applicable.

        Notes
        -----
        From the settings profile this method reads:

        - ``ms2_min_fe_score`` — minimum entropy score to keep a hit
        - ``include_fragment_types`` — lipid-style fragment-type metrics
        - default peak separation when *peak_sep_da* is None

        FlashEntropy **library** build fields (``max_ms2_tolerance_in_da``, etc.)
        should already have been applied when constructing *fe_lib*. Configure
        them on the same ``spectral_similarity_search`` instance before library
        generation.

        For an :class:`~corems.mass_spectra.factory.lc_class.LCMSCollection`,
        sample parameters are equal at construction; build one shared library
        with ``settings_from_lcms_collection(collection)``, then call
        ``fe_search`` on each sample (or via collection pipeline operations).

        Examples
        --------
        Single LCMSBase object (library + search share profile ``"ms2"``)::

            from corems.encapsulation.factory.parameters import settings_from_lcms

            settings = settings_from_lcms(lcms_obj, profile="ms2")
            fe_lib, _meta = msp.get_metabolomics_spectra_library(
                polarity="positive",
                format="flashentropy",
                settings=settings,
            )
            lcms_obj.fe_search(scan_list=ms2_scans, fe_lib=fe_lib)
            # peak_sep / score gate from mass_spectrum["ms2"].spectral_similarity_search

        Alternate MS2 profile (e.g. CID bag)::

            lcms_obj.fe_search(
                scan_list=cid_scans,
                fe_lib=fe_lib_cid,
                ms_params_key="ms2_cid",
            )
        """
        # Annotation knobs from nested SpectralSimilaritySearchSettings (default profile "ms2")
        if ms_params_key not in self.parameters.mass_spectrum:
            raise KeyError(
                f"ms_params_key={ms_params_key!r} not in parameters.mass_spectrum "
                f"(keys={list(self.parameters.mass_spectrum)})"
            )
        ms2_p = self.parameters.mass_spectrum[ms_params_key].spectral_similarity_search
        include_fragment_types = ms2_p.include_fragment_types
        min_match_score = ms2_p.ms2_min_fe_score
        if peak_sep_da is None:
            peak_sep_da = ms2_p.resolved_peak_sep_da

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
                        # Add database name, if present
                        db_name = [
                            fe_lib[x].get("database_name") for x in match_inds
                        ]
                        if db_name is not None:
                            overall_results_dict[scan_oi][precursor_mz].update(
                                {"database_name": db_name}
                            )
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
        if accumulate_results:
            # Merge results with existing spectral_search_results
            for scan_id, precursor_dict in overall_results_dict.items():
                if scan_id in self.spectral_search_results:
                    # Scan already has results, merge precursor_mz dictionaries
                    self.spectral_search_results[scan_id].update(precursor_dict)
                else:
                    # New scan, add entire dictionary
                    self.spectral_search_results[scan_id] = precursor_dict
        else:
            # Replace existing results (original behavior)
            self.spectral_search_results.update(overall_results_dict)

        # If there are mass features, associate the results with each mass feature
        if len(self.mass_features) > 0:
            # Determine which results to associate with mass features
            if accumulate_results:
                # When accumulating, only associate new results from this search
                # to avoid duplicating previously associated results
                results_to_associate = overall_results_dict
            else:
                # When not accumulating, clear existing associations and re-associate all results
                for mass_feature_id in self.mass_features.keys():
                    self.mass_features[mass_feature_id].ms2_similarity_results = []
                results_to_associate = self.spectral_search_results
            
            for mass_feature_id, mass_feature in self.mass_features.items():
                scan_ids = mass_feature.ms2_scan_numbers
                for ms2_scan_id in scan_ids:
                    precursor_mz = mass_feature.mz
                    try:
                        results_to_associate[ms2_scan_id][precursor_mz]
                    except KeyError:
                        pass
                    else:
                        self.mass_features[
                            mass_feature_id
                        ].ms2_similarity_results.append(
                            results_to_associate[ms2_scan_id][precursor_mz]
                        )
