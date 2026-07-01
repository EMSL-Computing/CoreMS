from corems.molecular_id.calc.ClusterFilter import ClusteringFilter


class MolecularFormulaSearchFilters:
    """Class containing static methods for filtering molecular formulas in a mass spectrum.

    Methods
    -------
    * filter_kendrick(ms_peak_indexes, mass_spectrum_obj).
        Apply Kendrick filter to the mass spectrum.
    * check_min_peaks(ms_peak_indexes, mass_spectrum_obj).
        Check if the number of peaks per class meets the minimum requirement.
    * filter_isotopologue(ms_peak_indexes, mass_spectrum_obj).
        Apply isotopologue filter to the mass spectrum.

    """

    @staticmethod
    def filter_kendrick(ms_peak_indexes, mass_spectrum_obj):
        """Apply Kendrick filter to the mass spectrum.

        Parameters
        ----------
        ms_peak_indexes : list
            List of peak indexes and their associated molecular formula objects.
        mass_spectrum_obj : MassSpectrum
            The mass spectrum object.

        Returns
        -------
        filtered_ms_peak_indexes : list
            List of peak indexes and their associated molecular formula objects after applying the Kendrick filter.
        """
        index_to_remove = []

        if mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter:
            index_to_remove = ClusteringFilter().filter_kendrick_by_index(
                ms_peak_indexes, mass_spectrum_obj
            )

            # for index in noise_indexes: self.mass_spectrum_obj[index].clear_molecular_formulas()

        all_index_to_remove = []

        for peak_index, mf_obj in index_to_remove:
            ms_peak_indexes.remove((peak_index, mf_obj))

            all_index_to_remove.extend(mf_obj.mspeak_mf_isotopologues_indexes)

        all_index_to_remove = list(set(all_index_to_remove + index_to_remove))

        for peak_index, mf_obj in all_index_to_remove:
            mass_spectrum_obj[peak_index].remove_molecular_formula(mf_obj)

        return ms_peak_indexes

    @staticmethod
    def check_min_peaks(ms_peak_indexes, mass_spectrum):
        """Check if the number of peaks per class meets the minimum requirement.

        Parameters
        --------
        ms_peak_indexes : list
            List of peak indexes and their associated molecular formula objects.
        mass_spectrum_obj : MassSpectrum
            The mass spectrum object.

        """
        if mass_spectrum.molecular_search_settings.use_min_peaks_filter:
            if (
                not len(ms_peak_indexes)
                >= mass_spectrum.molecular_search_settings.min_peaks_per_class
            ):
                for peak_index, mf_obj in ms_peak_indexes:
                    mass_spectrum[peak_index].remove_molecular_formula(mf_obj)

    @staticmethod
    def filter_isotopologue(ms_peak_indexes, mass_spectrum):
        """Apply isotopologue filter to the mass spectrum.

        Parameters
        --------
        ms_peak_indexes : list
            List of peak indexes and their associated molecular formula objects.
        mass_spectrum_obj : MassSpectrum
            The mass spectrum object.

        Returns
        ------------
        filtered_ms_peak_indexes : list
            List of peak indexes and their associated molecular formula objects after applying the isotopologue filter.
        """
        index_to_remove = []
        # print(len(ms_peak_indexes))
        if mass_spectrum.molecular_search_settings.use_isotopologue_filter:
            atoms_iso_filter = (
                mass_spectrum.molecular_search_settings.isotopologue_filter_atoms
            )

            isotopologue_count_threshold = (
                mass_spectrum.molecular_search_settings.isotopologue_filter_threshold
            )

            for mspeak_index, mf_obj in ms_peak_indexes:
                if mf_obj.isotopologue_count_percentile < isotopologue_count_threshold:
                    if set(mf_obj.atoms).intersection(atoms_iso_filter):
                        # removes tuple obj from initial list to be used on next filter steps
                        ms_peak_indexes.remove((mspeak_index, mf_obj))

                        # current mf_obj
                        index_to_remove.append((mspeak_index, mf_obj))
                        # all other associated isotopolgues
                        index_to_remove.extend(mf_obj.mspeak_mf_isotopologues_indexes)

        # iterate over all indexes to be remove and remove the mf from the mspeak

        # print(len(ms_peak_indexes))
        for peak_index, mf_obj in index_to_remove:
            # print(peak_index, mf_obj)
            mass_spectrum[peak_index].remove_molecular_formula(mf_obj)

        return ms_peak_indexes
