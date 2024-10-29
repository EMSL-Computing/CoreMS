class KendrickGrouping:
    """Class for Kendrick grouping of mass spectra.

    Methods
    -------
    * mz_odd_even_index_lists().
        Get odd and even indexes lists.
    * calc_error(current, test).
        Calculate the error between two values.
    * populate_kendrick_index_dict_error(list_indexes, sort=True).
        Populate the Kendrick index dictionary based on error.
    * populate_kendrick_index_dict_rounding(list_indexes, sort=True).
        Populate the Kendrick index dictionary based on rounding.
    * sort_abundance_kendrick_dict(even_kendrick_group_index, odd_kendrick_group_index).
        Sort the Kendrick index dictionary based on abundance.
    * kendrick_groups_indexes(sort=True).
        Get the Kendrick groups indexes dictionary.

    """

    def mz_odd_even_index_lists(self):
        """Get odd and even indexes lists.

        Returns
        -------
        tuple
            A tuple containing the lists of even and odd indexes.

        """
        even_idx = []
        odd_idx = []

        for i, mspeak in enumerate(self.mspeaks):
            if mspeak.nominal_mz_exp % 2 == 0:
                even_idx.append(i)
            else:
                odd_idx.append(i)

        return even_idx, odd_idx

    def calc_error(self, current: float, test: float):
        """Calculate the error between two values.

        Parameters
        ----------
        current : float
            The current value.
        test : float
            The test value.

        Returns
        -------
        float
            The calculated error.

        """
        return ((test - current) / current) * 1e6

    def populate_kendrick_index_dict_error(self, list_indexes: list, sort: bool = True):
        """Populate the Kendrick index dictionary based on error.

        Parameters
        ----------
        list_indexes : list
            The list of indexes.
        sort : bool, optional
            Whether to sort the dictionary by abundance (default is True).

        Returns
        -------
        dict
            The Kendrick index dictionary.

        """

        def error():
            return abs(current_kmd_reference - next_mspeak.kmd)

        already_found = []

        all_results = []

        for i in list_indexes:
            result_indexes = []

            mspeak = self.mspeaks[i]

            current_kmd_reference = mspeak.kmd

            for j in list_indexes:
                if j not in already_found and j != i:
                    next_mspeak = self.mspeaks[j]

                    if error() <= 0.001:
                        result_indexes.append(j)
                        already_found.append(j)

                        current_kmd_reference = next_mspeak.kmd

            if result_indexes and len(result_indexes) > 3:
                already_found.append(i)

                result_indexes.insert(0, i)

                all_results.append(result_indexes)
            else:
                for w in result_indexes:
                    already_found.remove(w)

        kendrick_group_index = {
            i: indexes_list for i, indexes_list in enumerate(all_results)
        }

        # return dictionary with the keys sorted by sum of the abundances
        if sort:
            return dict(
                sorted(
                    kendrick_group_index.items(),
                    key=lambda it: sum([self.mspeaks[i].abundance for i in it[1]]),
                    reverse=False,
                )
            )

        else:
            return kendrick_group_index

    def populate_kendrick_index_dict_rounding(
        self, list_indexes: list, sort: bool = True
    ):
        """Populate the Kendrick index dictionary based on rounding.

        Parameters
        ----------
        list_indexes : list
            The list of indexes.
        sort : bool, optional
            Whether to sort the dictionary by abundance (default is True).

        Returns
        -------
        dict
            The Kendrick index dictionary.

        """
        kendrick_group_index = {}

        for i in list_indexes:
            mspeak = self.mspeaks[i]

            group = round(mspeak.kmd * 100)

            if group not in kendrick_group_index:
                kendrick_group_index[group] = [i]

            else:
                last_index = kendrick_group_index[group][-1]

                if self.parameters.mass_spectrum.verbose_processing:
                    print(abs(mspeak.kmd - self.mspeaks[last_index].kmd))

                if abs(mspeak.kmd - self.mspeaks[last_index].kmd) < 0.001:
                    kendrick_group_index[group].append(i)

            # return dictionary with the keys sorted by sum of the abundances
        if sort:
            return dict(
                sorted(
                    kendrick_group_index.items(),
                    key=lambda it: sum([self.mspeaks[i].abundance for i in it[1]]),
                    reverse=True,
                )
            )

        else:
            return kendrick_group_index

    def sort_abundance_kendrick_dict(
        self, even_kendrick_group_index: dict, odd_kendrick_group_index: dict
    ):
        """Sort the Kendrick index dictionary based on abundance.

        Parameters
        ----------
        even_kendrick_group_index : dict
            The Kendrick index dictionary for even indexes.
        odd_kendrick_group_index : dict
            The Kendrick index dictionary for odd indexes.

        Returns
        -------
        dict
            The sorted Kendrick index dictionary.

        """
        all_even_indexes = [i for v in even_kendrick_group_index.values() for i in v]

        all_odd_indexes = [i for v in odd_kendrick_group_index.values() for i in v]

        sum_even = sum([self.mspeaks[i].abundance for i in all_even_indexes])

        sum_odd = sum([self.mspeaks[i].abundance for i in all_odd_indexes])

        if sum_even >= sum_odd:
            even_kendrick_group_index.update(odd_kendrick_group_index)

            return even_kendrick_group_index

        else:
            odd_kendrick_group_index.update(even_kendrick_group_index)

            return odd_kendrick_group_index

    def kendrick_groups_indexes(self, sort: bool = True):
        """Get the Kendrick groups indexes dictionary.

        Parameters
        ----------
        sort : bool, optional
            Whether to sort the dictionary by abundance (default is True).

        Returns
        -------
        dict
            The Kendrick groups indexes dictionary.

        """
        even_idx, odd_idx = self.mz_odd_even_index_lists()

        even_kendrick_group_index = self.populate_kendrick_index_dict_error(
            even_idx, sort=sort
        )

        odd_kendrick_group_index = self.populate_kendrick_index_dict_error(
            odd_idx, sort=sort
        )

        return self.sort_abundance_kendrick_dict(
            even_kendrick_group_index, odd_kendrick_group_index
        )
