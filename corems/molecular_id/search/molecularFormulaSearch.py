__author__ = "Yuri E. Corilo"
__date__ = "Jul 29, 2019"


from typing import List

import tqdm

from corems import chunks, timeit
from corems.encapsulation.constant import Labels
from corems.molecular_formula.factory.MolecularFormulaFactory import (
    LCMSLibRefMolecularFormula,
    MolecularFormula,
)
from corems.molecular_id.factory.MolecularLookupTable import MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
from corems.ms_peak.factory.MSPeakClasses import _MSPeak

last_error = 0
last_dif = 0
closest_error = 0
error_average = 0
nbValues = 0


class SearchMolecularFormulas:
    """Class for searching molecular formulas in a mass spectrum.

    Parameters
    ----------
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    sql_db : MolForm_SQL, optional
        The SQL database object, by default None.
    first_hit : bool, optional
        Flag to indicate whether to skip peaks that already have a molecular formula assigned, by default False.
    find_isotopologues : bool, optional
        Flag to indicate whether to find isotopologues, by default True.

    Attributes
    ----------
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    sql_db : MolForm_SQL
        The SQL database object.
    first_hit : bool
        Flag to indicate whether to skip peaks that already have a molecular formula assigned.
    find_isotopologues : bool
        Flag to indicate whether to find isotopologues.


    Methods
    -------
    * run_search().
        Run the molecular formula search.
    * run_worker_mass_spectrum().
        Run the molecular formula search on the mass spectrum object.
    * run_worker_ms_peaks().
        Run the molecular formula search on the given list of mass spectrum peaks.
    * database_to_dict().
        Convert the database results to a dictionary.
    * run_molecular_formula().
        Run the molecular formula search on the given list of mass spectrum peaks.
    * search_mol_formulas().
        Search for molecular formulas in the mass spectrum.

    """

    def __init__(
        self,
        mass_spectrum_obj,
        sql_db=None,
        first_hit: bool = False,
        find_isotopologues: bool = True,
    ):
        self.first_hit = first_hit

        self.find_isotopologues = find_isotopologues

        self.mass_spectrum_obj = mass_spectrum_obj

        if not sql_db:
            self.sql_db = MolForm_SQL(
                url=mass_spectrum_obj.molecular_search_settings.url_database
            )

        else:
            self.sql_db = sql_db

    def __enter__(self):
        """Open the SQL database connection."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close the SQL database connection."""
        self.sql_db.close()

        return False

    def run_search(
        self,
        mspeaks: list,
        query: dict,
        min_abundance: float,
        ion_type: str,
        ion_charge: int,
        adduct_atom=None,
    ):
        """Run the molecular formula search.

        Parameters
        ----------
        mspeaks : list of MSPeak
            The list of mass spectrum peaks.
        query : dict
            The query dictionary containing the possible molecular formulas.
        min_abundance : float
            The minimum abundance threshold.
        ion_type : str
            The ion type.
        ion_charge : int
            The ion charge.
        adduct_atom : str, optional
            The adduct atom, by default None.
        """

        def get_formulas(nominal_overlay: float = 0.1):
            """
            Get the list of formulas based on the nominal overlay.

            Parameters
            ----------
            nominal_overlay : float, optional
                The nominal overlay, by default 0.1.

            Returns
            -------
            list
                The list of formulas.
            """
            nominal_mz = ms_peak.nominal_mz_exp

            defect_mass = ms_peak.mz_exp - nominal_mz
            nominal_masses = [nominal_mz]

            if (defect_mass) >= 1 - nominal_overlay:
                nominal_masses.append(nominal_mz + 1)
            elif (defect_mass) <= nominal_overlay:
                nominal_masses.append(nominal_mz - 1)

            list_formulas_candidates = []

            for nominal_mass in nominal_masses:
                if nominal_mass in query.keys():
                    list_formulas_candidates.extend(query.get(nominal_mass))

            return list_formulas_candidates

        all_assigned_indexes = list()

        # molecular_search_settings = self.mass_spectrum_obj.molecular_search_settings

        search_molfrom = SearchMolecularFormulaWorker(
            find_isotopologues=self.find_isotopologues
        )

        for ms_peak in mspeaks:
            # already assigned a molecular formula
            if self.first_hit:
                if ms_peak.is_assigned:
                    continue

            ms_peak_indexes = search_molfrom.find_formulas(
                get_formulas(),
                min_abundance,
                self.mass_spectrum_obj,
                ms_peak,
                ion_type,
                ion_charge,
                adduct_atom,
            )

            all_assigned_indexes.extend(ms_peak_indexes)

        # all_assigned_indexes = MolecularFormulaSearchFilters().filter_isotopologue(all_assigned_indexes, self.mass_spectrum_obj)

        # all_assigned_indexes = MolecularFormulaSearchFilters().filter_kendrick(all_assigned_indexes, self.mass_spectrum_obj)

        # MolecularFormulaSearchFilters().check_min_peaks(all_assigned_indexes, self.mass_spectrum_obj)
        # filter per min peaks per mono isotopic class

    def run_worker_mass_spectrum(self):
        """Run the molecular formula search on the mass spectrum object."""
        self.run_molecular_formula(
            self.mass_spectrum_obj.sort_by_abundance(), 
            print_time=self.mass_spectrum_obj.molecular_search_settings.verbose_processing
            )

    def run_worker_ms_peaks(self, ms_peaks):
        """Run the molecular formula search on the given list of mass spectrum peaks.

        Parameters
        ----------
        ms_peaks : list of MSPeak
            The list of mass spectrum peaks.
        """
        self.run_molecular_formula(
            ms_peaks,
            print_time=self.mass_spectrum_obj.molecular_search_settings.verbose_processing
        )

    @staticmethod
    def database_to_dict(classe_str_list, nominal_mzs, mf_search_settings, ion_charge):
        """Convert the database results to a dictionary.

        Parameters
        ----------
        classe_str_list : list
            The list of class strings.
        nominal_mzs : list
            The list of nominal m/z values.
        mf_search_settings : MolecularFormulaSearchSettings
            The molecular formula search settings.
        ion_charge : int
            The ion charge.

        Returns
        -------
        dict
            The dictionary containing the database results.
        """
        sql_db = MolForm_SQL(url=mf_search_settings.url_database)

        dict_res = {}

        if mf_search_settings.isProtonated:
            dict_res[Labels.protonated_de_ion] = sql_db.get_dict_by_classes(
                classe_str_list,
                Labels.protonated_de_ion,
                nominal_mzs,
                ion_charge,
                mf_search_settings,
            )

        if mf_search_settings.isRadical:
            dict_res[Labels.radical_ion] = sql_db.get_dict_by_classes(
                classe_str_list,
                Labels.radical_ion,
                nominal_mzs,
                ion_charge,
                mf_search_settings,
            )

        if mf_search_settings.isAdduct:
            adduct_list = (
                mf_search_settings.adduct_atoms_neg
                if ion_charge < 0
                else mf_search_settings.adduct_atoms_pos
            )
            dict_res[Labels.adduct_ion] = sql_db.get_dict_by_classes(
                classe_str_list,
                Labels.adduct_ion,
                nominal_mzs,
                ion_charge,
                mf_search_settings,
                adducts=adduct_list,
            )

        return dict_res

    @timeit(print_time=True)
    def run_molecular_formula(self, ms_peaks, **kwargs):
        """Run the molecular formula search on the given list of mass spectrum peaks.

        Parameters
        ----------
        ms_peaks : list of MSPeak
            The list of mass spectrum peaks.
        **kwargs
            Additional keyword arguments. 
            Most notably, print_time, which is a boolean flag to indicate whether to print the time 
            and passed to the timeit decorator.
        """
        ion_charge = self.mass_spectrum_obj.polarity
        min_abundance = self.mass_spectrum_obj.min_abundance
        nominal_mzs = self.mass_spectrum_obj.nominal_mz

        verbose = self.mass_spectrum_obj.molecular_search_settings.verbose_processing
        # reset average error, only relevant is average mass error method is being used
        SearchMolecularFormulaWorker(
            find_isotopologues=self.find_isotopologues
        ).reset_error(self.mass_spectrum_obj)

        # check database for all possible molecular formula combinations based on the setting passed to self.mass_spectrum_obj.molecular_search_settings
        classes = MolecularCombinations(self.sql_db).runworker(
            self.mass_spectrum_obj.molecular_search_settings,
            print_time=self.mass_spectrum_obj.molecular_search_settings.verbose_processing
        )

        # split the database load to not blowout the memory
        # TODO add to the settings
        for classe_chunk in chunks(
            classes, self.mass_spectrum_obj.molecular_search_settings.db_chunk_size
        ):
            classes_str_list = [class_tuple[0] for class_tuple in classe_chunk]

            # load the molecular formula objs binned by ion type and heteroatoms classes, {ion type:{classe:[list_formula]}}
            # for adduct ion type a third key is added {atoms:{ion type:{classe:[list_formula]}}}
            dict_res = self.database_to_dict(
                classes_str_list,
                nominal_mzs,
                self.mass_spectrum_obj.molecular_search_settings,
                ion_charge,
            )
            pbar = tqdm.tqdm(classe_chunk, disable = not verbose)
            for classe_tuple in pbar:
                # class string is a json serialized dict
                classe_str = classe_tuple[0]
                classe_dict = classe_tuple[1]

                if self.mass_spectrum_obj.molecular_search_settings.isProtonated:
                    ion_type = Labels.protonated_de_ion
                    if verbose:
                        pbar.set_description_str(
                            desc="Started molecular formula search for class %s, (de)protonated "
                            % classe_str,
                            refresh=True,
                        )

                    candidate_formulas = dict_res.get(ion_type).get(classe_str)

                    if candidate_formulas:
                        self.run_search(
                            ms_peaks,
                            candidate_formulas,
                            min_abundance,
                            ion_type,
                            ion_charge,
                        )

                if self.mass_spectrum_obj.molecular_search_settings.isRadical:
                    if verbose:
                        pbar.set_description_str(
                            desc="Started molecular formula search for class %s, radical "
                            % classe_str,
                            refresh=True,
                        )

                    ion_type = Labels.radical_ion

                    candidate_formulas = dict_res.get(ion_type).get(classe_str)

                    if candidate_formulas:
                        self.run_search(
                            ms_peaks,
                            candidate_formulas,
                            min_abundance,
                            ion_type,
                            ion_charge,
                        )
                # looks for adduct, used_atom_valences should be 0
                # this code does not support H exchance by halogen atoms
                if self.mass_spectrum_obj.molecular_search_settings.isAdduct:
                    if verbose:
                        pbar.set_description_str(
                            desc="Started molecular formula search for class %s, adduct "
                            % classe_str,
                            refresh=True,
                        )

                    ion_type = Labels.adduct_ion
                    dict_atoms_formulas = dict_res.get(ion_type)

                    for adduct_atom, dict_by_class in dict_atoms_formulas.items():
                        candidate_formulas = dict_by_class.get(classe_str)

                        if candidate_formulas:
                            self.run_search(
                                ms_peaks,
                                candidate_formulas,
                                min_abundance,
                                ion_type,
                                ion_charge,
                                adduct_atom=adduct_atom,
                            )
        self.sql_db.close()

    def search_mol_formulas(
        self,
        possible_formulas_list: List[MolecularFormula],
        ion_type: str,
        neutral_molform=True,
        find_isotopologues=True,
        adduct_atom=None,
    ) -> List[_MSPeak]:
        """Search for molecular formulas in the mass spectrum.

        Parameters
        ----------
        possible_formulas_list : list of MolecularFormula
            The list of possible molecular formulas.
        ion_type : str
            The ion type.
        neutral_molform : bool, optional
            Flag to indicate whether the molecular formulas are neutral, by default True.
        find_isotopologues : bool, optional
            Flag to indicate whether to find isotopologues, by default True.
        adduct_atom : str, optional
            The adduct atom, by default None.

        Returns
        -------
        list of MSPeak
            The list of mass spectrum peaks with assigned molecular formulas.
        """
        # neutral_molform: some reference files already present the formula on ion mode, for instance, bruker reference files
        #    if that is the case than turn neutral_molform off

        SearchMolecularFormulaWorker(find_isotopologues=find_isotopologues).reset_error(
            self.mass_spectrum_obj
        )

        initial_min_peak_bool = (
            self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter
        )
        initial_runtime_kendrick_filter = (
            self.mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter
        )

        # Are the following 3 lines redundant?
        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = (
            False  # TODO check this line
        )
        self.mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter = (
            False
        )

        possible_formulas_dict_nm = {}

        for mf in possible_formulas_list:
            if neutral_molform:
                nm = int(mf.protonated_mz)
            else:
                nm = int(mf.mz_nominal_calc)

            if nm in possible_formulas_dict_nm.keys():
                possible_formulas_dict_nm[nm].append(mf)

            else:
                possible_formulas_dict_nm[nm] = [mf]

        min_abundance = self.mass_spectrum_obj.min_abundance

        ion_type = ion_type

        self.run_search(
            self.mass_spectrum_obj,
            possible_formulas_dict_nm,
            min_abundance,
            ion_type,
            self.mass_spectrum_obj.polarity,
            adduct_atom=adduct_atom,
        )

        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = (
            initial_min_peak_bool
        )
        self.mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter = (
            initial_runtime_kendrick_filter
        )

        mspeaks = [mspeak for mspeak in self.mass_spectrum_obj if mspeak.is_assigned]

        self.sql_db.close()

        return mspeaks


class SearchMolecularFormulaWorker:
    """Class for searching molecular formulas in a mass spectrum.

    Parameters
    ----------
    find_isotopologues : bool, optional
        Flag to indicate whether to find isotopologues, by default True.

    Attributes
    ----------
    find_isotopologues : bool
        Flag to indicate whether to find isotopologues.

    Methods
    -------
    * reset_error().
        Reset the error variables.
    * set_last_error().
        Set the last error.
    * find_formulas().
        Find the formulas.
    * calc_error().
        Calculate the error.
    """

    # TODO add reset error function
    # needs this wraper to pass the class to multiprocessing

    def __init__(self, find_isotopologues=True):
        self.find_isotopologues = find_isotopologues

    def __call__(self, args):
        """Call the find formulas function.

        Parameters
        ----------
        args : tuple
            The arguments.

        Returns
        -------
        list
            The list of mass spectrum peaks with assigned molecular formulas.
        """
        return self.find_formulas(*args)  # ,args[1]

    def reset_error(self, mass_spectrum_obj):
        """Reset the error variables.

        Parameters
        ----------
        mass_spectrum_obj : MassSpectrum
            The mass spectrum object.

        Notes
        -----
        This function resets the error variables for the given mass spectrum object.
        """
        global last_error, last_dif, closest_error, error_average, nbValues
        last_error, last_dif, closest_error, nbValues = 0.0, 0.0, 0.0, 0.0

    def set_last_error(self, error, mass_spectrum_obj):
        """Set the last error.

        Parameters
        ----------
        error : float
            The error.
        mass_spectrum_obj : MassSpectrum
            The mass spectrum object.
        """
        # set the changes to the global variables, not internal ones
        global last_error, last_dif, closest_error, error_average, nbValues

        if mass_spectrum_obj.molecular_search_settings.error_method == "distance":
            dif = error - last_error
            if dif < last_dif:
                last_dif = dif
                closest_error = error
                mass_spectrum_obj.molecular_search_settings.min_ppm_error = (
                    closest_error
                    - mass_spectrum_obj.molecular_search_settings.mz_error_range
                )
                mass_spectrum_obj.molecular_search_settings.max_ppm_error = (
                    closest_error
                    + mass_spectrum_obj.molecular_search_settings.mz_error_range
                )

        elif mass_spectrum_obj.molecular_search_settings.error_method == "lowest":
            if error < last_error:
                mass_spectrum_obj.molecular_search_settings.min_ppm_error = (
                    error - mass_spectrum_obj.molecular_search_settings.mz_error_range
                )
                mass_spectrum_obj.molecular_search_settings.max_ppm_error = (
                    error + mass_spectrum_obj.molecular_search_settings.mz_error_range
                )
                last_error = error

        elif mass_spectrum_obj.molecular_search_settings.error_method == "symmetrical":
            mass_spectrum_obj.molecular_search_settings.min_ppm_error = (
                mass_spectrum_obj.molecular_search_settings.mz_error_average
                - mass_spectrum_obj.molecular_search_settings.mz_error_range
            )
            mass_spectrum_obj.molecular_search_settings.max_ppm_error = (
                mass_spectrum_obj.molecular_search_settings.mz_error_average
                + mass_spectrum_obj.molecular_search_settings.mz_error_range
            )

        elif mass_spectrum_obj.molecular_search_settings.error_method == "average":
            nbValues += 1
            error_average = error_average + ((error - error_average) / nbValues)
            mass_spectrum_obj.molecular_search_settings.min_ppm_error = (
                error_average
                - mass_spectrum_obj.molecular_search_settings.mz_error_range
            )
            mass_spectrum_obj.molecular_search_settings.max_ppm_error = (
                error_average
                + mass_spectrum_obj.molecular_search_settings.mz_error_range
            )

        else:
            # using set mass_spectrum_obj.molecular_search_settings.min_ppm_error  and max_ppm_error range
            pass

        # returns the error based on the selected method at mass_spectrum_obj.molecular_search_settings.method

    @staticmethod
    def calc_error(mz_exp, mz_calc, method="ppm"):
        """Calculate the error.

        Parameters
        ----------
        mz_exp : float
            The experimental m/z value.
        mz_calc : float
            The calculated m/z value.
        method : str, optional
            The method, by default 'ppm'.

        Raises
        -------
        Exception
            If the method is not ppm or ppb.

        Returns
        -------
        float
            The error.
        """

        if method == "ppm":
            multi_factor = 1_000_000

        elif method == "ppb":
            multi_factor = 1_000_000_000

        elif method == "perc":
            multi_factor = 100

        else:
            raise Exception(
                "method needs to be ppm or ppb, you have entered %s" % method
            )

        if mz_exp:
            return ((mz_exp - mz_calc) / mz_calc) * multi_factor

        else:
            raise Exception("Please set mz_calc first")

    def find_formulas(
        self,
        formulas,
        min_abundance,
        mass_spectrum_obj,
        ms_peak,
        ion_type,
        ion_charge,
        adduct_atom=None,
    ):
        """Find the formulas.

        Parameters
        ----------
        formulas : list of MolecularFormula
            The list of molecular formulas.
        min_abundance : float
            The minimum abundance threshold.
        mass_spectrum_obj : MassSpectrum
            The mass spectrum object.
        ms_peak : MSPeak
            The mass spectrum peak.
        ion_type : str
            The ion type.
        ion_charge : int
            The ion charge.
        adduct_atom : str, optional
            The adduct atom, by default None.

        Returns
        -------
        list of MSPeak
            The list of mass spectrum peaks with assigned molecular formulas.

        Notes
        -----
        Uses the closest error the next search (this is not ideal, it needs to use confidence
        metric to choose the right candidate then propagate the error using the error from the best candidate).
        It needs to add s/n to the equation.
        It need optimization to define the mz_error_range within a m/z unit since it is directly proportional
        with the mass, and inversely proportional to the rp. It's not linear, i.e., sigma mass.
        The idea it to correlate sigma to resolving power, signal to noise and sample complexity per mz unit.
        Method='distance'
        """
        mspeak_assigned_index = list()

        min_ppm_error = mass_spectrum_obj.molecular_search_settings.min_ppm_error
        max_ppm_error = mass_spectrum_obj.molecular_search_settings.max_ppm_error

        min_abun_error = mass_spectrum_obj.molecular_search_settings.min_abun_error
        max_abun_error = mass_spectrum_obj.molecular_search_settings.max_abun_error

        # f = open("abundance_error.txt", "a+")
        ms_peak_mz_exp, ms_peak_abundance = ms_peak.mz_exp, ms_peak.abundance
        # min_error = min([pmf.mz_error for pmf in possible_formulas])

        def mass_by_ion_type(possible_formula_obj):
            if ion_type == Labels.protonated_de_ion:
                return possible_formula_obj._protonated_mz(ion_charge)

            elif ion_type == Labels.radical_ion:
                return possible_formula_obj._radical_mz(ion_charge)

            elif ion_type == Labels.adduct_ion and adduct_atom:
                return possible_formula._adduct_mz(ion_charge, adduct_atom)

            else:
                # will return externally calculated mz if is set, #use on Bruker Reference list import
                # if the ion type is known the ion mass based on molecular formula ion type
                # if ion type is unknow will return neutral mass
                return possible_formula.mz_calc

        if formulas:
            if isinstance(formulas[0], LCMSLibRefMolecularFormula):
                possible_mf_class = True

            else:
                possible_mf_class = False

        for possible_formula in formulas:
            if possible_formula:
                error = self.calc_error(
                    ms_peak_mz_exp, mass_by_ion_type(possible_formula)
                )

                # error = possible_formula.mz_error

                if min_ppm_error <= error <= max_ppm_error:
                    # update the error

                    self.set_last_error(error, mass_spectrum_obj)

                    # add molecular formula match to ms_peak

                    # get molecular formula dict from sql obj
                    # formula_dict = pickle.loads(possible_formula.mol_formula)
                    # if possible_mf_class:

                    #    molecular_formula = deepcopy(possible_formula)

                    # else:

                    formula_dict = possible_formula.to_dict()
                    # create the molecular formula obj to be stored
                    if possible_mf_class:
                        molecular_formula = LCMSLibRefMolecularFormula(
                            formula_dict,
                            ion_charge,
                            ion_type=ion_type,
                            adduct_atom=adduct_atom,
                        )

                        molecular_formula.name = possible_formula.name
                        molecular_formula.kegg_id = possible_formula.kegg_id
                        molecular_formula.cas = possible_formula.cas

                    else:
                        molecular_formula = MolecularFormula(
                            formula_dict,
                            ion_charge,
                            ion_type=ion_type,
                            adduct_atom=adduct_atom,
                        )
                    # add the molecular formula obj to the mspeak obj
                    # add the mspeak obj and it's index for tracking next assignment step

                    if self.find_isotopologues:
                        # calculates isotopologues
                        isotopologues = molecular_formula.isotopologues(
                            min_abundance,
                            ms_peak_abundance,
                            mass_spectrum_obj.dynamic_range,
                        )

                        # search for isotopologues
                        for isotopologue_formula in isotopologues:
                            molecular_formula.expected_isotopologues.append(
                                isotopologue_formula
                            )
                            # move this outside to improve preformace
                            # we need to increase the search space to -+1 m_z
                            first_index, last_index = (
                                mass_spectrum_obj.get_nominal_mz_first_last_indexes(
                                    isotopologue_formula.mz_nominal_calc
                                )
                            )

                            for ms_peak_iso in mass_spectrum_obj[
                                first_index:last_index
                            ]:
                                error = self.calc_error(
                                    ms_peak_iso.mz_exp, isotopologue_formula.mz_calc
                                )

                                if min_ppm_error <= error <= max_ppm_error:
                                    # need to define error distribution for abundance measurements

                                    # if mass_spectrum_obj.is_centroid:

                                    abundance_error = self.calc_error(
                                        isotopologue_formula.abundance_calc,
                                        ms_peak_iso.abundance,
                                        method="perc",
                                    )

                                    # area_error = self.calc_error(ms_peak.area, ms_peak_iso.area, method='perc')

                                    # margin of error was set empirically/ needs statistical calculation
                                    #  of margin of error for the measurement of the abundances
                                    if (
                                        min_abun_error
                                        <= abundance_error
                                        <= max_abun_error
                                    ):
                                        # update the error

                                        self.set_last_error(error, mass_spectrum_obj)

                                        # isotopologue_formula.mz_error = error

                                        # isotopologue_formula.area_error = area_error

                                        # isotopologue_formula.abundance_error = abundance_error

                                        isotopologue_formula.mspeak_index_mono_isotopic = ms_peak.index

                                        mono_isotopic_formula_index = len(ms_peak)

                                        isotopologue_formula.mspeak_index_mono_isotopic = ms_peak.index

                                        isotopologue_formula.mono_isotopic_formula_index = mono_isotopic_formula_index

                                        # add mspeaks isotopologue index to the mono isotopic MolecularFormula obj and the respective formula position

                                        # add molecular formula match to ms_peak
                                        x = ms_peak_iso.add_molecular_formula(
                                            isotopologue_formula
                                        )

                                        molecular_formula.mspeak_mf_isotopologues_indexes.append(
                                            (ms_peak_iso.index, x)
                                        )
                                        # add mspeaks mono isotopic index to the isotopologue MolecularFormula obj

                    y = ms_peak.add_molecular_formula(molecular_formula)

                    mspeak_assigned_index.append((ms_peak.index, y))

        return mspeak_assigned_index


class SearchMolecularFormulasLC:
    """Class for searching molecular formulas in a LC object.

    Parameters
    ----------
    lcms_obj : LCMSBase
        The LCMSBase object.
    sql_db : MolForm_SQL, optional
        The SQL database object, by default None.
    first_hit : bool, optional
        Flag to indicate whether to skip peaks that already have a molecular formula assigned, by default False.
    find_isotopologues : bool, optional
        Flag to indicate whether to find isotopologues, by default True.

    Methods
    -------

    * search_spectra_against_candidates().
        Search a list of mass spectra against a list of candidate formulas with a given ion type and charge.
    * bulk_run_molecular_formula_search().
        Run the molecular formula search on the given list of mass spectra.
        Pulls the settings from the LCMSBase object to set ion type and charge to search for. 
    * run_mass_feature_search().
        Run the molecular formula search on mass features.
        Calls bulk_run_molecular_formula_search() with specified mass spectra and mass peaks.
    * run_untargeted_worker_ms1().
        Run untargeted molecular formula search on the ms1 mass spectrum.
        DEPRECATED: use run_mass_feature_search() or bulk_run_molecular_formula_search() instead.
    * run_target_worker_ms1().
        Run targeted molecular formula search on the ms1 mass spectrum.
        DEPRECATED: use run_mass_feature_search() or bulk_run_molecular_formula_search() instead.
    """

    def __init__(self, lcms_obj, sql_db=None, first_hit=False, find_isotopologues=True):
        self.first_hit = first_hit

        self.find_isotopologues = find_isotopologues

        self.lcms_obj = lcms_obj

        if not sql_db:
            self.sql_db = MolForm_SQL(
                url=self.lcms_obj.parameters.mass_spectrum['ms1'].molecular_search.url_database
            )

        else:
            self.sql_db = sql_db

    def search_spectra_against_candidates(self, mass_spectrum_list, ms_peaks_list, candidate_formulas, ion_type, ion_charge):
        """Search a list of mass spectra against a list of candidate formulas with a given ion type and charge.

        Parameters
        ----------
        mass_spectrum_list : list of MassSpectrum
            The list of mass spectra to perform the search on.
        ms_peaks_list : list of lists of MSPeak objects
            The list of mass spectrum peaks to search within each mass spectrum.
        candidate_formulas : dict
            The candidate formulas.
        ion_type : str
            The ion type.
        ion_charge : int
            The ion charge, either 1 or -1.

        Notes
        -----
        This function is designed to be used with the bulk_run_molecular_formula_search function.
        """
        for mass_spectrum, ms_peaks in zip(mass_spectrum_list, ms_peaks_list):
            single_ms_search = SearchMolecularFormulas(
                mass_spectrum,
                sql_db=self.sql_db,
                first_hit=self.first_hit,
                find_isotopologues=self.find_isotopologues,
            )
            single_ms_search.run_search(
                ms_peaks,
                candidate_formulas,
                mass_spectrum.min_abundance,
                ion_type,
                ion_charge,
            )

    def bulk_run_molecular_formula_search(self, mass_spectrum_list, ms_peaks_list, mass_spectrum_setting_key='ms1'):
        """Run the molecular formula search on the given list of mass spectra

        Parameters
        ----------
        mass_spectrum_list : list of MassSpectrum
            The list of mass spectra to search.
        ms_peaks_list : list of lists of MSPeak objects 
            The mass peaks to perform molecular formula search within each mass spectrum
        mass_spectrum_setting_key : str, optional
            The mass spectrum setting key, by default 'ms1'.
            This is used to get the appropriate molecular search settings from the LCMSBase object
        """
        # Set min_abundance and nominal_mzs
        if self.lcms_obj.polarity == "positive":
            ion_charge = 1
        elif self.lcms_obj.polarity == "negative":
            ion_charge = -1
        else:
            raise ValueError("Polarity must be either 'positive' or 'negative'")
        
        # Check that the length of the mass spectrum list and the ms_peaks list are the same
        if len(mass_spectrum_list) != len(ms_peaks_list):
            raise ValueError("The length of the mass spectrum list and the ms_peaks list must be the same")
        
        nominal_mzs = [x.nominal_mz for x in mass_spectrum_list]
        nominal_mzs = list(set([item for sublist in nominal_mzs for item in sublist]))
        verbose = self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search.verbose_processing 

        # reset average error, only relevant if average mass error method is being used
        SearchMolecularFormulaWorker(
            find_isotopologues=self.find_isotopologues
        ).reset_error(mass_spectrum_list[0])

        # check database for all possible molecular formula combinations based on the setting passed to self.mass_spectrum_obj.molecular_search_settings
        classes = MolecularCombinations(self.sql_db).runworker(
            self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search,
            print_time=self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search.verbose_processing
        )
        
        # split the database load to not blowout the memory
        for classe_chunk in chunks(
            classes, self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search.db_chunk_size
        ):
            classes_str_list = [class_tuple[0] for class_tuple in classe_chunk]

            # load the molecular formula objs binned by ion type and heteroatoms classes, {ion type:{classe:[list_formula]}}
            # for adduct ion type a third key is added {atoms:{ion type:{classe:[list_formula]}}}
            dict_res = SearchMolecularFormulas.database_to_dict(
                classes_str_list,
                nominal_mzs,
                self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search,
                ion_charge,
            )

            pbar = tqdm.tqdm(classe_chunk, disable = not verbose)
            for classe_tuple in pbar:
                # class string is a json serialized dict
                classe_str = classe_tuple[0]

                # Perform search for (de)protonated ion type
                if self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search.isProtonated:
                    ion_type = Labels.protonated_de_ion

                    pbar.set_description_str(
                        desc="Started molecular formula search for class %s, (de)protonated "
                        % classe_str,
                        refresh=True,
                    )

                    candidate_formulas = dict_res.get(ion_type).get(classe_str)

                    if candidate_formulas:
                        self.search_spectra_against_candidates(
                            mass_spectrum_list=mass_spectrum_list,
                            ms_peaks_list=ms_peaks_list,
                            candidate_formulas=candidate_formulas,
                            ion_type=ion_type,
                            ion_charge=ion_charge
                        )


                # Perform search for radical ion type
                if self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search.isRadical:
                    pbar.set_description_str(
                        desc="Started molecular formula search for class %s, radical "
                        % classe_str,
                        refresh=True,
                    )

                    ion_type = Labels.radical_ion

                    candidate_formulas = dict_res.get(ion_type).get(classe_str)

                    if candidate_formulas:
                        self.search_spectra_against_candidates(
                            mass_spectrum_list=mass_spectrum_list,
                            ms_peaks_list=ms_peaks_list,
                            candidate_formulas=candidate_formulas,
                            ion_type=ion_type,
                            ion_charge=ion_charge
                        )

                # Perform search for adduct ion type
                # looks for adduct, used_atom_valences should be 0
                # this code does not support H exchance by halogen atoms
                if self.lcms_obj.parameters.mass_spectrum[mass_spectrum_setting_key].molecular_search.isAdduct:
                    pbar.set_description_str(
                        desc="Started molecular formula search for class %s, adduct "
                        % classe_str,
                        refresh=True,
                    )

                    ion_type = Labels.adduct_ion
                    dict_atoms_formulas = dict_res.get(ion_type)

                    for adduct_atom, dict_by_class in dict_atoms_formulas.items():
                        candidate_formulas = dict_by_class.get(classe_str)

                        if candidate_formulas:
                            self.search_spectra_against_candidates(
                                mass_spectrum_list=mass_spectrum_list,
                                ms_peaks_list=ms_peaks_list,
                                candidate_formulas=candidate_formulas,
                                ion_type=ion_type,
                                ion_charge=ion_charge
                            )
        self.sql_db.close()
        
    def run_mass_feature_search(self):
        """Run the molecular formula search on the mass features.
        
        Calls bulk_run_molecular_formula_search() with specified mass spectra and mass peaks.
        """
        mass_features_df = self.lcms_obj.mass_features_to_df()

        # Get the list of mass spectrum (and peaks to search with each mass spectrum) for all mass features
        scan_list = mass_features_df.apex_scan.unique()
        mass_spectrum_list = [self.lcms_obj._ms[x] for x in scan_list]
        ms_peaks = []
        for scan in scan_list:
            mf_df_scan = mass_features_df[mass_features_df.apex_scan == scan]
            peaks_to_search = [
                self.lcms_obj.mass_features[x].ms1_peak for x in mf_df_scan.index.tolist()
            ]
            ms_peaks.append(peaks_to_search)
        
        # Run the molecular formula search
        self.bulk_run_molecular_formula_search(mass_spectrum_list, ms_peaks)
    
    def run_untargeted_worker_ms1(self):
        """Run untargeted molecular formula search on the ms1 mass spectrum."""
        raise NotImplementedError("run_untargeted_worker_ms1 search is not implemented in CoreMS 3.0 and greater")

    def run_target_worker_ms1(self):
        """Run targeted molecular formula search on the ms1 mass spectrum."""
        raise NotImplementedError("run_target_worker_ms1 formula search is not yet implemented in CoreMS 3.0 and greater")
