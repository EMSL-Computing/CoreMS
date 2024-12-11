import os
import sys

from copy import deepcopy
from threading import Thread
from itertools import product

import tqdm

from corems.encapsulation.constant import Labels, Atoms
from corems.molecular_id.calc.MolecularFilter import MolecularFormulaSearchFilters
from corems.molecular_id.search.findOxygenPeaks import FindOxygenPeaks
from corems.molecular_id.search.molecularFormulaSearch import (
    SearchMolecularFormulaWorker,
)
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
from corems.molecular_id.calc.ClusterFilter import ClusteringFilter


class OxygenPriorityAssignment(Thread):
    """A class for assigning priority to oxygen classes in a molecular search.

    Parameters
    ----------
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    sql_db : bool, optional
        Whether to use an SQL database. The default is False.

    Attributes
    ----------
    mass_spectrum_obj : MassSpectrum
        The mass spectrum object.
    sql_db : MolForm_SQL
        The SQL database object.

    Methods
    -------
    * run().
        Run the priority assignment process.
    * create_data_base().
        Create the molecular database for the specified heteroatomic classes.
    * run_worker_mass_spectrum(assign_classes_order_tuples).
        Run the molecular formula search for each class in the specified order.
    * get_dict_molecular_database(classe_str_list).
        Get the molecular database as a dictionary.
    * ox_classes_and_peaks_in_order_().
        Get the oxygen classes and associated peaks in order.
    * get_classes_in_order(dict_ox_class_and_ms_peak)
        Get the classes in order.
    """

    def __init__(self, mass_spectrum_obj, sql_db=False):
        # TODO:- add support for other atoms and adducts: Done
        #        - add dbe range on search runtime : Done
        #        - add docs
        #        - improve performace : Done

        Thread.__init__(self)
        self.mass_spectrum_obj = mass_spectrum_obj
        #  initiated at create_molecular_database()
        # self.dict_molecular_lookup_table = None

        if not sql_db:
            self.sql_db = MolForm_SQL(
                url=mass_spectrum_obj.molecular_search_settings.url_database
            )

        else:
            self.sql_db = sql_db

    def run(self):
        """Run the priority assignment process."""
        # get Oxygen classes dict and the associate mspeak class
        # list_of_classes_min_max_dbe = self.class_and_dbes_in_order()
        # create database separated to give the user the chance to use mass spec filters

        assign_classes_order_str_dict_tuple_list = self.create_data_base()

        if assign_classes_order_str_dict_tuple_list:
            self.run_worker_mass_spectrum(assign_classes_order_str_dict_tuple_list)

        else:
            raise RuntimeError("call create_data_base() first")

        self.sql_db.close()

    def create_data_base(self):
        """Create the molecular database for the specified heteroatomic classes.

        Returns
        -------
        assign_classes_order_str_dict_tuple_ : list
            A list of tuples containing the class names and dictionaries of class attributes.
        """

        def create_molecular_database():
            """Checks and creates the database entries for the specified heteroatomic classes."""
            min_o = min(self.mass_spectrum_obj, key=lambda msp: msp[0]["O"])[0]["O"] - 2

            if min_o <= 0:
                min_o = 1

            max_o = max(self.mass_spectrum_obj, key=lambda msp: msp[0]["O"])[0]["O"] + 2

            # min_dbe = min(self.mass_spectrum_obj, key=lambda msp: msp[0].dbe)[0].dbe

            # max_dbe = max(self.mass_spectrum_obj, key=lambda msp: msp[0].dbe)[0].dbe

            # self.lookupTableSettings.use_pah_line_rule = False

            # self.lookupTableSettings.min_dbe = min_dbe/2#min_dbe - 7 if  (min_dbe - 7) > 0 else 0

            # self.lookupTableSettings.max_dbe = max_dbe * 2 #max_dbe + 7

            self.mass_spectrum_obj.reset_indexes()

            self.mass_spectrum_obj.filter_by_noise_threshold()

            # initial_ox = deepcopy(self.mass_spectrum_obj.molecular_search_settings.usedAtoms)

            self.mass_spectrum_obj.molecular_search_settings.usedAtoms["O"] = (
                min_o,
                max_o,
            )

            self.nominal_mzs = self.mass_spectrum_obj.nominal_mz

        # get the most abundant peak and them every 14Da, only allow Ox and its derivatives
        if self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Getting Oxygen Series")
        find_formula_thread = FindOxygenPeaks(self.mass_spectrum_obj, self.sql_db)
        find_formula_thread.run()

        # mass spec obj indexes are set to interate over only the peaks with a molecular formula candidate
        if self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Getting Oxygen Series")
        find_formula_thread.set_mass_spec_indexes_by_found_peaks()

        # get the Ox class and the DBE for the lowest error molecular formula candidate
        dict_ox_class_and_ms_peak = self.ox_classes_and_peaks_in_order_()

        # sort the classes by abundance
        if self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Getting Oxygen Series Order")
        assign_classes_order_str_dict_tuple_list = self.get_classes_in_order(
            dict_ox_class_and_ms_peak
        )

        create_molecular_database()

        return assign_classes_order_str_dict_tuple_list

    def run_worker_mass_spectrum(self, assign_classes_order_tuples):
        """Run the molecular formula search for each class in the specified order.

        Parameters
        ----------
        assign_classes_order_tuples : list
            A list of tuples containing the class names and dictionaries of class attributes.
        """

        def check_adduct_class(classe_dict):
            """Check if the class contains any adduct atoms.

            Parameters
            ----------
            classe_dict : dict
                The dictionary of class attributes.

            Returns
            -------
            bool
                True if the class contains adduct atoms, False otherwise.
            """
            return any(
                [
                    key in classe_dict.keys()
                    for key in self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_neg
                ]
            )

        def set_min_max_dbe_by_oxygen(classe_dict):
            """Calculate the minimum and maximum DBE based on the number of oxygen atoms.

            Parameters
            ----------
            classe_dict : dict
                The dictionary of class attributes.
            """
            # calculates min and max DBE based on the Oxygen number
            # ref :https://pubs.acs.org/doi/full/10.1021/ac200464q
            # if class does not has O it use the pha rule
            # ref : Vlad Lobodin manuscript to be include here

            # atoms_exchanges = ['N']
            # if 'O' in classe_dict.keys():
            #
            #    Oxygen_number = classe_dict.get("O")
            #    for atom in atoms_exchanges:
            #        if atom in classe_dict.keys():
            #            Oxygen_number += classe_dict.get(atom)
            #
            #    self.mass_spectrum_obj.molecular_search_settings.min_dbe = (Oxygen_number/3) - 0.5
            #    self.mass_spectrum_obj.molecular_search_settings.max_dbe = Oxygen_number*3 + 0.5 + 2
            #
            # else:

            self.mass_spectrum_obj.molecular_search_settings.use_pah_line_rule = True

        def run_search(possible_formulas_dict, mass_spectrum_obj, min_abundance):
            """Run the molecular formula search for each mass spectrum peak.

            Parameters
            ----------
            possible_formulas_dict : dict
                A dictionary of possible molecular formulas.
            mass_spectrum_obj : MassSpectrum
                The mass spectrum object.
            min_abundance : float
                The minimum abundance threshold.

            Returns
            -------
            list
                A list of assigned peak indexes.
            """
            all_assigned_indexes = list()

            for ms_peak in mass_spectrum_obj.sort_by_abundance():
                if ms_peak:
                    continue
                # already assigned a molecular formula

                nominal_mz = ms_peak.nominal_mz_exp

                # get mono isotopic peaks that was added a molecular formula obj
                # TODO update error variables

                possible_formulas_nominal = possible_formulas_dict.get(nominal_mz)

                if possible_formulas_nominal:
                    ms_peak_indexes = SearchMolecularFormulaWorker().find_formulas(
                        possible_formulas_nominal,
                        min_abundance,
                        mass_spectrum_obj,
                        ms_peak,
                    )

                    all_assigned_indexes.extend(ms_peak_indexes)

            # filter peaks by percentile threshold of found isotopologues
            all_assigned_indexes = MolecularFormulaSearchFilters().filter_isotopologue(
                all_assigned_indexes, mass_spectrum_obj
            )

            # filter noise by kendrick density
            all_assigned_indexes = MolecularFormulaSearchFilters().filter_kendrick(
                all_assigned_indexes, mass_spectrum_obj
            )

            # filter per min peaks per mono isotopic class
            # this function should always be the last filter,
            # thefore no need to return remaining indexes
            MolecularFormulaSearchFilters().check_min_peaks(
                all_assigned_indexes, mass_spectrum_obj
            )

        # error_average = self.mass_spectrum_obj.molecular_search_settings.mz_error_average

        kmd_base = self.mass_spectrum_obj.mspeaks_settings.kendrick_base

        self.mass_spectrum_obj.change_kendrick_base_all_mspeaks(kmd_base)

        ClusteringFilter().filter_kendrick(self.mass_spectrum_obj)

        min_abundance = self.mass_spectrum_obj.min_abundance

        list_classes_str = [i[0] for i in assign_classes_order_tuples]
        verbose = self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing
        pbar = tqdm.tqdm(assign_classes_order_tuples, disable= not verbose)
        dict_molecular_lookup_table = self.get_dict_molecular_database(list_classes_str)

        for classe_tuple in pbar:
            classe_str = classe_tuple[0]
            classe_dict = classe_tuple[1]

            set_min_max_dbe_by_oxygen(classe_dict)

            # if len(classe_dict.keys()) == 2:
            #    if classe_dict.get('S') == 1:
            #       continue
            # limits the dbe by the Ox class most abundant,
            # need to add other atoms contribution to be more accurate
            # but +-7 should be sufficient to cover the range

            if self.mass_spectrum_obj.molecular_search_settings.isProtonated:
                # tqdm.set_description_str(desc=None, refresh=True)
                if verbose:
                    pbar.set_description_str(
                        desc="Started molecular formula search for class %s, (de)protonated "
                        % classe_str,
                        refresh=True,
                    )

                ion_type = Labels.protonated_de_ion

                possible_formulas_dict = dict_molecular_lookup_table.get(ion_type).get(
                    classe_str
                )

                if possible_formulas_dict:
                    run_search(
                        possible_formulas_dict, self.mass_spectrum_obj, min_abundance
                    )

            if self.mass_spectrum_obj.molecular_search_settings.isRadical:
                # print("Started molecular formula search for class %s,  radical" % classe_str)
                if verbose:
                    pbar.set_description_str(
                        desc="Started molecular formula search for class %s, radical"
                        % classe_str,
                        refresh=True,
                    )

                ion_type = Labels.radical_ion

                possible_formulas_dict = dict_molecular_lookup_table.get(ion_type).get(
                    classe_str
                )

                if possible_formulas_dict:
                    run_search(
                        possible_formulas_dict, self.mass_spectrum_obj, min_abundance
                    )

            # looks for adduct, used_atom_valences should be 0
            # this code does not support H exchance by halogen atoms
            if self.mass_spectrum_obj.molecular_search_settings.isAdduct:
                if verbose:
                    pbar.set_description_str(
                        desc="Started molecular formula search for class %s, adduct"
                        % classe_str,
                        refresh=True,
                    )
                # print("Started molecular formula search for class %s, adduct" % classe_str)

                ion_type = Labels.radical_ion

                possible_formulas_dict = dict_molecular_lookup_table.get(ion_type).get(
                    classe_str
                )

                """ commenting  unfinished code for release 2.0, see end of file for details"""
                # possible_formulas_adduct =self.add_adducts(possible_formulas_dict)

                # if possible_formulas_adduct:

                run_search(
                    possible_formulas_dict, self.mass_spectrum_obj, min_abundance
                )

    def get_dict_molecular_database(self, classe_str_list):
        """Get the molecular database as a dictionary.

        Parameters
        ----------
        classe_str_list : list
            A list of class names.

        Returns
        -------
        dict
            A dictionary containing the molecular database.
        """
        nominal_mzs = self.nominal_mzs
        mf_search_settings = self.mass_spectrum_obj.molecular_search_settings
        ion_charge = self.mass_spectrum_obj.polarity

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

    def ox_classes_and_peaks_in_order_(self) -> dict:
        """Get the oxygen classes and associated peaks in order.

        Returns
        -------
        dict
            A dictionary containing the oxygen classes and associated peaks.
        """
        # order is only valid in python 3.4 and above
        # change to OrderedDict if your version is lower
        dict_ox_class_and_ms_peak = dict()

        for mspeak in self.mass_spectrum_obj.sort_by_abundance(reverse=True):
            # change this filter to cia filter, give more option here, confidence, number of isotopologue found etc

            ox_classe = mspeak.best_molecular_formula_candidate.class_label

            if ox_classe in dict_ox_class_and_ms_peak.keys():
                # get the most abundant of the same ox class
                if mspeak.abundance > dict_ox_class_and_ms_peak[ox_classe].abundance:
                    dict_ox_class_and_ms_peak[ox_classe] = mspeak
            else:
                dict_ox_class_and_ms_peak[ox_classe] = mspeak

        return dict_ox_class_and_ms_peak

    def get_classes_in_order(self, dict_ox_class_and_ms_peak) -> [(str, dict)]:
        """Get the classes in order.

        Parameters
        ----------
        dict_ox_class_and_ms_peak : dict
            A dictionary containing the oxygen classes and associated peaks.

        Returns
        -------
        list
            A list of tuples containing the class names and dictionaries of class attributes.

        Notes
        -----
        structure is
            ('HC', {'HC': 1})
        """

        usedAtoms = deepcopy(self.mass_spectrum_obj.molecular_search_settings.usedAtoms)

        usedAtoms.pop("C")
        usedAtoms.pop("H")
        usedAtoms.pop("O")

        min_n, max_n = usedAtoms.get("N") if usedAtoms.get("N") else (0, 0)
        min_s, max_s = usedAtoms.get("S") if usedAtoms.get("S") else (0, 0)
        min_p, max_p = usedAtoms.get("P") if usedAtoms.get("P") else (0, 0)

        possible_n = [n for n in range(min_n, max_n + 1)]
        possible_s = [s for s in range(min_s, max_s + 1)]
        possible_p = [p for p in range(min_p, max_p + 1)]

        # used to enforce order for commum atoms
        # and track the atom index in on the tuple in all_atoms_tuples
        atoms_in_order = ["N", "S", "P"]

        # do number atoms prodcut and remove then from the usedAtoms dict
        all_atoms_tuples = product(possible_n, possible_s, possible_p)
        for atom in atoms_in_order:
            usedAtoms.pop(atom, None)

        # iterate over other atoms besides C,H, N, O, S and P

        for selected_atom_label, min_max_tuple in usedAtoms.items():
            min_x = min_max_tuple[0]
            max_x = min_max_tuple[1]

            possible_x = [x for x in range(min_x, max_x + 1)]
            all_atoms_tuples = product(all_atoms_tuples, possible_x)

            # merge tuples
            all_atoms_tuples = [
                all_atoms_combined[0] + (all_atoms_combined[1],)
                for all_atoms_combined in all_atoms_tuples
            ]

            # add atom label to the atoms_in_order list

            # important to index where the atom position is in on the tuple in all_atoms_tuples
            atoms_in_order.append(selected_atom_label)

        classes_strings_dict_tuples, hc_class = self.get_class_strings_dict(
            all_atoms_tuples, atoms_in_order
        )

        combined_classes = self.combine_ox_class_with_other(
            atoms_in_order, classes_strings_dict_tuples, dict_ox_class_and_ms_peak
        )

        combination_classes_ordered = self.sort_classes(
            atoms_in_order, combined_classes
        )

        oxygen_class_str_dict_tuple = [
            (ox_class, mspeak[0].class_dict)
            for ox_class, mspeak in dict_ox_class_and_ms_peak.items()
        ]

        ## add classes together and ignores classes selected from the main series
        for class_tuple in combination_classes_ordered:
            if class_tuple not in oxygen_class_str_dict_tuple:
                oxygen_class_str_dict_tuple.append(class_tuple)

        return oxygen_class_str_dict_tuple

    @staticmethod
    def get_class_strings_dict(all_atoms_tuples, atoms_in_order) -> [(str, dict)]:
        """Get the class strings and dictionaries.

        Parameters
        ----------
        all_atoms_tuples : tuple
            A tuple containing the atoms.
        atoms_in_order : list
            A list of atoms in order.

        Returns
        --------
        list
            A list of tuples containing the class strings and dictionaries.

        """
        classe_list = []
        hc_class = []

        for all_atoms_tuple in all_atoms_tuples:
            classe_str = ""
            classe_dict = dict()

            for each_atoms_index, atoms_number in enumerate(all_atoms_tuple):
                if atoms_number != 0:
                    classe_str = (
                        classe_str
                        + atoms_in_order[each_atoms_index]
                        + str(atoms_number)
                        + " "
                    )

                    classe_dict[atoms_in_order[each_atoms_index]] = atoms_number

            classe_str = classe_str.strip()

            if len(classe_str) > 0:
                classe_list.append((classe_str, classe_dict))

            elif len(classe_str) == 0:
                hc_class.append(("HC", {"HC": 1}))

        return classe_list, hc_class

    @staticmethod
    def combine_ox_class_with_other(
        atoms_in_order, classes_strings_dict_tuples, dict_ox_class_and_ms_peak
    ) -> [dict]:
        """Combine the oxygen classes with other classes.

        Parameters
        ----------
        atoms_in_order : list
            A list of atoms in order.
        classes_strings_dict_tuples : list

        dict_ox_class_and_ms_peak : dict
            A dictionary containing the oxygen classes and associated peaks.

        Returns
        -------
        list
            A list of dictionaries.
        """
        # sort methods that uses the key of classes dictionary and the atoms_in_order as reference
        # c_tuple[1] = class_dict, because is one key:value map we loop through keys and get the first item only
        # sort by len first then sort based on the atoms_in_order list
        atoms_in_order = Atoms.atoms_order

        Oxygen_mfs = dict_ox_class_and_ms_peak.values()

        # sort_method = lambda word: (len(word[0]), [atoms_in_order.index(atom) for atom in list( word[1].keys())])

        # print(classes_strings_dict_tuples)
        # classe_in_order = sorted(classes_strings_dict_tuples, key = sort_method)
        # print(classe_in_order)

        combination = []

        # _ ignoring the class_str
        for _, other_classe_dict in classes_strings_dict_tuples:
            # combination.extend([[other_classe_str + ' ' + Oxygen_mf[0].class_label , {**other_classe_dict, **Oxygen_mf[0].class_dict}] for Oxygen_mf in Oxygen_mfs])
            combination.extend(
                [
                    {**other_classe_dict, **Oxygen_mf[0].class_dict}
                    for Oxygen_mf in Oxygen_mfs
                ]
            )

        return combination

    @staticmethod
    def sort_classes(atoms_in_order, combination_tuples) -> [(str, dict)]:
        """Sort the classes.

        Parameters
        ----------
        atoms_in_order : list
            A list of atoms in order.
        combination_tuples : list

        Returns
        -------
        list
            A list of tuples containing the class strings and dictionaries.
        """
        join_list_of_list_classes = list()
        atoms_in_order = ["N", "S", "P", "O"] + atoms_in_order[3:]

        sort_method = (
            lambda atoms_keys: [atoms_in_order.index(atoms_keys)]
        )  # (len(word[0]), print(word[1]))#[atoms_in_order.index(atom) for atom in list( word[1].keys())])
        for class_dict in combination_tuples:
            sorted_dict_keys = sorted(class_dict, key=sort_method)
            class_str = " ".join(
                [atom + str(class_dict[atom]) for atom in sorted_dict_keys]
            )
            new_class_dict = {atom: class_dict[atom] for atom in sorted_dict_keys}
            join_list_of_list_classes.append((class_str, new_class_dict))

        return join_list_of_list_classes

    '''
    The code bellow is unfinished, might be added to next release, 2.1
    def add_adducts(self, possible_formulas):
        """ Add adducts to the molecular formula candidates.

        Parameters
        ----------
        possible_formulas : dict
            A dictionary of possible molecular formulas.
        
        Returns
        -------
        dict 
            A dictionary of possible molecular formulas with adducts.
        
        """
        ion_type = Labels.adduct_ion

        if self.mass_spectrum_obj.polarity < 0:
            adduct_atoms = self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_neg
            molform_model = MolecularFormulaDict
        else:
            adduct_atoms = self.mass_spectrum_obj.molecular_search_settings.adduct_atoms_pos
            molform_model = MolecularFormulaTablePos

        new_dict = {}
        
        for nominal_mz, list_formulas in possible_formulas.items():
            
            for adduct_atom in adduct_atoms:
                
                adduct_atom_mass= Atoms.atomic_masses.get(adduct_atom) 

                for molecularFormulaTable in  list_formulas:
                    
                    formula_dict = json.loads(molecularFormulaTable.mol_formula)
                    
                    if adduct_atom in formula_dict.keys():
                        formula_dict[adduct_atom] += 1  
                    else:
                        formula_dict[adduct_atom] = 1      
                    
                    mz = adduct_atom_mass + molecularFormulaTable.mz
                    nm = int(mz)
                    
                    new_formul_obj = molform_model( **{"mol_formula" : json.dumps(formula_dict),
                                            "mz" : mz,
                                            "ion_type" : ion_type,
                                            "nominal_mz" : nm,
                                            "ion_charge" : molecularFormulaTable.ion_charge,
                                            "classe" : molecularFormulaTable.classe,
                                            "C" : molecularFormulaTable.C,
                                            "H" : molecularFormulaTable.H,
                                            "N" : molecularFormulaTable.N,
                                            "O" : molecularFormulaTable.O,
                                            "S" : molecularFormulaTable.S,
                                            "P" : molecularFormulaTable.P,
                                            "H_C" : molecularFormulaTable.H_C,
                                            "O_C" : molecularFormulaTable.O_C,
                                            "DBE" : molecularFormulaTable.DBE,
                                            })
                    if nm in new_dict:
                        new_dict[nm].append(new_formul_obj)
                    
                    else:
                        new_dict[nm]= [new_formul_obj]
                    
        return new_dict

    '''
