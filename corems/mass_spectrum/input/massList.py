__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

import warnings

from corems.encapsulation.constant import Atoms, Labels
from corems.mass_spectrum.factory.MassSpectrumClasses import (
    MassSpecCentroid,
    MassSpecProfile,
)
from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula


class ReadCoremsMasslist(MassListBaseClass):
    """
    The ReadCoremsMasslist object reads processed mass list data types
    and returns the mass spectrum obj with the molecular formula obj

    **Only available for centroid mass spectrum type:** it will ignore the parameter **isCentroid**
    Please see MassListBaseClass for more details

    """

    def get_mass_spectrum(self, loadSettings: bool = True) -> MassSpecCentroid:
        """
        Get the mass spectrum object from the processed mass list data.

        Parameters
        ----------
        loadSettings : bool, optional
            Whether to load the settings for the mass spectrum. Default is True.

        Returns
        -------
        MassSpecCentroid
            The mass spectrum object.

        Raises
        ------
        ValueError
            If the input file is not a valid CoreMS file.
        """

        dataframe = self.get_dataframe()

        if not set(
            ["H/C", "O/C", "Heteroatom Class", "Ion Type", "Is Isotopologue"]
        ).issubset(dataframe.columns):
            raise ValueError(
                "%s it is not a valid CoreMS file" % str(self.file_location)
            )

        self.check_columns(dataframe.columns)

        dataframe.rename(columns=self.parameters.header_translate, inplace=True)

        polarity = dataframe["Ion Charge"].values[0]

        output_parameters = self.get_output_parameters(polarity)

        mass_spec_obj = MassSpecCentroid(
            dataframe.to_dict(orient="list"), output_parameters
        )

        if loadSettings is True:
            self.load_settings(mass_spec_obj, output_parameters)

        self.add_molecular_formula(mass_spec_obj, dataframe)

        return mass_spec_obj

    def add_molecular_formula(self, mass_spec_obj, dataframe):
        """
        Add molecular formula information to the mass spectrum object.

        Parameters
        ----------
        mass_spec_obj : MassSpecCentroid
            The mass spectrum object to add the molecular formula to.
        dataframe : pandas.DataFrame
            The processed mass list data.

        """

        # check if is coreMS file
        if "Is Isotopologue" in dataframe:
            # Reindex dataframe to row index to avoid issues with duplicated indexes (e.g. when multiple formula map to single mz_exp)
            dataframe = dataframe.reset_index(drop=True)

            mz_exp_df = dataframe[Labels.mz].astype(float)
            formula_df = dataframe[
                dataframe.columns.intersection(Atoms.atoms_order)
            ].copy()
            formula_df.fillna(0, inplace=True)
            formula_df.replace(b"nan", 0, inplace=True)

            ion_type_df = dataframe["Ion Type"]
            ion_charge_df = dataframe["Ion Charge"]
            is_isotopologue_df = dataframe["Is Isotopologue"]
            if "Adduct" in dataframe:
                adduct_df = dataframe["Adduct"]
            else:
                adduct_df = None

        mass_spec_mz_exp_list = mass_spec_obj.mz_exp

        for df_index, mz_exp in enumerate(mz_exp_df):
            bad_mf = False
            counts = 0

            ms_peak_index = list(mass_spec_mz_exp_list).index(float(mz_exp))

            if "Is Isotopologue" in dataframe:
                atoms = list(formula_df.columns.astype(str))
                counts = list(formula_df.iloc[df_index].astype(int))

                formula_dict = dict(zip(atoms, counts))

                # Drop any atoms with 0 counts
                formula_dict = {
                    atom: formula_dict[atom]
                    for atom in formula_dict
                    if formula_dict[atom] > 0
                }

            if sum(counts) > 0:
                ion_type = str(Labels.ion_type_translate.get(ion_type_df[df_index]))
                if adduct_df is not None:
                    adduct_atom = str(adduct_df[df_index])
                    if adduct_atom == "None":
                        adduct_atom = None
                else:
                    adduct_atom = None

                # If not isotopologue, cast as MolecularFormula
                if not bool(int(is_isotopologue_df[df_index])):
                    mfobj = MolecularFormula(
                        formula_dict,
                        int(ion_charge_df[df_index]),
                        mspeak_parent=mass_spec_obj[ms_peak_index],
                        ion_type=ion_type,
                        adduct_atom=adduct_atom,
                    )

                # if is isotopologue, recast as MolecularFormulaIsotopologue
                if bool(int(is_isotopologue_df[df_index])):
                    # First make a MolecularFormula object for the parent so we can get probabilities etc
                    formula_list_parent = {}
                    for atom in formula_dict:
                        if atom in Atoms.isotopes.keys():
                            formula_list_parent[atom] = formula_dict[atom]
                        else:
                            # remove any numbers from the atom name to cast as a mono-isotopic atom
                            atom_mono = atom.strip("0123456789")
                            if (
                                atom_mono in Atoms.isotopes.keys()
                                and atom_mono in formula_list_parent.keys()
                            ):
                                formula_list_parent[atom_mono] = (
                                    formula_list_parent[atom_mono] + formula_dict[atom]
                                )
                            elif atom_mono in Atoms.isotopes.keys():
                                formula_list_parent[atom_mono] = formula_dict[atom]
                            else:
                                warnings.warn(f"Atom {atom} not in Atoms.atoms_order")
                    mono_index = int(dataframe.iloc[df_index]["Mono Isotopic Index"])
                    mono_mfobj = MolecularFormula(
                        formula_list_parent,
                        int(ion_charge_df[df_index]),
                        mspeak_parent=mass_spec_obj[mono_index],
                        ion_type=ion_type,
                        adduct_atom=adduct_atom,
                    )

                    # Next, generate isotopologues from the parent
                    isos = list(
                        mono_mfobj.isotopologues(
                            min_abundance=mass_spec_obj.abundance.min()*0.01,
                            current_mono_abundance=mass_spec_obj[mono_index].abundance,
                            dynamic_range=mass_spec_obj.dynamic_range,
                        )
                    )

                    # Finally, find the isotopologue that matches the formula_dict
                    matched_isos = []
                    for iso in isos:
                        # If match was already found, exit the loop
                        if len(matched_isos) > 0:
                            break
                        else:
                            # Check the atoms match
                            if set(iso.atoms) == set(formula_dict.keys()):
                                # Check the values of the atoms match
                                if all(
                                    [
                                        iso[atom] == formula_dict[atom]
                                        for atom in formula_dict
                                    ]
                                ):
                                    matched_isos.append(iso)

                    if len(matched_isos) == 0:
                        #FIXME: This should not occur see https://code.emsl.pnl.gov/mass-spectrometry/corems/-/issues/190
                        warnings.warn(f"No isotopologue matched the formula_dict: {formula_dict}")
                        bad_mf = True
                    else:
                        bad_mf = False                   
                        mfobj = matched_isos[0]

                        # Add the mono isotopic index, confidence score and isotopologue similarity
                        mfobj.mspeak_index_mono_isotopic = int(
                            dataframe.iloc[df_index]["Mono Isotopic Index"]
                        )
                if not bad_mf:
                    # Add the confidence score and isotopologue similarity and average MZ error score
                    if "m/z Error Score" in dataframe:
                        mfobj._mass_error_average_score = float(
                            dataframe.iloc[df_index]["m/z Error Score"]
                        )
                    if "Confidence Score" in dataframe:
                        mfobj._confidence_score = float(
                            dataframe.iloc[df_index]["Confidence Score"]
                        )
                    if "Isotopologue Similarity" in dataframe:
                        mfobj._isotopologue_similarity = float(
                            dataframe.iloc[df_index]["Isotopologue Similarity"]
                        )
                    mass_spec_obj[ms_peak_index].add_molecular_formula(mfobj)


class ReadMassList(MassListBaseClass):
    """
    The ReadMassList object reads unprocessed mass list data types
    and returns the mass spectrum object.

    Parameters
    ----------
    MassListBaseClass : class
        The base class for reading mass list data types.

    Methods
    -------
    * get_mass_spectrum(polarity, scan=0, auto_process=True, loadSettings=True). Reads mass list data types and returns the mass spectrum object.

    """

    def get_mass_spectrum(
        self,
        polarity: int,
        scan: int = 0,
        auto_process: bool = True,
        loadSettings: bool = True,
    ):
        """
        Reads mass list data types and returns the mass spectrum object.

        Parameters
        ----------
        polarity : int
            The polarity of the mass spectrum (+1 or -1).
        scan : int, optional
            The scan number of the mass spectrum (default is 0).
        auto_process : bool, optional
            Flag indicating whether to automatically process the mass spectrum (default is True).
        loadSettings : bool, optional
            Flag indicating whether to load settings for the mass spectrum (default is True).

        Returns
        -------
        mass_spec : MassSpecCentroid or MassSpecProfile
            The mass spectrum object.

        """

        # delimiter = "  " or " " or  "," or "\t" etc

        if self.isCentroid:
            dataframe = self.get_dataframe()

            self.check_columns(dataframe.columns)

            self.clean_data_frame(dataframe)

            dataframe.rename(columns=self.parameters.header_translate, inplace=True)

            output_parameters = self.get_output_parameters(polarity)

            mass_spec = MassSpecCentroid(
                dataframe.to_dict(orient="list"),
                output_parameters,
                auto_process=auto_process,
            )

            if loadSettings:
                self.load_settings(mass_spec, output_parameters)

            return mass_spec

        else:
            dataframe = self.get_dataframe()

            self.check_columns(dataframe.columns)

            output_parameters = self.get_output_parameters(polarity)

            self.clean_data_frame(dataframe)

            dataframe.rename(columns=self.parameters.header_translate, inplace=True)

            mass_spec = MassSpecProfile(
                dataframe.to_dict(orient="list"),
                output_parameters,
                auto_process=auto_process,
            )

            if loadSettings:
                self.load_settings(mass_spec, output_parameters)

            return mass_spec


class ReadBrukerXMLList(MassListBaseClass):
    """
    The ReadBrukerXMLList object reads Bruker XML objects
    and returns the mass spectrum object.
    See MassListBaseClass for details

    Parameters
    ----------
    MassListBaseClass : class
        The base class for reading mass list data types and returning the mass spectrum object.

    Methods
    -------
    * get_mass_spectrum(polarity: bool = None, scan: int = 0, auto_process: bool = True, loadSettings: bool = True). Reads mass list data types and returns the mass spectrum object.

    """

    def get_mass_spectrum(
        self,
        polarity: bool = None,
        scan: int = 0,
        auto_process: bool = True,
        loadSettings: bool = True,
    ):
        """
        Reads mass list data types and returns the mass spectrum object.

        Parameters
        ----------
        polarity : bool, optional
            The polarity of the mass spectrum. Can be +1 or -1. If not provided, it will be determined from the XML file.
        scan : int, optional
            The scan number of the mass spectrum. Default is 0.
        auto_process : bool, optional
            Whether to automatically process the mass spectrum. Default is True.
        loadSettings : bool, optional
            Whether to load the settings for the mass spectrum. Default is True.

        Returns
        -------
        mass_spec : MassSpecCentroid
            The mass spectrum object representing the centroided mass spectrum.
        """
        # delimiter = "  " or " " or  "," or "\t" etc

        if polarity == None:
            polarity = self.get_xml_polarity()
        dataframe = self.get_dataframe()

        self.check_columns(dataframe.columns)

        self.clean_data_frame(dataframe)

        dataframe.rename(columns=self.parameters.header_translate, inplace=True)

        output_parameters = self.get_output_parameters(polarity)

        mass_spec = MassSpecCentroid(
            dataframe.to_dict(orient="list"),
            output_parameters,
            auto_process=auto_process,
        )

        if loadSettings:
            self.load_settings(mass_spec, output_parameters)

        return mass_spec
