import re

from corems.encapsulation.constant import Atoms, Labels
from corems.molecular_formula.calc.MolecularFormulaCalc import MolecularFormulaCalc

__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"


class MolecularFormulaBase(MolecularFormulaCalc):
    """Base class for representing a molecular formula.

    Parameters
    ----------
    molecular_formula : dict, list, str
        The molecular formula.
    ion_charge : int
        The ion charge.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.
    mspeak_parent : _MSPeak, optional
        The parent mass spectrum peak object instance. Defaults to None.
    external_mz : float, optional
        The external m/z value. Defaults to None.

    Raises
    ------
    TypeError
        If the ion type is not 'DE_OR_PROTONATED', 'RADICAL' or  'ADDUCT'.

    Attributes
    ----------
    isotopologue_count_percentile : float
        The isotopologue count percentile.
    O_C : float
        The O/C ratio.
    H_C : float
        The H/C ratio.
    dbe : float
        The double bond equivalent.
    mz_nominal_calc : int
        The nominal m/z value.
    mz_error : float
        The m/z error.
    mz_calc : float
        The m/z value.
    protonated_mz : float
        The protonated or deprotonated m/z value.
    radical_mz : float
        The radical m/z value.
    neutral_mass : float
        The neutral mass.
    ion_type : str
        The ion type.
    ion_charge : int
        The ion charge.
    atoms : list
        The atoms in the molecular formula.
    confidence_score : float
        Composite assignment confidence score in about [0, 1] (higher is
        better). Weighted sum of ``average_mz_error_score`` and
        ``isotopologue_similarity`` using
        ``mz_error_score_weight`` / ``isotopologue_score_weight`` from the
        parent spectrum's molecular search settings. See the
        ``corems.molecular_formula`` package documentation for the full
        formulation and literature reference.
    isotopologue_similarity : float
        Similarity between expected and observed isotopologue abundance
        patterns (about [0, 1]). Zero when no isotopologues are expected or
        none were computed during search.
    average_mz_error_score : float
        Mean Gaussian mass-accuracy score over the monoisotopic peak and
        expected isotopologue peaks. Missing expected isotopologues
        contribute 0.
    mz_error_score : float
        Gaussian mass-accuracy score for this formula on its parent peak
        only (ppm error vs ``predicted_std`` or a 1.66 ppm fallback).
    kmd : float
        The Kendrick mass defect (KMD).
    kendrick_mass : float
        The Kendrick mass.
    knm : float
        The nominal Kendrick mass.
    string : str
        The molecular formula string.
    string_formated : str
        The molecular formula string formated with subscripts and superscripts.
    class_label : str
        The class label.
    class_dict : dict
        The class dictionary.

    Methods
    -------
    * change_kendrick_base(kendrick_dict_base).
        Change the Kendrick base.
    * isotopologues(min_abundance, current_mono_abundance, dynamic_range).
        Calculate the isotopologues.
    * atoms_qnt(atom).
        Get the atom quantity.
    * atoms_symbol(atom).
        Get the atom symbol without the mass number.
    * to_dict().
        Get the molecular formula as a dictionary.
    * to_list().
        Get the molecular formula as a list.
    """

    def __init__(
        self,
        molecular_formula,
        ion_charge,
        ion_type=None,
        adduct_atom=None,
        mspeak_parent=None,
        external_mz=None,
    ):
        # clear dictionary of atoms with 0 value
        if type(molecular_formula) is dict:
            self._from_dict(molecular_formula, ion_type, adduct_atom)

        elif type(molecular_formula) is list:
            self._from_list(molecular_formula, ion_type, adduct_atom)

        elif type(molecular_formula) is str:
            self._from_str(molecular_formula, ion_type, adduct_atom)

        self._ion_charge = ion_charge
        self._external_mz = external_mz
        self._confidence_score = None
        self._isotopologue_similarity = None
        self._mz_error_score = None
        self._mass_error_average_score = None

        self.is_isotopologue = False

        # parent mass spectrum peak obj instance
        self._mspeak_parent = mspeak_parent

        self.expected_isotopologues = []
        self.mspeak_mf_isotopologues_indexes = []

        if self._mspeak_parent:
            kendrick_dict_base = (
                self._mspeak_parent._ms_parent.mspeaks_settings.kendrick_base
            )
        else:
            kendrick_dict_base = {"C": 1, "H": 2}
        self._kmd, self._kendrick_mass, self._nominal_km = self._calc_kmd(
            kendrick_dict_base
        )

    def __repr__(self):
        return "MolecularFormula({0},{1},ion type = {2}".format(
            self._d_molecular_formula, self.ion_charge, self.ion_type
        )

    def __str__(self):
        return "MolecularFormula {0}, ion_charge:{1}, ion type:{2}, m/z:{3} ".format(
            self.string, self.ion_charge, self.ion_type, self.mz_calc
        )

    def __len__(self):
        # crash if keys are not ordered
        return len(self._d_molecular_formula.keys())

    def __getitem__(self, atom):
        # atom = list(self._d_molecular_formula.keys())[position]
        if atom in self._d_molecular_formula.keys():
            return self._d_molecular_formula[atom]
        else:
            return 0

    def get(self, atom):
        """Get the atom quantity of a specific atom.

        Parameters
        ----------
        atom : str
            The atom symbol.

        Returns
        -------
        int
            The atom quantity.
        """
        # atom = list(self._d_molecular_formula.keys())[position]
        if atom in self._d_molecular_formula.keys():
            return self._d_molecular_formula[atom]
        else:
            return 0

    def _from_dict(self, molecular_formula, ion_type, adduct_atom):
        self._d_molecular_formula = {
            key: val for key, val in molecular_formula.items() if val != 0
        }

        if ion_type is not None:
            self._d_molecular_formula[Labels.ion_type] = ion_type

        if adduct_atom:
            if adduct_atom in self._d_molecular_formula:
                self._d_molecular_formula[adduct_atom] += 1
            else:
                self._d_molecular_formula[adduct_atom] = 1
        self.adduct_atom = adduct_atom

    def _from_list(self, molecular_formula_list, ion_type, adduct_atom):
        # list has to be in the format
        # ['C', 10, 'H', 21, '13C', 1, 'Cl', 1, etc]
        self._d_molecular_formula = {}
        for each in range(0, len(molecular_formula_list), 2):
            atoms_label = molecular_formula_list[each]
            atoms_count = int(molecular_formula_list[each + 1])

            if atoms_count > 0:
                self._d_molecular_formula[atoms_label] = int(atoms_count)

        self._d_molecular_formula[Labels.ion_type] = ion_type
        if adduct_atom:
            self.adduct_atom = adduct_atom
            if adduct_atom in self._d_molecular_formula:
                self._d_molecular_formula[adduct_atom] += 1
            else:
                self._d_molecular_formula[adduct_atom] = 1
        else:
            self.adduct_atom = None

    def _from_str(self, molecular_formula_str, ion_type, adduct_atom):
        # string has to be in the format
        #'C10 H21 13C1 Cl1 37Cl1 etc'
        # Check if there are spaces in the string
        if " " not in molecular_formula_str:
            raise ValueError(
                "The molecular formula string should have spaces, input: %s"
                % molecular_formula_str
            )

        # Split the string by spaces
        # Grab the text before a digit for each element after splitting on spaces (atoms)
        elements = [re.sub(r"\d+$", "", x) for x in molecular_formula_str.split()]
        # Grab the digits at the end of each element after splitting on spaces (counts)
        counts = [re.findall(r"\d+$", x)[0] for x in molecular_formula_str.split()]
        # Check that the number of elements and counts are the same
        if len(elements) != len(counts):
            raise ValueError(
                "The number of elements and counts do not match, input: %s"
                % molecular_formula_str
            )

        # Create a dictionary from the elements and counts and add it to the molecular formula
        dict_ = dict(zip(elements, counts))
        # Cast counts to integers
        dict_ = {key: int(val) for key, val in dict_.items()}
        self._from_dict(dict_, ion_type, adduct_atom)

    def split(self, delimiters, string, maxsplit=0):  # pragma: no cover
        """Splits the molecular formula string.

        Parameters
        ----------
        delimiters : list
            The list of delimiters.
        string : str
            The molecular formula string.
        maxsplit : int, optional
            The maximum number of splits. Defaults to 0.

        Returns
        -------
        list
            The molecular formula list.

        Notes
        -----
        Does not work when formula has atoms with same characters in a row that below to different atoms, i.e. C10H21NNa.
        """
        regexPattern = "|".join(map(re.escape, delimiters))  # pragma: no cover
        isotopes = re.findall(regexPattern, string)  # pragma: no cover
        counts = re.split(regexPattern, string, maxsplit)  # pragma: no cover

        return [isotopes[0], int(counts[1])]

    @property
    def isotopologue_count_percentile(
        self,
    ):
        if not len(self.expected_isotopologues) == 0:
            return (
                len(self.mspeak_mf_isotopologues_indexes)
                / len(self.expected_isotopologues)
            ) * 100
        else:
            return 100

    @property
    def O_C(self):
        if "O" in self._d_molecular_formula.keys():
            # gather all the Os and Hs, regardless of the isotopic composition
            Os = sum(
                [
                    self._d_molecular_formula.get(key)
                    for key in ["O"] + Atoms.isotopes["O"][1]
                    if key in self._d_molecular_formula.keys()
                ]
            )
            Cs = sum(
                [
                    self._d_molecular_formula.get(key)
                    for key in ["C"] + Atoms.isotopes["C"][1]
                    if key in self._d_molecular_formula.keys()
                ]
            )
            return Os / Cs
        else:
            return 0

    @property
    def H_C(self):
        # gather all the Cs and Hs, regardless of the isotopic composition
        Cs = sum(
            [
                self._d_molecular_formula.get(key)
                for key in ["C"] + Atoms.isotopes["C"][1]
                if key in self._d_molecular_formula.keys()
            ]
        )
        Hs = sum(
            [
                self._d_molecular_formula.get(key)
                for key in ["H"] + Atoms.isotopes["H"][1]
                if key in self._d_molecular_formula.keys()
            ]
        )
        return Hs / Cs

    @property
    def A_I(self):
        """Aromaticity index"""
        return self._calc_aromaticity_index()

    @property
    def A_I_mod(self):
        """Modified aromaticity index"""
        return self._calc_aromaticity_index_mod()

    @property
    def nosc(self):
        """Nominal oxidation state of carbon"""
        return self._calc_nosc()

    @property
    def dbe(self):
        return self._calc_dbe()

    @property
    def mz_nominal_calc(self):
        return int(self._calc_mz())

    @property
    def mz_error(self):
        return self._calc_assignment_mass_error()

    @property
    def mz_calc(self):
        return self._calc_mz()

    @property
    def protonated_mz(self):
        return self._protonated_mz(self.ion_charge)

    @property
    def radical_mz(self):
        return self._radical_mz(self.ion_charge)

    @property
    def neutral_mass(self):
        return self._neutral_mass()

    def adduct_mz(self, adduct_atom):
        """Get m/z of an adducted ion version of the molecular formula.

        Parameters
        ----------
        adduct_atom : str
            The adduct atom.

        Returns
        -------
        float
            The m/z value of the adducted ion version of the molecular formula.
        """
        return self._adduct_mz(adduct_atom, self.ion_charge)

    @property
    def ion_type(self):
        ion_type = self._d_molecular_formula.get(Labels.ion_type)
        if ion_type == Labels.protonated_de_ion:
            if self.ion_charge > 0:
                return Labels.protonated
            else:
                return Labels.de_protonated
        else:
            return ion_type

    @ion_type.setter
    def ion_type(self, ion_type):
        if ion_type in [
            Labels.protonated_de_ion,
            Labels.adduct_ion,
            Labels.radical_ion,
        ]:
            self._d_molecular_formula[Labels.ion_type] = ion_type
        else:
            raise TypeError(
                "Ion type can only be: 'DE_OR_PROTONATED', 'RADICAL' or  'ADDUCT', not %s"
                % ion_type
            )

    @property
    def ion_charge(self):
        return self._ion_charge

    @property
    def atoms(self):
        """Get the atoms in the molecular formula."""
        # if there is an adduct_atom, them reduce it from the atoms list
        if self.adduct_atom is None:
            return [
                key
                for key in self._d_molecular_formula.keys()
                if key != Labels.ion_type
            ]
        else:
            temp_dict = self._d_molecular_formula.copy()
            temp_dict[self.adduct_atom] -= 1
            return [
                key
                for key, val in temp_dict.items()
                if key != Labels.ion_type and val > 0
            ]

    @property
    def confidence_score(self):
        """Composite molecular formula assignment confidence score.

        Combines formula-level mass accuracy and isotopologue pattern
        agreement:

        ```
        CS = (w_err * m_err) + (w_iso * m_iso)
        ```

        where ``m_err`` is :attr:`average_mz_error_score`, ``m_iso`` is
        :attr:`isotopologue_similarity`, and the weights are
        ``molecular_search_settings.mz_error_score_weight`` (default 0.6)
        and ``isotopologue_score_weight`` (default 0.4).

        Returns
        -------
        float
            Score typically in about [0, 1]; higher values indicate better
            agreement with mass error and (when available) isotopologue
            evidence. Used for ranking when ``score_method`` is
            ``\"prob_score\"``.

        Notes
        -----
        Requires an associated parent mass spectral peak. Isotopologue
        terms are only informative when search was run with isotopologue
        detection enabled. See the ``corems.molecular_formula`` package
        documentation and Dewey et al., *Anal. Chem.* 2025, 97, 13031–13039
        (https://doi.org/10.1021/acs.analchem.4c06826).
        """
        if not self._confidence_score:
            self._confidence_score = self._calc_confidence_score()

        return self._confidence_score

    @property
    def isotopologue_similarity(self):
        """Isotopologue abundance pattern similarity for this formula.

        Compares expected natural-abundance isotopologue intensities
        (from dynamic-range-limited expansion) to abundances of peaks
        assigned as those partners. Implemented as a normalized Manhattan
        distance mapped to a similarity in about [0, 1]
        (higher is better).

        Returns
        -------
        float
            Similarity score. Returns 0.0 when no isotopologues are
            expected for this monoisotopic formula (or none were computed).

        Notes
        -----
        This is the ``m_iso`` term in :attr:`confidence_score`. See the
        ``corems.molecular_formula`` package documentation for details.
        """
        if not self._isotopologue_similarity:
            self._isotopologue_similarity = self._calc_isotopologue_confidence()

        return self._isotopologue_similarity

    @property
    def average_mz_error_score(self):
        """Mean mass-accuracy score over mono and expected isotopologues.

        Averages :attr:`mz_error_score` for the monoisotopic assignment and
        for each expected isotopologue formula that has a parent peak.
        Expected isotopologues without a matched peak contribute 0.

        Returns
        -------
        float
            Mean Gaussian mass-error score in about [0, 1].

        Notes
        -----
        This is the ``m_err`` term in :attr:`confidence_score`.
        """
        # includes the isotopologues

        if not self._mass_error_average_score:
            self._mass_error_average_score = self._calc_average_mz_score()

        return self._mass_error_average_score

    @property
    def mz_error_score(self):
        """Gaussian mass-accuracy score for this formula on its parent peak.

        Uses the assignment mass error ``delta`` (ppm) and a width ``sigma``
        equal to the parent peak's ``predicted_std`` when set, otherwise a
        fallback of 1.66 ppm:

        ```
        m = exp( -((delta - mu) ** 2) / (2 * (sigma ** 2)) )
        ```

        with mean ``mu = 0`` (calibrated spectrum with near-zero mean error
        assumed).

        Returns
        -------
        float
            Score typically in about 0 to 1; exact mass match is 1,
            and larger ppm error yields values closer to 0.

        Notes
        -----
        Formula-level ranking uses :attr:`average_mz_error_score`, which
        folds in expected isotopologue peaks as well.
        """
        if not self._mz_error_score:
            self._mz_error_score = self._calc_mz_confidence()

        return self._mz_error_score

    @property
    def kmd(self):
        return self._kmd

    @property
    def kendrick_mass(self):
        return self._kendrick_mass

    @property
    def knm(self):
        return self._nominal_km

    def change_kendrick_base(self, kendrick_dict_base):
        """Change the Kendrick base.

        Parameters
        ----------
        kendrick_dict_base : dict
            The Kendrick base dictionary. Ex: {"C": 1, "H": 2}
        """
        self._kmd, self._kendrick_mass, self._nominal_km = self._calc_kmd(
            kendrick_dict_base
        )

    def isotopologues(self, min_abundance, current_mono_abundance, dynamic_range):
        """Calculate the isotopologues for a given molecular formula.

        Parameters
        ----------
        min_abundance : float
            The minimum abundance.
        current_mono_abundance : float
            The current monoisotopic abundance.
        dynamic_range : float
            The dynamic range.

        Yields
        ------
        MolecularFormulaIsotopologue
            The molecular formula isotopologue.

        Notes
        -----
        This calculation ignores the hydrogen isotopes.
        """
        isotopologues = []
        for mf in self._cal_isotopologues(
            self._d_molecular_formula,
            min_abundance,
            current_mono_abundance,
            dynamic_range,
        ):
            isotopologues.append(mf)

        # To account for differences in how the isotopologue outputs are sorted between IsoSpec versions.
        sorted_isotopologues = sorted(isotopologues, key=lambda mf: mf[1], reverse=True)

        for mf in sorted_isotopologues:
            yield MolecularFormulaIsotopologue(
                *mf,
                current_mono_abundance,
                self.ion_charge,
                ion_type=self.ion_type,
                adduct_atom=self.adduct_atom,
            )

    def atoms_qnt(self, atom):
        """Get the atom quantity of a specific atom in the molecular formula."""
        if atom in self._d_molecular_formula:
            return self._d_molecular_formula.get(atom)
        else:
            raise Warning(
                "Could not find %s in this Molecular Formula object" % str(atom)
            )

    def atoms_symbol(self, atom):
        """Get the atom symbol without the mass number."""
        return "".join([i for i in atom if not i.isdigit()])

    @property
    def string(self):
        """Returns the molecular formula as a string."""
        if self._d_molecular_formula:
            if self.adduct_atom is None:
                mol_form_dict = self._d_molecular_formula
            else:
                mol_form_dict = self._d_molecular_formula.copy()
                if self.adduct_atom not in mol_form_dict.keys():
                    raise Exception("Adduct atom not found in molecular formula dict")
                mol_form_dict[self.adduct_atom] -= 1
                mol_form_dict = {
                    key: val for key, val in mol_form_dict.items() if val != 0
                }
            formula_srt = ""
            for atom in Atoms.atoms_order:
                if atom in mol_form_dict.keys():
                    formula_srt += atom + str(int(mol_form_dict.get(atom))) + " "
            return formula_srt.strip()

        else:
            raise Exception("Molecular formula identification not performed yet")

    @property
    def string_formated(self):
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

        if self._d_molecular_formula:
            formula_srt = ""
            for atom in Atoms.atoms_order:
                if atom in self.to_dict().keys():
                    formula_srt += atom.translate(SUP) + str(
                        int(self.to_dict().get(atom))
                    ).translate(SUB)
            return formula_srt

        else:
            raise Exception("Molecular formula identification not performed yet")

    def to_dict(self):
        """Returns the molecular formula as a dictionary.

        Returns
        -------
        dict
            The molecular formula as a dictionary.
        """
        return self._d_molecular_formula

    def to_list(self):
        """Returns the molecular formula as a list.

        Returns
        -------
        list
            The molecular formula as a list.

        Raises
        ------
        Exception
            If the molecular formula identification was not performed yet.
        """
        # TODO ensure self._d_molecular_formula is a orderedDict

        if self._d_molecular_formula:
            formula_list = []

            for atom, atom_number in self._d_molecular_formula.items():
                if atom != Labels.ion_type:
                    formula_list.append(atom)
                    formula_list.append(atom_number)

            return formula_list
        else:
            raise Exception("Molecular formula identification not performed yet")

    @property
    def class_label(self):
        if self._d_molecular_formula:
            formulalist = self.to_list()
            classstring = ""

            for each in range(0, len(formulalist), 2):
                if (
                    formulalist[each] != "C"
                    and formulalist[each] != "H"
                    and formulalist[each] != "HC"
                ):
                    classstring = (
                        classstring
                        + str(formulalist[each])
                        + str(formulalist[each + 1])
                        + " "
                    )

            if classstring == "":
                classstring = "HC"

            classstring = classstring.strip()

            if self._d_molecular_formula.get(Labels.ion_type) == Labels.radical_ion:
                return classstring + " -R"

            # elif self._d_molecular_formula.get(Labels.ion_type) == Labels.adduct_ion:

            #    return classstring + ' -A'

            else:
                return classstring

            #'dict, tuple or string'

        else:
            raise Exception("Molecular formula identification not performed yet")

    @property
    def class_dict(self):
        if self._d_molecular_formula:
            class_dict = {}

            for atom, qnt in self._d_molecular_formula.items():
                if atom != Labels.ion_type and atom != "C" and atom != "H":
                    class_dict[atom] = qnt

            return class_dict

        raise Exception("Molecular formula identification not performed yet")


class MolecularFormulaIsotopologue(MolecularFormulaBase):
    """Class for representing a molecular formula isotopologue.

    Parameters
    ----------
    _d_molecular_formula : dict
        The molecular formula as a dictionary.
    prob_ratio : float
        The probability ratio.
    mono_abundance : float
        The monoisotopic abundance.
    ion_charge : int
        The ion charge.
    mspeak_parent : object, optional
        The parent mass spectrum peak object instance. Defaults to None.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.

    Attributes
    ----------
    prob_ratio : float
        The probability ratio.
    abundance_calc : float
        The calculated abundance.
    area_error : float
        The area error.
    abundance_error : float
        The abundance error.
    is_isotopologue : bool
        The isotopologue flag. Defaults to True.
    mspeak_index_mono_isotopic : int
        The index of the monoisotopic peak in the mass spectrum peak list. Defaults to None.
    mono_isotopic_formula_index : int
        The index of the monoisotopic formula in the molecular formula list. Defaults to None.
    """

    def __init__(
        self,
        _d_molecular_formula,
        prob_ratio,
        mono_abundance,
        ion_charge,
        mspeak_parent=None,
        ion_type=None,
        adduct_atom=None,
    ):
        if ion_type is None:
            # check if ion type or adduct_atom is in the molecular formula dict
            if Labels.ion_type in _d_molecular_formula:
                ion_type = _d_molecular_formula.get(Labels.ion_type)
            else:
                ion_type = None
        else:
            ion_type = Labels.ion_type_translate.get(ion_type)

        if ion_type == Labels.adduct_ion:
            adduct_atom_int = None
            if adduct_atom in _d_molecular_formula.keys():
                adduct_atom_int = adduct_atom
            else:
                # Check to see if adduct_atom should actually be an isotope of the adduct atom
                for adduct_iso in Atoms.isotopes.get(adduct_atom)[1]:
                    if adduct_iso in _d_molecular_formula.keys():
                        adduct_atom_int = adduct_iso
            adduct_atom = adduct_atom_int
            if adduct_atom is None:
                raise Exception("adduct_atom is required for adduct ion")
            _d_molecular_formula[adduct_atom] -= 1
            _d_molecular_formula = {
                key: val for key, val in _d_molecular_formula.items() if val != 0
            }

        super().__init__(
            molecular_formula=_d_molecular_formula,
            ion_charge=ion_charge,
            ion_type=ion_type,
            adduct_atom=adduct_atom,
        )
        # prob_ratio is relative to the monoisotopic peak p_isotopologue/p_mono_isotopic

        self.prob_ratio = prob_ratio

        self.abundance_calc = mono_abundance * prob_ratio

        self.is_isotopologue = True

        self.mspeak_index_mono_isotopic = None

        self.mono_isotopic_formula_index = None
        # parent mass spectrum peak obj instance
        self._mspeak_parent = mspeak_parent

    @property
    def area_error(self):
        return self._calc_area_error()

    @property
    def abundance_error(self):
        return self._calc_abundance_error()


class LCMSLibRefMolecularFormula:
    """Removed LC-MS library reference molecular formula type.

    This class is no longer supported. Canonical LC-MS MS1 formula assignment
    uses :class:`MolecularFormula` via combinatorial / SQL search. Construction
    raises :class:`NotImplementedError`. The name is retained as a stub until
    the next major release.

    Parameters
    ----------
    *args
        Ignored; retained only for call-site compatibility.
    **kwargs
        Ignored; retained only for call-site compatibility.
    """

    def __init__(self, *args, **kwargs) -> None:
        raise NotImplementedError(
            "LCMSLibRefMolecularFormula is no longer supported. "
            "Use MolecularFormula for formula objects. "
            "This stub will be removed in the next major release."
        )


class MolecularFormula(MolecularFormulaBase):
    """General class for representing a molecular formula.

    Parameters
    ----------
    molecular_formula : dict, list, str
        The molecular formula.
    ion_charge : int
        The ion charge.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.
    mspeak_parent : object, optional
        The parent mass spectrum peak object instance. Defaults to None.
    external_mz : float, optional
        The external m/z value. Defaults to False.
    """

    def __init__(
        self,
        molecular_formula,
        ion_charge,
        ion_type=None,
        adduct_atom=None,
        mspeak_parent=None,
        external_mz=False,
    ):
        super().__init__(
            molecular_formula,
            ion_charge,
            ion_type=ion_type,
            adduct_atom=adduct_atom,
            mspeak_parent=mspeak_parent,
            external_mz=external_mz,
        )
