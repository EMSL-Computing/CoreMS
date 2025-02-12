__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"

import sys

from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.encapsulation.constant import Labels


def test_molecular_formula_from_dict():
    """Test the parsing of a molecular formula string and the calculation of isotopologues"""
    formula_dict = {"C": 10, "H": 0, "O": 10, "Cl": 2, Labels.ion_type: "RADICAL"}
    ion_charge = 1
    formula_obj = MolecularFormula(
        molecular_formula=formula_dict, ion_charge=ion_charge
    )
    assert formula_obj.ion_type == "RADICAL"
    assert round(formula_obj.mz_calc, 2) == round(349.886303060457, 2)
    assert formula_obj.kmd == -50
    assert round(formula_obj.kendrick_mass, 2) == round(349.4956152007638, 2)
    assert formula_obj.knm == 349
    assert formula_obj.class_label == "O10 Cl2 -R"
    assert formula_obj.atoms_qnt("C") == 10
    assert formula_obj.atoms_symbol("13C") == "C"
    assert formula_obj.string == "C10 O10 Cl2"

    # Create isotopologues of the formula_obj
    isotopologues = list(formula_obj.isotopologues(0.01, 1, 500))
    assert round(isotopologues[0].mz_calc, 2) == round(351.883352980637, 2)
    assert round(isotopologues[0].prob_ratio, 2) == round(0.6399334750069298, 2)
    assert isotopologues[0].string == "C10 O10 Cl1 37Cl1"


def test_molecular_formula_from_string():
    """Test the parsing of a molecular formula string and the dealing with neutral mass"""
    ion_charge = 1
    formula_str = "C10 H21 N1"
    formula_obj = MolecularFormula(formula_str, ion_charge)
    assert formula_obj.string == "C10 H21 N1"
    assert formula_obj.ion_type is None
    # This returns a neutral mass since the ion type is not set and therefore the ion_type is interpreted as None - is that expected?
    assert round(formula_obj.mz_calc, 2) == round(formula_obj.neutral_mass, 2)


def test_molecular_formula_adducts():
    """Test the parsing of a molecular formula string with adducts and the calculation of isotopologues"""
    formula_obj = MolecularFormula(
        {"C": 6, "H": 10, "O": 6}, ion_charge=-1, ion_type="ADDUCT", adduct_atom="Cl"
    )

    isotopologues = list(formula_obj.isotopologues(0.05, 1, dynamic_range=1000))

    assert round(formula_obj.mz_calc, 2) == round(213.01713930162907, 2)
    assert round(isotopologues[0].mz_calc, 2) == round(215.01418922162907, 2)
    assert round(isotopologues[0].prob_ratio, 2) == round(0.3199577613516368, 2)
    assert isotopologues[0].string == "C6 H10 O6"
    assert isotopologues[0].adduct_atom == "37Cl"
