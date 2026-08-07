"""Tests for retired LC-MS library formula stubs."""

import pytest

from corems.molecular_formula.factory.MolecularFormulaFactory import (
    LCMSLibRefMolecularFormula,
)
from corems.molecular_formula.input.masslist_ref import ImportMassListRef


def test_lcms_lib_ref_molecular_formula_raises():
    with pytest.raises(NotImplementedError, match="no longer supported"):
        LCMSLibRefMolecularFormula({"C": 1, "H": 4}, ion_charge=1)


def test_from_lcms_lib_file_raises(tmp_path):
    # Constructor requires an existing path; body must still raise without reading.
    ref = tmp_path / "unused_lcms_lib.csv"
    ref.write_text("placeholder\n", encoding="utf-8")
    importer = ImportMassListRef(str(ref))
    with pytest.raises(NotImplementedError, match="from_lcms_lib_file"):
        importer.from_lcms_lib_file(ion_charge=1, ion_types=["[M+H]+"])
