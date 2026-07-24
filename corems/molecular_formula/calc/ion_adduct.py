"""Ion-type / adduct helpers for neutral formulas and precursor m/z."""

from __future__ import annotations

import re

from corems.encapsulation.constant import Atoms
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula

# adduct : [atoms to add, atoms to subtract when calculating formula of ion]
ion_type_dict = {
    "M+": [{}, {}],
    "[M]+": [{}, {}],
    "protonated": [{"H": 1}, {}],
    "[M+H]+": [{"H": 1}, {}],
    "[M+NH4]+": [{"N": 1, "H": 4}, {}],  # ammonium
    "[M+Na]+": [{"Na": 1}, {}],
    "[M+K]+": [{"K": 1}, {}],
    "[M+2Na+Cl]+": [{"Na": 2, "Cl": 1}, {}],
    "[M+2Na-H]+": [{"Na": 2}, {"H": 1}],
    "[M+C2H3Na2O2]+": [{"C": 2, "H": 3, "Na": 2, "O": 2}, {}],
    "[M+C4H10N3]+": [{"C": 4, "H": 10, "N": 3}, {}],
    "[M+NH4+ACN]+": [{"C": 2, "H": 7, "N": 2}, {}],
    "[M+H-H2O]+": [{}, {"H": 1, "O": 1}],
    "de-protonated": [{}, {"H": 1}],
    "[M-H]-": [{}, {"H": 1}],
    "[M+Cl]-": [{"Cl": 1}, {}],
    "[M+HCOO]-": [{"C": 1, "H": 1, "O": 2}, {}],  # formate
    "[M+CH3COO]-": [{"C": 2, "H": 3, "O": 2}, {}],  # acetate
    "[M+2NaAc+Cl]-": [{"Na": 2, "C": 2, "H": 3, "O": 2, "Cl": 1}, {}],
    "[M+K-2H]-": [{"K": 1}, {"H": 2}],
    "[M+Na-2H]-": [{"Na": 1}, {"H": 2}],
}

# Common alternate adduct strings → keys in ion_type_dict
ADDUCT_ALIASES = {
    "[M+HCOOH-H]-": "[M+HCOO]-",
    "[M+CH3COOH-H]-": "[M+CH3COO]-",
    "[M+FA-H]-": "[M+HCOO]-",
    "[M+AcOH-H]-": "[M+CH3COO]-",
}


def normalize_ion_type(ion_type: str) -> str:
    """Map alternate adduct labels onto ``ion_type_dict`` keys."""
    return ADDUCT_ALIASES.get(ion_type, ion_type)


def charge_from_adduct(adduct: str) -> int:
    """Parse integer charge from an adduct string (e.g. ``[M+H]+``, ``[M+2H]2+``)."""
    a = (adduct or "").strip()
    m = re.search(r"(\d+)([+-])$", a)
    if m:
        return int(m.group(1)) * (1 if m.group(2) == "+" else -1)
    if a.endswith("+"):
        return 1
    if a.endswith("-"):
        return -1
    raise ValueError(f"cannot parse charge from adduct {adduct!r}")


def get_ion_formula(neutral_formula, ion_type):
    """From a neutral formula and an ion type, return the formula of the ion.

    Parameters
    ----------
    neutral_formula : str
        Neutral formula as MolecularFormula-style (``'C2 H4 O2'``) or compact
        (``'C2H4O2'``; isotopes only handled in the spaced form).
    ion_type : str
        Ion type / adduct key (see ``ion_type_dict``), or a known alias.

    Returns
    -------
    str or None
        Ion formula string (e.g. ``'C2 H5 O2'``), or None if
        ``neutral_formula`` is not a string.
    """
    if not isinstance(neutral_formula, str):
        return None

    ion_type = normalize_ion_type(ion_type)
    if ion_type not in ion_type_dict:
        raise KeyError(f"unsupported ion type {ion_type!r}")

    if re.search(r"\s", neutral_formula):
        formula_obj = MolecularFormula(neutral_formula, ion_charge=0)
    else:
        form_pre = re.sub(r"([A-Z])", r" \1", neutral_formula)[1:]
        elements = [re.findall(r"[A-Z][a-z]*", x) for x in form_pre.split()]
        counts = [re.findall(r"\d+", x) for x in form_pre.split()]
        formula_obj = MolecularFormula(
            dict(
                zip(
                    [x[0] for x in elements],
                    [int(x[0]) if x else 1 for x in counts],
                )
            ),
            ion_charge=0,
        )
    neutral_formula_dict = formula_obj.to_dict().copy()

    adduct_add_dict = ion_type_dict[ion_type][0]
    for key in adduct_add_dict:
        if key in neutral_formula_dict:
            neutral_formula_dict[key] += adduct_add_dict[key]
        else:
            neutral_formula_dict[key] = adduct_add_dict[key]

    adduct_subtract = ion_type_dict[ion_type][1]
    for key in adduct_subtract:
        neutral_formula_dict[key] -= adduct_subtract[key]

    return MolecularFormula(neutral_formula_dict, ion_charge=0).string


def precursor_mz_from_formula(neutral_formula, ion_type, charge=None):
    """Calculated precursor m/z from neutral formula and adduct / ion type.

    Parameters
    ----------
    neutral_formula : str
        Neutral molecular formula.
    ion_type : str
        Adduct or ion type (``ion_type_dict`` key or alias).
    charge : int, optional
        Ion charge. If omitted, inferred from ``ion_type``.

    Returns
    -------
    float
        Precursor m/z of the ion.

    Raises
    ------
    KeyError
        If ``ion_type`` is not supported.
    ValueError
        If the ion formula cannot be built or charge cannot be inferred.
    """
    ion_form = get_ion_formula(neutral_formula, ion_type)
    if ion_form is None:
        raise ValueError(
            f"could not build ion formula for {neutral_formula!r} {ion_type!r}"
        )
    if charge is None:
        charge = charge_from_adduct(normalize_ion_type(ion_type))
    if charge == 0:
        raise ValueError("charge must be non-zero for precursor m/z")
    mf = MolecularFormula(ion_form, ion_charge=0)
    return float(
        (mf.neutral_mass + (charge * -1 * Atoms.electron_mass)) / abs(charge)
    )
