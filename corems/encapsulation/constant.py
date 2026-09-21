

from corems.encapsulation import nist_atoms


class Labels:  # pragma: no cover
    """Class for Labels used in CoreMS

    These labels are used to define:
    * types of columns in plaintext data inputs,
    * types of data/mass spectra
    * types of assignment for ions

    """

    mz = "m/z"
    abundance = "Peak Height"
    rp = "Resolving Power"
    s2n = "S/N"

    label = "label"
    bruker_profile = "Bruker_Profile"
    thermo_profile = "Thermo_Profile"
    simulated_profile = "Simulated Profile"
    booster_profile = "Booster Profile"
    bruker_frequency = "Bruker_Frequency"
    midas_frequency = "Midas_Frequency"
    thermo_centroid = "Thermo_Centroid"
    corems_centroid = "CoreMS_Centroid"
    gcms_centroid = "Thermo_Centroid"

    unassigned = "unassigned"

    radical_ion = "RADICAL"
    protonated_de_ion = "DE_OR_PROTONATED"
    protonated = "protonated"
    de_protonated = "de-protonated"
    adduct_ion = "ADDUCT"
    neutral = "neutral"
    ion_type = "IonType"

    ion_type_translate = {
        "protonated": "DE_OR_PROTONATED",
        "de-protonated": "DE_OR_PROTONATED",
        "radical": "RADICAL",
        "adduct": "ADDUCT",
        "ADDUCT": "ADDUCT",
    }


class Atoms:  # pragma: no cover
    """Class for Atoms in CoreMS

    Public API for exact masses, isotopic abundances, rare-isotope lists,
    formula-string order, English names, and covalence. Do not use
    ``nist_atoms.isotopes`` as ``Atoms.isotopes``: the generated module stores
    a list of rare keys (or ``[None]``); this class stores
    ``[English name, rare keys]``.

    IUPAC monoisotopic mass uses the most abundant isotope of each element.
    Bare element symbols are the most abundant isotope (C, Cl). Rare isotopes
    use mass-number keys (13C, 37Cl). Hydrogen uses H and D. Masses and
    abundances may also be looked up by nuclide key (12C, 1H). Formula strings
    and mass-list columns use ``atoms_order`` (canonical keys only).

    Masses, abundances, rare lists, and ``atoms_order`` come from a pinned NIST
    snapshot (``corems.encapsulation.nist_atoms``). English names and covalence
    are hand-maintained. ``Atoms.isotopes`` membership is not enough to search
    an element: ``usedAtoms`` still needs a valence in
    ``used_atom_valences`` (usually from ``atoms_covalence``). Pin deltas
    and breaking notes are in ``tools/nist_atoms/CHANGES.md``.

    References
    ----------

    1. NIST Atomic Weights and Isotopic Compositions (Coursey et al., version 4.1)
    https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses

    """

    electron_mass = 0.0005_485_799_090_65  # NIST value

    atomic_masses = nist_atoms.atomic_masses
    isotopic_abundance = nist_atoms.isotopic_abundance
    atoms_order = nist_atoms.atoms_order

    # Not a NIST field (the dump has no name column).
    element_names = {
        "H": "Hydrogen",
        "He": "Helium",
        "Li": "Lithium",
        "Be": "Beryllium",
        "B": "Boron",
        "C": "Carbon",
        "N": "Nitrogen",
        "O": "Oxygen",
        "F": "Fluorine",
        "Ne": "Neon",
        "Na": "Sodium",
        "Mg": "Magnesium",
        "Al": "Aluminum",
        "Si": "Silicon",
        "P": "Phosphorus",
        "S": "Sulfur",
        "Cl": "Chlorine",
        "Ar": "Argon",
        "K": "Potassium",
        "Ca": "Calcium",
        "Sc": "Scandium",
        "Ti": "Titanium",
        "V": "Vanadium",
        "Cr": "Chromium",
        "Mn": "Manganese",
        "Fe": "Iron",
        "Co": "Cobalt",
        "Ni": "Nickel",
        "Cu": "Copper",
        "Zn": "Zinc",
        "Ga": "Gallium",
        "Ge": "Germanium",
        "As": "Arsenic",
        "Se": "Selenium",
        "Br": "Bromine",
        "Kr": "Krypton",
        "Rb": "Rubidium",
        "Sr": "Strontium",
        "Y": "Yttrium",
        "Zr": "Zirconium",
        "Nb": "Niobium",
        "Mo": "Molybdenum",
        "Ru": "Ruthenium",
        "Rh": "Rhodium",
        "Pd": "Palladium",
        "Ag": "Silver",
        "Cd": "Cadmium",
        "In": "Indium",
        "Sn": "Tin",
        "Sb": "Antimony",
        "Te": "Tellurium",
        "I": "Iodine",
        "Xe": "Xenon",
        "Cs": "Cesium",
        "Ba": "Barium",
        "La": "Lanthanum",
        "Ce": "Cerium",
        "Pr": "Praseodymium",
        "Nd": "Neodymium",
        "Sm": "Samarium",
        "Eu": "Europium",
        "Gd": "Gadolinium",
        "Tb": "Terbium",
        "Dy": "Dysprosium",
        "Ho": "Holmium",
        "Er": "Erbium",
        "Tm": "Thulium",
        "Yb": "Ytterbium",
        "Lu": "Lutetium",
        "Hf": "Hafnium",
        "Ta": "Tantalum",
        "W": "Tungsten",
        "Re": "Rhenium",
        "Os": "Osmium",
        "Ir": "Iridium",
        "Pt": "Platinum",
        "Au": "Gold",
        "Hg": "Mercury",
        "Tl": "Thallium",
        "Pb": "Lead",
        "Bi": "Bismuth",
        "Th": "Thorium",
        "Pa": "Protactinium",
        "U": "Uranium",
    }

    isotopes = {}
    for symbol, rares in nist_atoms.isotopes.items():
        isotopes[symbol] = [element_names.get(symbol, symbol), list(rares)]
    del symbol, rares

    atoms_covalence = {
        "C": (4),
        "13C": (4),
        "N": (3),
        "O": (2),
        "S": (2),
        "H": (1),
        "F": (1, 0),
        "Cl": (1, 0),
        "Br": (1, 0),
        "I": (1, 0),
        "Li": (1, 0),
        "Na": (1, 0),
        "K": (1, 0),
        "Rb": (1),
        "Cs": (1),
        "B": (4, 3, 2, 1),
        "In": (3, 2, 1),
        "Al": (3, 1, 2),
        "P": (3, 5, 4, 2, 1),
        "Ga": (3, 1, 2),
        "Mg": (2, 1),
        "Be": (2, 1),
        "Ca": (2, 1),
        "Sr": (2, 1),
        "Ba": (2),
        "V": (5, 4, 3, 2, 1),
        "Fe": (3, 2, 4, 5, 6),
        "Si": (4, 3, 2),
        "Sc": (3, 2, 1),
        "Ti": (4, 3, 2, 1),
        "Cr": (1, 2, 3, 4, 5, 6),
        "Mn": (1, 2, 3, 4, 5, 6, 7),
        "Co": (1, 2, 3, 4, 5),
        "Ni": (1, 2, 3, 4),
        "Cu": (2, 1, 3, 4),
        "Zn": (2, 1),
        "Ge": (4, 3, 2, 1),
        "As": (5, 3, 2, 1),
        "Se": (6, 4, 2, 1),
        "Y": (3, 2, 1),
        "Zr": (4, 3, 2, 1),
        "Nb": (5, 4, 3, 2, 1),
        "Mo": (6, 5, 4, 3, 2, 1),
        "Ru": (8, 7, 6, 5, 4, 3, 2, 1),
        "Rh": (6, 5, 4, 3, 2, 1),
        "Pd": (4, 2, 1),
        "Ag": (0, 1, 2, 3, 4),
        "Cd": (2, 1),
        "Sn": (4, 2),
        "Sb": (5, 3),
        "Te": (6, 5, 4, 2),
        "La": (3, 2),
        "Hf": (4, 3, 2),
        "Ta": (5, 4, 3, 2),
        "W": (6, 5, 4, 3, 2, 1),
        "Re": (4, 7, 6, 5, 3, 2, 1),
        "Os": (4, 8, 7, 6, 5, 3, 2, 1),
        "Ir": (4, 8, 6, 5, 3, 2, 1),
        "Pt": (4, 6, 5, 3, 2, 1),
        "Au": (3, 5, 2, 1),
        "Hg": (1, 2, 4),
        "Tl": (3, 1),
        "Pb": (4, 2),
        "Bi": (3, 1, 5),
    }

ION_TYPE_DICT = {
    'M+': {
        "add": {},
        "sub": {},
        "polarity": 'positive',
    },
    '[M]+': {
        "add": {},
        "sub": {},
        "polarity": 'positive',
    },
    'protonated': {
        "add": {'H': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+H]+': {
        "add": {'H': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+2H]2+': {
        "add": {'H': 2},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+3H]3+': {
        "add": {'H': 3},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+NH4]+': {
        "add": {'N': 1, 'H': 4},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+Na]+': {
        "add": {'Na': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+K]+': {
        "add": {'K': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+2Na]2+': {
        "add": {'Na': 2},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+H+Na]2+': {
        "add": {'H': 1, 'Na': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+H+K]2+': {
        "add": {'H': 1, 'K': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+2Na+Cl]+': {
        "add": {'Na': 2, 'Cl': 1},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+2Na-H]+': {
        "add": {'Na': 2},
        "sub": {'H': 1},
        "polarity": 'positive',
    },
    '[M-H+2Na]+': {
        "add": {'Na': 2},
        "sub": {'H': 1},
        "polarity": 'positive',
    },
    '[M+C2H3Na2O2]+': {
        "add": {'C': 2, 'H': 3, 'Na': 2, 'O': 2},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+C4H10N3]+': {
        "add": {'C': 4, 'H': 10, 'N': 3},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+NH4+ACN]+': {
        "add": {'C': 2, 'H': 7, 'N': 2},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+H-H2O]+': {
        "add": {},
        "sub": {'H': 1, 'O': 1},
        "polarity": 'positive',
    },
    '[M+H-2H2O]+': {
        "add": {},
        "sub": {'H': 3, 'O': 2},
        "polarity": 'positive',
    },
    '[M+H-NH3]+': {
        "add": {},
        "sub": {'N': 1, 'H': 2},
        "polarity": 'positive',
    },
    '[M+2H-NH3]2+': {
        "add": {},
        "sub": {'N': 1, 'H': 1},
        "polarity": 'positive',
    },
    '[M+2H-H2O]2+': {
        "add": {},
        "sub": {'O': 1},
        "polarity": 'positive',
    },
    '[M+NH4-H2O]+': {
        "add": {'N': 1, 'H': 2},
        "sub": {},
        "polarity": 'positive',
    },
    '[M+H+H2O]+': {
        "add": {'H': 3, 'O': 1},
        "sub": {},
        "polarity": 'positive',
    },
    'de-protonated': {
        "add": {},
        "sub": {'H': 1},
        "polarity": 'negative',
    },
    '[M-H]-': {
        "add": {},
        "sub": {'H': 1},
        "polarity": 'negative',
    },
    '[M-2H]2-': {
        "add": {},
        "sub": {'H': 2},
        "polarity": 'negative',
    },
    '[M-H-H2O]-': {
        "add": {},
        "sub": {'H': 3, 'O': 1},
        "polarity": 'negative',
    },
    '[M-H+H2O]-': {
        "add": {'H': 1, 'O': 1},
        "sub": {},
        "polarity": 'negative',
    },
    '[M+Cl]-': {
        "add": {'Cl': 1},
        "sub": {},
        "polarity": 'negative',
    },
    '[M+HCOO]-': {
        "add": {'C': 1, 'H': 1, 'O': 2},
        "sub": {},
        "polarity": 'negative',
    },
    '[M+CH3COO]-': {
        "add": {'C': 2, 'H': 3, 'O': 2},
        "sub": {},
        "polarity": 'negative',
    },
    '[M+2NaAc+Cl]-': {
        "add": {'Na': 2, 'C': 2, 'H': 3, 'O': 2, 'Cl': 1},
        "sub": {},
        "polarity": 'negative',
    },
    '[M+K-2H]-': {
        "add": {'K': 1},
        "sub": {'H': 2},
        "polarity": 'negative',
    },
    '[M+Na-2H]-': {
        "add": {'Na': 1},
        "sub": {'H': 2},
        "polarity": 'negative',
    },
}

ADDUCT_ALIASES = {
    '[M+HCOOH-H]-': '[M+HCOO]-',
    '[M+CH3COOH-H]-': '[M+CH3COO]-',
    '[M+FA-H]-': '[M+HCOO]-',
    '[M+AcOH-H]-': '[M+CH3COO]-',
}
