"""Tests for NIST-pinned Atoms tables (masses, abundances, isotope lists)."""

import importlib.util
import sys
from pathlib import Path

import pytest

from corems.encapsulation.constant import Atoms
from corems.encapsulation.factory.processingSetting import validate_used_atoms_keys


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATE_PY = REPO_ROOT / "tools" / "nist_atoms" / "generate.py"


def _load_generate_module():
    spec = importlib.util.spec_from_file_location("nist_atoms_generate", GENERATE_PY)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_pin_metadata_present():
    from corems.encapsulation import nist_atoms

    assert "4.1" in nist_atoms.NIST_TABLE_ID
    assert "nist.gov" in nist_atoms.NIST_TABLE_ID
    assert nist_atoms.NIST_DUMP_URL.startswith("https://physics.nist.gov/")
    assert nist_atoms.NIST_RETRIEVED


def test_dual_lookup_most_abundant_and_nuclide():
    assert Atoms.atomic_masses["C"] == Atoms.atomic_masses["12C"]
    assert Atoms.atomic_masses["C"] == pytest.approx(12.0)
    assert Atoms.atomic_masses["13C"] != Atoms.atomic_masses["C"]
    assert Atoms.atomic_masses["Cl"] == Atoms.atomic_masses["35Cl"]
    assert Atoms.atomic_masses["37Cl"] != Atoms.atomic_masses["Cl"]
    assert Atoms.isotopic_abundance["C"] == Atoms.isotopic_abundance["12C"]
    assert Atoms.isotopic_abundance["13C"] < Atoms.isotopic_abundance["C"]


def test_hydrogen_canonical_and_nuclide_aliases():
    assert Atoms.atomic_masses["H"] == Atoms.atomic_masses["1H"]
    assert Atoms.atomic_masses["D"] == Atoms.atomic_masses["2H"]
    assert "H" in Atoms.atoms_order
    assert "D" in Atoms.atoms_order
    assert "1H" not in Atoms.atoms_order
    assert "2H" not in Atoms.atoms_order
    assert "T" not in Atoms.atomic_masses
    assert "3H" not in Atoms.atomic_masses


def test_most_abundant_nuclide_keys_not_in_atoms_order():
    assert "12C" not in Atoms.atoms_order
    assert "35Cl" not in Atoms.atoms_order
    assert "C" in Atoms.atoms_order
    assert "13C" in Atoms.atoms_order
    assert "37Cl" in Atoms.atoms_order


def test_atoms_table_key_coverage():
    """Canonical tables share keys; nuclide aliases may exist only on masses."""
    masses = set(Atoms.atomic_masses)
    abund = set(Atoms.isotopic_abundance)
    order = set(Atoms.atoms_order)
    names = set(Atoms.element_names)
    iso = set(Atoms.isotopes)
    cov = set(Atoms.atoms_covalence)

    assert masses == abund, sorted(masses.symmetric_difference(abund))
    assert names == iso, sorted(names.symmetric_difference(iso))
    assert iso <= masses, sorted(iso - masses)
    assert order <= masses, sorted(order - masses)
    assert order <= abund, sorted(order - abund)
    assert cov <= masses, sorted(cov - masses)

    rares = set()
    for symbol, (_name, rare_keys) in Atoms.isotopes.items():
        assert symbol in order
        for rare in rare_keys:
            if rare is None:
                continue
            rares.add(rare)
            assert rare in masses, rare
            assert rare in order, rare
    assert iso | rares == order, sorted((iso | rares).symmetric_difference(order))


def test_boron_and_used_atoms_elements_are_assignable():
    assert "B" in Atoms.atomic_masses
    assert "B" in Atoms.isotopes
    validate_used_atoms_keys({"C": (1, 10), "H": (1, 20), "B": (0, 2), "Cl": (0, 4)})


def test_empty_isotopic_composition_excluded():
    assert "T" not in Atoms.atomic_masses
    assert "14C" not in Atoms.atomic_masses
    assert "14C" not in Atoms.isotopes["C"][1]
    assert "Tc" not in Atoms.isotopes
    assert "Og" not in Atoms.isotopes
    assert "Pm" not in Atoms.isotopes


def test_monoisotopic_rare_list_is_none_sentinel():
    assert Atoms.isotopes["F"][1] == [None]
    assert Atoms.isotopes["P"][1] == [None]
    assert Atoms.isotopes["Na"][1] == [None]


def test_carbon_rare_isotopes_complete_and_mass_ordered():
    rares = Atoms.isotopes["C"][1]
    assert rares == ["13C"]
    assert Atoms.isotopes["C"][0] == "Carbon"
    assert Atoms.element_names["C"] == "Carbon"


def test_nist_module_isotopes_are_rare_keys_only():
    from corems.encapsulation import nist_atoms

    assert nist_atoms.isotopes["C"] == ["13C"]
    assert nist_atoms.isotopes["F"] == [None]
    assert all(isinstance(v, list) for v in nist_atoms.isotopes.values())
    for rares in nist_atoms.isotopes.values():
        if rares == [None]:
            continue
        assert all(isinstance(item, str) and item != "Carbon" for item in rares)


def test_cadmium_bare_symbol_is_most_abundant_114cd():
    assert Atoms.atomic_masses["Cd"] == Atoms.atomic_masses["114Cd"]
    assert "112Cd" in Atoms.isotopes["Cd"][1]
    assert "114Cd" not in Atoms.isotopes["Cd"][1]
    assert "114Cd" not in Atoms.atoms_order


@pytest.mark.skipif(not GENERATE_PY.is_file(), reason="generator not in this install")
def test_nist_atoms_module_matches_vendored_dump():
    gen = _load_generate_module()
    gen.check_committed_module()


@pytest.mark.skipif(not GENERATE_PY.is_file(), reason="generator not in this install")
def test_records_body_strips_html_and_comments():
    gen = _load_generate_module()
    raw = (
        "<html><pre>\n"
        "title\n"
        "# comment\n"
        "Atomic Number = 1\n"
        "Atomic Symbol = H\n"
        "Notes = &nbsp;\n"
        "</pre></html>\n"
    )
    body = gen.records_body(raw)
    assert body.startswith("Atomic Number = 1\n")
    assert "&nbsp;" not in body
    assert "title" not in body
    assert "# comment" not in body


@pytest.mark.skipif(not GENERATE_PY.is_file(), reason="generator not in this install")
def test_generate_errors_when_nist_download_fails():
    gen = _load_generate_module()
    original = gen.fetch_nist_dump

    def boom():
        raise SystemExit("Failed to download NIST dump from example: timed out")

    gen.fetch_nist_dump = boom
    try:
        with pytest.raises(SystemExit, match="Failed to download NIST dump"):
            gen.generate()
    finally:
        gen.fetch_nist_dump = original


@pytest.mark.skipif(not GENERATE_PY.is_file(), reason="generator not in this install")
def test_generate_writes_nothing_when_dump_and_module_match():
    gen = _load_generate_module()
    original = gen.fetch_nist_dump
    vendored = gen.NIST_TXT.read_text(encoding="utf-8")
    module_before = gen.OUTPUT_PY.read_text(encoding="utf-8")
    changes_mtime = gen.CHANGES_MD.stat().st_mtime

    gen.fetch_nist_dump = lambda: vendored
    try:
        gen.generate()
    finally:
        gen.fetch_nist_dump = original

    assert gen.NIST_TXT.read_text(encoding="utf-8") == vendored
    assert gen.OUTPUT_PY.read_text(encoding="utf-8") == module_before
    assert gen.CHANGES_MD.stat().st_mtime == changes_mtime
