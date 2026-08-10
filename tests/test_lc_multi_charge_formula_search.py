"""Tests for multi-charge molecular formula search (min/max ion charge).

m/z inputs are hard-coded (not taken from MolecularFormula.mz_calc) so a bug in
formula mass math cannot make search tests pass circularly. Values use the same
element masses as Atoms (C, H, N, Na, 13C, electron) with:

  neutral mass M = sum n_i * mass_i
  radical  [M].z+   : (M - z*e) / |z|
  protonated [M+zH]z+ : (M + z*H - z*e) / |z|
  adduct   [M+Na]+  : (M + Na - e) / 1   (adduct search is single-charge only)
"""

import sys

sys.path.append(".")

import pytest

from corems.encapsulation.factory.processingSetting import MolecularFormulaSearchSettings
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas

# C6 H6 radical, z=2  -> (6*C + 6*H - 2*e) / 2
MZ_C6H6_RADICAL_Z2 = 39.0229265168

# C56 H73 N1 and mono-13C partners (independent atomic-mass arithmetic)
# [M+H]+ / [M+H]+ 13C
MZ_C56H73N1_MH_Z1 = 760.5815778095
MZ_C55H73N1_13C_MH_Z1 = 761.5849326446
# [M+Na]+ / [M+Na]+ 13C (single-charge adduct only)
MZ_C56H73N1_MNA_Z1 = 782.5635220593
MZ_C55H73N1_13C_MNA_Z1 = 783.5668768944
# [M+2H]2+ / [M+2H]2+ 13C
MZ_C56H73N1_MH_Z2 = 380.7944271309
MZ_C55H73N1_13C_MH_Z2 = 381.2961045485


def test_ion_charges_for_search():
    """Signed charge tuples from polarity and absolute min/max bounds."""
    assert SearchMolecularFormulas.ion_charges_for_search(1) == (1,)
    assert SearchMolecularFormulas.ion_charges_for_search(1, 1, 1) == (1,)
    assert SearchMolecularFormulas.ion_charges_for_search(-1, 1, 2) == (-1, -2)
    assert SearchMolecularFormulas.ion_charges_for_search("positive", 1, 3) == (1, 2, 3)
    assert SearchMolecularFormulas.ion_charges_for_search("negative", 1, 2) == (-1, -2)
    # lower |z| first (relevant when first_hit=True)
    assert SearchMolecularFormulas.ion_charges_for_search(1, 1, 3) == (1, 2, 3)


def test_ion_charges_for_search_invalid():
    """Invalid polarity or charge bounds raise ValueError."""
    with pytest.raises(ValueError, match="min_ion_charge"):
        SearchMolecularFormulas.ion_charges_for_search(1, 0, 1)
    with pytest.raises(ValueError, match="max_ion_charge"):
        SearchMolecularFormulas.ion_charges_for_search(1, 2, 1)
    with pytest.raises(ValueError, match="positive"):
        SearchMolecularFormulas.ion_charges_for_search("both", 1, 1)
    with pytest.raises(ValueError, match="non-zero"):
        SearchMolecularFormulas.ion_charges_for_search(0, 1, 1)


def test_min_max_ion_charge_defaults():
    """Defaults preserve single-charge search (min=max=1)."""
    s = MolecularFormulaSearchSettings()
    assert s.min_ion_charge == 1
    assert s.max_ion_charge == 1


def test_legacy_mfss_ion_charge_accepted_but_unused():
    """Legacy MFSS.ion_charge must not break construct/load or control search.

    Older parameter files and constructors still pass ion_charge (often -1).
    Search polarity comes from spectrum/LC data; multi-z from min/max only.
    """
    # Constructor BC (YAML/TOML-style kwargs)
    s = MolecularFormulaSearchSettings(ion_charge=-1)
    assert s.ion_charge == -1
    assert s.min_ion_charge == 1
    assert s.max_ion_charge == 1
    s.validate_ion_charge_settings()

    # Positive legacy value must not be treated as multi-z or polarity
    s_pos = MolecularFormulaSearchSettings(
        ion_charge=1, min_ion_charge=1, max_ion_charge=2, isAdduct=False
    )
    assert s_pos.ion_charge == 1
    assert s_pos.max_ion_charge == 2
    s_pos.validate_ion_charge_settings()

    # setattr as parameter loaders do
    s_loaded = MolecularFormulaSearchSettings()
    setattr(s_loaded, "ion_charge", -1)
    s_loaded.validate_ion_charge_settings()
    assert s_loaded.ion_charge == -1

    # Charge list ignores settings.ion_charge; uses polarity + min/max only
    assert SearchMolecularFormulas.ion_charges_for_search(
        1, s_pos.min_ion_charge, s_pos.max_ion_charge
    ) == (1, 2)
    assert SearchMolecularFormulas.ion_charges_for_search(
        -1, s_pos.min_ion_charge, s_pos.max_ion_charge
    ) == (-1, -2)


def test_isAdduct_incompatible_with_multi_charge():
    """isAdduct=True with max_ion_charge > 1 is rejected at construction."""
    with pytest.raises(ValueError, match="isAdduct=True is incompatible"):
        MolecularFormulaSearchSettings(isAdduct=True, max_ion_charge=2)

    s = MolecularFormulaSearchSettings(isAdduct=True, max_ion_charge=1)
    s.max_ion_charge = 2
    with pytest.raises(ValueError, match="isAdduct=True is incompatible"):
        s.validate_ion_charge_settings()


def test_export_unassigned_ion_charge_is_polarity():
    """Unassigned export 'Ion Charge' is peak polarity (±1)."""
    mz = [100.0]
    abundance = [1.0]
    rp, s2n = [[1000.0], [10.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "unassigned export", polarity=-1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0
    mass_spectrum_obj.process_mass_spec()

    assert mass_spectrum_obj[0].polarity == -1
    assert not mass_spectrum_obj[0].is_assigned
    assert list(mass_spectrum_obj.to_dataframe()["Ion Charge"]) == [-1]


@pytest.mark.molecular_db
def test_export_assigned_ion_charge_is_formula_charge(postgres_database):
    """Assigned export 'Ion Charge' is formula ion_charge (may be multi-charge)."""
    mz = [MZ_C6H6_RADICAL_Z2]
    abundance = [1.0]
    rp, s2n = [[1000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "assigned export z2", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -10
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 10
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = False
    mass_spectrum_obj.molecular_search_settings.isRadical = True
    mass_spectrum_obj.molecular_search_settings.isAdduct = False
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.min_ion_charge = 1
    mass_spectrum_obj.molecular_search_settings.max_ion_charge = 2
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (6, 6),
        "H": (6, 6),
        "O": (0, 0),
        "N": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=False
    ).run_worker_ms_peaks([mass_spectrum_obj[0]])

    assert mass_spectrum_obj[0].is_assigned
    assert mass_spectrum_obj[0].polarity == 1
    assert 2 in set(mass_spectrum_obj.to_dataframe()["Ion Charge"])


@pytest.mark.molecular_db
def test_multi_charge_radical_assigns_formula_charge(postgres_database):
    """max_ion_charge=2 can assign a z=2 radical; peak polarity stays ±1."""
    mz = [MZ_C6H6_RADICAL_Z2]
    abundance = [1.0]
    rp, s2n = [[1000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "multi charge radical", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -10
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 10
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = False
    mass_spectrum_obj.molecular_search_settings.isRadical = True
    mass_spectrum_obj.molecular_search_settings.isAdduct = False
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.min_ion_charge = 1
    mass_spectrum_obj.molecular_search_settings.max_ion_charge = 2
    # Legacy parameter-file field must not break search or flip polarity
    mass_spectrum_obj.molecular_search_settings.ion_charge = -1
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (6, 6),
        "H": (6, 6),
        "O": (0, 0),
        "N": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    peak = mass_spectrum_obj[0]
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=False
    ).run_worker_ms_peaks([peak])

    assert peak.polarity == 1
    assert peak.ion_charge == peak.polarity  # deprecated alias
    assert peak.is_assigned
    assert peak[0].string == "C6 H6"
    assert peak[0].ion_charge == 2
    # settings.ion_charge=-1 did not force negative search (would miss z=+2 radical)
    assert mass_spectrum_obj.molecular_search_settings.ion_charge == -1


@pytest.mark.molecular_db
def test_default_max_ion_charge_skips_z2_only_peak(postgres_database):
    """With default max_ion_charge=1, a pure z=2 m/z does not get a z=2 formula."""
    mz = [MZ_C6H6_RADICAL_Z2]
    abundance = [1.0]
    rp, s2n = [[1000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "single charge only", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.isProtonated = False
    mass_spectrum_obj.molecular_search_settings.isRadical = True
    mass_spectrum_obj.molecular_search_settings.isAdduct = False
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    # defaults min=max=1
    assert mass_spectrum_obj.molecular_search_settings.min_ion_charge == 1
    assert mass_spectrum_obj.molecular_search_settings.max_ion_charge == 1
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (6, 6),
        "H": (6, 6),
        "O": (0, 0),
        "N": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    peak = mass_spectrum_obj[0]
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=False
    ).run_worker_ms_peaks([peak])

    if peak.is_assigned:
        assert all(mf.ion_charge == 1 for mf in peak)
    assert not any(mf.ion_charge == 2 for mf in peak)


@pytest.mark.molecular_db
def test_multi_charge_protonated_with_c13(postgres_database):
    """Recover [M+H]+ and [M+2H]2+ with 13C partners (adducts off).

    Peak layout (hard-coded m/z; mono then 13C for each species):
      0,1  [M+H]+ and 13C
      2,3  [M+2H]2+ and 13C
    """
    mz = [
        MZ_C56H73N1_MH_Z1,
        MZ_C55H73N1_13C_MH_Z1,
        MZ_C56H73N1_MH_Z2,
        MZ_C55H73N1_13C_MH_Z2,
    ]
    abundance = [1.0, 0.4, 1.0, 0.4]
    rp, s2n = [[10000.0] * 4, [100.0] * 4]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "mh multi-charge c13", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = True
    mass_spectrum_obj.molecular_search_settings.isRadical = False
    mass_spectrum_obj.molecular_search_settings.isAdduct = False
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False
    mass_spectrum_obj.molecular_search_settings.min_ion_charge = 1
    mass_spectrum_obj.molecular_search_settings.max_ion_charge = 2
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (56, 56),
        "H": (73, 73),
        "N": (1, 1),
        "O": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    mono_peaks = [mass_spectrum_obj[0], mass_spectrum_obj[2]]
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=True
    ).run_worker_ms_peaks(mono_peaks)

    # [M+H]+
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].ion_charge == 1
    assert mass_spectrum_obj[0].polarity == 1
    assert mass_spectrum_obj[1][0].string == "C55 H73 N1 13C1"
    assert mass_spectrum_obj[1][0].ion_charge == 1
    assert mass_spectrum_obj[1][0].is_isotopologue

    # [M+2H]2+
    assert mass_spectrum_obj[2][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[2][0].ion_charge == 2
    assert mass_spectrum_obj[2].polarity == 1
    assert mass_spectrum_obj[3][0].string == "C55 H73 N1 13C1"
    assert mass_spectrum_obj[3][0].ion_charge == 2
    assert mass_spectrum_obj[3][0].is_isotopologue


@pytest.mark.molecular_db
def test_single_charge_mna_with_c13(postgres_database):
    """Recover [M+Na]+ with 13C partner when isAdduct on and max_ion_charge=1."""
    mz = [MZ_C56H73N1_MNA_Z1, MZ_C55H73N1_13C_MNA_Z1]
    abundance = [1.0, 0.4]
    rp, s2n = [[10000.0, 10000.0], [100.0, 100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "mna single-charge c13", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = True
    mass_spectrum_obj.molecular_search_settings.isRadical = False
    mass_spectrum_obj.molecular_search_settings.isAdduct = True
    mass_spectrum_obj.molecular_search_settings.adduct_atoms_pos = ("Na",)
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False
    mass_spectrum_obj.molecular_search_settings.min_ion_charge = 1
    mass_spectrum_obj.molecular_search_settings.max_ion_charge = 1
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (56, 56),
        "H": (73, 73),
        "N": (1, 1),
        "O": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=True
    ).run_worker_ms_peaks([mass_spectrum_obj[0]])

    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].ion_charge == 1
    assert mass_spectrum_obj[0][0].adduct_atom == "Na"
    assert mass_spectrum_obj[1][0].string == "C55 H73 N1 13C1"
    assert mass_spectrum_obj[1][0].ion_charge == 1
    assert mass_spectrum_obj[1][0].adduct_atom == "Na"
    assert mass_spectrum_obj[1][0].is_isotopologue


@pytest.mark.molecular_db
def test_search_errors_if_adduct_and_multi_charge(postgres_database):
    """Search fails loudly if settings mix isAdduct with max_ion_charge > 1."""
    mz = [MZ_C56H73N1_MH_Z1]
    abundance = [1.0]
    rp, s2n = [[10000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "bad settings", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.isProtonated = True
    mass_spectrum_obj.molecular_search_settings.isAdduct = True
    mass_spectrum_obj.molecular_search_settings.max_ion_charge = 2
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (56, 56),
        "H": (73, 73),
        "N": (1, 1),
        "O": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    with pytest.raises(ValueError, match="isAdduct=True is incompatible"):
        SearchMolecularFormulas(
            mass_spectrum_obj, find_isotopologues=False
        ).run_worker_ms_peaks([mass_spectrum_obj[0]])
