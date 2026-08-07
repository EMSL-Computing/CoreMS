"""Tests for multi-charge molecular formula search (min/max ion charge).

No 13C charge determination: search expands absolute charge range and
assignment charge lives on MolecularFormula.ion_charge. Peak polarity
(±1) is separate; peak.ion_charge is a deprecated alias for polarity.
"""

import sys

sys.path.append(".")

import pytest

from corems.encapsulation.factory.processingSetting import MolecularFormulaSearchSettings
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


class TestIonChargesForSearch:
    def test_defaults_positive_int(self):
        assert SearchMolecularFormulas.ion_charges_for_search(1) == (1,)
        assert SearchMolecularFormulas.ion_charges_for_search(1, 1, 1) == (1,)

    def test_negative_int_range(self):
        assert SearchMolecularFormulas.ion_charges_for_search(-1, 1, 2) == (-1, -2)

    def test_positive_string_range(self):
        assert SearchMolecularFormulas.ion_charges_for_search("positive", 1, 3) == (
            1,
            2,
            3,
        )

    def test_negative_string_range(self):
        assert SearchMolecularFormulas.ion_charges_for_search("negative", 1, 2) == (
            -1,
            -2,
        )

    def test_sort_order_lower_abs_first(self):
        # first_hit prefers lower |z| when locking early
        charges = SearchMolecularFormulas.ion_charges_for_search(1, 1, 3)
        assert charges == (1, 2, 3)
        assert abs(charges[0]) <= abs(charges[-1])

    def test_invalid_min(self):
        with pytest.raises(ValueError, match="min_ion_charge"):
            SearchMolecularFormulas.ion_charges_for_search(1, 0, 1)

    def test_invalid_max_lt_min(self):
        with pytest.raises(ValueError, match="max_ion_charge"):
            SearchMolecularFormulas.ion_charges_for_search(1, 2, 1)

    def test_invalid_polarity_string(self):
        with pytest.raises(ValueError, match="positive"):
            SearchMolecularFormulas.ion_charges_for_search("both", 1, 1)

    def test_invalid_polarity_zero(self):
        with pytest.raises(ValueError, match="non-zero"):
            SearchMolecularFormulas.ion_charges_for_search(0, 1, 1)


class TestMolecularFormulaSearchSettingsChargeRange:
    def test_defaults_are_single_charge(self):
        s = MolecularFormulaSearchSettings()
        assert s.min_ion_charge == 1
        assert s.max_ion_charge == 1


class TestExportIonChargeSemantics:
    """Export 'Ion Charge': polarity when unassigned, formula charge when assigned."""

    def test_unassigned_export_uses_peak_polarity(self):
        mz = [100.0]
        abundance = [1.0]
        rp, s2n = [[1000.0], [10.0]]
        mass_spectrum_obj = ms_from_array_centroid(
            mz, abundance, rp, s2n, "unassigned export", polarity=-1, auto_process=False
        )
        mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
        mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0
        mass_spectrum_obj.process_mass_spec()
        peak = mass_spectrum_obj[0]
        assert peak.polarity == -1
        assert not peak.is_assigned
        df = mass_spectrum_obj.to_dataframe()
        assert list(df["Ion Charge"]) == [-1]

    @pytest.mark.molecular_db
    def test_assigned_export_uses_formula_ion_charge(self, postgres_database):
        formula = MolecularFormula({"C": 6, "H": 6}, ion_charge=2, ion_type="RADICAL")
        mz_z2 = formula.mz_calc
        mz = [mz_z2]
        abundance = [1.0]
        rp, s2n = [[1000.0], [100.0]]
        mass_spectrum_obj = ms_from_array_centroid(
            mz,
            abundance,
            rp,
            s2n,
            "assigned export z2",
            polarity=1,
            auto_process=False,
        )
        mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
        mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

        settings = mass_spectrum_obj.molecular_search_settings
        settings.url_database = postgres_database
        settings.error_method = "None"
        settings.min_ppm_error = -10
        settings.max_ppm_error = 10
        settings.mz_error_range = 1
        settings.isProtonated = False
        settings.isRadical = True
        settings.isAdduct = False
        settings.use_min_peaks_filter = False
        settings.min_ion_charge = 1
        settings.max_ion_charge = 2
        settings.usedAtoms = {
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

        assert peak.is_assigned
        assert peak.polarity == 1
        df = mass_spectrum_obj.to_dataframe()
        # Assigned rows export formula charge (z=2), not peak polarity
        assert 2 in set(df["Ion Charge"])
        assert all(z in (1, 2) for z in df["Ion Charge"])


@pytest.mark.molecular_db
def test_di_multi_charge_search_assigns_formula_charge(postgres_database):
    """DI search at max_ion_charge=2 can assign a z=2 radical formula.

    Peak polarity stays ±1; formula ion_charge is the search charge.
    """
    # C6H6 radical at z=2: m/z = (neutral - 2*e) / 2
    formula = MolecularFormula({"C": 6, "H": 6}, ion_charge=2, ion_type="RADICAL")
    mz_z2 = formula.mz_calc
    # Approximate 13C partner abundance not required for assignment of mono
    mz = [mz_z2]
    abundance = [1.0]
    rp, s2n = [[1000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "multi charge radical", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    settings = mass_spectrum_obj.molecular_search_settings
    settings.url_database = postgres_database
    settings.error_method = "None"
    settings.min_ppm_error = -10
    settings.max_ppm_error = 10
    settings.mz_error_range = 1
    settings.isProtonated = False
    settings.isRadical = True
    settings.isAdduct = False
    settings.use_min_peaks_filter = False
    settings.min_ion_charge = 1
    settings.max_ion_charge = 2
    settings.usedAtoms = {
        "C": (6, 6),
        "H": (6, 6),
        "O": (0, 0),
        "N": (0, 0),
    }

    mass_spectrum_obj.process_mass_spec()
    peak = mass_spectrum_obj[0]
    peak_polarity_before = peak.polarity

    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=False
    ).run_worker_ms_peaks([peak])

    # Peak carries polarity only; ion_charge is deprecated alias
    assert peak.polarity == peak_polarity_before
    assert peak.polarity == 1
    assert peak.ion_charge == peak.polarity

    assert peak.is_assigned
    formula_charges = {mf.ion_charge for mf in peak}
    assert 2 in formula_charges
    z2_hits = [mf for mf in peak if mf.ion_charge == 2]
    assert any("C6" in mf.string and "H6" in mf.string for mf in z2_hits)


@pytest.mark.molecular_db
def test_di_default_max_one_skips_z2_only_peak(postgres_database):
    """With default max_ion_charge=1, a pure z=2 m/z need not get a z=2 formula."""
    formula = MolecularFormula({"C": 6, "H": 6}, ion_charge=2, ion_type="RADICAL")
    mz_z2 = formula.mz_calc
    mz = [mz_z2]
    abundance = [1.0]
    rp, s2n = [[1000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "single charge only", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    settings = mass_spectrum_obj.molecular_search_settings
    settings.url_database = postgres_database
    settings.error_method = "None"
    settings.min_ppm_error = -5
    settings.max_ppm_error = 5
    settings.isProtonated = False
    settings.isRadical = True
    settings.isAdduct = False
    settings.use_min_peaks_filter = False
    # defaults min=max=1
    assert settings.min_ion_charge == 1
    assert settings.max_ion_charge == 1
    settings.usedAtoms = {
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
    # z=2 formula must not appear when max is 1
    assert not any(getattr(mf, "ion_charge", None) == 2 for mf in peak)
