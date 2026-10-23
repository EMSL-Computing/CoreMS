import sys

import pytest

from corems.encapsulation.constant import Labels
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.molecular_formula.input.masslist_ref import ImportMassListRef
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


@pytest.mark.molecular_db
def test_search_imported_ref_files(mass_spectrum_ftms, ref_file_location, postgres_database):
    mass_spectrum_obj = mass_spectrum_ftms
    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mf_references_list = ImportMassListRef(ref_file_location).from_bruker_ref_file()
    assert len(mf_references_list) == 60
    assert round(mf_references_list[0].mz_calc, 2) == 149.06
    assert mf_references_list[0].class_label == "O2"

    ion_type = "unknown"

    ms_peaks_assigned = SearchMolecularFormulas(mass_spectrum_obj).search_mol_formulas(
        mf_references_list, ion_type, neutral_molform=False, find_isotopologues=False
    )

    assert (len(ms_peaks_assigned)) > 10


# ---------------------------------------------------------------------------
# Regression tests for ion-type-aware nominal-mass binning in
# SearchMolecularFormulas.search_mol_formulas().
#
# Neutral candidates are binned by the nominal m/z of the *requested* ion
# hypothesis. Before the fix, every neutral candidate was binned by its
# protonated m/z regardless of ion_type, so radical and adduct candidates
# landed in the wrong bin and were never retrieved for the matching peak.
#
# Every experimental peak below is placed at the m/z CoreMS itself calculates
# for the candidate/ion hypothesis (candidate.protonated_mz / .radical_mz /
# .adduct_mz(atom)); no mass constants are duplicated here. Formula mass math
# itself is covered by tests/test_molecular_formula.py.
# ---------------------------------------------------------------------------

# A neutral composition and a clearly different decoy composition. The decoy
# sits in a different nominal bin, so a correct search never assigns it.
_TRUE_COMPOSITION = {"C": 56, "H": 73, "N": 1}
_DECOY_COMPOSITION = {"C": 40, "H": 30, "O": 15}


def _neutral_candidate(composition, ion_charge):
    return MolecularFormula(dict(composition), ion_charge=ion_charge)


def _run_targeted_search(
    target_mz, ion_type, ion_charge, polarity, adduct_atom, postgres_database
):
    """Build a two-peak spectrum (one true peak, one far-off decoy peak) and
    search a candidate list of [true, decoy] neutral formulas against it."""
    decoy_peak_mz = target_mz + 50.0
    mz = [target_mz, decoy_peak_mz]
    abundance = [1.0, 1.0]
    rp, s2n = [[10000.0, 10000.0], [100.0, 100.0]]

    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "targeted mf search", polarity=polarity, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1

    mass_spectrum_obj.process_mass_spec()

    candidates = [
        _neutral_candidate(_TRUE_COMPOSITION, ion_charge),
        _neutral_candidate(_DECOY_COMPOSITION, ion_charge),
    ]

    assigned = SearchMolecularFormulas(mass_spectrum_obj).search_mol_formulas(
        candidates,
        ion_type,
        neutral_molform=True,
        find_isotopologues=False,
        adduct_atom=adduct_atom,
    )

    return mass_spectrum_obj, assigned


@pytest.mark.molecular_db
def test_search_mol_formulas_protonated_positive(postgres_database):
    """[M+H]+: binning and matching both use the protonated m/z (unchanged behavior)."""
    candidate = _neutral_candidate(_TRUE_COMPOSITION, ion_charge=1)
    target_mz = candidate.protonated_mz

    mass_spectrum_obj, assigned = _run_targeted_search(
        target_mz, Labels.protonated_de_ion, 1, 1, None, postgres_database
    )

    assert len(assigned) == 1
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].ion_charge == 1
    assert not mass_spectrum_obj[1].is_assigned


@pytest.mark.molecular_db
def test_search_mol_formulas_deprotonated_negative(postgres_database):
    """[M-H]-: negative-charge protonated_mz resolves to the deprotonated mass."""
    candidate = _neutral_candidate(_TRUE_COMPOSITION, ion_charge=-1)
    target_mz = candidate.protonated_mz

    mass_spectrum_obj, assigned = _run_targeted_search(
        target_mz, Labels.protonated_de_ion, -1, -1, None, postgres_database
    )

    assert len(assigned) == 1
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].ion_charge == -1
    assert not mass_spectrum_obj[1].is_assigned


@pytest.mark.molecular_db
def test_search_mol_formulas_radical_cation(postgres_database):
    """[M].+ radical: the primary regression case.

    The radical m/z is ~1 Da (a proton) below the protonated m/z, so before the
    fix the candidate was binned under int(protonated_mz) while the peak's
    nominal bin is int(radical_mz) -- the candidate was never retrieved and the
    peak went unassigned. After the fix, binning uses radical_mz and it matches.
    """
    candidate = _neutral_candidate(_TRUE_COMPOSITION, ion_charge=1)
    target_mz = candidate.radical_mz

    mass_spectrum_obj, assigned = _run_targeted_search(
        target_mz, Labels.radical_ion, 1, 1, None, postgres_database
    )

    assert len(assigned) == 1
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].ion_type == Labels.radical_ion
    assert mass_spectrum_obj[0][0].ion_charge == 1
    assert not mass_spectrum_obj[1].is_assigned


@pytest.mark.molecular_db
def test_search_mol_formulas_adduct_sodium(postgres_database):
    """[M+Na]+ adduct: the adduct m/z is ~22 Da above the protonated m/z, an
    unambiguous nominal-bin displacement that the old protonated-only binning
    got wrong. After the fix, binning uses adduct_mz('Na') and it matches."""
    candidate = _neutral_candidate(_TRUE_COMPOSITION, ion_charge=1)
    target_mz = candidate.adduct_mz("Na")

    mass_spectrum_obj, assigned = _run_targeted_search(
        target_mz, Labels.adduct_ion, 1, 1, "Na", postgres_database
    )

    assert len(assigned) == 1
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].ion_type == Labels.adduct_ion
    assert mass_spectrum_obj[0][0].adduct_atom == "Na"
    assert mass_spectrum_obj[0][0].ion_charge == 1
    assert not mass_spectrum_obj[1].is_assigned
