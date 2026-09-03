"""Unit tests for shared DI/LC formula search job orchestration."""

from types import SimpleNamespace

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.processingSetting import (
    MolecularFormulaSearchSettings,
)
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.molecular_id.search.molecularFormulaSearch import (
    FormulaSearchJob,
    SearchMolecularFormulas,
    SearchMolecularFormulasLC,
)


def _settings(**kwargs):
    s = MolecularFormulaSearchSettings()
    s.isProtonated = kwargs.get("isProtonated", True)
    s.isRadical = kwargs.get("isRadical", False)
    s.isAdduct = kwargs.get("isAdduct", False)
    s.verbose_processing = False
    s.db_chunk_size = kwargs.get("db_chunk_size", 300)
    return s


def _dict_res():
    return {
        Labels.protonated_de_ion: {"C1H1": {100: ["prot"]}},
        Labels.radical_ion: {"C1H1": {101: ["rad"]}},
        Labels.adduct_ion: {
            "Na": {"C1H1": {102: ["na"]}},
            "K": {"C1H1": {}},
        },
    }


def test_iter_ion_type_jobs_protonated_only():
    jobs = list(
        SearchMolecularFormulas.iter_ion_type_jobs(
            _settings(isProtonated=True, isRadical=False, isAdduct=False),
            _dict_res(),
            "C1H1",
            1,
        )
    )
    assert len(jobs) == 1
    job = jobs[0]
    assert isinstance(job, FormulaSearchJob)
    assert job.ion_charge == 1
    assert job.ion_type == Labels.protonated_de_ion
    assert job.adduct_atom is None
    assert job.candidate_formulas == {100: ["prot"]}
    assert "protonated" in job.progress
    assert "C1H1" in job.progress
    assert "z=1" in job.progress


def test_iter_ion_type_jobs_all_types_at_z1():
    jobs = list(
        SearchMolecularFormulas.iter_ion_type_jobs(
            _settings(isProtonated=True, isRadical=True, isAdduct=True),
            _dict_res(),
            "C1H1",
            1,
        )
    )
    assert [j.ion_type for j in jobs] == [
        Labels.protonated_de_ion,
        Labels.radical_ion,
        Labels.adduct_ion,
    ]
    assert [j.adduct_atom for j in jobs] == [None, None, "Na"]
    # Empty K candidate list is skipped
    assert all(j.adduct_atom != "K" for j in jobs)


def test_iter_ion_type_jobs_skips_adduct_when_abs_z_not_1():
    jobs = list(
        SearchMolecularFormulas.iter_ion_type_jobs(
            _settings(isProtonated=True, isRadical=False, isAdduct=True),
            _dict_res(),
            "C1H1",
            2,
        )
    )
    assert [j.ion_type for j in jobs] == [Labels.protonated_de_ion]
    jobs_neg = list(
        SearchMolecularFormulas.iter_ion_type_jobs(
            _settings(isProtonated=False, isRadical=False, isAdduct=True),
            _dict_res(),
            "C1H1",
            -2,
        )
    )
    assert jobs_neg == []


def test_iter_ion_type_jobs_missing_class_or_ion_key_is_empty():
    jobs = list(
        SearchMolecularFormulas.iter_ion_type_jobs(
            _settings(isProtonated=True, isRadical=True, isAdduct=True),
            {Labels.protonated_de_ion: {}},
            "C1H1",
            1,
        )
    )
    assert jobs == []


# Independent arithmetic: [M+Na]+ of C56 H73 N1 (same constant as multi-charge tests)
MZ_C56H73N1_MNA_Z1 = 782.5635220593


def test_run_formula_search_jobs_one_db_load_per_charge_per_chunk(monkeypatch):
    db_calls = []

    def fake_database_to_dict(
        classe_str_list, nominal_mzs, mf_search_settings, ion_charge, sql_db=None
    ):
        db_calls.append((tuple(classe_str_list), ion_charge))
        return {
            Labels.protonated_de_ion: {
                "C1": {10: ["c1"]},
                "C2": {20: ["c2"]},
            }
        }

    monkeypatch.setattr(
        SearchMolecularFormulas, "database_to_dict", fake_database_to_dict
    )

    settings = _settings(isProtonated=True, isRadical=False, isAdduct=False)
    settings.db_chunk_size = 1
    jobs = []
    SearchMolecularFormulas.run_formula_search_jobs(
        classes=[("C1", {}), ("C2", {})],
        nominal_mzs=[10, 20],
        mf_search_settings=settings,
        search_charges=(1, 2),
        sql_db=None,
        apply_fn=jobs.append,
    )

    assert db_calls == [
        (("C1",), 1),
        (("C1",), 2),
        (("C2",), 1),
        (("C2",), 2),
    ]
    assert [j.ion_charge for j in jobs] == [1, 2, 1, 2]
    assert [j.candidate_formulas for j in jobs] == [
        {10: ["c1"]},
        {10: ["c1"]},
        {20: ["c2"]},
        {20: ["c2"]},
    ]
    assert all(j.ion_type == Labels.protonated_de_ion for j in jobs)


class _DummySQL:
    def close(self):
        pass


def _fake_lcms(mol_search):
    return SimpleNamespace(
        polarity="positive",
        parameters=SimpleNamespace(
            mass_spectrum={"ms1": SimpleNamespace(molecular_search=mol_search)}
        ),
    )


def test_lc_search_spectra_stores_adduct_atom():
    mz = [MZ_C56H73N1_MNA_Z1]
    abundance = [1.0]
    rp, s2n = [[10000.0], [100.0]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "lc adduct atom", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False
    mass_spectrum_obj.process_mass_spec()

    mf = MolecularFormula({"C": 56, "H": 73, "N": 1}, ion_charge=1)
    query = {int(MZ_C56H73N1_MNA_Z1): [mf]}

    lc = SearchMolecularFormulasLC(
        _fake_lcms(mass_spectrum_obj.molecular_search_settings),
        sql_db=_DummySQL(),
        find_isotopologues=False,
    )
    lc.search_spectra_against_candidates(
        mass_spectrum_list=[mass_spectrum_obj],
        ms_peaks_list=[[mass_spectrum_obj[0]]],
        candidate_formulas=query,
        ion_type=Labels.adduct_ion,
        ion_charge=1,
        adduct_atom="Na",
    )

    assert mass_spectrum_obj[0].is_assigned
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].adduct_atom == "Na"
    assert mass_spectrum_obj[0][0].ion_charge == 1
