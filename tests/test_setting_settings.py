import sys
import os

import pytest

from corems.encapsulation.output import parameter_to_json, parameter_to_dict
from corems.encapsulation.factory.processingSetting import (
    MolecularLookupDictSettings,
    MolecularFormulaSearchSettings,
    validate_used_atoms_keys,
)
from corems.encapsulation.factory.parameters import hush_output, reset_ms_parameters, reset_gcms_parameters, reset_lcms_parameters, LCMSParameters, MSParameters, GCMSParameters

def test_toml():
      
    parameter_to_json.dump_all_settings_toml()
    assert os.path.exists('SettingsCoreMS.toml')
    os.remove('SettingsCoreMS.toml')

    parameter_to_json.dump_gcms_settings_toml()
    assert os.path.exists('SettingsCoreMS.toml')
    os.remove('SettingsCoreMS.toml')

    parameter_to_json.dump_ms_settings_toml()
    assert os.path.exists('SettingsCoreMS.toml')
    os.remove('SettingsCoreMS.toml')


def test_toml_omits_null_settings(tmp_path):
    """Unset numeric settings are left out of the TOML file.

    TOML has no null. Those fields stay at their in-memory defaults and are
    not written as values.
    """
    import tomlkit

    path = tmp_path / "SettingsCoreMS.toml"
    parameter_to_json.dump_all_settings_toml(file_path=path)
    with path.open(encoding="utf-8") as handle:
        data = tomlkit.load(handle)

    section = data["SpectralSimilaritySearch"]
    assert "precursor_ions_removal_da" not in section
    assert "peak_sep_da" not in section
    assert section["max_ms2_tolerance_in_da"] == 0.01
    assert section["search_type"] == "open"


def test_toml_reloads_used_atoms_as_tuples(tmp_path):
    """Atom-count ranges stay tuples after a TOML round trip.

    Each range is stored as an array. Loading puts a tuple back on the
    parameter object so the reloaded settings compare equal to the original.
    """
    from corems.encapsulation.input.parameter_from_json import (
        load_and_set_toml_parameters_lcms,
    )

    class Holder:
        def __init__(self):
            self.parameters = LCMSParameters(use_defaults=True)

    original = Holder()
    original.parameters.mass_spectrum["ms1"].molecular_search.usedAtoms = {
        "C": (5, 30),
        "H": (18, 200),
    }
    path = tmp_path / "SettingsCoreMS.toml"
    parameter_to_json.dump_lcms_settings_toml(file_path=path, lcms_obj=original)

    loaded = Holder()
    load_and_set_toml_parameters_lcms(loaded, parameters_path=str(path))
    used = loaded.parameters.mass_spectrum["ms1"].molecular_search.usedAtoms
    assert used["C"] == (5, 30)
    assert isinstance(used["C"], tuple)
    assert used["H"] == (18, 200)
    assert isinstance(used["H"], tuple)
    assert loaded.parameters == original.parameters


def test_json():
      
    parameter_to_json.dump_all_settings_json()
    assert os.path.exists('SettingsCoreMS.json')
    os.remove('SettingsCoreMS.json')

    parameter_to_json.dump_gcms_settings_json()
    assert os.path.exists('SettingsCoreMS.json')
    os.remove('SettingsCoreMS.json')

    parameter_to_json.dump_ms_settings_json()
    assert os.path.exists('SettingsCoreMS.json')
    os.remove('SettingsCoreMS.json')
   
def test_data():
    
    param_dict = parameter_to_dict.get_dict_ms_default_data()
    assert len(param_dict) > 4
    param_dict = parameter_to_dict.get_dict_gcms_default_data()
    assert  len(param_dict) > 1
    

def test_settings_search():

    test = MolecularLookupDictSettings().__dict__
    assert len(test) > 10
    assert "usedAtoms" in test
    assert "url_database" in test


def test_validate_used_atoms_keys_accepts_elements():
    """Element symbols (mono codes) are valid usedAtoms keys."""
    validate_used_atoms_keys({"C": (1, 10), "H": (1, 20), "Fe": (0, 1), "Cl": (0, 2)})
    settings = MolecularFormulaSearchSettings(
        usedAtoms={"C": (1, 20), "H": (4, 40), "O": (0, 5), "Fe": (0, 1)}
    )
    assert "Fe" in settings.usedAtoms
    # full reassignment with element keys
    settings.usedAtoms = {"C": (1, 10), "H": (1, 20), "N": (0, 2)}
    assert settings.usedAtoms["N"] == (0, 2)


def test_validate_used_atoms_keys_rejects_rare_isotopes():
    """Specific isotope labels are not valid usedAtoms keys."""
    with pytest.raises(ValueError, match="54Fe"):
        validate_used_atoms_keys({"C": (1, 10), "H": (1, 20), "54Fe": (0, 1)})

    with pytest.raises(ValueError, match="13C"):
        MolecularFormulaSearchSettings(
            usedAtoms={"C": (1, 10), "H": (1, 20), "13C": (0, 1)}
        )

    settings = MolecularFormulaSearchSettings()
    with pytest.raises(ValueError, match="37Cl"):
        settings.usedAtoms = {"C": (1, 10), "H": (1, 20), "37Cl": (0, 1)}


def test_search_entry_rejects_in_place_rare_used_atoms(mass_spectrum_ftms):
    """In-place usedAtoms mutation is caught when search is constructed."""
    from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas

    # Bypass assignment validation by mutating the dict in place
    mass_spectrum_ftms.molecular_search_settings.usedAtoms["54Fe"] = (0, 1)
    with pytest.raises(ValueError, match="54Fe"):
        SearchMolecularFormulas(mass_spectrum_ftms, find_isotopologues=True)


def test_hush_setting():
    LCMSParameters.lc_ms.eic_tolerance_ppm = 10 # set to 10

    assert LCMSParameters.lc_ms.eic_tolerance_ppm == 10

    reset_ms_parameters()
    reset_gcms_parameters()
    reset_lcms_parameters()

    assert LCMSParameters.lc_ms.eic_tolerance_ppm == 5 # default value
    assert LCMSParameters.lc_ms.verbose_processing
    assert MSParameters.mass_spectrum.verbose_processing
    assert MSParameters.molecular_search.verbose_processing
    assert GCMSParameters.gc_ms.verbose_processing

    hush_output()

    assert not LCMSParameters.lc_ms.verbose_processing
    assert not MSParameters.mass_spectrum.verbose_processing
    assert not MSParameters.molecular_search.verbose_processing
    assert not GCMSParameters.gc_ms.verbose_processing

    reset_ms_parameters()
    reset_gcms_parameters()
    reset_lcms_parameters()
