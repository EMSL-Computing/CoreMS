import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from corems.molecular_id.factory.classification import  HeteroatomsClassification, Labels
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


@pytest.fixture
def classified_mass_spectrum(mass_spectrum_ftms, postgres_database):
    """Run molecular formula search and return a HeteroatomsClassification object."""
    mass_spectrum_ftms.molecular_search_settings.url_database = postgres_database
    mass_spectrum_ftms.molecular_search_settings.error_method = 'None'
    mass_spectrum_ftms.molecular_search_settings.min_ppm_error = -10
    mass_spectrum_ftms.molecular_search_settings.max_ppm_error = 10
    mass_spectrum_ftms.molecular_search_settings.mz_error_range = 1
    mass_spectrum_ftms.molecular_search_settings.isProtonated = True
    mass_spectrum_ftms.molecular_search_settings.isRadical = False
    mass_spectrum_ftms.molecular_search_settings.isAdduct = False
    usedAtoms = {'C': (1, 100), 'H': (4, 200), 'O': (1, 18)}
    mass_spectrum_ftms.molecular_search_settings.usedAtoms = usedAtoms

    assert mass_spectrum_ftms.percentile_assigned()[2] == 0
    SearchMolecularFormulas(mass_spectrum_ftms).run_worker_mass_spectrum()
    assert mass_spectrum_ftms.percentile_assigned()[2] > 0

    return HeteroatomsClassification(mass_spectrum_ftms)


def test_heteroatoms_classification(classified_mass_spectrum):
    mass_spectrum_by_classes = classified_mass_spectrum

    # Check that the plot is created
    mass_spectrum_by_classes.plot_ms_assigned_unassigned()

    # Check that ratios, DBE, carbon number, abundance and mz_exp are calculated

    assert len(mass_spectrum_by_classes.atoms_ratio_all("H", "C")) > 0
    assert len(mass_spectrum_by_classes.dbe_all()) > 0
    assert len(mass_spectrum_by_classes.abundance_assigned()) > 0
    assert len(mass_spectrum_by_classes.mz_exp_assigned()) > 0
    assert mass_spectrum_by_classes.abundance_count_percentile(Labels.unassigned) > 0
    assert mass_spectrum_by_classes.peaks_count_percentile(Labels.unassigned) > 0


def test_plot_van_krevelen_single_class(classified_mass_spectrum):
    """Test Van Krevelen plot for a single heteroatom class (backward compat)."""
    mass_spectrum = classified_mass_spectrum
    # Get a valid class name from the assigned classes
    assigned_classes = [c for c in mass_spectrum.get_classes() if c != Labels.unassigned]
    assert len(assigned_classes) > 0, "No assigned classes found"

    plt.figure()
    result = mass_spectrum.plot_van_krevelen(assigned_classes[0])
    ax, abun_perc = result
    assert ax is not None
    assert abun_perc > 0
    assert assigned_classes[0] in ax.get_title()
    assert ax.get_xlabel() == "O/C"
    assert ax.get_ylabel() == "H/C"
    plt.close()


def test_plot_van_krevelen_all_classes(classified_mass_spectrum):
    """Test Van Krevelen plot for all assigned classes (new functionality)."""
    plt.figure()
    ax = classified_mass_spectrum.plot_van_krevelen()
    assert ax is not None
    assert "All Assigned Classes" in ax.get_title()
    assert ax.get_xlabel() == "O/C"
    assert ax.get_ylabel() == "H/C"
    plt.close()


def test_plot_van_krevelen_log_abundance(classified_mass_spectrum):
    """Test Van Krevelen plot with log10 abundance scaling."""
    plt.figure()
    ax = classified_mass_spectrum.plot_van_krevelen(log_abundance=True)
    assert ax is not None
    # Check that colorbar has log label
    cbar = ax.figure.axes[-1]  # colorbar is the last axes
    assert "log" in cbar.get_ylabel().lower() or "log" in cbar.get_ylabel()
    plt.close()


def test_plot_dbe_vs_carbon_number_single_class(classified_mass_spectrum):
    """Test DBE vs Carbon Number plot for a single class (backward compat)."""
    mass_spectrum = classified_mass_spectrum
    assigned_classes = [c for c in mass_spectrum.get_classes() if c != Labels.unassigned]
    assert len(assigned_classes) > 0

    plt.figure()
    result = mass_spectrum.plot_dbe_vs_carbon_number(assigned_classes[0])
    ax, abun_perc = result
    assert ax is not None
    assert abun_perc > 0
    assert assigned_classes[0] in ax.get_title()
    assert ax.get_xlabel() == "Carbon number"
    assert ax.get_ylabel() == "DBE"
    plt.close()


def test_plot_dbe_vs_carbon_number_all_classes(classified_mass_spectrum):
    """Test DBE vs Carbon Number plot for all assigned classes (new functionality)."""
    plt.figure()
    ax = classified_mass_spectrum.plot_dbe_vs_carbon_number()
    assert ax is not None
    assert "All Assigned Classes" in ax.get_title()
    assert ax.get_xlabel() == "Carbon number"
    assert ax.get_ylabel() == "DBE"
    plt.close()


def test_plot_dbe_vs_carbon_number_log_abundance(classified_mass_spectrum):
    """Test DBE vs Carbon Number plot with log10 abundance scaling."""
    plt.figure()
    ax = classified_mass_spectrum.plot_dbe_vs_carbon_number(log_abundance=True)
    assert ax is not None
    cbar = ax.figure.axes[-1]
    assert "log" in cbar.get_ylabel().lower() or "log" in cbar.get_ylabel()
    plt.close()